#include <algorithm>
#include <bitset>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <random>
#include <sstream>


#include "functions_for_python.h"
#include "input_structure.h"
#include "output_structure.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>

std::mutex file_mutex;
std::queue<std::string> output_queue;
std::mutex queue_mutex;
std::condition_variable cv;
bool finished = false;

void print_output(std::ostringstream& ss, const output_structure& output) {
    ss << "[";
    if (!output.convergence.empty()) {
        std::copy(output.convergence.begin(), output.convergence.end() - 1,
                  std::ostream_iterator<double>(ss, ","));
        ss << output.convergence.back();
    }
    ss << "],[";
    if (!output.result.empty()) {
        std::copy(output.result.begin(), output.result.end() - 1,
                  std::ostream_iterator<uint32_t>(ss, ","));
        ss << static_cast<uint32_t>(output.result.back());
    }
    ss << "]";
}

void output_thread(std::ofstream& file) {
    while (true) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        cv.wait(lock, [] { return !output_queue.empty() || finished; });

        if (finished && output_queue.empty()) {
            break;
        }

        if (!output_queue.empty()) {
            std::string output = output_queue.front();
            output_queue.pop();
            lock.unlock();

            std::lock_guard<std::mutex> file_lock(file_mutex);
            file << output << std::endl;
            file.flush(); // Zajistí okamžitý zápis do souboru
        }
    }
}

void process_iteration(uint64_t iteration, int64_t dim, int64_t i, const input_structure& input_structure) {
    std::ostringstream ss;
    ss << iteration << "," << dim << "," << i << ",";

    output_structure output;

    output = GeneticAlgorithm_MaxOne(input_structure);
    print_output(ss, output);
    ss << ",";

    output = GeneticAlgorithm_LeadingOne(input_structure);
    print_output(ss, output);
    ss << ",";

    output = GeneticAlgorithm_TradingOne(input_structure);
    print_output(ss, output);
    ss << ",";

    output = GeneticAlgorithm_MaxZero(input_structure);
    print_output(ss, output);

    std::unique_lock<std::mutex> lock(queue_mutex);
    output_queue.push(ss.str());
    lock.unlock();
    cv.notify_one();
}

int main() {
    input_structure input_structure{0.02, 0.005, 0.085, .1, .5, .75};
    int64_t dimensions[] = {10, 30, 100, 500};

    std::ofstream result_file("result.csv", std::ios::out | std::ios::trunc);
    if (!result_file.is_open()) {
        std::cerr << "Nelze otevřít soubor result.csv pro zápis." << std::endl;
        return 1;
    }

    result_file << "id,dimension,iteration,convergence_MaxOne,bit_result_MaxOne,"
                << "convergence_LeadingOne,bit_result_LeadingOne,convergence_TradingOne,"
                << "bit_result_TradingOne,convergence_MaxZero,bit_result_MaxZero" << std::endl;

    std::thread output_t(output_thread, std::ref(result_file));

    uint64_t iteration = 0;
    std::vector<std::thread> threads;

    for (auto & dim : dimensions) {
        input_structure.initial_population_size = 50 * 2 * dim/10;
        input_structure.population_size = 50 * dim/10;
        input_structure.problem_dimension = dim;
        input_structure.fitness_rating_count = dim * 100;

        for (int64_t i = 0; i < 100; i++) {
            threads.emplace_back(process_iteration, iteration, dim, i, input_structure);
            iteration++;

            if (threads.size() >= std::thread::hardware_concurrency()) {
                for (auto& t : threads) {
                    t.join();
                }
                threads.clear();
            }
        }
    }

    for (auto& t : threads) {
        t.join();
    }

    finished = true;
    cv.notify_one();
    output_t.join();

    result_file.close();

    return 0;
}