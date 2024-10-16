#include <array>
#include <bitset>
#include <cstdint>
#include <iostream>

#include "gen_algo.h"

#define MUTATION_RATE 0.15
#define ELITE_ROULETTE 0.1
#define ELITE_SELECT 0.3
#define CROSSING_ROULETTE 0.1
#define FIRST_PARENT 0.5
#define SECOND_PARENT 0.7
#define INIT_POPULATION_SIZE 25
#define GENERATION_SIZE 15
#define DIMENSION 100
#define STOP_SIZE (DIMENSION * 100)

int main() {
  constexpr std::size_t CHROMOSOME_SIZE = (DIMENSION + 63) / 64;
#ifdef DEBUG
  std::cout << "CHROMOSOME_SIZE = " << CHROMOSOME_SIZE << std::endl;
#endif
  auto fitness_function =
      [](const std::array<std::uint64_t, CHROMOSOME_SIZE>& chromosome){
    std::uint64_t result = 0;
    for (const auto& gene : chromosome) {
      result += std::bitset<64>(gene).count();
    }
    return result;
  };
  for (std::size_t i = 0; i < 100; i++) {
    auto ga = GeneticAlgorithm<CHROMOSOME_SIZE>(MUTATION_RATE, ELITE_ROULETTE,
        ELITE_SELECT, CROSSING_ROULETTE, FIRST_PARENT, SECOND_PARENT,
        INIT_POPULATION_SIZE, GENERATION_SIZE, DIMENSION, STOP_SIZE,
        fitness_function);

    ga.Compute();

    auto convergence = ga.GetConvergence();


    std::cout << convergence.back() << ", ";
  }

  return 0;
}
