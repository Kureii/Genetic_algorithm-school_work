#include <bitset>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <ranges>

#include "gen_algo.h"

#define MUTATION_RATE 0.005
#define ELITE_ROULETTE 0.02
#define ELITE_SELECT 0.2
#define CROSSING_ROULETTE 0.0
#define FIRST_PARENT 0.5
#define SECOND_PARENT 0.75
#define INIT_POPULATION_SIZE 10
#define GENERATION_SIZE 50
#define DIMENSION 10
#define STOP_SIZE (DIMENSION * 10)

int main() {
  constexpr std::size_t CHROMOSOME_SIZE = (DIMENSION + 63) / 64;
#ifdef DEBUG
  std::cout << "CHROMOSOME_SIZE = " << CHROMOSOME_SIZE << std::endl;
#endif
  auto fitness_function =
      [](const myArrayEncapsulation<uint64_t>&chromosome){
    std::uint64_t result = 0;
    for (const auto& gene : chromosome) {
      result += std::bitset<64>(gene).count();
    }
    return result;
  };
  for (std::size_t i = 0; i < 1; i++) {
    auto ga = GeneticAlgorithm(CHROMOSOME_SIZE, MUTATION_RATE, ELITE_ROULETTE,
        ELITE_SELECT, CROSSING_ROULETTE, FIRST_PARENT, SECOND_PARENT,
        INIT_POPULATION_SIZE, GENERATION_SIZE, DIMENSION, STOP_SIZE,
        fitness_function);

    ga.Compute();

    auto convergence = ga.GetConvergence();
    auto result = ga.GetResult();

    std::ranges::for_each(result, [&](auto gene) {
      std::cout << static_cast<uint32_t>(gene) << ", ";
    });

    std::cout << convergence.back() << ", ";
  }

  return 0;
}
