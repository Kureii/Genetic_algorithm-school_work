/*
 * Created by kureii on 10/18/24.
*/

#include <bitset>
#include <cstdint>
#include <random>

#include "../gen_algo.h"
#include "functions_for_python.h"

output_structure_t GeneticAlgorithm_MaxOne(input_structure_t input_structure) {
  const std::size_t chromozome_size = (input_structure.problem_dimension + 63) / 64;
  auto max_one_fitness_function =
      [](const myArrayEncapsulation<uint64_t>&chromosome){
        std::uint64_t result = 0;
        for (const auto& gene : chromosome) {
          result += std::bitset<64>(gene).count();
        }
        return result;
  };

  auto ga = GeneticAlgorithm(input_structure, chromozome_size, max_one_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_LeadingOne(input_structure_t input_structure) {
  const std::size_t chromozome_size = (input_structure.problem_dimension + 63) / 64;
  auto leading_ones_fitness_function =
      [](const myArrayEncapsulation<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        bool found_zero = false;

        for (const auto& gene : chromosome) {
          for (int i = 0; i < 64; ++i) {
            if (!found_zero) {
              if ((gene >> i) & 0b1) {
                ++result;
              } else {
                found_zero = true;
                break;
              }
            }
          }
          if (found_zero) break;
        }

        return result;
  };
  auto ga = GeneticAlgorithm(input_structure, chromozome_size, leading_ones_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_TradingOne(input_structure_t input_structure) {
  const std::size_t chromosome_size = (input_structure.problem_dimension + 63) / 64;
  auto trailing_ones_fitness_function =
      [](const myArrayEncapsulation<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        bool found_zero = false;

        for (auto it = chromosome.rbegin(); it != chromosome.rend(); ++it) {
          const auto& gene = *it;
          for (int i = 63; i >= 0; --i) {
            if (!found_zero) {
              if ((gene >> i) & 0b1) {
                ++result;
              } else {
                found_zero = true;
                break;
              }
            }
          }
          if (found_zero) break;
        }

        return result;
  };

  auto ga = GeneticAlgorithm(input_structure, chromosome_size, trailing_ones_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_MaxZero(input_structure_t input_structure) {
  const std::size_t chromosome_size = (input_structure.problem_dimension + 63) / 64;
  auto max_zero_fitness_function =
      [](const myArrayEncapsulation<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        for (const auto& gene : chromosome) {
          result += 64 - std::bitset<64>(gene).count();
        }
        return result;
  };

  auto ga = GeneticAlgorithm(input_structure, chromosome_size, max_zero_fitness_function);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();
  return output;
}
