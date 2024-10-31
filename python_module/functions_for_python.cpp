/*
 * Created by kureii on 10/18/24.
 */

#include "functions_for_python.h"

#include <bitset>
#include <cstdint>
#include <random>

#include "../gen_algo.h"
#include "../real_convertor.h"

output_structure_t GeneticAlgorithm_MaxOne(input_structure_t input_structure) {
  const std::size_t chromozome_size =
      (input_structure.problem_dimension + 63) / 64;
  auto max_one_fitness_function =
      [](const ChromosomeArray<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        for (size_t i = 0; i < chromosome.size(); i++) {
          uint64_t bits;
          auto chromosome_part = chromosome[i];
          if (i != chromosome.size() - 1) {
            bits = 64;
          } else {
            bits = chromosome.chromosome_length() % 64;
          }
          for (size_t j = 0; j < bits; ++j) {
            result += chromosome_part & 0b1;
            chromosome_part >>= 1;
          }
        }
        return result;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromozome_size, max_one_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_LeadingOne(
    input_structure_t input_structure) {
  const std::size_t chromozome_size =
      (input_structure.problem_dimension + 63) / 64;
  auto leading_ones_fitness_function =
      [](const ChromosomeArray<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        for (size_t i = chromosome.size(); i > 0; --i) {
          int bits;
          auto chromosome_part = chromosome[i-1];
          if (i != chromosome.size()) {
            bits = 64;
          } else {
            bits = chromosome.chromosome_length() % 64;
          }
          for (uint64_t j = bits-1; j > 0; --j) {
            uint64_t mask = static_cast<uint64_t>(0b1) << j;
            if (chromosome_part & mask) {
              result++;
            }else {
              return result;
            }
          }
        }
        return result;
      };
  auto ga = GeneticAlgorithm(
      input_structure, chromozome_size, leading_ones_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_TradingOne(
    input_structure_t input_structure) {
  const std::size_t chromosome_size =
      (input_structure.problem_dimension + 63) / 64;
  auto trailing_ones_fitness_function =
      [](const ChromosomeArray<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        for (size_t i = 0; i < chromosome.size(); i++) {
          uint64_t bits;
          auto chromosome_part = chromosome[i];
          if (i != chromosome.size() - 1) {
            bits = 64;
          } else {
            bits = chromosome.chromosome_length() % 64;
          }
          for (size_t j = 0; j < bits; ++j) {
            if(chromosome_part & 0b1) {
              result++;
            }else {
              return result;
            }

            chromosome_part >>= 1;
          }
        }
        return result;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, trailing_ones_fitness_function);
  ga.Compute();

  output_structure_t output;

  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_MaxZero(input_structure_t input_structure) {
  const std::size_t chromosome_size =
      (input_structure.problem_dimension + 63) / 64;
  auto max_zero_fitness_function =
      [](const ChromosomeArray<uint64_t>& chromosome) {
        std::uint64_t result = 0;
        for (size_t i = 0; i < chromosome.size(); i++) {
          uint64_t bits;
          auto chromosome_part = chromosome[i];
          if (i != chromosome.size() - 1) {
            bits = 64;
          } else {
            bits = chromosome.chromosome_length() % 64;
          }
          for (size_t j = 0; j < bits; ++j) {
            result += (chromosome_part & 0b1) ^ 0b1;
            chromosome_part >>= 1;
          }
        }
        return result;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, max_zero_fitness_function);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();
  return output;
}
