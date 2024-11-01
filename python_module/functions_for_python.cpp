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

output_structure_t GeneticAlgorithm_Sphere(
    input_structure_t input_structure, mapping_structure_t mapping_structure) {
  const std::size_t chromosome_size = input_structure.problem_dimension;
  const std::size_t inpu = 64 * input_structure.problem_dimension;

  auto sphere_fitness_function =
      [mapping_structure](const ChromosomeArray<uint64_t>& chromosome) {
        uint64_t sum_squares = 0;

        switch (mapping_structure.mapping_method) {
          using enum mapping_method_t;
          case BIT_AS_DOUBLE: {
            for (const auto& gene : chromosome) {
              uint64_t x = gene;
              sum_squares += x * x;
            }
            break;
          }
          case FIXED_POINT: {
            const uint64_t scaling_factor = 1ULL << mapping_structure.fractional_bits;
            const uint64_t mask = ~0ULL;

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene & mask;
              int64_t int_value = static_cast<int64_t>(bits_value);

              if (bits_value & (1ULL << 63)) {
                // Negative number, sign extension
                int_value |= ~((1ULL << 64) - 1);
              }

              int64_t x = int_value / scaling_factor;
              sum_squares += static_cast<uint64_t>(x * x);
            }
            break;
          }
          case BINARY_CODED_DECIMAL: {
            int num_digits = 16; // 64 bits / 4 bits per digit

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene;
              uint64_t x = 0;

              for (int d = 0; d < num_digits; ++d) {
                uint64_t digit_bits = (bits_value >> (4 * d)) & 0xF;
                int digit = static_cast<int>(digit_bits);
                if (digit > 9) digit = 9; // Limit to digits 0-9
                x += digit * static_cast<uint64_t>(pow(10, num_digits - d - 1));
              }

              sum_squares += x * x;
            }
            break;
          }
          case MAPPED_RANGE: {
              double min_value = mapping_structure.min_value;
              double max_value = mapping_structure.max_value;
              double range = max_value - min_value;
              uint64_t max_int_value = ~0ULL;

              for (const auto& gene : chromosome) {
                uint64_t bits_value = gene;
                double x = min_value + (static_cast<double>(bits_value) / max_int_value) * range;
                sum_squares += x * x;
              }
              break;
          }
        }

        // Since we are maximizing in the GA, invert the sum_squares
        uint64_t fitness = std::numeric_limits<uint64_t>::max() - sum_squares;
        return fitness;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, sphere_fitness_function, true);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_Schwefel(
    input_structure_t input_structure, mapping_structure_t mapping_structure) {
  const std::size_t chromosome_size = input_structure.problem_dimension;

  auto schwefel_fitness_function =
      [input_structure, mapping_structure](const ChromosomeArray<uint64_t>& chromosome) {
        double sum = 0.0;

        switch (mapping_structure.mapping_method) {
          using enum mapping_method_t;
          case BIT_AS_DOUBLE: {
            for (const auto& gene : chromosome) {
              double x = std::bit_cast<double>(gene);
              sum += x * sin(sqrt(fabs(x)));
            }
            break;
          }
          case FIXED_POINT: {
            const auto scaling_factor =
                static_cast<double>(1ULL << mapping_structure.fractional_bits);
            const uint64_t mask = ~0ULL;

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene & mask;
              int64_t int_value = static_cast<int64_t>(bits_value);

              if (bits_value & (1ULL << 63)) {
                // Negative number, sign extension
                int_value |= ~((1ULL << 64) - 1);
              }

              double x = static_cast<double>(int_value) / scaling_factor;
              sum += x * sin(sqrt(fabs(x)));
            }
            break;
          }
          case BINARY_CODED_DECIMAL: {
            int num_digits = 16; // 64 bits / 4 bits per digit

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene;
              double x = 0.0;

              for (int d = 0; d < num_digits; ++d) {
                uint64_t digit_bits = (bits_value >> (4 * d)) & 0xF;
                int digit = static_cast<int>(digit_bits);
                if (digit > 9) digit = 9; // Limit to digits 0-9
                x += digit * pow(10, num_digits - d - 1);
              }

              sum += x * sin(sqrt(fabs(x)));
            }
            break;
          }
          case MAPPED_RANGE: {
            double min_value = mapping_structure.min_value;
            double max_value = mapping_structure.max_value;
            double range = max_value - min_value;
            uint64_t max_int_value = ~0ULL;

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene;
              double x = min_value + (static_cast<double>(bits_value) / max_int_value) * range;
              sum += x * sin(sqrt(fabs(x)));
            }
            break;
          }
        }

        double fitness_double = 418.9829 * input_structure.problem_dimension - sum;

        // Scale the fitness value to preserve precision and convert to uint64_t
        uint64_t fitness_uint64 = static_cast<uint64_t>(fitness_double * 1e6);

        return fitness_uint64;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, schwefel_fitness_function, true);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}


output_structure_t GeneticAlgorithm_Rosenbrock(
    input_structure_t input_structure, mapping_structure_t mapping_structure) {
  const std::size_t chromosome_size = input_structure.problem_dimension;
  const std::size_t total_bits = 64 * input_structure.problem_dimension;

  auto rosenbrock_fitness_function =
      [input_structure, mapping_structure](const ChromosomeArray<uint64_t>& chromosome) {
        uint64_t sum = 0;
        std::vector<int64_t> variables;

        switch (mapping_structure.mapping_method) {
          using enum mapping_method_t;
          case BIT_AS_DOUBLE: {
            for (const auto& gene : chromosome) {
              int64_t x = static_cast<int64_t>(gene);
              variables.push_back(x);
            }
            break;
          }
          case FIXED_POINT: {
            const uint64_t scaling_factor = 1ULL << mapping_structure.fractional_bits;
            const uint64_t mask = ~0ULL;

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene & mask;
              int64_t int_value = static_cast<int64_t>(bits_value);

              if (bits_value & (1ULL << 63)) {
                // Negative number, sign extension
                int_value |= ~((1ULL << 64) - 1);
              }

              int64_t x = int_value / scaling_factor;
              variables.push_back(x);
            }
            break;
          }
          case BINARY_CODED_DECIMAL: {
            int num_digits = 16; // 64 bits / 4 bits per digit

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene;
              int64_t x = 0;

              for (int d = 0; d < num_digits; ++d) {
                uint64_t digit_bits = (bits_value >> (4 * d)) & 0xF;
                int digit = static_cast<int>(digit_bits);
                if (digit > 9) digit = 9; // Limit to digits 0-9
                x += digit * static_cast<int64_t>(pow(10, num_digits - d - 1));
              }

              variables.push_back(x);
            }
            break;
          }
          case MAPPED_RANGE: {
            double min_value = mapping_structure.min_value;
            double max_value = mapping_structure.max_value;
            double range = max_value - min_value;
            uint64_t max_int_value = ~0ULL;

            for (const auto& gene : chromosome) {
              uint64_t bits_value = gene;
              double x = min_value + (static_cast<double>(bits_value) / max_int_value) * range;
              variables.push_back(x);
            }
            break;
          }
        }

        for (std::size_t i = 0; i < variables.size() - 1; ++i) {
          int64_t xi = variables[i];
          int64_t x_next = variables[i + 1];
          sum += 100 * (x_next - xi * xi) * (x_next - xi * xi) + (xi - 1) * (xi - 1);
        }

        uint64_t fitness = std::numeric_limits<uint64_t>::max() - sum;
        return fitness;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, rosenbrock_fitness_function, true);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}
