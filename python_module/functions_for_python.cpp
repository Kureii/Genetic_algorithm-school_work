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
  std::size_t bits_per_variable = mapping_structure.bits_per_variable;
  std::size_t total_bits =
      bits_per_variable * input_structure.problem_dimension;
  const std::size_t chromosome_size = (total_bits + 63) / 64;

  auto sphere_fitness_function =
      [](
          const ChromosomeArray<uint64_t>& chromosome) {

        double sum_squares = 0.0;
        auto mapping = chromosome.mapping_structure();
        switch (chromosome.mapping_structure().mapping_method) {
          using enum mapping_method_t;
          case BIT_AS_DOUBLE:
              sum_squares += std::bit_cast<double>(chromosome[0]) * std::bit_cast<double>(chromosome[0]);
            break;
          case FIXED_POINT: {
            const auto scaling_factor =
                static_cast<double>(1ULL << mapping.fractional_bits);

            const uint64_t mask = (mapping.bits_per_variable < 64) ? ((1ULL << mapping.bits_per_variable) - 1) : ~0ULL;

            for (const auto& x_tmp : chromosome) {
                uint64_t bits_value = x_tmp & mask;

                auto int_value = static_cast<int64_t>(bits_value);

                if (mapping.bits_per_variable < 64) {
                  if (uint64_t sign_bit =
                          bits_value &
                          (1ULL << (mapping.bits_per_variable - 1))) {
                    // Záporné číslo, rozšíření znaménka
                    uint64_t sign_extension =
                        ~((1ULL << mapping.bits_per_variable) - 1);
                    int_value |= sign_extension;
                  }
                }

                const double x = static_cast<double>(int_value) / scaling_factor;

                sum_squares += x * x;
            }
            break;
          }
          case BINARY_CODED_DECIMAL: {
            int num_digits = mapping.bits_per_variable / 4;

            for (const auto& x_tmp : chromosome) {
                uint64_t bits_value = x_tmp;

                double x = 0.0;
                for (int d = 0; d < num_digits; ++d) {
                    uint64_t digit_bits = (bits_value >> (4 * d)) & 0xF;
                    int digit = static_cast<int>(digit_bits);
                    if (digit > 9) digit = 9; // Omezení na číslice 0-9
                    x += digit * pow(10, num_digits - d - 1); // Nejvýznamnější číslice první
                }

                sum_squares += x * x;
            }
            break;
          }
          case MAPPED_RANGE: {
            int bits_per_variable = mapping.bits_per_variable;
            double min_value = mapping.min_value;
            double max_value = mapping.max_value;
            double range = max_value - min_value;
            uint64_t max_int_value = (1ULL << bits_per_variable) - 1;

            for (const auto& x_tmp : chromosome) {
                uint64_t bits_value = x_tmp & max_int_value;

                double x = min_value + (static_cast<double>(bits_value) / max_int_value) * range;

                sum_squares += x * x;
            }
            break;
          }
        }

        return -sum_squares;
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
  std::size_t bits_per_variable = mapping_structure.bits_per_variable;
  std::size_t total_bits =
      bits_per_variable * input_structure.problem_dimension;
  const std::size_t chromosome_size = (total_bits + 63) / 64;

  auto schwefel_fitness_function =
      [input_structure, mapping_structure](
          const ChromosomeArray<uint64_t>& chromosome) {
        std::vector<double> variables = RealConvertor::ConvertChromosomeToReal(
            chromosome, mapping_structure, input_structure.problem_dimension);
        double sum = 0.0;
        for (double x : variables) {
          sum += x * sin(sqrt(fabs(x)));
        }
        double fitness = 418.9829 * input_structure.problem_dimension - sum;
        return fitness;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, schwefel_fitness_function);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}

output_structure_t GeneticAlgorithm_Rosenbrock(
    input_structure_t input_structure, mapping_structure_t mapping_structure) {
  std::size_t bits_per_variable = mapping_structure.bits_per_variable;
  std::size_t total_bits =
      bits_per_variable * input_structure.problem_dimension;
  const std::size_t chromosome_size = (total_bits + 63) / 64;

  auto rosenbrock_fitness_function =
      [input_structure, mapping_structure](
          const ChromosomeArray<uint64_t>& chromosome) {
        std::vector<double> variables = RealConvertor::ConvertChromosomeToReal(
            chromosome, mapping_structure, input_structure.problem_dimension);
        double sum = 0.0;
        for (std::size_t i = 0; i < variables.size() - 1; ++i) {
          double xi = variables[i];
          double x_next = variables[i + 1];
          sum += 100.0 * pow((x_next - xi * xi), 2) + pow((xi - 1), 2);
        }
        // Protože maximalizujeme, vracíme zápornou hodnotu
        return -sum;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, rosenbrock_fitness_function);
  ga.Compute();

  output_structure_t output;
  output.convergence = ga.GetConvergence();
  output.result = ga.GetResult();

  return output;
}