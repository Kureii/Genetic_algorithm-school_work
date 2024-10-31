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
      [input_structure, mapping_structure](
          const ChromosomeArray<uint64_t>& chromosome) {
        std::vector<double> variables = RealConvertor::ConvertChromosomeToReal(
            chromosome, mapping_structure, input_structure.problem_dimension);
        double sum_squares = 0.0;
        for (double x : variables) {
          sum_squares += x * x;
        }
        // Protože maximalizujeme, vracíme zápornou hodnotu
        return -sum_squares;
      };

  auto ga = GeneticAlgorithm(
      input_structure, chromosome_size, sphere_fitness_function);
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
        // Protože maximalizujeme, vracíme fitness přímo
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