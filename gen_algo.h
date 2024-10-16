/*
 * Created by kureii on 10/16/24.
 */

#pragma once
#include <cstdint>
#include <functional>
#include <random>

#include "gen_structure.h"

#define MINIMUM_PARENT_COUNT 2

template <std::size_t SIZE>
class GeneticAlgorithm {
 public:
  using Chromosome = std::array<std::uint64_t, SIZE>;
  using FitnessFunction = std::function<std::uint64_t(const Chromosome&)>;
  /**
   *
   * @param mutation_probability Probability of changing one bit in the child.
   * The range is 0-1
   * @param elitism_roulette_percent Percentage of random selection for elitism.
   * The range is 0-1
   * @param elitism_selection_percent Percentage of the best selection for
   * elitism. The range is 0-1
   * @param crossing_roulette_probability The probability that a parent will be
   * chosen completely at random. The range is 0-1
   * @param first_parent_range The percentage of the best result from which the
   * first parent is randomly selected. If the value is too low (lower than two
   * elements), the selected element is always in random range. Range 0-1.
   * @param second_parent_range The percentage of the best result from which the
   * second parent is randomly selected. If the value is too low (lower than two
   * elements), the selected element is always in random range. Range 0-1.
   * @param init_population_size First generation size.
   * @param generation_size Second to last generation size.
   * @param dimension_size Dimensionality of the problem.
   * @param stop_limit Limit the number of times a fitness function can be
   * performed.
   * @param fitness_function Function for optimization, input and output of
   * uin64_t.
   */
  GeneticAlgorithm(double mutation_probability, double elitism_roulette_percent,
      double elitism_selection_percent, double crossing_roulette_probability,
      double first_parent_range, double second_parent_range,
      std::uint64_t init_population_size, std::uint64_t generation_size,
      std::uint32_t dimension_size, std::uint32_t stop_limit,
      const FitnessFunction& fitness_function);
  void Compute();
  [[nodiscard]] std::vector<uint64_t> GetConvergence() const;

 private:
  void GenerateInitialPopulation();
  void Crossing();
  std::vector<FitnessFunctionResults_t<SIZE>> MakeChildren(Chromosome A, Chromosome B);

  double mutation_probability_;
  double elitism_roulette_percent_;
  double elitism_selection_percent_;
  double crossing_roulette_probability_;
  double first_parent_range_;
  double second_parent_range_;
  uint64_t init_population_size_;
  uint64_t generation_size_;
  uint64_t dimension_size_;
  uint64_t stop_limit_;
  uint64_t actual_fitness_count_ = 0;
  size_t array_size_;

  std::vector<FitnessFunctionResults_t<SIZE>> old_fitness_function_results_;
  std::vector<FitnessFunctionResults_t<SIZE>> new_fitness_function_results_;

  std::vector<uint64_t> convergence_;
  std::mt19937_64 generator_;
  FitnessFunction fitness_function_;
};

#include "gen_algo.tpp"