/*
 * Created by kureii on 10/16/24.
 */

#pragma once
#include <cstdint>
#include <functional>
#include <random>

#include "fitness_function_results_struct.h"
#include "input_structure.h"
#define MINIMUM_PARENT_COUNT 2

class GeneticAlgorithm {
 public:
  using Chromosome = myArrayEncapsulation<uint64_t>;
  using FitnessFunction = std::function<uint64_t(const Chromosome&)>;
  /**
   * @param chromosome_size Size of chromosome
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
  GeneticAlgorithm(const input_structure_t& input_structure, uint64_t chromosome_size,
      const FitnessFunction& fitness_function);
  void Compute();
  [[nodiscard]] std::vector<uint64_t> GetConvergence() const;
  [[nodiscard]] std::vector<uint8_t> GetResult() const;

 private:
  void GenerateInitialPopulation();
  void Crossing();
  void Elitism();
  void ParentSelection(Chromosome& first_parent, Chromosome& second_parent);
  void Mutation (FitnessFunctionResults& first_child, FitnessFunctionResults& second_child);
  std::vector<FitnessFunctionResults> MakeChildren(Chromosome& A, Chromosome& B);



  double mutation_probability_;
  double elitism_roulette_percent_;
  double elitism_selection_percent_;
  double crossing_roulette_probability_;
  double first_parent_range_;
  double second_parent_range_;
  uint64_t chromosome_size_;
  uint64_t init_population_size_;
  uint64_t generation_size_;
  uint64_t dimension_size_;
  uint64_t stop_limit_;
  uint64_t actual_fitness_count_ = 0;
  Chromosome best_from_generation_;

  std::vector<FitnessFunctionResults> old_fitness_function_results_;
  std::vector<FitnessFunctionResults> new_fitness_function_results_;

  std::vector<uint64_t> convergence_;
  std::mt19937_64 generator_;
  std::uniform_real_distribution<> zero_one_distribution_ = std::uniform_real_distribution(0.0, 1.0);
  std::uniform_int_distribution<uint64_t> roulette_distribution_;
  FitnessFunction fitness_function_;
};