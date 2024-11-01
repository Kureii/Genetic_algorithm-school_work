/*
 * Created by kureii on 10/16/24.
 */

#pragma once
#include <cstdint>
#include <functional>
#include <random>

#include "fitness_function_results_struct.h"
#include "input_structure.h"
#include "chromosome_array.h"
#define MINIMUM_PARENT_COUNT 2

class GeneticAlgorithm {
 public:
  using Chromosome = ChromosomeArray<uint64_t>;
  using FitnessFunction = std::function<uint64_t(const Chromosome&)>;
  /**
   * @param input_structure
   * @param chromosome_size Size of chromosome
   * @param fitness_function Function for optimization, input and output of
   * uin64_t.
   * @param real_numbers
   */
  GeneticAlgorithm(const input_structure_t& input_structure, uint64_t chromosome_size,
      FitnessFunction  fitness_function,
      bool real_numbers = false );
  void Compute();
  [[nodiscard]] std::vector<uint64_t> GetConvergence() const;
  [[nodiscard]] std::vector<uint64_t> GetResult() const;

 private:
  void CheckInit() const;
  void GenerateInitialPopulation();
  void GenerateInitialPopulationBit(uint64_t chromosome_array_size, uint64_t all_bit_size);
  void GenerateInitialPopulationRealIEEE();
  void CrossingBit(uint64_t chromosome_array_size, uint64_t all_bit_size);
  bool IsFixedPointGenValid(uint64_t gen);
  void MakeValidChromosomeBinaryCodedDecimal(Chromosome& A);
  void CrossingReal();
  void Elitism();
  void ParentSelection(Chromosome& first_parent, Chromosome& second_parent);
  std::pair<ChromosomeArray<double>, ChromosomeArray<double>>ParentSelectionReal();
  void Mutation (FitnessFunctionResults& first_child, FitnessFunctionResults& second_child);
  void MutationReal(FitnessFunctionResults& first_child, FitnessFunctionResults& second_child);
  void ComputeFitnessFunction();
  std::vector<FitnessFunctionResults> MakeChildren(Chromosome& A, Chromosome& B);
  std::vector<FitnessFunctionResults> MakeChildrenReal();


  bool real_numbers_;

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
  ChromosomeArray<double> best_from_generation_double_;
  std::optional<mapping_structure_t> mapping_structure_;

  std::vector<FitnessFunctionResults> old_fitness_function_results_;
  std::vector<FitnessFunctionResults> new_fitness_function_results_;

  std::vector<uint64_t> convergence_;
  std::mt19937_64 generator_;
  std::uniform_real_distribution<> zero_one_distribution_ = std::uniform_real_distribution(0.0, 1.0);
  std::uniform_int_distribution<> mutate_full_range_distribution_ = std::uniform_int_distribution(0, 64);
  std::uniform_real_distribution<> ieee_distribution_;
  std::uniform_int_distribution<uint64_t> roulette_distribution_;
  FitnessFunction fitness_function_;
};