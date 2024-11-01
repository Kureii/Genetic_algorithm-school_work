/*
 * Created by kureii on 10/18/24.
*/

#pragma once
#include <cstdint>
#include <optional>

#include "real_structures.h"

struct input_structure {
  double mutation_probability;
  double elitism_random_select;
  double elitism_best_select;
  double crossing_full_random_probability;
  double first_parent_select_range;
  double second_parent_select_range;
  uint64_t initial_population_size;
  uint64_t population_size;
  uint64_t problem_dimension;
  uint64_t fitness_rating_count;
  std::optional<mapping_structure_t> mapping_structure;

};

using input_structure_t = input_structure;