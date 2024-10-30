/*
 * Created by kureii on 10/18/24.
*/

#pragma once
#include <cstdint>

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
};

using input_structure_t = struct input_structure;