/*
 * Created by kureii on 10/16/24.
 */
#pragma once

#include <random>
#include <stdexcept>

#include "gen_structure.h"

template <std::size_t SIZE>
GeneticAlgorithm<SIZE>::GeneticAlgorithm(double mutation_probability,
    double elitism_roulette_percent, double elitism_selection_percent,
    double crossing_roulette_probability, double first_parent_range,
    double second_parent_range, std::uint64_t init_population_size,
    std::uint64_t generation_size, std::uint32_t dimension_size,
    std::uint32_t stop_limit, const FitnessFunction& fitness_function)
    : mutation_probability_(mutation_probability),
      elitism_roulette_percent_(elitism_roulette_percent),
      elitism_selection_percent_(elitism_selection_percent),
      crossing_roulette_probability_(crossing_roulette_probability),
      first_parent_range_(first_parent_range),
      second_parent_range_(second_parent_range),
      init_population_size_(init_population_size),
      generation_size_(generation_size),
      dimension_size_(dimension_size),
      stop_limit_(stop_limit),
      fitness_function_(fitness_function) {
  if (SIZE == 0) {
    throw std::invalid_argument("Chromosome size cannot be zero");
  }

  if (dimension_size_ > SIZE * 64) {
    throw std::invalid_argument(
        "Dimension size is too large for the given chromosome size");
  }

  if (elitism_roulette_percent_ + elitism_selection_percent_ > 0.75) {
    throw std::invalid_argument("Elitism percent cannot be more than 0.75");
  }

  if (mutation_probability_ > 1.0 || elitism_roulette_percent_ > 1.0 ||
      elitism_selection_percent_ > 1.0 ||
      crossing_roulette_probability_ > 1.0 || first_parent_range_ > 1.0 ||
      second_parent_range_ > 1.0 || mutation_probability_ < 0.0 ||
      elitism_roulette_percent_ < 0.0 || elitism_selection_percent_ < 0.0 ||
      crossing_roulette_probability_ < 0.0 || first_parent_range_ < 0.0 ||
      second_parent_range_ < 0.0) {
    throw std::invalid_argument("All values of type double must be in the range 0-1");
  }

  old_fitness_function_results_.resize(init_population_size_);
  new_fitness_function_results_.resize(init_population_size_);

  std::random_device random_device;
  generator_ = std::mt19937_64(random_device());
}

template <std::size_t SIZE>
void GeneticAlgorithm<SIZE>::GenerateInitialPopulation() {
  old_fitness_function_results_.clear();
  std::vector<Chromosome> old_chromosomes;
  old_chromosomes.resize(init_population_size_);
  for (std::size_t i = 0; i < init_population_size_; i++) {
    old_chromosomes[i] = std::array<std::uint64_t, SIZE>();
  }
  std::uniform_int_distribution<uint64_t> distribution(0, UINT64_MAX);
#ifdef DEBUG
  std::cout << "Initial generation" << std::endl;
#endif
  for (std::size_t i = 0; i < init_population_size_; ++i) {
    for (std::size_t j = 0; j < SIZE; ++j) {
      old_chromosomes[i][j] = distribution(generator_);
    }
    if (SIZE > 1) {
      uint64_t iterator = dimension_size_ - (SIZE - 1) * 64;
      uint64_t mask = 0;
      for (std::size_t j = 0; j < iterator; ++j) {
        mask <<= 1;
        mask++;
      }
      old_chromosomes[i][SIZE - 1] &= mask;
    } else {
      uint64_t mask = 0;
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mask <<= 1;
        mask++;
      }
      old_chromosomes[i][0] &= mask;
    }
#ifdef DEBUG
    for (std::size_t j = 0; j < SIZE; ++j) {
      std::cout << std::bitset<64>(old_chromosomes[i][SIZE - 1 - j]) << " ";
    }
    std::cout << std::endl;
#endif
    if (!old_fitness_function_results_.empty()) {
      old_fitness_function_results_.clear();
    }
  }
  for (uint64_t i = 0; i < old_chromosomes.size(); i++) {
    old_fitness_function_results_.emplace_back(
        FitnessFunctionResults_t<SIZE>(old_chromosomes[i]));
  }
}

template <std::size_t SIZE>
void GeneticAlgorithm<SIZE>::Compute() {
  convergence_.clear();
  GenerateInitialPopulation();

  while (actual_fitness_count_ < stop_limit_) {
    auto emergency_brake = true;
    for (std::size_t i = 0; i < old_fitness_function_results_.size(); ++i) {
      if (!old_fitness_function_results_[i].fitness.has_value()) {
        emergency_brake = false;
        break;
      }
    }
    if (emergency_brake) {
#ifdef DEBUG
        std::cout << "emergency break" << std::endl;
#endif
      break;
    }
#ifdef DEBUG
    for (auto &input : old_fitness_function_results_) {
      if (input.input[0] > 0b111) {
        std::cout << std::bitset<64>(input.input[0]) << std::endl;
      }
    }
#endif
    new_fitness_function_results_.clear();
    for (std::size_t i = 0; i < old_fitness_function_results_.size(); ++i) {
      if (!old_fitness_function_results_[i].fitness.has_value()) {
        old_fitness_function_results_[i].fitness =
            fitness_function_(old_fitness_function_results_[i].input);
        actual_fitness_count_++;
      }
      new_fitness_function_results_.push_back(old_fitness_function_results_[i]);
    }
    old_fitness_function_results_.clear();

#ifdef DEBUG
    for (auto &input : new_fitness_function_results_) {
      if (input.input[0] > 0b111) {
        std::cout << std::bitset<64>(input.input[0]) << std::endl;
      }
    }
#endif
    {
      auto elit_selection_vector = new_fitness_function_results_;
      std::sort(elit_selection_vector.begin(), elit_selection_vector.end());
      convergence_.push_back(
          elit_selection_vector[elit_selection_vector.size() - 1]
              .fitness.value_or(0));
      // elitismus selection
      uint64_t number_of_selections =
          elitism_selection_percent_ * generation_size_;
      for (uint64_t i = 0; i < number_of_selections; i++) {
        auto index = elit_selection_vector.size() - 1;
        old_fitness_function_results_.push_back(elit_selection_vector[index]);
        elit_selection_vector.pop_back();
      }

      // elitismus roulette
      uint64_t number_of_roulete_elites =
          elitism_roulette_percent_ * generation_size_;
      for (int i = 0; i < number_of_roulete_elites; ++i) {
        std::uniform_int_distribution<uint64_t> elit_distribution(
            0, elit_selection_vector.size()-1);
        auto to_delete = elit_distribution(generator_);
        auto new_element = elit_selection_vector[to_delete];
        old_fitness_function_results_.push_back(new_element);
        elit_selection_vector.erase(elit_selection_vector.begin() + to_delete);
      }
    }
#ifdef DEBUG
    for (auto &input : old_fitness_function_results_) {
      if (input.input[0] > 0b111) {
        std::cout << std::bitset<64>(input.input[0]) << std::endl;
      }
    }
#endif

#ifdef DEBUG
    std::cout << std::endl << "Fitness function results:" << std::endl;
    for (auto tzname : new_fitness_function_results_) {
      std::cout << tzname.fitness.value_or(0) << std::endl;
    }
#endif

    Crossing();
  }
}

template <std::size_t SIZE>
void GeneticAlgorithm<SIZE>::Crossing() {

  std::uniform_int_distribution<uint64_t> roulette_distribution(0, new_fitness_function_results_.size()-1);
  std::uniform_real_distribution<double> use_roulette_distribution(0, 1);

  size_t number_of_first_parent = first_parent_range_ * new_fitness_function_results_.size();
  size_t number_of_second_parent = first_parent_range_ * new_fitness_function_results_.size();
  if (number_of_first_parent < MINIMUM_PARENT_COUNT) {
    number_of_first_parent = roulette_distribution(generator_);
    if (number_of_first_parent < MINIMUM_PARENT_COUNT) {
      number_of_first_parent = MINIMUM_PARENT_COUNT;
    }
  }

  if (number_of_second_parent < MINIMUM_PARENT_COUNT) {
    number_of_second_parent = roulette_distribution(generator_);
    if (number_of_second_parent < MINIMUM_PARENT_COUNT) {
      number_of_second_parent = MINIMUM_PARENT_COUNT;
    }
  }

  auto first_parent = Chromosome();
  auto second_parent = Chromosome();

  while (old_fitness_function_results_.size() < generation_size_) {
#ifdef DEBUG
    for (auto &input : old_fitness_function_results_) {
      if (input.input[0] > 0b111) {
        std::cout << std::bitset<64>(input.input[0]) << std::endl;
      }
    }
#endif

    std::uniform_int_distribution<uint64_t> parent_distribution(0, number_of_first_parent);
    auto first_index = parent_distribution(generator_);
    auto second_index = parent_distribution(generator_);
    if (use_roulette_distribution(generator_) < crossing_roulette_probability_) {
      first_parent = new_fitness_function_results_[roulette_distribution(generator_)].input;
      second_parent = new_fitness_function_results_[roulette_distribution(generator_)].input;
#ifdef DEBUG
        if (first_parent[0] > 0b111 || second_parent[0] > 0b111) {
          std::cout << std::bitset<64>(first_parent[0]) << " "<< std::bitset<64>(second_parent[0]) << std::endl;
        }
#endif
    } else {
      first_parent = new_fitness_function_results_[first_index].input;
      second_parent = new_fitness_function_results_[second_index].input;
#ifdef DEBUG
      if (first_parent[0] > 0b111 || second_parent[0] > 0b111) {
        std::cout << std::bitset<64>(first_parent[0]) << " "<< std::bitset<64>(second_parent[0]) << std::endl;
      }
#endif
    }
    auto children = MakeChildren(first_parent, second_parent);
    old_fitness_function_results_.push_back(children[0]);
    if (old_fitness_function_results_.size() == generation_size_) {
      break;
    }
    old_fitness_function_results_.push_back(children[1]);
  }
#ifdef DEBUG
  for (auto &input : old_fitness_function_results_) {
    if (input.input[0] > 0b111) {
      std::cout << std::bitset<64>(input.input[0]) << std::endl;
    }
  }
#endif
}

template <std::size_t SIZE>
std::vector<FitnessFunctionResults_t<SIZE>> GeneticAlgorithm<SIZE>::MakeChildren(
    Chromosome A, Chromosome B) {
  FitnessFunctionResults_t<SIZE> first_child;
  FitnessFunctionResults_t<SIZE> second_child;
  std::uniform_int_distribution<uint64_t> mask_distribution(0, UINT64_MAX);
  std::uniform_real_distribution<double> mutate_distribution(0.0, 1.0);

  Chromosome first_parent_first_half;
  Chromosome first_parent_second_half;
  Chromosome second_parent_first_half;
  Chromosome second_parent_second_half;
  auto mask = Chromosome();
  auto mask_inverted = Chromosome();
  for (std::size_t i = 0; i < SIZE; ++i) {
    mask[i] = mask_distribution(generator_);
    mask_inverted[i] = ~mask[i];

    first_parent_first_half[i] = mask[i] & A[i];
    first_parent_second_half[i] &= mask_inverted[i] & A[i];
    second_parent_first_half[i] &= mask_inverted[i] & B[i];
    second_parent_second_half[i] &= mask[i] & B[i];

    first_child.input[i] = first_parent_first_half[i] | second_parent_first_half[i];
    second_child.input[i] = first_parent_second_half[i] | second_parent_second_half[i];
  }

  uint64_t mutate_mask = 0;
  if (mutate_distribution(generator_) < mutation_probability_) {
    std::uniform_int_distribution<uint64_t> full_position_distribution(0, 63);
    if (A.size() == 1) {
      first_child.input[0] ^= 1 << full_position_distribution(generator_);
      second_child.input[0] ^= 1 << full_position_distribution(generator_);
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      first_child.input[0] &= mutate_mask;
      second_child.input[0] &= mutate_mask;
    } else {
      std::uniform_int_distribution<uint64_t> part_position_distribution(0, SIZE -1);
      first_child.input[part_position_distribution(generator_)] ^= 1 << full_position_distribution(generator_);
      second_child.input[part_position_distribution(generator_)] ^= 1 << full_position_distribution(generator_);
      uint64_t iterator = dimension_size_ - (SIZE - 1) * 64;
      for (std::size_t j = 0; j < iterator; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      first_child.input[SIZE - 1] &= mutate_mask;
      second_child.input[SIZE - 1] &= mutate_mask;
    }

  }
#ifdef DEBUG
  if (first_child.input[0] > 0b111 || second_child.input[0] > 0b111) {
    std::cout << first_child.input[0] << " " << second_child.input[0] << std::endl;
  }
#endif
  return std::vector<FitnessFunctionResults_t<SIZE>>{
    first_child, second_child
  };
}

template <std::size_t SIZE>
std::vector<uint64_t> GeneticAlgorithm<SIZE>::GetConvergence() const {
  return convergence_;
}