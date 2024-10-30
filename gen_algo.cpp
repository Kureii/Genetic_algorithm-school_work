/*
 * Created by kureii on 10/16/24.
 */

#include "gen_algo.h"
#include <bitset>
#include <iostream>
#include <random>
#include <algorithm>
#include <ranges>
#include <stdexcept>

#include "my_array_encapsulation.h"

GeneticAlgorithm::GeneticAlgorithm(const input_structure_t& input_structure, uint64_t chromosome_size,
      const FitnessFunction& fitness_function)
    : chromosome_size_(chromosome_size),
mutation_probability_(input_structure.mutation_probability),
      elitism_roulette_percent_(input_structure.elitism_random_select),
      elitism_selection_percent_(input_structure.elitism_best_select),
      crossing_roulette_probability_(input_structure.crossing_full_random_probability),
      first_parent_range_(input_structure.first_parent_select_range),
      second_parent_range_(input_structure.second_parent_select_range),
      init_population_size_(input_structure.initial_population_size),
      generation_size_(input_structure.population_size),
      dimension_size_(input_structure.problem_dimension),
      stop_limit_(input_structure.fitness_rating_count),
      fitness_function_(fitness_function) {
  if (chromosome_size_ == 0) {
    throw std::invalid_argument("Chromosome size cannot be zero");
  }

  if (dimension_size_ > chromosome_size_ * 64) {
    throw std::invalid_argument(
        "Dimension size is too large for the given chromosome size");
  }

  if (elitism_roulette_percent_ + elitism_selection_percent_ > 0.75) {
    throw std::invalid_argument("Elitism percent cannot be more than 0.75");
  }

  if (init_population_size_ <= 2 || generation_size_ <= 2) {
    throw std::invalid_argument("Population must be bigger than 2");
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


void GeneticAlgorithm::GenerateInitialPopulation() {
  old_fitness_function_results_.clear();
  std::vector<Chromosome> old_chromosomes;
  old_chromosomes.resize(init_population_size_);
  for (std::size_t i = 0; i < init_population_size_; i++) {
    old_chromosomes[i] = myArrayEncapsulation<uint64_t>(chromosome_size_);
  }
  std::uniform_int_distribution<uint64_t> distribution(0, UINT64_MAX);
#ifdef DEBUG
  std::cout << "Initial generation" << std::endl;
#endif
  for (std::size_t i = 0; i < init_population_size_; ++i) {
    for (std::size_t j = 0; j < chromosome_size_; ++j) {
      old_chromosomes[i][j] = distribution(generator_);
    }
    if (chromosome_size_ > 1) {
      uint64_t iterator = dimension_size_ - (chromosome_size_ - 1) * 64;
      uint64_t mask = 0;
      for (std::size_t j = 0; j < iterator; ++j) {
        mask <<= 1;
        mask++;
      }
      old_chromosomes[i][chromosome_size_ - 1] &= mask;
    } else {
      uint64_t mask = 0;
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mask <<= 1;
        mask++;
      }
      old_chromosomes[i][0] &= mask;
    }
#ifdef DEBUG
    for (std::size_t j = 0; j < chromosome_size_; ++j) {
      std::cout << std::bitset<64>(old_chromosomes[i][chromosome_size_ - 1 - j]) << " ";
    }
    std::cout << std::endl;
#endif
    if (!old_fitness_function_results_.empty()) {
      old_fitness_function_results_.clear();
    }
  }
  for (uint64_t i = 0; i < old_chromosomes.size(); i++) {
    old_fitness_function_results_.emplace_back(
        old_chromosomes[i]);
  }
}


void GeneticAlgorithm::Compute() {
  convergence_.clear();
  GenerateInitialPopulation();

  while (actual_fitness_count_ < stop_limit_) {
    auto emergency_brake = !std::ranges::any_of(old_fitness_function_results_, [](const auto& result) {
    return !result.fitness.has_value();
});
    if (emergency_brake) {
#ifdef DEBUG
        std::cout << "emergency break" << std::endl;
#endif
      break;
    }
    new_fitness_function_results_.clear();
    std::ranges::for_each(old_fitness_function_results_, [&](auto& old_fitness_function_result) {
    if (!old_fitness_function_result.fitness.has_value()) {
        old_fitness_function_result.fitness = fitness_function_(old_fitness_function_result.input);
        actual_fitness_count_++;
    }
    new_fitness_function_results_.push_back(old_fitness_function_result);
});

    old_fitness_function_results_.clear();

    Elitism();

#ifdef DEBUG
    std::cout << std::endl << "Fitness function results:" << std::endl;
    for (auto const& f_f_result : new_fitness_function_results_) {
      std::cout << f_f_result.fitness.value_or(0) << ", ";
    }
    std::cout << std::endl;
#endif

    Crossing();
  }
  std::ranges::sort(old_fitness_function_results_, std::greater<>());
  convergence_.emplace_back(old_fitness_function_results_[0].fitness.value_or(0));
}


void GeneticAlgorithm::Crossing() {

  roulette_distribution_ = std::uniform_int_distribution<uint64_t>(0, new_fitness_function_results_.size()-1);

  auto first_parent = Chromosome(chromosome_size_);
  auto second_parent = Chromosome(chromosome_size_);

  std::ranges::sort(new_fitness_function_results_, std::greater<>());
#ifdef DEBUG
  std::cout << std::endl << "Select parent array:" << std::endl;
  for (std::size_t i = 0; i < new_fitness_function_results_.size(); ++i) {
    std::cout << new_fitness_function_results_[i].fitness.value_or(-1) << ", ";
  }
  std::cout << std::endl;

#endif

  while (old_fitness_function_results_.size() < generation_size_) {
    first_parent.reset();
    second_parent.reset();
    ParentSelection(first_parent,second_parent);

    auto children = MakeChildren(first_parent, second_parent);
    old_fitness_function_results_.push_back(children[0]);
    if (old_fitness_function_results_.size() == generation_size_) {
      break;
    }
    old_fitness_function_results_.push_back(children[1]);
  }
}


void GeneticAlgorithm::Elitism() {
  auto elit_selection_vector = new_fitness_function_results_;
  std::ranges::sort(elit_selection_vector);
  convergence_.push_back(
      elit_selection_vector[elit_selection_vector.size() - 1]
          .fitness.value_or(0));
  best_from_generation_ =
      elit_selection_vector[elit_selection_vector.size() - 1].input;


  auto multiplier = std::min(generation_size_, init_population_size_);
  // elitismus selection
  uint64_t number_of_selections =
      elitism_selection_percent_ * multiplier;
  for (uint64_t i = 0; i < number_of_selections; i++) {
    auto index = elit_selection_vector.size() - 1;
    old_fitness_function_results_.push_back(elit_selection_vector[index]);
    elit_selection_vector.pop_back();
  }

  // elitismus roulette
  uint64_t number_of_roulete_elites =
      elitism_roulette_percent_ * multiplier;
  for (int i = 0; i < number_of_roulete_elites; ++i) {
    std::uniform_int_distribution<uint64_t> elit_distribution(
        0, elit_selection_vector.size()-1);
    auto to_delete = elit_distribution(generator_);
    auto new_element = elit_selection_vector[to_delete];
    old_fitness_function_results_.push_back(new_element);
    elit_selection_vector.erase(elit_selection_vector.begin() + to_delete);
  }
}

void GeneticAlgorithm::ParentSelection(
    Chromosome& first_parent, Chromosome& second_parent) {
  if (zero_one_distribution_(generator_) < crossing_roulette_probability_) {
    first_parent = new_fitness_function_results_[roulette_distribution_(generator_)].input;
    second_parent = new_fitness_function_results_[roulette_distribution_(generator_)].input;
#ifdef DEBUG
    if(first_parent.data() == nullptr || second_parent.data() == nullptr) {
      std::cerr << "Problem with random select: " << std::endl;
    }
#endif
  } else {
    size_t number_of_first_parent = first_parent_range_ * new_fitness_function_results_.size();
    size_t number_of_second_parent = first_parent_range_ * new_fitness_function_results_.size();
    while (number_of_first_parent < MINIMUM_PARENT_COUNT) {
      number_of_first_parent = roulette_distribution_(generator_);
    }

    while (number_of_second_parent < MINIMUM_PARENT_COUNT) {
      number_of_second_parent = roulette_distribution_(generator_);
    }

    std::vector<Chromosome> select_first_parent_vector;
    for (std::size_t i = 0; i < new_fitness_function_results_.size() - number_of_first_parent; ++i) {
      uint64_t iterator = new_fitness_function_results_[i].fitness.value_or(0) + 1;
      for (std::size_t j = 0; j < (iterator); ++j) {
        select_first_parent_vector.emplace_back(new_fitness_function_results_[i].input);
      }
    }

    std::vector<Chromosome> select_second_parent_vector;
    for (std::size_t i = 0; i < new_fitness_function_results_.size() - number_of_second_parent; ++i) {
      uint64_t iterator = new_fitness_function_results_[i].fitness.value_or(0) + 1;
      for (std::size_t j = 0; j < (iterator); ++j) {
        select_second_parent_vector.emplace_back(new_fitness_function_results_[i].input);
      }
    }
    std::uniform_int_distribution<uint64_t> first_parent_distribution(0, select_first_parent_vector.size()-1);
    std::uniform_int_distribution<uint64_t> second_parent_distribution(0, select_second_parent_vector.size()-1);

    auto first_index = first_parent_distribution(generator_);
    auto second_index = second_parent_distribution(generator_);
    first_parent = select_first_parent_vector[first_index];
    second_parent = select_second_parent_vector[second_index];

#ifdef DEBUG
    if(first_parent.data() == nullptr || second_parent.data() == nullptr) {
      std::cerr << "Problem with roulette select: " << std::endl;
    }
#endif
  }

#ifdef DEBUG
  if(first_parent.data() == nullptr || second_parent.data() == nullptr) {
    std::cerr << "Some parent is null" << std::endl;
  }
#endif
}

void GeneticAlgorithm::Mutation(
    FitnessFunctionResults& first_child, FitnessFunctionResults& second_child) {
  std::uniform_int_distribution<uint64_t> full_position_distribution(0, 63);
  std::uniform_int_distribution<uint64_t> part_position_distribution(0, chromosome_size_ -1);
  uint64_t mutate_mask = 0;
  // first child
  if (zero_one_distribution_(generator_) < mutation_probability_) {
    if (first_child.input.size() == 1) {
      first_child.input[0] ^= 1 << full_position_distribution(generator_);
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      first_child.input[0] &= mutate_mask;
    } else {
      first_child.input[part_position_distribution(generator_)] ^= 1 << full_position_distribution(generator_);
      uint64_t iterator = dimension_size_ - (chromosome_size_ - 1) * 64;
      for (std::size_t j = 0; j < iterator; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      first_child.input[chromosome_size_ - 1] ^= mutate_mask;
    }
  }

  // second child
  if (zero_one_distribution_(generator_) < mutation_probability_) {
    if (second_child.input.size() == 1) {
      second_child.input[0] ^= 1 << full_position_distribution(generator_);
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      second_child.input[0] &= mutate_mask;
    } else {
      second_child.input[part_position_distribution(generator_)] ^= 1 << full_position_distribution(generator_);
      uint64_t iterator = dimension_size_ - (chromosome_size_ - 1) * 64;
      for (std::size_t j = 0; j < iterator; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      second_child.input[chromosome_size_ - 1] ^= mutate_mask;
    }
  }
}


std::vector<FitnessFunctionResults> GeneticAlgorithm::MakeChildren(
    Chromosome& A, Chromosome& B) {
  FitnessFunctionResults first_child {myArrayEncapsulation<uint64_t>(chromosome_size_)};
  FitnessFunctionResults second_child {myArrayEncapsulation<uint64_t>(chromosome_size_)};
  std::uniform_int_distribution<uint64_t> mask_distribution(0, UINT64_MAX);
#ifdef DEBUG
  uint8_t bit_len_normal, bit_len = 64;
  uint8_t bit_len_sort = dimension_size_ - (chromosome_size_ - 1) * 64;
  std::cout << "Crossing:" << std::endl << "A: " << fitness_function_(A) << "; ";
  for (int64_t j = static_cast<int64_t>(chromosome_size_)-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(A[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << " B: " << fitness_function_(B) << "; ";
  for (int64_t j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(B[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
#endif

  Chromosome first_parent_first_half(chromosome_size_);
  Chromosome first_parent_second_half(chromosome_size_);
  Chromosome second_parent_first_half(chromosome_size_);
  Chromosome second_parent_second_half(chromosome_size_);
  Chromosome mask(chromosome_size_);
  Chromosome mask_inverted(chromosome_size_);
  for (std::size_t i = 0; i < chromosome_size_; ++i) {
    mask[i] = mask_distribution(generator_);
    mask_inverted[i] = ~mask[i];

    first_parent_first_half[i] = mask[i] & A[i];
    first_parent_second_half[i] = mask_inverted[i] & A[i];
    second_parent_first_half[i] = mask_inverted[i] & B[i];
    second_parent_second_half[i] = mask[i] & B[i];

    first_child.input[i] = first_parent_first_half[i] | second_parent_first_half[i];
    second_child.input[i] = first_parent_second_half[i] | second_parent_second_half[i];
  }
#ifdef DEBUG

  std::cout << std::endl << "mask: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(mask[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "mask_inverted: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
        std::bitset<64> bits(mask_inverted[j]);
std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_parent_first_half: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(first_parent_first_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_parent_second_half: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(first_parent_second_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_parent_first_half: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(second_parent_first_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_parent_second_half: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(second_parent_second_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_child: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(first_child.input[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_child: ";
  for (int j = chromosome_size_-1; 0 <= j; --j) {
    bit_len = j==0? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(second_child.input[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "--------------------------------" << std::endl;
#endif
  Mutation(first_child, second_child);

  return std::vector{
    first_child, second_child
  };
}


std::vector<uint64_t> GeneticAlgorithm::GetConvergence() const {
  return convergence_;
}

std::vector<uint8_t> GeneticAlgorithm::GetResult() const {
  std::vector<uint8_t> final_vector;
  std::ranges::for_each(best_from_generation_, [&](auto bit) {
    for (auto i = 0; i < 64 ; ++i) {
      final_vector.push_back((bit >> i) & 0b1 );
    }
  } );
  final_vector.resize(dimension_size_);
  std::ranges::reverse(final_vector);

  return final_vector;
}
