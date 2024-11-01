/*
 * Created by kureii on 10/16/24.
 */

#include "gen_algo.h"

#include <algorithm>
#include <bitset>
#include <complex>
#include <iostream>
#include <random>
#include <ranges>
#include <stdexcept>
#include <utility>

#include "chromosome_array.h"

GeneticAlgorithm::GeneticAlgorithm(const input_structure_t& input_structure,
    uint64_t chromosome_size, FitnessFunction fitness_function,
    bool real_numbers)
    : real_numbers_(real_numbers),
      mutation_probability_(input_structure.mutation_probability),
      elitism_roulette_percent_(input_structure.elitism_random_select),
      elitism_selection_percent_(input_structure.elitism_best_select),
      crossing_roulette_probability_(
          input_structure.crossing_full_random_probability),
      first_parent_range_(input_structure.first_parent_select_range),
      second_parent_range_(input_structure.second_parent_select_range),
      chromosome_size_(chromosome_size),
      init_population_size_(input_structure.initial_population_size),
      generation_size_(input_structure.population_size),
      dimension_size_(input_structure.problem_dimension),
      stop_limit_(input_structure.fitness_rating_count),
      mapping_structure_(input_structure.mapping_structure),
      fitness_function_(std::move(fitness_function)) {
  CheckInit();

  old_fitness_function_results_.resize(init_population_size_);
  new_fitness_function_results_.resize(init_population_size_);

  std::random_device random_device;
  generator_ = std::mt19937_64(random_device());
  if (real_numbers_) {
    ieee_distribution_ =
        std::uniform_real_distribution(mapping_structure_.value().min_value,
            mapping_structure_.value().max_value);
  }
}

void GeneticAlgorithm::CheckInit() const {
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
    throw std::invalid_argument(
        "All values of type double must be in the range 0-1");
  }

  if (real_numbers_ && !mapping_structure_.has_value()) {
    throw std::invalid_argument("mapping structure does not exist");
  }
}

void GeneticAlgorithm::GenerateInitialPopulation() {
  old_fitness_function_results_.clear();
  if (real_numbers_ && mapping_structure_.value().mapping_method ==
                           mapping_method_t::BIT_AS_DOUBLE) {
    GenerateInitialPopulationRealIEEE();
  } else if (real_numbers_) {
    GenerateInitialPopulationBit(dimension_size_, dimension_size_ * 64);
  } else {
    GenerateInitialPopulationBit(chromosome_size_, dimension_size_);
  }
}
void GeneticAlgorithm::GenerateInitialPopulationBit(
    uint64_t chromosome_array_size, uint64_t all_bit_size) {
  std::vector<Chromosome> old_chromosomes(init_population_size_);
  for (std::size_t i = 0; i < init_population_size_; i++) {
    old_chromosomes[i] =
        ChromosomeArray<uint64_t>(chromosome_array_size, all_bit_size);
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
      std::cout << std::bitset<64>(old_chromosomes[i][chromosome_size_ - 1 - j])
                << " ";
    }
    std::cout << std::endl;
#endif
  }
  if (!old_fitness_function_results_.empty()) {
    old_fitness_function_results_.clear();
  }
  for (auto& old_chromosome : old_chromosomes) {
    old_fitness_function_results_.emplace_back(old_chromosome);
  }
}

void GeneticAlgorithm::GenerateInitialPopulationRealIEEE() {
  for (uint64_t i = 0; i < init_population_size_; ++i) {
    ChromosomeArray<uint64_t> chromosome(dimension_size_, dimension_size_ * 64);
    for (uint64_t j = 0; j < dimension_size_; ++j) {
      chromosome[j] = std::bit_cast<uint64_t>(ieee_distribution_(generator_));
    }
    old_fitness_function_results_.emplace_back(chromosome);
  }
}

void GeneticAlgorithm::Compute() {
  convergence_.clear();
  GenerateInitialPopulation();

  while (actual_fitness_count_ < stop_limit_) {
    auto emergency_brake = !std::ranges::any_of(
        old_fitness_function_results_, [](const auto& result) {
          return std::holds_alternative<std::monostate>(result.fitness);
        });

    if (emergency_brake) {
#ifdef DEBUG
      std::cout << "emergency break" << std::endl;
#endif
      break;
    }
    new_fitness_function_results_.clear();
    ComputeFitnessFunction();

    old_fitness_function_results_.clear();

    Elitism();

#ifdef DEBUG
    std::cout << std::endl << "Fitness function results:" << std::endl;
    for (auto const& f_f_result : new_fitness_function_results_) {
      std::cout << f_f_result.GetFitnessAsUInt64() << ", ";
    }
    std::cout << std::endl;
#endif
  }
  std::ranges::sort(old_fitness_function_results_, std::greater());
  if (real_numbers_) {
    if (mapping_structure_.value().mapping_method ==
        mapping_method_t::BIT_AS_DOUBLE) {
      CrossingReal();
    } else {
      CrossingBit(dimension_size_, dimension_size_ * 64);
    }
    convergence_.emplace_back(
        old_fitness_function_results_[0].GetFitnessAsDouble());
  } else {
    CrossingBit(chromosome_size_, dimension_size_);
    convergence_.emplace_back(
        old_fitness_function_results_[0].GetFitnessAsUInt64());
  }
}

void GeneticAlgorithm::ComputeFitnessFunction() {
  auto compute_f_f = [this]<typename T>(T&& input) -> uint64_t {
    using InputType = std::decay_t<std::remove_reference_t<T>>;
    static_assert(
        std::is_same_v<InputType, ChromosomeArray<uint64_t>> ||
        std::is_same_v<InputType, ChromosomeArray<double>>,
        "Unexpected input type"
    );
    if constexpr (std::is_same_v<InputType, ChromosomeArray<uint64_t>>) {
      return this->fitness_function_(input);
    } else {
      ChromosomeArray<uint64_t> converted_input;
      for (size_t i = 0; i < input.size(); ++i) {
        converted_input[i] = static_cast<uint64_t>(input[i]);
      }
      return this->fitness_function_(converted_input);
    }
  };
  std::ranges::for_each(
      old_fitness_function_results_, [&](auto& old_fitness_function_result) {
        if (!old_fitness_function_result.HasValidFitness()) {
          old_fitness_function_result.fitness = std::visit(
              compute_f_f, old_fitness_function_result.input);
          actual_fitness_count_++;
        }
        new_fitness_function_results_.push_back(old_fitness_function_result);
      });
}


void GeneticAlgorithm::CrossingBit(
    uint64_t chromosome_array_size, uint64_t all_bit_size) {
  roulette_distribution_ = std::uniform_int_distribution<uint64_t>(
      0, new_fitness_function_results_.size() - 1);

  auto first_parent = Chromosome(chromosome_array_size, all_bit_size);
  auto second_parent = Chromosome(chromosome_array_size, all_bit_size);

  std::ranges::sort(new_fitness_function_results_, std::greater<>());
#ifdef DEBUG
  std::cout << std::endl << "Select parent array:" << std::endl;
  for (std::size_t i = 0; i < new_fitness_function_results_.size(); ++i) {
    std::cout << new_fitness_function_results_[i].GetFitnessAsUInt64() << ", ";
  }
  std::cout << std::endl;

#endif

  while (old_fitness_function_results_.size() < generation_size_) {
    first_parent.reset();
    second_parent.reset();
    ParentSelection(first_parent, second_parent);

    auto children = MakeChildren(first_parent, second_parent);
    old_fitness_function_results_.push_back(children[0]);
    if (old_fitness_function_results_.size() == generation_size_) {
      break;
    }
    old_fitness_function_results_.push_back(children[1]);
  }
}

void GeneticAlgorithm::CrossingReal() {
  roulette_distribution_ = std::uniform_int_distribution<uint64_t>(
          0, new_fitness_function_results_.size() - 1);

  std::ranges::sort(new_fitness_function_results_, std::greater<>());

  while (old_fitness_function_results_.size() < generation_size_) {
    auto children = MakeChildrenReal();
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
  convergence_.push_back(elit_selection_vector[elit_selection_vector.size() - 1]
                             .GetFitnessAsUInt64());
  if (real_numbers_ && mapping_structure_.value().mapping_method == mapping_method_t::BIT_AS_DOUBLE) {
    best_from_generation_double_ =
        std::get<ChromosomeArray<double>>(elit_selection_vector[elit_selection_vector.size() - 1].input);
  } else {
    best_from_generation_ =
        std::get<ChromosomeArray<uint64_t>>(elit_selection_vector[elit_selection_vector.size() - 1].input);
  }
  auto multiplier = std::min(generation_size_, init_population_size_);
  // elitismus selection
  uint64_t number_of_selections = elitism_selection_percent_ * multiplier;
  for (uint64_t i = 0; i < number_of_selections; i++) {
    auto index = elit_selection_vector.size() - 1;
    old_fitness_function_results_.push_back(elit_selection_vector[index]);
    elit_selection_vector.pop_back();
  }

  // elitismus roulette
  uint64_t number_of_roulete_elites = elitism_roulette_percent_ * multiplier;
  for (int i = 0; i < number_of_roulete_elites; ++i) {
    std::uniform_int_distribution<uint64_t> elit_distribution(
        0, elit_selection_vector.size() - 1);
    auto to_delete = elit_distribution(generator_);
    auto new_element = elit_selection_vector[to_delete];
    old_fitness_function_results_.push_back(new_element);
    elit_selection_vector.erase(elit_selection_vector.begin() + to_delete);
  }
}

void GeneticAlgorithm::ParentSelection(
    Chromosome& first_parent, Chromosome& second_parent) {
  if (zero_one_distribution_(generator_) < crossing_roulette_probability_) {
    first_parent =
        std::get<ChromosomeArray<uint64_t>>(new_fitness_function_results_[roulette_distribution_(generator_)].input);
    second_parent =
        std::get<ChromosomeArray<uint64_t>>(new_fitness_function_results_[roulette_distribution_(generator_)].input);
#ifdef DEBUG
    if (first_parent.data() == nullptr || second_parent.data() == nullptr) {
      std::cerr << "Problem with random select: " << std::endl;
    }
#endif
  } else {
    size_t number_of_first_parent =
        first_parent_range_ * new_fitness_function_results_.size();
    size_t number_of_second_parent =
        first_parent_range_ * new_fitness_function_results_.size();
    while (number_of_first_parent < MINIMUM_PARENT_COUNT) {
      number_of_first_parent = roulette_distribution_(generator_);
    }

    while (number_of_second_parent < MINIMUM_PARENT_COUNT) {
      number_of_second_parent = roulette_distribution_(generator_);
    }

    std::vector<Chromosome> select_first_parent_vector;
    for (std::size_t i = 0;
         i < new_fitness_function_results_.size() - number_of_first_parent;
         ++i) {
      uint64_t iterator =
          new_fitness_function_results_[i].GetFitnessAsUInt64() + 1;
      for (std::size_t j = 0; j < (iterator); ++j) {
        select_first_parent_vector.emplace_back(
            std::get<ChromosomeArray<uint64_t>>(new_fitness_function_results_[i].input));
      }
    }

    std::vector<Chromosome> select_second_parent_vector;
    for (std::size_t i = 0;
         i < new_fitness_function_results_.size() - number_of_second_parent;
         ++i) {
      uint64_t iterator =
          new_fitness_function_results_[i].GetFitnessAsUInt64() + 1;
      for (std::size_t j = 0; j < (iterator); ++j) {
        select_second_parent_vector.emplace_back(
            std::get<ChromosomeArray<uint64_t>>(new_fitness_function_results_[i].input));
      }
    }
    std::uniform_int_distribution<uint64_t> first_parent_distribution(
        0, select_first_parent_vector.size() - 1);
    std::uniform_int_distribution<uint64_t> second_parent_distribution(
        0, select_second_parent_vector.size() - 1);

    auto first_index = first_parent_distribution(generator_);
    auto second_index = second_parent_distribution(generator_);
    first_parent = select_first_parent_vector[first_index];
    second_parent = select_second_parent_vector[second_index];

#ifdef DEBUG
    if (first_parent.data() == nullptr || second_parent.data() == nullptr) {
      std::cerr << "Problem with roulette select: " << std::endl;
    }
#endif
  }

#ifdef DEBUG
  if (first_parent.data() == nullptr || second_parent.data() == nullptr) {
    std::cerr << "Some parent is null" << std::endl;
  }
#endif
}

std::pair<ChromosomeArray<double>, ChromosomeArray<double>>
GeneticAlgorithm::ParentSelectionReal() {
  auto select_parent = [this]() -> ChromosomeArray<double> {
    if (zero_one_distribution_(generator_) < crossing_roulette_probability_) {
      return std::get<ChromosomeArray<double>>(
          new_fitness_function_results_[roulette_distribution_(generator_)].input);
    } else {
      // Turnajová selekce
      std::uniform_int_distribution<size_t> tournament_dist(0,
          new_fitness_function_results_.size() - 1);

      const size_t tournament_size = 3;
      size_t best_idx = tournament_dist(generator_);
      uint64_t best_fitness = new_fitness_function_results_[best_idx].GetFitnessAsUInt64();

      for (size_t i = 1; i < tournament_size; ++i) {
        size_t idx = tournament_dist(generator_);
        uint64_t fitness = new_fitness_function_results_[idx].GetFitnessAsUInt64();
        if (fitness > best_fitness) {
          best_idx = idx;
          best_fitness = fitness;
        }
      }

      return std::get<ChromosomeArray<double>>(
          new_fitness_function_results_[best_idx].input);
    }
  };

  return {select_parent(), select_parent()};
}

void GeneticAlgorithm::Mutation(
    FitnessFunctionResults& first_child, FitnessFunctionResults& second_child) {
  std::uniform_int_distribution<uint64_t> full_position_distribution(0, 63);
  std::uniform_int_distribution<uint64_t> part_position_distribution(
      0, chromosome_size_ - 1);
  uint64_t mutate_mask = 0;
  // first child
  if (zero_one_distribution_(generator_) < mutation_probability_) {
    if (std::get<ChromosomeArray<uint64_t>>(first_child.input).size() == 1) {
      std::get<ChromosomeArray<uint64_t>>(first_child.input)[0] ^= 1 << full_position_distribution(generator_);
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      std::get<ChromosomeArray<uint64_t>>(first_child.input)[0] &= mutate_mask;
    } else {
      std::get<ChromosomeArray<uint64_t>>(first_child.input)[part_position_distribution(generator_)] ^=
          1 << full_position_distribution(generator_);
      uint64_t iterator = dimension_size_ - (chromosome_size_ - 1) * 64;
      for (std::size_t j = 0; j < iterator; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      std::get<ChromosomeArray<uint64_t>>(first_child.input)[chromosome_size_ - 1] ^= mutate_mask;
    }
  }

  // second child
  if (zero_one_distribution_(generator_) < mutation_probability_) {
    if (std::get<ChromosomeArray<uint64_t>>(second_child.input).size() == 1) {
      std::get<ChromosomeArray<uint64_t>>(second_child.input)[0] ^= 1 << full_position_distribution(generator_);
      for (std::size_t j = 0; j < dimension_size_; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      std::get<ChromosomeArray<uint64_t>>(second_child.input)[0] &= mutate_mask;
    } else {
      std::get<ChromosomeArray<uint64_t>>(second_child.input)[part_position_distribution(generator_)] ^=
          1 << full_position_distribution(generator_);
      uint64_t iterator = dimension_size_ - (chromosome_size_ - 1) * 64;
      for (std::size_t j = 0; j < iterator; ++j) {
        mutate_mask <<= 1;
        mutate_mask++;
      }
      std::get<ChromosomeArray<uint64_t>>(second_child.input)[chromosome_size_ - 1] ^= mutate_mask;
    }
  }
}

void GeneticAlgorithm::MutationReal(FitnessFunctionResults& first_child,
                                   FitnessFunctionResults& second_child) {
  std::uniform_real_distribution<double> mutation_strength(0.8, 1.2);  // ±20% změna

  auto mutate_chromosome = [&](auto& child) {
    if (zero_one_distribution_(generator_) < mutation_probability_) {
      std::uniform_int_distribution<size_t> pos_dist(0, std::get<ChromosomeArray<double>>(child.input).size() - 1);
      size_t pos = pos_dist(generator_);

      auto& chromosome = std::get<ChromosomeArray<double>>(child.input);
      chromosome[pos] *= mutation_strength(generator_);

      // Omezení na platný rozsah
      if (chromosome[pos] < mapping_structure_.value().min_value) {
        chromosome[pos] = mapping_structure_.value().min_value;
      }
      if (chromosome[pos] > mapping_structure_.value().max_value) {
        chromosome[pos] = mapping_structure_.value().max_value;
      }
    }
  };

  mutate_chromosome(first_child);
  mutate_chromosome(second_child);
}

std::vector<FitnessFunctionResults> GeneticAlgorithm::MakeChildren(
    Chromosome& A, Chromosome& B) {
  FitnessFunctionResults first_child{
      ChromosomeArray<uint64_t>(A.size(), A.chromosome_length())};
  FitnessFunctionResults second_child{
      ChromosomeArray<uint64_t>(A.size(), A.chromosome_length())};
  std::uniform_int_distribution<uint64_t> mask_distribution(0, UINT64_MAX);
#ifdef DEBUG
  uint8_t bit_len_normal, bit_len = 64;
  uint8_t bit_len_sort = dimension_size_ - (chromosome_size_ - 1) * 64;
  std::cout << "Crossing:" << std::endl
            << "A: " << fitness_function_(A) << "; ";
  for (int64_t j = static_cast<int64_t>(chromosome_size_) - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(A[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << " B: " << fitness_function_(B) << "; ";
  for (int64_t j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(B[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
#endif

  Chromosome first_parent_first_half(A.size(), A.chromosome_length());
  Chromosome first_parent_second_half(A.size(), A.chromosome_length());
  Chromosome second_parent_first_half(A.size(), A.chromosome_length());
  Chromosome second_parent_second_half(A.size(), A.chromosome_length());
  Chromosome mask(A.size(), A.chromosome_length());
  Chromosome mask_inverted(A.size(), A.chromosome_length());
  for (std::size_t i = 0; i < chromosome_size_; ++i) {
    mask[i] = mask_distribution(generator_);
    mask_inverted[i] = ~mask[i];

    first_parent_first_half[i] = mask[i] & A[i];
    first_parent_second_half[i] = mask_inverted[i] & A[i];
    second_parent_first_half[i] = mask_inverted[i] & B[i];
    second_parent_second_half[i] = mask[i] & B[i];

    std::get<ChromosomeArray<uint64_t>>(first_child.input)[i] =
        first_parent_first_half[i] | second_parent_first_half[i];
    std::get<ChromosomeArray<uint64_t>>(second_child.input)[i] =
        first_parent_second_half[i] | second_parent_second_half[i];
  }
  if (real_numbers_) {
    switch (mapping_structure_.value().mapping_method) {
      using enum mapping_method_t;
      case FIXED_POINT: {
        auto make_valid = [this](auto& input) {
          std::ranges::for_each(std::get<ChromosomeArray<uint64_t>>(input), [this](auto& x) {
            constexpr int MAX_ATTEMPTS = 100;
            int attempts = 0;

            while (!IsFixedPointGenValid(x) && attempts < MAX_ATTEMPTS) {
              x = (x * 9) / 1;
              attempts++;
            }

            if (attempts >= MAX_ATTEMPTS) {
              x = mapping_structure_.value().min_value;
            }
          });
        };

        make_valid(first_child.input);
        make_valid(second_child.input);

        break;
      }
      case MAPPED_RANGE: {
        auto make_valid = [this](auto& input) {
          MakeValidChromosomeBinaryCodedDecimal(std::get<ChromosomeArray<uint64_t>>(input));
        };

        make_valid(first_child.input);
        make_valid(second_child.input);
        break;
      }
      case BINARY_CODED_DECIMAL: {
        break;
      }
      default:;
    }
  }
#ifdef DEBUG

  std::cout << std::endl << "mask: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(mask[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "mask_inverted: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(mask_inverted[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_parent_first_half: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(first_parent_first_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_parent_second_half: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(first_parent_second_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_parent_first_half: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(second_parent_first_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_parent_second_half: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(second_parent_second_half[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "first_child: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(std::get<ChromosomeArray<uint64_t>>(first_child.input)[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "second_child: ";
  for (int j = chromosome_size_ - 1; 0 <= j; --j) {
    bit_len = j == 0 ? bit_len_sort : bit_len_normal;
    std::bitset<64> bits(std::get<ChromosomeArray<uint64_t>>(second_child.input)[j]);
    std::cout << bits.to_string().substr(64 - bit_len);
  }
  std::cout << std::endl << "--------------------------------" << std::endl;
#endif
  Mutation(first_child, second_child);

  return std::vector{first_child, second_child};
}

std::vector<FitnessFunctionResults> GeneticAlgorithm::MakeChildrenReal() {
  auto [first_parent, second_parent] = ParentSelectionReal();

  // Vytvoříme FFR s ChromosomeArray<double>
  FitnessFunctionResults first_child{ChromosomeArray<double>(first_parent.size(),dimension_size_*64)};
  FitnessFunctionResults second_child{ChromosomeArray<double>(second_parent.size(),dimension_size_*64)};

  // Získáme reference na ChromosomeArray<double>
  auto& first_child_array = std::get<ChromosomeArray<double>>(first_child.input);
  auto& second_child_array = std::get<ChromosomeArray<double>>(second_child.input);

  // Interpolační křížení
  std::uniform_real_distribution<double> alpha_dist(0.0, 1.0);
  double alpha = alpha_dist(generator_);

  for (size_t i = 0; i < first_parent.size(); ++i) {
    first_child_array[i] = alpha * first_parent[i] + (1 - alpha) * second_parent[i];
    second_child_array[i] = (1 - alpha) * first_parent[i] + alpha * second_parent[i];
  }

  MutationReal(first_child, second_child);
  return std::vector{first_child, second_child};
}

bool GeneticAlgorithm::IsFixedPointGenValid(uint64_t gen) {
  bool is_valid = true;
  const auto reduced_gen = gen >> mapping_structure_.value().fractional_bits;
  const bool is_reduced_in_range =
      (reduced_gen > mapping_structure_.value().min_value &&
          reduced_gen < mapping_structure_.value().max_value);
  is_valid &= is_reduced_in_range;
  if (!is_reduced_in_range) {
    const bool is_reduced_equal =
        reduced_gen == mapping_structure_.value().min_value ||
        reduced_gen == mapping_structure_.value().max_value;
    is_valid &= is_reduced_equal;
    if (is_reduced_equal) {
      uint64_t mask = 0;
      for (uint64_t j = 0; j < mapping_structure_.value().fractional_bits;
           ++j) {
        mask <<= 1;
        mask |= 1;
      }
      double max_fractional, min_fractional;
      modf(mapping_structure_.value().max_value, &max_fractional);
      modf(mapping_structure_.value().min_value, &min_fractional);

      // Převod desetinných částí na fixed-point reprezentaci
      auto max_fract_fixed = static_cast<uint64_t>(
          max_fractional *
          (1ULL << mapping_structure_.value().fractional_bits));
      auto min_fract_fixed = static_cast<uint64_t>(
          min_fractional *
          (1ULL << mapping_structure_.value().fractional_bits));

      uint64_t gen_fract = gen & mask;
      is_valid &= (gen_fract > max_fract_fixed || gen_fract < min_fract_fixed);
    }
  }

  return is_valid;
}

void GeneticAlgorithm::MakeValidChromosomeBinaryCodedDecimal(Chromosome& A) {
  constexpr int BITS_PER_NIBBLE = 4;
  constexpr int NIBBLES_IN_UINT64 = 16;
  constexpr uint64_t NIBBLE_MASK = 0b1111;
  constexpr uint64_t MAX_BCD_VALUE = 9;

  std::ranges::for_each(A, [](auto& gen) {
    for (int i = 0; i < NIBBLES_IN_UINT64; ++i) {  // 64/4
      uint64_t mask = NIBBLE_MASK << (i * BITS_PER_NIBBLE);
      if (((gen & mask) >> (i * BITS_PER_NIBBLE)) > MAX_BCD_VALUE) {
        gen &= ~mask;
        gen |= MAX_BCD_VALUE << (i * BITS_PER_NIBBLE);
      }
    }
  });
}

std::vector<uint64_t> GeneticAlgorithm::GetConvergence() const {
  return convergence_;
}

std::vector<uint8_t> GeneticAlgorithm::GetResult() const {
  std::vector<uint8_t> final_vector;
  std::ranges::for_each(best_from_generation_, [&](auto bit) {
    for (auto i = 0; i < 64; ++i) {
      final_vector.push_back((bit >> i) & 0b1);
    }
  });
  final_vector.resize(dimension_size_);
  std::ranges::reverse(final_vector);

  return final_vector;
}
