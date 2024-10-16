/*
 * Created by kureii on 10/16/24.
*/

#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>

template <std::size_t SIZE>
struct FitnessFunctionResults {
  std::optional<uint64_t> fitness;
  std::array<uint64_t, SIZE> input;

  FitnessFunctionResults() = default;

  explicit FitnessFunctionResults(const std::array<uint64_t, SIZE>& inputArray, uint64_t fitnessValue)
      : fitness(fitnessValue), input(inputArray) {}

  explicit FitnessFunctionResults( uint64_t fitnessValue)
      : fitness(fitnessValue) {}

  explicit FitnessFunctionResults(const std::array<uint64_t, SIZE>& inputArray)
      : input(inputArray) {
    fitness = std::nullopt;
  }

  bool operator<(const FitnessFunctionResults& other) const {
    if (fitness.has_value() && other.fitness.has_value()) {
      return fitness.value() < other.fitness.value();
    }
    if (fitness.has_value() && !other.fitness.has_value()) {
      return false;
    }
    if (!fitness.has_value() && other.fitness.has_value()) {
      return true;
    }
    return false;
  }
};

template <std::size_t SIZE>
using FitnessFunctionResults_t = FitnessFunctionResults<SIZE>;
