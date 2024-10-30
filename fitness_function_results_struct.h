/*
 * Created by kureii on 10/16/24.
 */

#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>

#include "my_array_encapsulation.h"

struct FitnessFunctionResults {
  std::optional<uint64_t> fitness;
  myArrayEncapsulation<uint64_t> input{};

  FitnessFunctionResults() = default;

  explicit FitnessFunctionResults(
      const myArrayEncapsulation<uint64_t>& inputArray, uint64_t fitnessValue)
      : fitness(fitnessValue), input(inputArray) {}

  explicit FitnessFunctionResults(uint64_t fitnessValue)
      : fitness(fitnessValue) {}

  explicit FitnessFunctionResults(
      const myArrayEncapsulation<uint64_t>& inputArray)
      : input(inputArray) {
    fitness = std::nullopt;
  }

  std::strong_ordering operator<=>(const FitnessFunctionResults& other) const {
    if (fitness.has_value() && other.fitness.has_value()) {
      return fitness.value() <=> other.fitness.value();
    }
    if (fitness.has_value() && !other.fitness.has_value()) {
      return std::strong_ordering::greater;
    }
    if (!fitness.has_value() && other.fitness.has_value()) {
      return std::strong_ordering::less;
    }
    return std::strong_ordering::equal;
  }

  bool operator<=(const FitnessFunctionResults& other) const {
    return (*this <=> other) != std::strong_ordering::greater;
  }

  bool operator>(const FitnessFunctionResults& other) const {
    return (*this <=> other) == std::strong_ordering::greater;
  }

  bool operator>=(const FitnessFunctionResults& other) const {
    return (*this <=> other) != std::strong_ordering::less;
  }

  bool operator==(const FitnessFunctionResults& other) const {
    return (*this <=> other) == std::strong_ordering::equal;
  }
};
