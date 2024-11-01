/*
 * Created by kureii on 10/16/24.
 */

#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>

#include "chromosome_array.h"

#include <variant>
#include <iostream>
#include <compare>

struct FitnessFunctionResults {
    std::variant<std::monostate, uint64_t, double> fitness;
    std::variant<ChromosomeArray<uint64_t>,ChromosomeArray<double>> input{};

    FitnessFunctionResults() = default;

    FitnessFunctionResults(const ChromosomeArray<uint64_t>& inputArray, uint64_t fitnessValue)
        : fitness(fitnessValue), input(inputArray) {}

    explicit FitnessFunctionResults(uint64_t fitnessValue)
        : fitness(fitnessValue) {}

    explicit FitnessFunctionResults(const ChromosomeArray<uint64_t>& inputArray)
        : input(inputArray) {
        fitness = std::monostate{};
    }

  explicit FitnessFunctionResults(const ChromosomeArray<double>& inputArray)
        : input(inputArray) {
      fitness = std::monostate{};
    }

  FitnessFunctionResults(const ChromosomeArray<double>& inputArray, double fitnessValue)
      : fitness(fitnessValue), input(inputArray) {}

  [[nodiscard]] bool HasValidFitness() const {
      return !std::holds_alternative<std::monostate>(fitness);
    }

  [[nodiscard]] double GetFitnessAsDouble(double invalid_return = 0.0) const {
      return std::visit([invalid_return](const auto& arg) -> double {
          using T = std::decay_t<decltype(arg)>;
          if constexpr (std::is_same_v<T, std::monostate>) {
              return invalid_return;
          } else if constexpr (std::is_same_v<T, uint64_t>) {
              return static_cast<double>(arg);
          } else {
              return arg;
          }
      }, fitness);
    }

  [[nodiscard]] uint64_t GetFitnessAsUInt64(uint64_t invalid_return = 0) const {
      return std::visit([invalid_return](const auto& arg) -> uint64_t {
          using T = std::decay_t<decltype(arg)>;
          if constexpr (std::is_same_v<T, std::monostate>) {
              return invalid_return;
          } else if constexpr (std::is_same_v<T, double>) {
              return static_cast<uint64_t>(arg);
          } else {
              return arg;
          }
      }, fitness);
    }

    // Přetížení spaceship operátoru <=> (std::strong_ordering)
  std::partial_ordering operator<=>(const FitnessFunctionResults& other) const {
      return std::visit(
          [](const auto& lhs, const auto& rhs) -> std::partial_ordering {
              using LhsT = std::decay_t<decltype(lhs)>;
              using RhsT = std::decay_t<decltype(rhs)>;

              if constexpr (std::is_same_v<LhsT, std::monostate> && std::is_same_v<RhsT, std::monostate>) {
                  return std::partial_ordering::equivalent;
              } else if constexpr (std::is_same_v<LhsT, std::monostate>) {
                  return std::partial_ordering::less;
              } else if constexpr (std::is_same_v<RhsT, std::monostate>) {
                  return std::partial_ordering::greater;
              } else if constexpr (std::is_same_v<LhsT, RhsT>) {
                  if constexpr (std::is_same_v<LhsT, double>) {
                      return (lhs == rhs) ? std::partial_ordering::equivalent : (lhs < rhs ? std::partial_ordering::less : std::partial_ordering::greater);
                  } else {
                      return lhs <=> rhs; // Porovnání stejného typu (uint64_t vs uint64_t)
                  }
              } else {
                  // Různé typy - považuj double za větší než uint64_t
                  return std::is_same_v<LhsT, double> ? std::partial_ordering::greater : std::partial_ordering::less;
              }
          },
          fitness, other.fitness);
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
