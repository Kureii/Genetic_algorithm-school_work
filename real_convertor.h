/*
 * Created by kureii on 10/30/24.
 */

#pragma once
#include <cstdint>
#include <iostream>
#include <vector>

#include "chromosome_array.h"
#include "real_structures.h"

class RealConvertor {
 public:
  RealConvertor() = default;
  static uint64_t ExtractBits(const ChromosomeArray<uint64_t>& chromosome,
      std::size_t bit_index, int bits_per_variable);
  static std::vector<double> ConvertChromosomeToReal(
      const ChromosomeArray<uint64_t>& chromosome,
      const mapping_structure_t& mapping_structure,
      std::size_t problem_dimension);
};
