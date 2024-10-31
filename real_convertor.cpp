/*
 * Created by kureii on 10/30/24.
 */

#include "real_convertor.h"

#include <cmath>

#include "real_structures.h"

uint64_t RealConvertor::ExtractBits(
    const ChromosomeArray<uint64_t>& chromosome, std::size_t bit_index,
    int bits_per_variable) {
  std::size_t start_word = bit_index / 64;
  std::size_t end_word = (bit_index + bits_per_variable - 1) / 64;
  int start_offset = bit_index % 64;

  uint64_t value = 0;

  if (start_word == end_word) {
    // Všechny bity jsou v jednom slově
    uint64_t word = chromosome[start_word];
    uint64_t mask = (1ULL << bits_per_variable) - 1;
    value = (word >> start_offset) & mask;
  } else {
    // Bity se rozkládají přes dvě slova
    uint64_t low_part = chromosome[start_word] >> start_offset;
    int bits_in_low_part = 64 - start_offset;
    uint64_t high_part = chromosome[end_word] & ((1ULL << (bits_per_variable - bits_in_low_part)) - 1);
    value = (high_part << bits_in_low_part) | low_part;
  }

  return value;
}

std::vector<double> RealConvertor::ConvertChromosomeToReal(
    const ChromosomeArray<uint64_t>& chromosome,
    const mapping_structure_t& mapping_structure,
    std::size_t problem_dimension) {
  std::vector<double> variables(problem_dimension);
    int bits_per_variable = mapping_structure.bits_per_variable;

    std::size_t total_bits = bits_per_variable * problem_dimension;
    std::size_t chromosome_bits = chromosome.size() * 64;
    if (total_bits > chromosome_bits) {
        throw std::invalid_argument("Chromozom nemá dostatek bitů");
    }

    std::size_t bit_index = 0;

    for (std::size_t i = 0; i < problem_dimension; ++i) {
        uint64_t bits_value = ExtractBits(chromosome, bit_index, bits_per_variable);
        double variable = 0.0;

        switch (mapping_structure.mapping_method) {
          case mapping_method_t::BIT_AS_DOUBLE: {
                union {
                    uint64_t u64;
                    double d;
                } u;
                u.u64 = bits_value;
                variable = u.d;
                break;
            }
            case mapping_method_t::FIXED_POINT: {
                // Pevná desetinná čárka
                int64_t int_value = static_cast<int64_t>(bits_value);
                double scaling_factor = 1e6; // Upravit dle potřeby
                variable = static_cast<double>(int_value) / scaling_factor;
                break;
            }
            case mapping_method_t::BINARY_CODED_DECIMAL: {
                int num_digits = bits_per_variable / 4;
                int64_t int_value = 0;
                for (int d = 0; d < num_digits; ++d) {
                    int digit = (bits_value >> (4 * d)) & 0xF;
                    if (digit > 9) digit = 9; // Omezení na desetinné číslice
                    int_value += digit * static_cast<int64_t>(std::pow(10, d));
                }
                variable = static_cast<double>(int_value);
                break;
            }
            case mapping_method_t::MAPPED_RANGE: {
                uint64_t max_int_value = (1ULL << bits_per_variable) - 1;
                variable = mapping_structure.min_value + (static_cast<double>(bits_value) / max_int_value) * (mapping_structure.max_value - mapping_structure.min_value);
                break;
            }
            default:
                throw std::invalid_argument("Neplatná metoda mapování");
        }

        variables[i] = variable;
        bit_index += bits_per_variable;
    }

    return variables;
}

