/*
 * Created by kureii on 10/30/24.
*/

#pragma once

enum class mapping_method_t {
  BIT_AS_DOUBLE,
  FIXED_POINT,
  BINARY_CODED_DECIMAL,
  MAPPED_RANGE
};

struct mapping_structure_t {
  mapping_method_t mapping_method;
  double min_value;
  double max_value;
  int bits_per_variable;
};

using mapping_structure_t = struct mapping_structure_t;
