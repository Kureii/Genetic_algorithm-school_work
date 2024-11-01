/*
 * Created by kureii on 10/18/24.
*/

#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

using output_structure_t = struct output_structure {
  std::vector<uint64_t> convergence;
  std::vector<uint64_t> result;
};
