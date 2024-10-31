#include <format>
#include <vector>

#include "exporter.h"

int main() {
  input_structure input_struct{0.02, 0.005, 0.085, 0.1, 0.5, 0.75, 15, 10};

  input_struct.fitness_rating_count = 30;

  const std::vector<uint64_t> dimensions{10, 30, 100};
  constexpr size_t iterations = 10;

  Exporter::ExportBitInputs(
      input_struct, "BitInputProblemsStatistics", iterations, dimensions);

  return 0;
}
