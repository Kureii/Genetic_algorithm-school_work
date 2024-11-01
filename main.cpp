#include <format>
#include <vector>

#include "exporter.h"
#include "functions_for_python.h"

int main() {
  input_structure input_struct{0.02, 0.005, 0.085, 0.1, 0.5, 0.75, 15, 10};

  input_struct.fitness_rating_count = 10000;

  const std::vector<uint64_t> dimensions{10};
  constexpr size_t iterations = 10;

  /*Exporter::ExportBitInputs(
      input_struct, "BitInputProblemsStatistics.md", iterations, dimensions);*/

  input_struct.mapping_structure = mapping_structure_t {
    mapping_method_t::BIT_AS_DOUBLE,
    -500,500, 64,54
  };
  // input_struct.problem_dimension = 10;
  //
  // GeneticAlgorithm_Sphere(input_struct, input_struct.mapping_structure.value());

  Exporter::ExportRealInputs(input_struct, "RealInputProblemsStatisticsBitAsDouble.md", iterations, dimensions);
  input_struct.mapping_structure.value().mapping_method = mapping_method_t::FIXED_POINT;
  Exporter::ExportRealInputs(input_struct, "RealInputProblemsStatisticsFixPoint.md", iterations, dimensions);
  input_struct.mapping_structure.value().mapping_method = mapping_method_t::MAPPED_RANGE;
  Exporter::ExportRealInputs(input_struct, "RealInputProblemsStatisticsMap.md", iterations, dimensions);
  input_struct.mapping_structure.value().mapping_method = mapping_method_t::BINARY_CODED_DECIMAL;
  Exporter::ExportRealInputs(input_struct, "RealInputProblemsStatisticsBinDec.md", iterations, dimensions);

  return 0;
}
