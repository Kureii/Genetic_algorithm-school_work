/*
 * Created by kureii on 10/31/24.
 */

#pragma once

#include <mutex>
#include <vector>

#include "input_structure.h"
#include "output_structure.h"

class Exporter {
 public:
  Exporter() = delete;
  ~Exporter() = delete;
  Exporter(const Exporter&) = delete;
  Exporter& operator=(const Exporter&) = delete;

  static void ProcessDimensionsBitInputs(input_structure input_struct,
      size_t iterations,
      std::vector<std::vector<output_structure>>& output_structures_MaxOne,
      std::vector<std::vector<output_structure>>& output_structures_LeadingOne,
      std::vector<std::vector<output_structure>>& output_structures_TradingOne,
      std::vector<std::vector<output_structure>>& output_structures_MaxZero);

  static void ExportBitInputs(const input_structure& input_struct,
      const std::string& output_file, uint64_t iterations,
      const std::vector<uint64_t>& dimensions);

  static double compute_average(const std::vector<double>& data);

  static void WriteToFileBitInputs(const std::string& output_file,
  const std::vector<uint64_t>& dimensions,
      std::vector<std::vector<output_structure>>& output_structures_MaxOne,
      std::vector<std::vector<output_structure>>& output_structures_LeadingOne,
      std::vector<std::vector<output_structure>>& output_structures_TradingOne,
      std::vector<std::vector<output_structure>>& output_structures_MaxZero);

  static std::string CreateConvergencePlot(std::vector<double>& avg_convergence,
      const std::vector<output_structure>& dim_outputs,
      const std::string& fitness_name, uint64_t dim);

  static std::string PrintMaxFitness(
      const std::vector<output_structure>& dim_outputs,
      std::vector<double>& fitness_values);
  static std::string PrintMinFitness(
      const std::vector<output_structure>& dim_outputs,
      std::vector<double>& fitness_values);

  static void ProcessDimensionsRealInputs(input_structure_t input_struct,
    mapping_structure_t mapping_struct,
    size_t iterations,
    std::vector<std::vector<output_structure_t>>& output_structures_Sphere,
    std::vector<std::vector<output_structure_t>>& output_structures_Schwefel,
    std::vector<std::vector<output_structure_t>>& output_structures_Rosenbrock);

  static void ExportRealInputs(const input_structure_t& input_struct,
    const std::string& output_file, uint64_t iterations,
    const std::vector<uint64_t>& dimensions);

  static void WriteToFileRealInputs(const std::string& output_file,
    const std::vector<uint64_t>& dimensions,
    std::vector<std::vector<output_structure_t>>& output_structures_Sphere,
    std::vector<std::vector<output_structure_t>>& output_structures_Schwefel,
    std::vector<std::vector<output_structure_t>>& output_structures_Rosenbrock,
    mapping_method_t mapping_method);

static std::string CreateConvergencePlot(
    const std::vector<double>& avg_convergence,
    const std::string& fitness_name, const uint64_t dim, const std::string& color = "purple");

  static std::string PrintMaxFitness(
    const std::vector<output_structure_t>& dim_outputs,
    const std::vector<double>& fitness_values);

  static std::string PrintMinFitness(
    const std::vector<output_structure_t>& dim_outputs,
    const std::vector<double>& fitness_values);
};

