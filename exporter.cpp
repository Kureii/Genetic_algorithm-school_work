/*
 * Created by kureii on 10/31/24.
 */

#include "exporter.h"

#include <format>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <thread>
#include <Python.h>

#include "functions_for_python.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::mutex output_mutex_;


double Exporter::compute_average(const std::vector<double>& data) {
  if (data.empty()) return 0.0;
  double sum = std::accumulate(data.begin(), data.end(), 0.0);
  return sum / static_cast<double>(data.size());
}

void Exporter::ProcessDimensionsBitInputs(input_structure input_struct,
    size_t iterations,
    std::vector<std::vector<output_structure>>& output_structures_MaxOne,
    std::vector<std::vector<output_structure>>& output_structures_LeadingOne,
    std::vector<std::vector<output_structure>>& output_structures_TradingOne,
    std::vector<std::vector<output_structure>>& output_structures_MaxZero) {
  std::vector<output_structure> dim_outputs_MaxOne(iterations);
  std::vector<output_structure> dim_outputs_LeadingOne(iterations);
  std::vector<output_structure> dim_outputs_TradingOne(iterations);
  std::vector<output_structure> dim_outputs_MaxZero(iterations);

#pragma omp parallel for
  for (size_t i = 0; i < iterations; ++i) {
    dim_outputs_MaxOne[i] = GeneticAlgorithm_MaxOne(input_struct);
  }
  std::cout << std::format(
      "MaxOne in {} dimensions finished\n", input_struct.problem_dimension);

#pragma omp parallel for
  for (size_t i = 0; i < iterations; ++i) {
    dim_outputs_LeadingOne[i] = GeneticAlgorithm_LeadingOne(input_struct);
  }
  std::cout << std::format(
      "LeadingOne in {} dimensions finished\n", input_struct.problem_dimension);

#pragma omp parallel for
  for (size_t i = 0; i < iterations; ++i) {
    dim_outputs_TradingOne[i] = GeneticAlgorithm_TradingOne(input_struct);
  }
  std::cout << std::format(
      "TradingOne in {} dimensions finished\n", input_struct.problem_dimension);

#pragma omp parallel for
  for (size_t i = 0; i < iterations; ++i) {
    dim_outputs_MaxZero[i] = GeneticAlgorithm_MaxZero(input_struct);
  }
  std::cout << std::format(
      "MaxZero in {} dimensions finished\n", input_struct.problem_dimension);

  {
    std::lock_guard lock(output_mutex_);
    output_structures_MaxOne.push_back(dim_outputs_MaxOne);
    output_structures_LeadingOne.push_back(dim_outputs_LeadingOne);
    output_structures_TradingOne.push_back(dim_outputs_TradingOne);
    output_structures_MaxZero.push_back(dim_outputs_MaxZero);
  }
}

void Exporter::ExportBitInputs(const input_structure& input_struct,
    const std::string& output_file, uint64_t iterations,
    const std::vector<uint64_t>& dimensions) {
  std::vector<std::vector<output_structure>> output_structures_MaxOne;
  std::vector<std::vector<output_structure>> output_structures_LeadingOne;
  std::vector<std::vector<output_structure>> output_structures_TradingOne;
  std::vector<std::vector<output_structure>> output_structures_MaxZero;

  {
    std::vector<std::jthread> threads;

    for (auto& dim : dimensions) {
      input_structure thread_input_struct = input_struct;
      thread_input_struct.initial_population_size *= dim / 3;
      thread_input_struct.population_size *= dim / 3;
      thread_input_struct.fitness_rating_count *= dim / 3;
      thread_input_struct.problem_dimension = dim;

      threads.emplace_back(ProcessDimensionsBitInputs, thread_input_struct,
          iterations, std::ref(output_structures_MaxOne),
          std::ref(output_structures_LeadingOne),
          std::ref(output_structures_TradingOne),
          std::ref(output_structures_MaxZero));
    }
  }

  WriteToFileBitInputs(output_file, dimensions, output_structures_MaxOne,
      output_structures_LeadingOne, output_structures_TradingOne,
      output_structures_MaxZero);
}

void Exporter::WriteToFileBitInputs(const std::string& output_file,
    const std::vector<uint64_t>& dimensions,
    std::vector<std::vector<output_structure>>& output_structures_MaxOne,
    std::vector<std::vector<output_structure>>& output_structures_LeadingOne,
    std::vector<std::vector<output_structure>>& output_structures_TradingOne,
    std::vector<std::vector<output_structure>>& output_structures_MaxZero) {
  struct FitnessFunction {
    const char* name;
    std::vector<std::vector<output_structure>>& outputs;
  };

  std::vector<FitnessFunction> fitness_functions = {
      {"MaxOne", output_structures_MaxOne},
      {"LeadingOne", output_structures_LeadingOne},
      {"TradingOne", output_structures_TradingOne},
      {"MaxZero", output_structures_MaxZero}};

  std::ofstream markdown_file(output_file.c_str());

  for (const auto& [fitness_name, outputs] : fitness_functions) {
    markdown_file << "# " << fitness_name << " Problem\n\n";

    for (size_t dim_index = 0; dim_index < dimensions.size(); ++dim_index) {
      auto dim = static_cast<int64_t>(dimensions[dim_index]);

      const std::vector<output_structure>& dim_outputs = outputs[dim_index];

      size_t convergence_length = 0;
      for (const auto& [convergence, result] : dim_outputs) {
        convergence_length = convergence.size() > convergence_length
                                 ? convergence.size()
                                 : convergence_length;
      }

      std::vector avg_convergence(convergence_length, 0.0);

      markdown_file << std::format(
          "![]({})\n\n", CreateConvergencePlot(
                             avg_convergence, dim_outputs, fitness_name, dim));

      std::vector<double> fitness_values;
      for (const auto& [convergence, result] : dim_outputs) {
        fitness_values.push_back(static_cast<double>(convergence.back()));
      }

      double best_fitness = *std::ranges::max_element(fitness_values);
      double worst_fitness = *std::ranges::min_element(fitness_values);
      double average_fitness =
          std::accumulate(fitness_values.begin(), fitness_values.end(), 0.0) /
          static_cast<double>(fitness_values.size());
      std::vector<double> sorted_fitness = fitness_values;
      std::ranges::sort(sorted_fitness);
      double median_fitness = sorted_fitness[sorted_fitness.size() / 2];
      double variance =
          std::accumulate(sorted_fitness.begin(), sorted_fitness.end(), 0.0,
              [average_fitness](const double sum, const double val) {
                return sum + (val - average_fitness) * (val - average_fitness);
              }) /
          static_cast<double>(fitness_values.size());
      double stddev_fitness = std::sqrt(variance);

      markdown_file << std::format("- Dimension: {}D\n", dim);
      markdown_file << std::format("- Best Fitness: {}\n", best_fitness);
      markdown_file << std::format("- Worst Fitness: {}\n", worst_fitness);
      markdown_file << std::format(
          "- Average Fitness: {:.2f}\n", average_fitness);
      markdown_file << std::format("- Median Fitness: {}\n", median_fitness);
      markdown_file << std::format(
          "- Standard Deviation: {}\n\n", stddev_fitness);
      markdown_file << PrintMaxFitness(dim_outputs, fitness_values);
      markdown_file << "\n\n";
      markdown_file << PrintMinFitness(dim_outputs, fitness_values);
      markdown_file << "`\n\n---\n\n";
    }
  }

  markdown_file.close();
}

std::string Exporter::CreateConvergencePlot(
    std::vector<double>& avg_convergence,
    const std::vector<output_structure>& dim_outputs,
    const std::string& fitness_name, const uint64_t dim) {
  for (const auto& [convergence, result] : dim_outputs) {
    for (size_t i = 0; i < convergence.size(); ++i) {
      avg_convergence[i] += static_cast<double>(convergence[i]);
    }
  }

  const size_t num_runs = dim_outputs.size();
  for (double& i : avg_convergence) {
    i /= static_cast<double>(num_runs);
  }

  std::vector<double> x(avg_convergence.size());
  std::vector<double> y(avg_convergence.size());

  for (size_t i = 0; i < avg_convergence.size(); ++i) {
    x[i] = static_cast<double>(i + 1);
    y[i] = avg_convergence[i];
  }

  plt::figure();
  plt::plot(x, y);
  const std::string title =
      std::format("Average Convergence for {} (Dimension {} )",
          std::string(fitness_name), std::to_string(dim));
  plt::title(title);
  plt::xlabel("Iteration");
  plt::ylabel("Fitness");
  std::string filename = std::format("convergence_{}_dim{}.png",
      std::string(fitness_name), std::to_string(dim));
  plt::save(filename);
  plt::close();

  return filename;
}

std::string Exporter::PrintMaxFitness(
    const std::vector<output_structure>& dim_outputs,
    std::vector<double>& fitness_values) {
  std::string result;
  if (!dim_outputs[std::distance(fitness_values.begin(),
                       std::ranges::max_element(fitness_values))]
           .result.empty()) {
    for (const auto& val :
        dim_outputs[std::distance(fitness_values.begin(),
                        std::ranges::max_element(fitness_values))]
            .result) {
      result = std::format("Best Result Vector:\n\n`{}`", val);
    }
  }
  return result;
}

std::string Exporter::PrintMinFitness(
    const std::vector<output_structure>& dim_outputs,
    std::vector<double>& fitness_values) {
  std::string result;
  if (!dim_outputs[std::distance(fitness_values.begin(),
                       std::ranges::min_element(fitness_values))]
           .result.empty()) {
    for (const auto& val :
        dim_outputs[std::distance(fitness_values.begin(),
                        std::ranges::min_element(fitness_values))]
            .result) {
      result = std::format("Best Result Vector:\n\n`{}`", val);
    }
  }
  return result;
}

void Exporter::ProcessDimensionsRealInputs(input_structure_t input_struct,
    mapping_structure_t mapping_struct,
    size_t iterations,
    std::vector<std::vector<output_structure_t>>& output_structures_Sphere,
    std::vector<std::vector<output_structure_t>>& output_structures_Schwefel,
    std::vector<std::vector<output_structure_t>>& output_structures_Rosenbrock) {
      std::vector<output_structure_t> dim_outputs_Sphere(iterations);
      std::vector<output_structure_t> dim_outputs_Schwefel(iterations);
      std::vector<output_structure_t> dim_outputs_Rosenbrock(iterations);

#pragma omp parallel for
      for (size_t i = 0; i < iterations; ++i) {
        dim_outputs_Sphere[i] = GeneticAlgorithm_Sphere(input_struct, mapping_struct);
      }
      std::cout << std::format(
          "Sphere in {} dimensions finished\n", input_struct.problem_dimension);

#pragma omp parallel for
      for (size_t i = 0; i < iterations; ++i) {
        dim_outputs_Schwefel[i] = GeneticAlgorithm_Schwefel(input_struct, mapping_struct);
      }
      std::cout << std::format(
          "Schwefel in {} dimensions finished\n", input_struct.problem_dimension);

#pragma omp parallel for
      for (size_t i = 0; i < iterations; ++i) {
        dim_outputs_Rosenbrock[i] = GeneticAlgorithm_Rosenbrock(input_struct, mapping_struct);
      }
      std::cout << std::format(
          "Rosenbrock in {} dimensions finished\n", input_struct.problem_dimension);

        {
        std::lock_guard<std::mutex> lock(output_mutex_);
        output_structures_Sphere.push_back(dim_outputs_Sphere);
        output_structures_Schwefel.push_back(dim_outputs_Schwefel);
        output_structures_Rosenbrock.push_back(dim_outputs_Rosenbrock);
        }
}


void Exporter::ExportRealInputs(const input_structure_t& input_struct,
    const std::string& output_file, uint64_t iterations,
    const std::vector<uint64_t>& dimensions) {
  std::vector<std::vector<output_structure_t>> output_structures_Sphere;
  std::vector<std::vector<output_structure_t>> output_structures_Schwefel;
  std::vector<std::vector<output_structure_t>> output_structures_Rosenbrock;

  {
    std::vector<std::jthread> threads;

    for (auto& dim : dimensions) {
      input_structure_t thread_input_struct = input_struct;
      thread_input_struct.problem_dimension = dim;
      thread_input_struct.fitness_rating_count*= dim * 2;

      mapping_structure_t thread_mapping_struct = input_struct.mapping_structure.value();
      // Adjust mapping_struct per dimension if needed

      threads.emplace_back(ProcessDimensionsRealInputs, thread_input_struct,
                thread_mapping_struct, iterations,
                std::ref(output_structures_Sphere),
                std::ref(output_structures_Schwefel),
                std::ref(output_structures_Rosenbrock));
    }
  }

  WriteToFileRealInputs(output_file, dimensions, output_structures_Sphere,
      output_structures_Schwefel, output_structures_Rosenbrock, input_struct.mapping_structure.value().mapping_method);
}


void Exporter::WriteToFileRealInputs(const std::string& output_file,
    const std::vector<uint64_t>& dimensions,
    std::vector<std::vector<output_structure_t>>& output_structures_Sphere,
    std::vector<std::vector<output_structure_t>>& output_structures_Schwefel,
    std::vector<std::vector<output_structure_t>>& output_structures_Rosenbrock,
    mapping_method_t mapping_method) {
  struct FitnessFunction {
    std::string name;
    std::vector<std::vector<output_structure_t>>& outputs;
    std::string color;
  };
  std::string way;

  switch (mapping_method) {
    using enum mapping_method_t;
    case FIXED_POINT: {
      way = "FIXED_POINT";
      break;
    }
    case BIT_AS_DOUBLE: {
      way = "BIT_AS_DOUBLE";
      break;
    }
    case MAPPED_RANGE: {
      way = "MAPPED_RANGE";
      break;
    }
    case BINARY_CODED_DECIMAL: {
      way = "BINARY_CODED_DECIMAL";
      break;
    }

  }

  std::vector<FitnessFunction> fitness_functions = {
      {std::format("Sphere_{}", way).c_str(), output_structures_Sphere,"purple"},
      {std::format("Schwefel_{}", way).c_str(), output_structures_Schwefel,"green"},
      {std::format("Rosenbrock_{}", way).c_str(), output_structures_Rosenbrock,"red"}};

  std::ofstream markdown_file(output_file.c_str());

  for (const auto& [fitness_name, outputs, color] : fitness_functions) {
    markdown_file << "# " << fitness_name << " Problem\n\n";

    for (size_t dim_index = 0; dim_index < dimensions.size(); ++dim_index) {
      auto dim = static_cast<int64_t>(dimensions[dim_index]);

      const std::vector<output_structure_t>& dim_outputs = outputs[dim_index];

      size_t convergence_length = 0;
      for (const auto& output : dim_outputs) {
        convergence_length = output.convergence.size() > convergence_length
                                 ? output.convergence.size()
                                 : convergence_length;
      }

      std::vector avg_convergence(convergence_length, 0.0);

      for (const auto& output : dim_outputs) {
        for (size_t i = 0; i < output.convergence.size(); ++i) {
          avg_convergence[i] += static_cast<double>(output.convergence[i]);
        }
      }

      const size_t num_runs = dim_outputs.size();
      for (double& val : avg_convergence) {
        val /= static_cast<double>(num_runs);
      }

      markdown_file << std::format(
          "![]({})\n\n", CreateConvergencePlot(
                             avg_convergence, fitness_name, dim, color));

      std::vector<double> fitness_values;
      for (const auto& output : dim_outputs) {
        fitness_values.push_back(static_cast<double>(output.convergence.back()));
      }

      double best_fitness = *std::ranges::max_element(fitness_values);
      double worst_fitness = *std::ranges::min_element(fitness_values);
      double average_fitness =
          std::accumulate(fitness_values.begin(), fitness_values.end(), 0.0) /
          static_cast<double>(fitness_values.size());
      std::vector<double> sorted_fitness = fitness_values;
      std::ranges::sort(sorted_fitness);
      double median_fitness = sorted_fitness[sorted_fitness.size() / 2];
      double variance =
          std::accumulate(sorted_fitness.begin(), sorted_fitness.end(), 0.0,
              [average_fitness](const double sum, const double val) {
                return sum + (val - average_fitness) * (val - average_fitness);
              }) /
          static_cast<double>(fitness_values.size());
      double stddev_fitness = std::sqrt(variance);

      markdown_file << std::format("- Dimension: {}D\n", dim);
      markdown_file << std::format("- Best Fitness: {}\n", best_fitness);
      markdown_file << std::format("- Worst Fitness: {}\n", worst_fitness);
      markdown_file << std::format(
          "- Average Fitness: {:.2f}\n", average_fitness);
      markdown_file << std::format("- Median Fitness: {}\n", median_fitness);
      markdown_file << std::format(
          "- Standard Deviation: {}\n\n", stddev_fitness);
      markdown_file << PrintMaxFitness(dim_outputs, fitness_values);
      markdown_file << "\n\n";
      markdown_file << PrintMinFitness(dim_outputs, fitness_values);
      markdown_file << "\n\n---\n\n";
    }
  }

  markdown_file.close();
}

std::string Exporter::CreateConvergencePlot(
    const std::vector<double>& avg_convergence,
    const std::string& fitness_name,
    const uint64_t dim,
    const std::string& color) {  // Přidán parametr pro barvu

  std::vector<double> x(avg_convergence.size());
  std::vector<double> y(avg_convergence.size());

  for (size_t i = 0; i < avg_convergence.size(); ++i) {
    x[i] = static_cast<double>(i + 1);
    y[i] = avg_convergence[i];
  }

  plt::figure();
  plt::plot(x, y);
  const std::string title =
      std::format("Average Convergence for {} (Dimension {} )",
          fitness_name, dim);
  plt::title(title);
  plt::xlabel("Generation");
  plt::ylabel("Fitness");
  std::string filename = std::format("convergence_{}_dim{}.png",
      fitness_name, dim);
  plt::save(filename);
  plt::close();

  return filename;
}

std::string Exporter::PrintMaxFitness(
    const std::vector<output_structure_t>& dim_outputs,
    const std::vector<double>& fitness_values) {
  auto max_iter = std::ranges::max_element(fitness_values);
  size_t index = std::distance(fitness_values.begin(), max_iter);
  const auto& best_output = dim_outputs[index];

  std::ostringstream oss;
  oss << "Best Result Vector:\n\n`";
  for (const auto& val : best_output.result) {
    oss << val << ", ";
  }
  oss << "`";
  return oss.str();
}

std::string Exporter::PrintMinFitness(
    const std::vector<output_structure_t>& dim_outputs,
    const std::vector<double>& fitness_values) {
  auto min_iter = std::ranges::min_element(fitness_values);
  size_t index = std::distance(fitness_values.begin(), min_iter);
  const auto& worst_output = dim_outputs[index];

  std::ostringstream oss;
  oss << "Worst Result Vector:\n\n`";
  for (const auto& val : worst_output.result) {
    oss << val << ", ";
  }
  oss << "`";
  return oss.str();
}

