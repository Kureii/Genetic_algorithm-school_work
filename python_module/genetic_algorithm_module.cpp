/*
 * Created by kureii on 10/18/24.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "functions_for_python.h"
#include "../my_array_encapsulation.h"
#include "../input_structure.h"
#include "../output_structure.h"

namespace py = pybind11;

template <typename T>
void bind_myArrayEncapsulation(py::module &m, const std::string &name) {
    py::class_<myArrayEncapsulation<T>>(m, name.c_str())
        .def(py::init<>())
        .def(py::init<uint64_t>())
        .def("size", &myArrayEncapsulation<T>::size)
        .def("__getitem__", [](const myArrayEncapsulation<T> &a, size_t i) {
            if (i >= a.size()) throw py::index_error();
            return a[i];
        })
        .def("__setitem__", [](myArrayEncapsulation<T> &a, size_t i, T v) {
            if (i >= a.size()) throw py::index_error();
            a[i] = v;
        })
        .def("__len__", &myArrayEncapsulation<T>::size);
}

PYBIND11_MODULE(genetic_algorithm_module, m) {
    // Expose input_structure_t
    py::class_<input_structure_t>(m, "InputStructure")
        .def(py::init<>())
        .def_readwrite("mutation_probability", &input_structure_t::mutation_probability)
        .def_readwrite("elitism_random_select", &input_structure_t::elitism_random_select)
        .def_readwrite("elitism_best_select", &input_structure_t::elitism_best_select)
        .def_readwrite("crossing_full_random_probability", &input_structure_t::crossing_full_random_probability)
        .def_readwrite("first_parent_select_range", &input_structure_t::first_parent_select_range)
        .def_readwrite("second_parent_select_range", &input_structure_t::second_parent_select_range)
        .def_readwrite("initial_population_size", &input_structure_t::initial_population_size)
        .def_readwrite("population_size", &input_structure_t::population_size)
        .def_readwrite("problem_dimension", &input_structure_t::problem_dimension)
        .def_readwrite("fitness_rating_count", &input_structure_t::fitness_rating_count);

    // Expose output_structure_t
    py::class_<output_structure_t>(m, "OutputStructure")
        .def(py::init<>())
        .def_readwrite("convergence", &output_structure_t::convergence)
        .def_readwrite("result", &output_structure_t::result);

    // Expose myArrayEncapsulation if required
    bind_myArrayEncapsulation<uint64_t>(m, "MyArrayEncapsulationUint64");

    // Expose GeneticAlgorithm functions
    m.def("GeneticAlgorithm_MaxOne", &GeneticAlgorithm_MaxOne, "Run the MaxOne genetic algorithm", py::arg("input_structure"));
    m.def("GeneticAlgorithm_LeadingOne", &GeneticAlgorithm_LeadingOne, "Run the LeadingOne genetic algorithm", py::arg("input_structure"));
    m.def("GeneticAlgorithm_TradingOne", &GeneticAlgorithm_TradingOne, "Run the TradingOne genetic algorithm", py::arg("input_structure"));

    // Extend with additional binary optimization functions
    m.def("GeneticAlgorithm_MaxZero", &GeneticAlgorithm_MaxZero, "Run the MaxZero genetic algorithm", py::arg("input_structure"));
}
