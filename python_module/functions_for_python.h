/*
 * Created by kureii on 10/18/24.
*/

#pragma once

#include "../input_structure.h"
#include "../output_structure.h"
#include "../real_structures.h"

extern output_structure_t GeneticAlgorithm_MaxOne(input_structure_t input_structure);

extern output_structure_t GeneticAlgorithm_LeadingOne(input_structure_t input_structure);

extern output_structure_t GeneticAlgorithm_TradingOne(input_structure_t input_structure);

extern output_structure_t GeneticAlgorithm_MaxZero(input_structure_t input_structure);

extern output_structure_t GeneticAlgorithm_Sphere(input_structure_t input_structure, mapping_structure_t mapping_structure);

extern output_structure_t GeneticAlgorithm_Schwefel(input_structure_t input_structure, mapping_structure_t mapping_structure);

extern output_structure_t GeneticAlgorithm_Rosenbrock(input_structure_t input_structure, mapping_structure_t mapping_structure);
