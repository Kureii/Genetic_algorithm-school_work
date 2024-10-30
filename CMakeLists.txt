cmake_minimum_required(VERSION 3.28)
project(genetic_algorithm_module LANGUAGES CXX)

#add_subdirectory(pybind11)

set(CMAKE_CXX_STANDARD 23
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

set( CMAKE_CXX_FLAGS  "-O3")
set(sources
        gen_algo.cpp
        gen_algo.h
        fitness_function_results_struct.h
        my_array_encapsulation.tpp
        my_array_encapsulation.h
        input_structure.h
        output_structure.h
        python_module/functions_for_python.h
        python_module/functions_for_python.cpp
#        python_module/genetic_algorithm_module.cpp
)
#pybind11_add_module(${PROJECT_NAME} ${sources})

#target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} python_module)

add_executable(${PROJECT_NAME}
        main.cpp
        ${sources}
)
target_include_directories(${PROJECT_NAME} PRIVATE python_module)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG)
endif()

if(WIN32)
    # Windows-specific settings
    target_compile_definitions(genetic_algorithm_module PRIVATE NOMINMAX)
endif()
