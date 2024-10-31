cmake_minimum_required(VERSION 3.28)
project(genetic_algorithm_module LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 23
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)


set(MATPLOTLIBCPP_DIR lib/matplotlibcpp)
include_directories(${MATPLOTLIBCPP_DIR})


set(sources
        gen_algo.cpp
        gen_algo.h
        fitness_function_results_struct.h
        chromosome_array.h
        chromosome_array.tpp
        input_structure.h
        output_structure.h
        python_module/functions_for_python.h
        python_module/functions_for_python.cpp
        real_structures.h
        real_convertor.cpp
        real_convertor.h
        exporter.cpp
        exporter.h
#        python_module/genetic_algorithm_module.cpp
)
#pybind11_add_module(${PROJECT_NAME} ${sources})

#target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} python_module)

add_executable(${PROJECT_NAME}
        main.cpp
        ${sources}
)

find_package(Python3 COMPONENTS Development NumPy REQUIRED)

target_include_directories(${PROJECT_NAME} PRIVATE
        python_module
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG)
else ()
    set( CMAKE_CXX_FLAGS  "-O3")
endif()

if(WIN32)
    # Windows-specific settings
    target_compile_definitions(${PROJECT_NAME} PRIVATE NOMINMAX)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})
