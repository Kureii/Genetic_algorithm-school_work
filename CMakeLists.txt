cmake_minimum_required(VERSION 3.29)
project(gan_algo)

set(CMAKE_CXX_STANDARD 23)

add_executable(${PROJECT_NAME}
        main.cpp
        gen_algo.cpp
        gen_algo.h
        fitness_function_results_struct.h
        my_array_encapsulation.tpp
        my_array_encapsulation.h
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG)
endif()
