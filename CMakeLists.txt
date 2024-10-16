cmake_minimum_required(VERSION 3.29)
project(gan_algo)

set(CMAKE_CXX_STANDARD 23)

add_executable(${PROJECT_NAME}
        main.cpp
        gen_algo.tpp
        gen_algo.h
        gen_structure.h
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG)
endif()
