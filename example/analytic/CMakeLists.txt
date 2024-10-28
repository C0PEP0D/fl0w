cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(analytic)

# Build target module
add_executable(analytic main.cpp)
target_include_directories(analytic PRIVATE .)

# Link s0s
target_include_directories(analytic PRIVATE ../../include)
target_include_directories(analytic PRIVATE ../../../../thirdparty/eigen)
