cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(kinematic)

# Build target module
add_executable(kinematic main.cpp)
target_include_directories(kinematic PRIVATE .)

# Link s0s
target_include_directories(kinematic PRIVATE ../../include)
target_include_directories(kinematic PRIVATE ../../../../thirdparty/eigen)
