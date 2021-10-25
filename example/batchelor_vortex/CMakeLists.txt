cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(batchelorvortex)

# Build target module
add_executable(batchelorvortex main.cpp)
target_include_directories(batchelorvortex PRIVATE .)

# Link s0s
target_include_directories(batchelorvortex PRIVATE ../../include)
target_include_directories(batchelorvortex PRIVATE ../../../../thirdparty/eigen)
