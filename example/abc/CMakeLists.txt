cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(abc)

# Build target module
add_executable(abc main.cpp)
target_include_directories(abc PRIVATE .)

# Link s0s
target_include_directories(abc PRIVATE ../../include)
target_include_directories(abc PRIVATE ../../../../thirdparty/eigen)
