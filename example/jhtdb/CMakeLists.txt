cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(test)

# build thirdparty library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../thirdparty/getjhtdb ${CMAKE_CURRENT_BINARY_DIR}/getjhtdb)

# build target module
add_executable(test main.cpp)
# include
target_include_directories(test PRIVATE .)
target_include_directories(test PRIVATE ../../include)
target_include_directories(test PRIVATE ../../thirdparty)
target_include_directories(test PRIVATE ../../../../thirdparty/eigen)
# link libraries
target_link_libraries(test getjhtdb)
