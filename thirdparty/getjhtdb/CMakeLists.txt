cmake_minimum_required(VERSION 3.11)

project(getjhtdb)

# build library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cjson ${CMAKE_CURRENT_BINARY_DIR}/cjson)

# build target
add_library(getjhtdb getData.c)
target_include_directories(getjhtdb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
# link libraries
target_link_libraries(getjhtdb cjson curl m)
