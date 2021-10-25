cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(jhtdb)

list(APPEND THIRDPARTIES "../../thirdparty/turblib/soapC.c" "../../thirdparty/turblib/soapClient.c" "../../thirdparty/turblib/stdsoap2.c" "../../thirdparty/turblib/turblib.c")

# Build target module
add_executable(jhtdb main.cpp ${THIRDPARTIES})
target_include_directories(jhtdb PRIVATE .)

# Link s0s
target_include_directories(jhtdb PRIVATE ../../include)
target_include_directories(jhtdb PRIVATE ../../../../thirdparty/eigen)
target_include_directories(jhtdb PRIVATE ../../thirdparty/turblib)
