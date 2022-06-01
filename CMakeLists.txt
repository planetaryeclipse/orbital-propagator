cmake_minimum_required(VERSION 3.22)
project(universal_propagator C)

set(CMAKE_C_STANDARD 99)

find_package(Catch2 3 REQUIRED)

add_library(universal_propagator include/stumpff.h src/stumpff.c include/anomaly.h src/anomaly.c)
target_include_directories(universal_propagator PUBLIC "include/")

project(universal_propagator CXX)

add_executable(test test.cpp)
target_include_directories(test PUBLIC "include/")
target_link_libraries(test PUBLIC universal_propagator Catch2::Catch2WithMain)