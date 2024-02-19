cmake_minimum_required(VERSION 3.16)
project(quantum)

find_package(OpenCL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

# Here is the path to the hpp
include_directories(ADD_POCL_INCLUDES)

add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_HPP_ENABLE_EXCEPTIONS")

add_executable(print.out main.cpp)
target_link_libraries(print.out OpenCL::OpenCL)
