cmake_minimum_required(VERSION 3.16)
project(quantum)

find_package(OpenCL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

#include_directories("/mnt/netapp1/Optcesga_FT2_RHEL7/2020/software/Core/pocl/3.1/include")
include_directories("/home/usc/ec/jvp/Proyectos/OpenQCL/src/PoCLQPU/POCL/pocl-3.0/include/hpp")

add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_HPP_ENABLE_EXCEPTIONS")

add_executable(print.out main.cpp)
target_link_libraries(print.out OpenCL::OpenCL)