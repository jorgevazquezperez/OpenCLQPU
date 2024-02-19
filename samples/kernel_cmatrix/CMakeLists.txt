cmake_minimum_required(VERSION 3.16)
project(quantum)

find_package(OpenCL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

set(POCL_INCLUDES ADD_POCL_INCLUDES)
set(KMA_INCLUDES ADD_SOURCE_DIRECTORY)
set(HEADERS ${POCL_INCLUDES} ${KMA_INCLUDES})

add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_HPP_ENABLE_EXCEPTIONS")

add_executable(matrix.out main.cpp)
target_include_directories(matrix.out PRIVATE ${HEADERS})
target_link_libraries(matrix.out OpenCL::OpenCL)

FILE(GLOB MyCSources *.c)
FILE(GLOB MyHeaders *.h)
set_source_files_properties(${MyCSources} ${MyHeaders} PROPERTIES LANGUAGE C)
target_sources(matrix.out PRIVATE malloc/kma.c)
