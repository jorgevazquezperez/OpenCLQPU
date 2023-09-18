cmake_minimum_required(VERSION 3.16)
project(quantum)

find_package(OpenCL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

set(POCL_INCLUDES /mnt/netapp1/Optcesga_FT2_RHEL7/2020/software/Core/pocl/3.1/include)
set(SRC_INCLUDES /home/usc/ec/jvp/Proyectos/OpenQCL/src)
set(HEADERS ${POCL_INCLUDES} ${KMA_INCLUDES})

add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_HPP_ENABLE_EXCEPTIONS")

add_executable(${MY_VARIABLE} main.cpp ${SRC_INCLUDES}/utils/utils.cpp)
target_include_directories(${MY_VARIABLE} PRIVATE ${HEADERS})
target_link_libraries(${MY_VARIABLE} OpenCL::OpenCL)

FILE(GLOB MyCSources *.c)
FILE(GLOB MyHeaders *.h)
set_source_files_properties(${MyCSources} ${MyHeaders} PROPERTIES LANGUAGE C)
target_sources(${MY_VARIABLE} PRIVATE ${SRC_INCLUDES}/malloc/kma.c)

