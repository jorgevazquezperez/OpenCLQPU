cmake_minimum_required(VERSION 3.16)
project(quantum)

find_package(OpenCL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

#set(POCL_INCLUDES /home/usc/ec/jvp/Proyectos/OpenQCL/src/PoCLQPU/POCL/pocl-3.0/include/hpp)
set(POCL_INCLUDES /mnt/netapp1/Optcesga_FT2_RHEL7/2020/software/Core/pocl/3.1/include)
set(KMA_INCLUDES /home/usc/ec/jvp/Proyectos/OpenQCL/samples/kernel_malloc/malloc)
set(HEADERS ${POCL_INCLUDES} ${KMA_INCLUDES})

add_definitions("-DCL_HPP_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_TARGET_OPENCL_VERSION=210")
add_definitions("-DCL_HPP_ENABLE_EXCEPTIONS")

add_executable(malloc.out main.cpp)
target_include_directories(malloc.out PRIVATE ${HEADERS})
target_link_libraries(malloc.out OpenCL::OpenCL)

FILE(GLOB MyCSources *.c)
FILE(GLOB MyHeaders *.h)
set_source_files_properties(${MyCSources} ${MyHeaders} PROPERTIES LANGUAGE C)
target_sources(malloc.out PRIVATE malloc/kma.c)

