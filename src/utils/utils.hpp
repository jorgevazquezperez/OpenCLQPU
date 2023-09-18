#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <CL/opencl.hpp>

using Sources = std::vector<std::string>;

using UINT = unsigned int;

typedef struct{
    UINT qubit_value[3];
    double prob;
} COUNT;

void add_kernel_files(Sources &src, std::string relative_path, std::string main_cl);
void print_probabilities(COUNT *counts, size_t size);
cl::Device search_device();