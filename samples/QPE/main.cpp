#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include <CL/opencl.hpp>

#include "../../src/utils/utils.hpp"

extern "C" {
   #include "../../src/malloc/kma.h"
}

int main(int argc, char *argv[]) {
    cl::Device device = search_device();
    cl::Context context(device);
    cl::CommandQueue queue(context, device);
    
    // Create the heap
    Sources src;
    add_kernel_files(src, "../../../", "../qpe.cl");
	cl::Program program(context, src);

    try {
        const char *includes = "-I../../../src/malloc";
        program.build(includes);
    } catch (cl::BuildError &err) {
        std::cerr << err.getBuildLog()[0].second << '\n';
        exit(-1);
    }

    try {
        size_t n_counts = 3;
        size_t count_dim = pow(2, n_counts); 
        COUNT counts[count_dim];

        // I create the heap
        cl::Buffer heap(kma_create(device.get(), context.get(), queue.get(), program.get(), 2048));
        cl::Buffer countsbuf(context, CL_MEM_READ_WRITE, sizeof(COUNT) * count_dim);
        
        cl::Kernel qpe(program, "qpe", nullptr);
        qpe.setArg(0, heap);
        qpe.setArg(1, countsbuf);
        qpe.setArg(2, &n_counts);

        queue.enqueueWriteBuffer(countsbuf, CL_TRUE, 0, sizeof(COUNT) * count_dim, counts);
        queue.enqueueNDRangeKernel(qpe, cl::NDRange(1), cl::NDRange(1));
        queue.enqueueReadBuffer(countsbuf, CL_TRUE, 0, sizeof(COUNT) * count_dim, counts);
        queue.finish();

        print_probabilities(counts, count_dim);
    } catch (cl::Error &err) {
        std::cerr << err.what() << '\n';
        exit(-1);
    }
    return 0;
}