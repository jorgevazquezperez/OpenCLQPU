#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#define VECTOR_SIZE 1024

#include "kma.h"
#include "test/cl.h"

cl_device_id get_device(cl_platform_id platforms){
    //Get the devices list and choose the device you want to run on
    cl_uint num_devices;
    cl_int clStatus = clGetDeviceIDs( platforms, CL_DEVICE_TYPE_CPU, 0, NULL, &num_devices);
    cl_device_id device_list[num_devices];
    clStatus = clGetDeviceIDs(platforms, CL_DEVICE_TYPE_CPU, num_devices, device_list, NULL);

    return device_list[0];
}

cl_platform_id get_platform(){
    cl_uint num_platforms;
    cl_int clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
    cl_platform_id platforms[num_platforms];
    clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);

    return platforms[0];
}

int main(void) {
    // Setting up the platform layer
    cl_int clStatus;

    cl_platform_id platform = get_platform();    
    cl_device_id device = get_device(platform);
    cl_context context = clCreateContext( NULL, 1, &device, NULL, NULL, &clStatus);
    cl_command_queue command_queue = clCreateCommandQueueWithProperties(context, device, 0, &clStatus);

    // Create the heap
    char *src[2];
    src[0] = kernel_read("clIndexedQueue.cl");
	src[1] = kernel_read("kernel_malloc.cl");
	cl_program program = program_compile(platform, context, &device, 2, src);
	if(program < 0)
		return -1;

    cl_mem heap = kma_create(device, context, command_queue, program, 127);

    cl_kernel kernel = clCreateKernel(program, "test_malloc_outside", &clStatus);
    if(clStatus != CL_SUCCESS) {
		printf("KMA_test: Could not create size test kernel: %i\n", clStatus);
		return -1;
	}

    unsigned int iter = 5;
	clSetKernelArg(kernel, 0, sizeof(cl_mem), &heap);
    clSetKernelArg(kernel, 1, sizeof(unsigned int), &iter);

    clStatus = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, NULL, NULL, 0, NULL, NULL);

    // Clean up and wait for all the comands to complete.
    //clStatus = clFlush(command_queue);
    clStatus = clFinish(command_queue);

    // Finally release all OpenCL allocated objects and host buffers.
    clStatus = clReleaseKernel(kernel);
    clStatus = clReleaseProgram(program);
    clStatus = clReleaseCommandQueue(command_queue);
    clStatus = clReleaseContext(context);
    return 0;
}