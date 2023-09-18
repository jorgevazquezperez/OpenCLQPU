#include <iostream>
#include <fstream>

#include <CL/opencl.hpp>

cl::Device search_device() {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty()) {
        std::cerr << "No platforms found!" << '\n';
        exit(-1);
    }

    std::vector<cl::Device> devices;
    cl::Device first_device;
    bool set = false;
    std::cerr << "Devices:"<< '\n';
    for (auto &platform: platforms) {
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
        for (auto &device: devices) {
            std::cerr << "  " << device.getInfo<CL_DEVICE_NAME>() << "(" << device.getInfo<CL_DEVICE_VERSION>() << ")"
                      << '\n';
            if (!set) {
                set = true;
                first_device = device;
            }
        }
    }
    if (!set) {
        std::cerr << "No devices found!" << '\n';
        exit(-1);
    }
    return first_device;
}

std::string getSource(const std::string &name) {

    std::ifstream file(name);
    if (!file.is_open()) {
        std::cerr << name + " no found!" << '\n';
        exit(-1);
    }

    return {std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>())};
}


int main(int argc, char *argv[]) {
    cl::Device device = search_device();
    std::string src = getSource("../kernel.cl");

    cl::Context context(device);
    cl::Program program(context, src);

    try {
        program.build();
    } catch (cl::BuildError &err) {
        std::cerr << err.getBuildLog()[0].second << '\n';
        exit(-1);
    }

    try {
        cl::CommandQueue queue(context, device);
        cl::Kernel kernel(program, "PrintingOnKernel", nullptr);

        queue.enqueueNDRangeKernel(kernel, cl::NDRange(1), cl::NDRange(1));
        queue.finish();
    } catch (cl::Error &err) {
        std::cerr << err.what() << '\n';
        exit(-1);
    }
    return 0;
}
