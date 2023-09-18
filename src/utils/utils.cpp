#include "utils.hpp"

std::string get_source(const std::string &name) {

    std::ifstream file(name);
    if (!file.is_open()) {
        std::cerr << name + " no found!" << '\n';
        exit(-1);
    }

    return {std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>())};
}

void add_kernel_files(Sources &src, std::string relative_path, std::string main_cl){
    src.push_back(get_source(relative_path + "/src/malloc/clIndexedQueue.cl"));
    src.push_back(get_source(relative_path + "/src/malloc/kma.cl"));
    src.push_back(get_source(relative_path + "/src/quantum/qulacs.cl"));
    src.push_back(get_source(relative_path + "/src/utils/utils.cl"));
    src.push_back(get_source(main_cl));
}

void print_probabilities(COUNT *counts, size_t size){
    std::cout << "Probabilities:" << "\n";
    for(int i = 0; i < size; i++){
        for (int j = 2; j >= 0; j--){
            std::cout << counts[i].qubit_value[j] << " ";
        }
        std::cout << ": " << counts[i].prob << "\n";
    }
}

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
            std::cerr << "  " << device.getInfo<CL_DEVICE_EXTENSIONS>() << "(" << device.getInfo<CL_DEVICE_VERSION>() << ")"
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