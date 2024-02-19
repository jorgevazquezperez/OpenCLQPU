# OpenCL QPU
Project which aims to add a QPU as part of OpenCL standard.

## Steps to test the program
In the `samples` folder there is the QPE and Shor examples. The only thing needed to do is to **modify the CMakeLists.txt** and change the `ADD_POCL_INCLUDES` for the path of the includes of PoCL in your device and `ADD_INCLUDE_DIRECTORY` for the path to the `src` file of this repository. Then run `run.sh` and get the results.

## Current status
As for now, QPU in OpenCL is just added simulated using CPU as a device (or every other classical device desired, but in this case CPU was chosen for the execution of the examples). A limited version of QPE and Shor algorithms can be found in `samples` folder.

## Following steps
Create a QPU driver to be able to execute code using a QPU as a device.
