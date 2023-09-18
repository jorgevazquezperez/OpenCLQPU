rm -r build/
cmake . -B build
cd build
make
./malloc.out