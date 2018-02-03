cd qhull

For GCC builds:

cd build
cmake ..

FOR CLANG builds:

mkdir build_clang
cd build_clang
CC=/usr/local/bin/clang CXX=/usr/local/bin/clang++ cmake ..
