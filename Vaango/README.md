Vaango
======

Material Point Method simulation of solids and Implicit Compressible Eulerian finite volume simulation of fluids.  Also includes coupled MPM+ICE and a test Peridynamics implementation.

Numerous material models are available.  

The source source is in the `src` directory and detailed documentation is in the `docs` directory.

The master has now been updated and partially tested - but only on Ubuntu 24.04 with gcc 13.3.  It requires c++20.  No backward compatibility is maintained.  For older versions of operating systems, use the single release version  from 2022.

## Building

### Debug build
A debug build can be created using:
```
mkdir dbg
cd dbg
cmake -DCMAKE_BUILD_TYPE=Debug ../src
```

Several packages may need to be installed before the build completes.  Use `apt get` to install those after iteratively running cmake.  To see which pacakages are needed, please see the `CMakeLists.txt` file in `src`.

### Release build
The release build is created using
```
mkdir opt
cd opt
cmake -DCMAKE_BUILD_TYPE=Release ../src
```

### Petsc+Hypre build
For building with PetSc and Hypre, use
```
mkdir dbg
cd dbg
cmake -DCMAKE_BUILD_TYPE=Debug -DPETSC=On -DHYPRE=On ../src
```

We typically use the Petsc and Hypre that is prebuilt for Ubuntu.

### CUDA+Kokkos build
For CUDA functionality, we use Kokkos.

The Ubuntu Kokkos library does not have CUDA enabled, so we have to build it ourselves using
```
cd /path/to/kokkos
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DKokkos_ENABLE_THREADS=ON \
 -DKokkos_ENABLE_CUDA=ON \
 -DKokkos_ENABLE_SERIAL=ON \
 -DKokkos_ARCH_NATIVE=ON \
 -DCMAKE_CXX_STANDARD=20 \
 -DBUILD_SHARED_LIBS=ON
cmake --build build
cmake --install build --prefix /path/to/kokkos/install
```

Then we build Vaango using
```
mkdir dbg
cd dbg
cmake -DCMAKE_BUILD_TYPE=Debug -DCUDA=On -DKokkos_ROOT=/path/to/kokkos/install ../src
```

