# Vaango Build System Modernization and Malloc Alignment Update

This document summarizes the changes made to the Vaango build system and the core memory allocation library to support modern C++ standards and ensure compatibility with high-performance libraries like Eigen.

## 1. CMake Build System Modernization

The build system was overhauled to follow modern "Target-based CMake" practices (CMake 3.28+).

### Key Improvements:
- **Centralized Configuration**: Enhanced the `Vaango::Options` (alias `vaango_options`) interface library to manage project-wide include paths, compiler flags, and common dependencies.
- **Modern Target Linking**: Replaced manual `include_directories` and variable-based linking with idiomatic `target_link_libraries`.
    - **MPI**: Linked via `MPI::MPI_CXX`.
    - **LibXml2**: Linked via `LibXml2::LibXml2`.
    - **HDF5**: Transitioned to modern `hdf5::hdf5_cpp` / `HDF5::HDF5` targets.
    - **Eigen3**: Linked via `Eigen3::Eigen`.
    - **Threads**: Linked via `Threads::Threads`.
- **GTest Integration**: Updated unit tests to use the `GTest::gtest_main` target, removing redundant header paths and improving test discovery via `gtest_discover_tests`.
- **Consistency**: Standardized `if/else/endif` and `message` syntax to be lowercase and idiomatic.
- **Bug Fixes**: Corrected a bug where `DEF_NO_FORTRAN` was being overwritten, and removed brittle hardcoded compiler paths for Clang.

## 2. Core/Malloc Fixes and Alignment Updates

Critical fixes were applied to the custom memory allocator to resolve C++20 compatibility issues and memory alignment requirements.

### C++ Standard Compliance:
- **`new.cc` Update**: Removed dynamic exception specifications (`throw(std::bad_alloc)`) which were deprecated in C++11 and removed in C++17.
- **`noexcept` Usage**: Replaced `throw()` with `noexcept` for delete operators and nothrow new variants to align with the C++20 standard.

### 16-Byte Alignment for Eigen:
Eigen and other SIMD-optimized libraries require strict 16-byte (or higher) alignment. The previous allocator was returning 8-byte aligned pointers, causing assertion failures.
- **`AllocPriv.h` Changes**:
    - Applied `alignas(16)` to the `Tag` and `Sentinel` structures.
    - Added explicit padding members to ensure structure sizes are multiples of 16 bytes.
- **`Allocator.cc` Logic**:
    - Updated `fill_bin` to round up `tsize` (object size + overhead) to the nearest multiple of `ALIGN` (16).
    - Updated `alloc_big` to ensure big object allocations and their tail sentinels maintain 16-byte alignment.
- **Header Inclusion Fix**: Linked `Vaango_Core_Malloc` to `Vaango::Options` to ensure `sci_defs/bits_defs.h` (and other generated headers) are correctly found.

## 3. Library Refactoring
- **`MIGUtils.cc`**: Created a dedicated `MPM_ConstitutiveModel_FORTRAN_Utils` library for this source file. This ensures its dependencies (like `Core/Parallel`) are correctly managed and linked, resolving "Header Not Found" errors in tests that compile this file.

## 4. Verification Results
- Successfully configured and built the project using **Ninja** in the `dbg` directory.
- Verified that all unit tests pass, including those previously failing due to Eigen alignment assertions or missing HDF5 dependencies.
