# CMakeLists.txt for building the Crack fuzz test

# Minimum version of CMake required
cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

# Project name
project(CrackFuzzTest LANGUAGES CXX)

# Force Clang as the C++ compiler
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
#message(STATUS "Compile = ${CMAKE_CXX_COMPILER_ID}")

# Set C++ standard to C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Prefer strict standard compliance

# Find Google Test based on include and lib directories
# User should set GTEST_ROOT to the base directory of their GTest installation,
# e.g., cmake -DGTEST_ROOT=/home/banerjee/ParSim/Vaango/dbg_install ..
# If GTEST_ROOT is not set, CMake will search in standard locations.
# Try to find the GTest include directory
find_path(GTEST_INCLUDE_DIR
          NAMES gtest/gtest.h
          HINTS ${GTEST_ROOT}/include
          PATH_SUFFIXES include
          DOC "Path to Google Test include directory"
          REQUIRED
)

# Try to find the GTest libraries (gtest_main and gtest)
find_library(GTEST_MAIN_LIBRARY
             NAMES gtest_main gtest_maind # 'gtest_maind' for debug builds on some systems
             HINTS ${GTEST_ROOT}/lib
             PATH_SUFFIXES lib lib64
             DOC "Path to gtest_main library"
             REQUIRED
)

find_library(GTEST_LIBRARY
             NAMES gtest gtestd # 'gtestd' for debug builds on some systems
             HINTS ${GTEST_ROOT}/lib
             PATH_SUFFIXES lib lib64
             DOC "Path to gtest library"
             REQUIRED
)

# Check if both libraries were found
if (NOT GTEST_MAIN_LIBRARY OR NOT GTEST_LIBRARY)
    message(FATAL_ERROR "Could not find both gtest and gtest_main libraries. "
                        "Please set GTEST_ROOT or ensure libraries are in standard paths.")
endif()


# Check for the compiler and add appropriate sanitizer flags
# libFuzzer is primarily a Clang feature. AddressSanitizer (ASan) is widely supported.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Clang compiler detected. Enabling libFuzzer and AddressSanitizer.")
    add_compile_options(-fsanitize=fuzzer -fsanitize=address)
    add_link_options(-fsanitize=fuzzer -fsanitize=address)

    # Add libFuzzer specific include directory.
    # On many Linux systems, this might be /usr/lib/llvm-<version>/include or similar.
    # You might need to adjust this path based on your Clang/LLVM installation.
    # If using a custom LLVM build, this would be <llvm-build-dir>/include.
    #find_path(LIBFUZZER_INCLUDE_DIR
    #          NAMES fuzzing/libfuzzer/libfuzzer_macro.h
    #          HINTS /usr/lib/llvm-*-dev/include # Common path for Debian/Ubuntu
    #                /usr/local/include          # Common local install path
    #                /usr/include                # Standard system include
    #                /usr/lib/llvm-<version>/include # Specific LLVM version include
    #          PATH_SUFFIXES llvm/include
    #          DOC "Path to libFuzzer include directory"
    #          REQUIRED
    #)
    #target_include_directories(testReadEllipticCracks PRIVATE ${LIBFUZZER_INCLUDE_DIR})
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(WARNING "GNU C++ compiler (g++) detected. -fsanitize=fuzzer is typically a Clang-specific flag.")
    message(WARNING "Attempting to build with -fsanitize=address, but libFuzzer functionality might be limited or unavailable.")
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
else()
    message(WARNING "Unknown compiler detected. Sanitizers and fuzzer may not be enabled.")
endif()

# Add the executable for the fuzz test
# The source file is the one containing LLVMFuzzerTestOneInput and your Google Test cases.
add_executable(testReadEllipticCracks testReadEllipticCracks.cc)
target_include_directories(testReadEllipticCracks PRIVATE ${GTEST_INCLUDE_DIR})
target_link_libraries(testReadEllipticCracks PRIVATE ${GTEST_MAIN_LIBRARY} ${GTEST_LIBRARY} pthread)

add_executable(testDiscretizeEllipticCracks testDiscretizeEllipticCracks.cc)
target_include_directories(testDiscretizeEllipticCracks PRIVATE ${GTEST_INCLUDE_DIR})
target_link_libraries(testDiscretizeEllipticCracks PRIVATE ${GTEST_MAIN_LIBRARY} ${GTEST_LIBRARY} pthread)

add_executable(testCrackFrontNodeIndices testCrackFrontNodeIndices.cc)
target_include_directories(testCrackFrontNodeIndices PRIVATE ${GTEST_INCLUDE_DIR})
target_link_libraries(testCrackFrontNodeIndices PRIVATE ${GTEST_MAIN_LIBRARY} ${GTEST_LIBRARY} pthread)

# Example of how to add a custom command to run the fuzzer
# This is optional, you can also just run ./crack_fuzz_test from the build directory
add_custom_target(run_fuzzer
    COMMAND $<TARGET_FILE:testReadEllipticCracks>
    WORKING_DIRECTORY $<TARGET_FILE_DIR:testReadEllipticCracks>
    COMMENT "Running the fuzzer. Press Ctrl+C to stop."
)

