# cmake/CheckCppFeatures.cmake
# This file contains functions to check for various system and compiler features.
# It is designed to be included by the main CMakeLists.txt.

# Specify the minimum required CMake version for this file.
#cmake_minimum_required(VERSION 3.28)
include_guard(GLOBAL)

# Function to check for various project features.
# This function sets cache variables and global compile definitions based on feature availability.
# It is designed to be called once from the top-level CMakeLists.txt.
function(vaango_check_features)
    message(STATUS "Performing feature checks for Vaango Project...")

    # Include CMake's built-in modules for feature checks.
    # CheckIncludeFiles is included for check_include_file_cxx.
    include(CheckIncludeFiles)
    # CheckCXXSourceCompiles is traditionally included, but direct try_compile is often preferred.
    # It's kept here as it was in the original, but try_compile is used explicitly for more complex checks.
    include(CheckCXXSourceCompiles)

    # Base directory for temporary CMake test projects.
    # Using a unique subdirectory for each try_compile to prevent "nested TRY_COMPILE" errors.
    set(VAANGO_TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/__vaango_cmake_tmp")
    file(MAKE_DIRECTORY ${VAANGO_TMP_DIR}) # Ensure the base temporary directory exists

    #-----------------------------------------------------------------------------
    # Check availability of semaphore.h
    # Using check_include_file_cxx, which is appropriate for C++ projects.
    #-----------------------------------------------------------------------------
    message(STATUS "Checking for semaphore.h...")
    check_include_file_cxx(semaphore.h HAVE_SEMAPHORE_H)
    if(HAVE_SEMAPHORE_H)
        message(STATUS "Found semaphore.h.")
    else()
        message(STATUS "semaphore.h not found.")
    endif()

    #-----------------------------------------------------------------------------
    # Check atomic builtins for threads support in GCC/MSVC
    # Using try_compile for a direct source compilation test.
    # This test verifies if compiler intrinsics for atomic operations are available.
    #-----------------------------------------------------------------------------
    message(STATUS "Checking for atomic builtins...")
    # Store the source code in a variable before passing it to try_compile.
    set(ATOMIC_CHECK_SOURCE "
#ifdef _MSC_VER
#include <windows.h> // Required for Windows-specific atomic functions
#endif
int main() {
#ifdef _MSC_VER
        volatile LONG val = 1;
        MemoryBarrier();              // Ensures memory ordering
        InterlockedCompareExchange(&val, 0, 1); // Atomic compare-and-swap
        InterlockedIncrement(&val);   // Atomic increment
        InterlockedDecrement(&val);   // Atomic decrement
#else
        volatile unsigned long val = 1;
        __sync_synchronize();         // Compiler barrier for memory ordering
        __sync_val_compare_and_swap(&val, 1, 0); // GCC/Clang intrinsic for compare-and-swap
        __sync_add_and_fetch(&val, 1);    // GCC/Clang intrinsic for atomic add
        __sync_sub_and_fetch(&val, 1);    // GCC/Clang intrinsic for atomic subtract
#endif
        return 0;
      }
    ")
    try_compile(
        GCC_HAS_ATOMICS               # CMake variable to store the result (TRUE if compiles, FALSE otherwise)
        SOURCE_FROM_VAR "atomic_check.cxx" ATOMIC_CHECK_SOURCE
        CMAKE_FLAGS "-DLINK_LIBRARIES=" # Prevents linking issues if libs are implicitly tried
        OUTPUT_VARIABLE ATOMIC_CHECK_OUTPUT # Stores the compiler output for debugging
    )

    if(NOT GCC_HAS_ATOMICS)
        message(STATUS "Warning: Atomic builtins are missing. Consider building thread-unsafe or providing alternative implementation.")
        # Set cache variables to indicate the chosen implementation files.
        set(VAANGO_REFCOUNT_IMPL "RefCounted_gcc.cc" CACHE INTERNAL "Implementation file for RefCounted based on atomic support.")
        set(VAANGO_ATOMIC_IMPL "AtomicCounter_gcc.cc" CACHE INTERNAL "Implementation file for AtomicCounter based on atomic support.")
    else()
        message(STATUS "Atomic builtins found.")
        set(VAANGO_REFCOUNT_IMPL "RefCounted.cc" CACHE INTERNAL "Implementation file for RefCounted based on atomic support.")
        set(VAANGO_ATOMIC_IMPL "AtomicCounter_default.cc" CACHE INTERNAL "Implementation file for AtomicCounter based on atomic support.")
    endif()

    # Determine time implementation based on operating system.
    if (WIN32)
        set(VAANGO_TIME_IMPL "Time_win32.cc" CACHE INTERNAL "Time implementation file for Windows.")
    else()
        set(VAANGO_TIME_IMPL "Time_unix.cc" CACHE INTERNAL "Time implementation file for Unix-like systems.")
    endif()
    message(STATUS "Using Time implementation: ${VAANGO_TIME_IMPL}")

    #-----------------------------------------------------------------------------
    # Check included header files and set configuration variables.
    # These variables are typically used to generate a configuration header file.
    #-----------------------------------------------------------------------------
    message(STATUS "Checking common headers...")

    # Check for <ext/algorithm> - a GNU extension often used in older GCC/Clang.
    set(EXT_ALGORITHM_CHECK_SOURCE "
#include <ext/algorithm>
#include <vector>
int main() {
        std::vector<int> vec; // Test a simple usage of STL with the header
        return 0;
      }
    ")
    try_compile(
        HAVE_EXT_ALGORITHM
        SOURCE_FROM_VAR "ext_algo_check.cxx" EXT_ALGORITHM_CHECK_SOURCE # Pass the VARIABLE NAME
        OUTPUT_VARIABLE EXT_ALGORITHM_CHECK_OUTPUT
    )
    if (HAVE_EXT_ALGORITHM)
        message(STATUS "ext/algorithm found.")
        # Add a global compile definition if this header is found.
        add_compile_definitions(HAVE_EXT_ALGORITHM)
    else()
        message(STATUS "ext/algorithm not found.")
    endif()

    # Check for standard C/C++ headers using check_include_file_cxx.
    # These are fundamental headers for many C++ applications.
    check_include_file_cxx(inttypes.h HAVE_INTTYPES_H)
    check_include_file_cxx(limits.h HAVE_LIMITS_H)
    check_include_file_cxx(memory.h HAVE_MEMORY_H)
    check_include_file_cxx(stdint.h HAVE_STDINT_H)
    check_include_file_cxx(stdlib.h HAVE_STDLIB_H)
    check_include_file_cxx(strings.h HAVE_STRINGS_H)
    check_include_file_cxx(string.h HAVE_STRING_H)
    check_include_file_cxx("sys/select.h" HAVE_SYS_SELECT_H)
    check_include_file_cxx("sys/stat.h" HAVE_SYS_STAT_H)
    check_include_file_cxx("sys/time.h" HAVE_SYS_TIME_H)
    check_include_file_cxx("sys/types.h" HAVE_SYS_TYPES_H)
    check_include_file_cxx(unistd.h HAVE_UNISTD_H)

    # STDC_HEADERS conditional logic.
    # This variable is often used in older codebases to indicate C standard header availability.
    if (${HAVE_STDLIB_H})
        set(STDC_HEADERS 1 CACHE INTERNAL "Set if standard C headers are available.")
        add_compile_definitions(STDC_HEADERS)
    else()
        set(STDC_HEADERS 0 CACHE INTERNAL "Set if standard C headers are available.")
    endif()
    message(STATUS "STDC_HEADERS: ${STDC_HEADERS}")

    #-----------------------------------------------------------------------------
    # Check for location of hash_map implementation (for hashmap_defs.h)
    # This section prioritizes modern C++11 <unordered_map> first, then falls back
    # to older or proprietary hash map implementations if needed for compatibility.
    # A single macro will be defined to indicate the found implementation.
    #-----------------------------------------------------------------------------
    message(STATUS "Looking for hash_map implementation...")
    # Initialize a variable to hold the name of the macro that will be defined.
    set(VAANGO_HASHMAP_IMPL_MACRO "" CACHE INTERNAL "Macro to define for the chosen hash map implementation.")

    # 1. Check for C++11 <unordered_map> (standard and preferred).
    set(C11_UNORDERED_MAP_CHECK_SOURCE "
#include <unordered_map>
int main() { std::unordered_map<int, int> xx; return 0; }
    ")
    try_compile(
        HAVE_C11_UNORDERED_MAP
        SOURCE_FROM_VAR "c11_unordered_map.cxx" C11_UNORDERED_MAP_CHECK_SOURCE # Pass the VARIABLE NAME
        OUTPUT_VARIABLE C11_UNORDERED_MAP_CHECK_OUTPUT
    )
    if(HAVE_C11_UNORDERED_MAP)
        message(STATUS "Found C++11 <unordered_map>.")
        set(VAANGO_HASHMAP_IMPL_MACRO "HAVE_C11_HASHMAP")
        add_compile_definitions(HAVE_C11_HASHMAP)
    endif()

    # 2. Check for TR1 <tr1/unordered_map> if C++11 <unordered_map> was not found.
    if(NOT HAVE_C11_UNORDERED_MAP)
        set(TR1_UNORDERED_MAP_CHECK_SOURCE "
#include <tr1/unordered_map>
int main() { std::tr1::unordered_map<int, int> xx; return 0; }
        ")
        try_compile(
            HAVE_TR1_UNORDERED_MAP
            "${VAANGO_TMP_DIR}/TR1UnorderedMapCheck" # Unique binary directory for this test
            CXX
            SOURCE_FROM_VAR
            TR1_UNORDERED_MAP_CHECK_SOURCE # Pass the VARIABLE NAME
            OUTPUT_VARIABLE TR1_UNORDERED_MAP_CHECK_OUTPUT
        )
        if(HAVE_TR1_UNORDERED_MAP)
            message(STATUS "Found <tr1/unordered_map>.")
            set(VAANGO_HASHMAP_IMPL_MACRO "HAVE_TR1_HASHMAP")
            add_compile_definitions(HAVE_TR1_HASHMAP)
        endif()
    endif()

    # 3. Check for GNU extension <ext/hash_map> with __gnu_cxx namespace if previous were not found.
    if(NOT HAVE_C11_UNORDERED_MAP AND NOT HAVE_TR1_UNORDERED_MAP)
        set(GNU_EXT_HASH_MAP_CHECK_SOURCE "
#include <ext/hash_map>
int main() { __gnu_cxx::hash_map<int, int> xx; return 0; }
        ")
        try_compile(
            HAVE_GNU_EXT_HASH_MAP
            "${VAANGO_TMP_DIR}/GNUHashMapCheck" # Unique binary directory for this test
            CXX
            SOURCE_FROM_VAR
            GNU_EXT_HASH_MAP_CHECK_SOURCE # Pass the VARIABLE NAME
            OUTPUT_VARIABLE GNU_EXT_HASH_MAP_CHECK_OUTPUT
        )
        if(HAVE_GNU_EXT_HASH_MAP)
            message(STATUS "Found __gnu_cxx::hash_map in <ext/hash_map>.")
            set(VAANGO_HASHMAP_IMPL_MACRO "HAVE_GNU_HASHMAP")
            add_compile_definitions(HAVE_GNU_HASHMAP)
        endif()
    endif()

    # 4. Check for <ext/hash_map> with std namespace if previous were not found.
    if(NOT HAVE_C11_UNORDERED_MAP AND NOT HAVE_TR1_UNORDERED_MAP AND NOT HAVE_GNU_EXT_HASH_MAP)
        set(EXT_HASH_MAP_CHECK_SOURCE "
#include <ext/hash_map>
int main() { std::hash_map<int, int> xx; return 0; }
        ")
        try_compile(
            HAVE_EXT_HASH_MAP
            "${VAANGO_TMP_DIR}/ExtHashMapCheck" # Unique binary directory for this test
            CXX
            SOURCE_FROM_VAR
            EXT_HASH_MAP_CHECK_SOURCE # Pass the VARIABLE NAME
            OUTPUT_VARIABLE EXT_HASH_MAP_CHECK_OUTPUT
        )
        if(HAVE_EXT_HASH_MAP)
            message(STATUS "Found std::hash_map in <ext/hash_map>.")
            set(VAANGO_HASHMAP_IMPL_MACRO "HAVE_EXT_HASHMAP")
            add_compile_definitions(HAVE_EXT_HASHMAP)
        endif()
    endif()

    # 5. Check for <hash_map> with stdext namespace (common in MSVC) if previous were not found.
    if(NOT HAVE_C11_UNORDERED_MAP AND NOT HAVE_TR1_UNORDERED_MAP AND NOT HAVE_GNU_EXT_HASH_MAP AND NOT HAVE_EXT_HASH_MAP)
        set(STDEXT_HASH_MAP_CHECK_SOURCE "
#include <hash_map>
int main() { stdext::hash_map<int, int> xx; return 0; }
        ")
        try_compile(
            HAVE_STDEXT_HASH_MAP
            "${VAANGO_TMP_DIR}/StdextHashMapCheck" # Unique binary directory for this test
            CXX
            SOURCE_FROM_VAR
            STDEXT_HASH_MAP_CHECK_SOURCE # Pass the VARIABLE NAME
            OUTPUT_VARIABLE STDEXT_HASH_MAP_CHECK_OUTPUT
        )
        if(HAVE_STDEXT_HASH_MAP)
            message(STATUS "Found stdext::hash_map in <hash_map>.")
            set(VAANGO_HASHMAP_IMPL_MACRO "HAVE_STDEXT_HASHMAP")
            add_compile_definitions(HAVE_STDEXT_HASHMAP)
        endif()
    endif()

    # Final message regarding hash map detection.
    if(NOT VAANGO_HASHMAP_IMPL_MACRO)
        message(STATUS "No suitable hash_map implementation found. Will use std::map instead or rely on default fallbacks.")
    else()
        message(STATUS "Selected hash_map implementation: ${VAANGO_HASHMAP_IMPL_MACRO}")
    endif()

    #-----------------------------------------------------------------------
    # Check for additional standard C/C++ headers.
    # These checks ensure the availability of common POSIX/system headers.
    #-----------------------------------------------------------------------
    message(STATUS "Checking more standard C/C++ headers...")
    check_include_file_cxx(dirent.h HAVE_DIRENT_H)
    check_include_file_cxx(sys/param.h HAVE_SYS_PARAM_H)
    check_include_file_cxx(values.h HAVE_VALUES_H) # Note: values.h is often deprecated or not available on modern systems.
    check_include_file_cxx(malloc.h HAVE_MALLOC_H)
    check_include_file_cxx(netdb.h HAVE_NETDB_H)
    check_include_file_cxx(sys/socket.h HAVE_SYS_SOCKET_H)
    check_include_file_cxx(sys/mman.h HAVE_SYS_MMAN_H)
    check_include_file_cxx(sys/ioctl.h HAVE_SYS_IOCTL_H)
    check_include_file_cxx(sys/resource.h HAVE_SYS_RESOURCE_H)
    check_include_file_cxx(sys/wait.h HAVE_SYS_WAIT_H)
    check_include_file_cxx(sys/utsname.h HAVE_SYS_UTSNAME_H)
    check_include_file_cxx(rpc/types.h HAVE_RPC_TYPES_H) # Often part of RPC libraries, may not be generally available.
    check_include_file_cxx(netinet/in.h HAVE_NETINET_IN_H)

    #-----------------------------------------------------------------------
    # Check for C++ standard template library headers.
    # This comprehensive check ensures fundamental STL headers compile together.
    #-----------------------------------------------------------------------
    message(STATUS "Checking C++ STL headers...")
    set(STL_HEADERS_CHECK_SOURCE "
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <string>
int main() {
        std::vector<int> vec; // Simple usage to ensure compilation
        return 0;
      }
    ")
    try_compile(
        HAVE_STL_HEADERS
        SOURCE_FROM_VAR "stl_headers_check.cxx" STL_HEADERS_CHECK_SOURCE # Pass the VARIABLE NAME
        OUTPUT_VARIABLE STL_HEADERS_CHECK_OUTPUT
    )
    if (HAVE_STL_HEADERS)
        message(STATUS "Standard template library headers found and compile.")
        add_compile_definitions(HAVE_STL_HEADERS) # Add as a global definition
    else()
        message(STATUS "Standard template library headers not found or failed to compile.")
    endif()

    # Store a list of all feature variables that should be exposed in a configuration header.
    # These are typically cache variables set by check_include_file_cxx and try_compile.
    set(VAANGO_CONFIG_H_VARS
        HAVE_SEMAPHORE_H
        GCC_HAS_ATOMICS
        HAVE_EXT_ALGORITHM
        HAVE_INTTYPES_H
        HAVE_LIMITS_H
        HAVE_MEMORY_H
        HAVE_STDINT_H
        HAVE_STDLIB_H
        HAVE_STRINGS_H
        HAVE_STRING_H
        HAVE_SYS_SELECT_H
        HAVE_SYS_STAT_H
        HAVE_SYS_TIME_H
        HAVE_SYS_TYPES_H
        HAVE_UNISTD_H
        STDC_HEADERS
        HAVE_C11_UNORDERED_MAP
        HAVE_TR1_UNORDERED_MAP
        HAVE_GNU_EXT_HASH_MAP
        HAVE_EXT_HASH_MAP
        HAVE_STDEXT_HASHMAP
        HAVE_DIRENT_H
        HAVE_SYS_PARAM_H
        HAVE_VALUES_H
        HAVE_MALLOC_H
        HAVE_NETDB_H
        HAVE_SYS_SOCKET_H
        HAVE_SYS_MMAN_H
        HAVE_SYS_IOCTL_H
        HAVE_SYS_RESOURCE_H
        HAVE_SYS_WAIT_H
        HAVE_SYS_UTSNAME_H
        HAVE_RPC_TYPES_H
        HAVE_NETINET_IN_H
        HAVE_STL_HEADERS
        VAANGO_HASHMAP_IMPL_MACRO # The chosen hashmap macro (e.g., HAVE_C11_HASHMAP)
        VAANGO_REFCOUNT_IMPL      # The chosen refcount implementation filename
        VAANGO_ATOMIC_IMPL        # The chosen atomic implementation filename
        VAANGO_TIME_IMPL          # The chosen time implementation filename
        # Add any other variables you want to expose in config_defs.h here.
    )
    # Store this list in the cache for later use (e.g., by configure_file).
    set(VAANGO_CONFIG_H_VARS "${VAANGO_CONFIG_H_VARS}" CACHE INTERNAL "List of variables for config.h generation.")

endfunction()
