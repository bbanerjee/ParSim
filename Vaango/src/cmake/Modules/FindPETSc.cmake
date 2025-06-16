#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

# - Try to find PETSc
# Once done this will define
#
# PETSc_FOUND - system has PETSc
# PETSc_INCLUDE_DIRS - the PETSc include directories
# PETSc_LIBRARIES - Link these to use PETSc
# PETSc_COMPILER - Compiler used by PETSc, helpful to find a compatible MPI
# PETSc_DEFINITIONS - Compiler switches for using PETSc
# PETSc_MPIEXEC - Executable for running MPI programs
# PETSc_VERSION - Version string (MAJOR.MINOR.SUBMINOR)
#
# An imported target `PETSc::PETSc` will also be created, which encapsulates
# all necessary include directories, libraries, and compile definitions.
#
# Usage:
# find_package(PETSc COMPONENTS CXX) - required if PETSc built with C++ (e.g., --with-clanguage=C++)
# find_package(PETSc COMPONENTS C) - standard behavior, checks build using a C compiler
# find_package(PETSc) - same as above, will default to C or CXX based on enabled languages
#
# Setting these changes the behavior of the search:
# PETSC_DIR - directory in which PETSc resides
# PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_minimum_required(VERSION 3.28...3.29)

# Attempt to find PETSc using its CMake config file first.
# Modern PETSc distributions (>= 3.14) typically provide a PETScConfig.cmake.
# This is the preferred method as it handles all complexities internally.
find_package(PETSc CONFIG QUIET
             PATHS "${PETSC_DIR}" # Look in user-provided PETSC_DIR first
             PATH_SUFFIXES "lib/cmake/petsc" # Standard installation path for config files
)

if (PETSc_FOUND)
    message(STATUS "Found PETSc (via config file) version ${PETSc_VERSION}")
    # The config file typically sets PETSc_INCLUDE_DIRS, PETSc_LIBRARIES, etc.
    # and creates an imported target like PETSc::PETSc.
    # We will map these to our expected variables for consistency.
    if (TARGET PETSc::PETSc)
        # Populate our output variables from the imported target properties
        get_target_property(PETSc_INCLUDE_DIRS PETSc::PETSc INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(PETSc_LIBRARIES PETSc::PETSc INTERFACE_LINK_LIBRARIES)
        get_target_property(PETSc_DEFINITIONS PETSc::PETSc INTERFACE_COMPILE_DEFINITIONS)

        # PETSc_COMPILER and PETSc_MPIEXEC are often not part of the config file.
        # If needed, these might still require probing or assuming system defaults.
        # For this module, we assume FindMPI handles MPIEXEC and the compiler is
        # compatible with the one used to build PETSc.
        set(PETSc_COMPILER "" CACHE FILEPATH "Compiler used by PETSc (Not available via config file)")
        set(PETSc_MPIEXEC "" CACHE FILEPATH "Executable for running MPI programs (Not available via config file)")

        # Finalize and report success using standard CMake args.
        include(FindPackageHandleStandardArgs)
        find_package_handle_standard_args(PETSc
            REQUIRED_VARS PETSc_INCLUDE_DIRS PETSc_LIBRARIES
            VERSION_VAR PETSc_VERSION
            FAIL_MESSAGE "PETSc could not be found via config file. Try setting PETSC_DIR or PETSC_ARCH."
        )
        return() # Exit the module if found via config file
    endif()
endif()

# Fallback to module mode if config file was not found or the target was not created.
message(STATUS "PETSc config file not found or target not created. Falling back to module mode search.")

# --- Component Handling (C or CXX language bindings) ---
set(PETSC_VALID_COMPONENTS C CXX)
if (NOT PETSc_FIND_COMPONENTS)
    # If no component is specified, default to C if enabled, else CXX.
    get_property(_enabled_langs GLOBAL PROPERTY ENABLED_LANGUAGES)
    if ("C" IN_LIST _enabled_langs)
        set(PETSC_LANGUAGE_BINDINGS "C")
    else ()
        set(PETSC_LANGUAGE_BINDINGS "CXX")
    endif ()
else()
    # Ensure only one valid component is specified.
    list(LENGTH PETSc_FIND_COMPONENTS components_length)
    if(${components_length} GREATER 1)
        message(FATAL_ERROR "Only one component for PETSc is allowed to be specified (C or CXX).")
    endif()
    if (NOT PETSc_FIND_COMPONENTS IN_LIST PETSC_VALID_COMPONENTS)
        message(FATAL_ERROR "Invalid PETSc component(s) specified: ${PETSc_FIND_COMPONENTS}. Only C or CXX is allowed.")
    endif()
    set(PETSC_LANGUAGE_BINDINGS ${PETSc_FIND_COMPONENTS})
endif()
message(STATUS "Finding PETSc for language binding: ${PETSC_LANGUAGE_BINDINGS}")

# --- Step 1: Find PETSC_DIR ---
find_path(PETSC_DIR include/petsc.h
    HINTS ENV PETSC_DIR # Check environment variable first
    PATHS
    /usr/lib/petsc       # Common Linux path
    /opt/petsc           # Common Linux/macOS path
    /opt/local/lib/petsc # MacPorts path
    $ENV{HOME}/petsc     # User home directory install
    DOC "Root directory of PETSc installation"
)

# --- Step 2: Find PETSC_ARCH if not explicitly set and PETSC_DIR is found ---
if (PETSC_DIR AND NOT PETSC_ARCH)
    set(_petsc_arches
        $ENV{PETSC_ARCH}
        linux-gnu-c-debug
        linux-gnu-c-opt
        x86_64-unknown-linux-gnu
        i386-unknown-linux-gnu
        # Add more common architectures if necessary (e.g., custom builds)
    )
    foreach (arch ${_petsc_arches})
        find_path(petscconf_h_probe petscconf.h
            HINTS "${PETSC_DIR}"
            PATH_SUFFIXES "${arch}/include" "bmake/${arch}" # Common locations for petscconf.h
            NO_DEFAULT_PATH
        )
        if (petscconf_h_probe)
            set(PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture" FORCE)
            break() # Found a suitable architecture, stop searching
        endif()
    endforeach()
endif()

message(STATUS "Determined PETSC_DIR = ${PETSC_DIR} ; PETSC_ARCH = ${PETSC_ARCH}")

# Only proceed if PETSC_DIR and PETSC_ARCH are identified
if (PETSC_DIR AND PETSC_ARCH)

    # --- Step 3: Define the petsc_get_version function ---
    function (petsc_get_version)
        if (EXISTS "${PETSC_DIR}/include/petscversion.h")
            file(STRINGS "${PETSC_DIR}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
            foreach (line ${vstrings})
                string(REGEX REPLACE " +" ";" fields ${line}) # Split line by spaces
                list(GET fields 1 var) # Get variable name (e.g., PETSC_VERSION_MAJOR)
                list(GET fields 2 val) # Get value
                set(${var} ${val} PARENT_SCOPE) # Set in parent scope for module to use
                set(${var} ${val}) # Also set in local scope for internal logic
            endforeach()

            # Construct the PETSc_VERSION string
            if (PETSC_VERSION_RELEASE) # If it's a release version
                if (${PETSC_VERSION_PATCH} GREATER 0)
                    set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" CACHE INTERNAL "PETSc version")
                else ()
                    set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" CACHE INTERNAL "PETSc version")
                endif ()
            else () # Development version
                # Make dev version compare higher than any patch level of a released version
                set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" CACHE INTERNAL "PETSc version")
            endif ()
        else ()
            message(SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist.")
        endif ()
    endfunction ()

    petsc_get_version() # Call to get the PETSc version information

    # --- Step 4: Determine PETSc configuration file paths ---
    # These files contain critical build information (compiler, flags, libs).
    set(petsc_conf_rules "")
    set(petsc_conf_variables "")
    set(petscconf_h_path "")

    find_path(petscconf_h_path petscconf.h
        HINTS "${PETSC_DIR}"
        PATH_SUFFIXES "${PETSC_ARCH}/include" "bmake/${PETSC_ARCH}" # Common locations for petscconf.h
        NO_DEFAULT_PATH
    )

    if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables") # PETSc >= 3.5 layout
        set(petsc_conf_rules "${PETSC_DIR}/lib/petsc/conf/rules")
        set(petsc_conf_variables "${PETSC_DIR}/lib/petsc/conf/variables")
    elseif (EXISTS "${PETSC_DIR}/conf/rules" AND EXISTS "${PETSC_DIR}/conf/variables" AND petscconf_h_path) # PETSc > 2.3.3 layout
        set(petsc_conf_rules "${PETSC_DIR}/conf/rules")
        set(petsc_conf_variables "${PETSC_DIR}/conf/variables")
    elseif (EXISTS "${PETSC_DIR}/bmake/common/rules" AND EXISTS "${PETSC_DIR}/bmake/common/variables" AND petscconf_h_path) # PETSc <= 2.3.3 layout
        set(petsc_conf_rules "${PETSC_DIR}/bmake/common/rules")
        set(petsc_conf_variables "${PETSC_DIR}/bmake/common/variables")
    else ()
        message(FATAL_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} does not specify a valid PETSc installation (configuration files not found).")
    endif ()

    message(STATUS "PETSc configuration rules: ${petsc_conf_rules}")
    message(STATUS "PETSc configuration variables: ${petsc_conf_variables}")

    # --- Step 5: Parse petscconf.h and petscvariables to extract build info ---
    file(READ "${petscconf_h_path}" petscconf_content)
    file(READ "${petsc_conf_variables}" petsc_variables_content)

    # Extract PETSC_INCLUDE from petscconf.h (old style)
    string(REGEX MATCH "define[ \t]+PETSC_INCLUDE[ \t]+\"([^\"]+)\"" _match "${petscconf_content}")
    set(petsc_include_from_header "${CMAKE_MATCH_1}")

    # Extract PETSC_CCPPFLAGS from petscconf.h
    string(REGEX MATCH "define[ \t]+PETSC_CCPPFLAGS[ \t]+\"([^\"]+)\"" _match "${petscconf_content}")
    set(petsc_cpp_line "${CMAKE_MATCH_1}")

    # Extract variables from petsc_conf_variables file (key=value format)
    string(REGEX MATCH "PETSC_LIB_DIR[ \t]*=[ \t]*(.*)" _match "${petsc_variables_content}")
    set(petsc_lib_dir "${CMAKE_MATCH_1}")

    string(REGEX MATCH "PETSC_EXTERNAL_LIB_BASIC[ \t]*=[ \t]*(.*)" _match "${petsc_variables_content}")
    set(petsc_libs_external "${CMAKE_MATCH_1}")

    string(REGEX MATCH "MPIEXEC[ \t]*=[ \t]*(.*)" _match "${petsc_variables_content}")
    set(petsc_mpiexec "${CMAKE_MATCH_1}")

    string(REGEX MATCH "PCC[ \t]*=[ \t]*(.*)" _match "${petsc_variables_content}")
    set(petsc_cc "${CMAKE_MATCH_1}")

    string(REGEX MATCH "PCC_FLAGS[ \t]*=[ \t]*(.*)" _match "${petsc_variables_content}")
    set(petsc_cc_flags "${CMAKE_MATCH_1}")

    message(STATUS "Extracted from PETSc configuration files:")
    message(STATUS "  PETSC_LIB_DIR: '${petsc_lib_dir}'")
    message(STATUS "  PETSC_EXTERNAL_LIB_BASIC: '${petsc_libs_external}'")
    message(STATUS "  MPIEXEC: '${petsc_mpiexec}'")
    message(STATUS "  PCC (Compiler): '${petsc_cc}'")
    message(STATUS "  PCC_FLAGS: '${petsc_cc_flags}'")

    # Process include paths from PETSC_CCPPFLAGS
    set(petsc_includes_all "")
    # Find all include paths specified with -I in PETSC_CCPPFLAGS
    string(REGEX MATCHALL "(-I[a-zA-Z0-9./_-]+)" _all_include_flags "${petsc_cpp_line}")
    foreach(inc_flag ${_all_include_flags})
        string(REPLACE "-I" "" _path "${inc_flag}")
        list(APPEND petsc_includes_all "${_path}")
    endforeach()

    # Add PETSc's standard include directories
    list(APPEND petsc_includes_all "${PETSC_DIR}/include")
    list(APPEND petsc_includes_all "${PETSC_DIR}/${PETSC_ARCH}/include")

    # Remove duplicate include paths
    list(REMOVE_DUPLICATES petsc_includes_all)
    message(STATUS "Discovered PETSc include paths: ${petsc_includes_all}")

    # --- Step 6: Find PETSc Libraries ---
    macro(PETSC_FIND_LIBRARY_HELPER suffix name)
        # Clear any stale value for the library variable
        set(PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
        # find_library automatically handles platform-specific prefixes/suffixes (e.g., lib/.).
        find_library(PETSC_LIBRARY_${suffix}
            NAMES ${name} "${name}_real" "${name}_complex" # Common PETSc library suffixes (for real/complex builds)
            HINTS "${petsc_lib_dir}" "${PETSC_DIR}/${PETSC_ARCH}/lib" # Search in specified PETSc lib directories
            NO_DEFAULT_PATH # Do not search system default paths (unless explicitly needed)
        )
        set(PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}")
        mark_as_advanced(PETSC_LIBRARY_${suffix})
    endmacro()

    set(PETSC_LIBRARIES_ALL "") # Initialize the list of all PETSc libraries

    # Determine if PETSc uses separate libraries for each package or a single monolithic library
    PETSC_FIND_LIBRARY_HELPER(VEC petscvec) # Check for petscvec as an indicator

    if (PETSC_LIBRARY_VEC)
        message(STATUS "PETSc installation uses separate libraries for each package (e.g., petscvec, petscmat).")
        PETSC_FIND_LIBRARY_HELPER(SYS petscsys petsc) # petscsys or older petsc
        PETSC_FIND_LIBRARY_HELPER(MAT petscmat)
        PETSC_FIND_LIBRARY_HELPER(DM petscdm)
        PETSC_FIND_LIBRARY_HELPER(KSP petscksp)
        PETSC_FIND_LIBRARY_HELPER(SNES petscsnes)
        PETSC_FIND_LIBRARY_HELPER(TS petscts)

        # Build the complete list of core PETSc libraries in dependency order
        set(PETSC_LIBRARIES_ALL
            ${PETSC_LIBRARIES_TS}
            ${PETSC_LIBRARIES_SNES}
            ${PETSC_LIBRARIES_KSP}
            ${PETSC_LIBRARIES_DM}
            ${PETSC_LIBRARIES_MAT}
            ${PETSC_LIBRARIES_VEC}
            ${PETSC_LIBRARIES_SYS}
        )
    else ()
        message(STATUS "PETSc installation uses a single monolithic library (e.g., petsc, petsc_real).")
        PETSC_FIND_LIBRARY_HELPER(SINGLE petsc)
        if (NOT PETSC_LIBRARY_SINGLE)
            PETSC_FIND_LIBRARY_HELPER(SINGLE petsc_real)
        endif()
        if (NOT PETSC_LIBRARY_SINGLE)
            PETSC_FIND_LIBRARY_HELPER(SINGLE petsc_complex)
        endif()
        if (NOT PETSC_LIBRARY_SINGLE)
            message(FATAL_ERROR "Could not find single PETSc library (petsc, petsc_real, or petsc_complex).")
        endif()
        set(PETSC_LIBRARIES_ALL ${PETSC_LIBRARY_SINGLE})
    endif ()

    # --- Step 7: Resolve external libraries (e.g., BLAS, LAPACK, etc.) ---
    set(petsc_resolved_external_libs "")
    # Split the PETSC_EXTERNAL_LIB_BASIC string by '-l' to get individual library names
    string(REPLACE "-l" ";" _external_libs_raw "${petsc_libs_external}")
    foreach(lib_name_raw ${_external_libs_raw})
        string(STRIP "${lib_name_raw}" lib_name) # Remove leading/trailing whitespace
        if (lib_name)
            # Try to find the full path for common external libraries
            find_library(FOUND_EXTERNAL_LIB NAMES ${lib_name})
            if (FOUND_EXTERNAL_LIB)
                list(APPEND petsc_resolved_external_libs ${FOUND_EXTERNAL_LIB})
            else()
                # If not found as a full path, treat it as a linker flag (-l<name>)
                list(APPEND petsc_resolved_external_libs "-l${lib_name}")
                message(WARNING "Could not find full path for external PETSc dependency: '${lib_name}'. Adding as linker flag. You might need to add system library paths to CMAKE_PREFIX_PATH.")
            endif()
        endif()
    endforeach()
    list(APPEND PETSC_LIBRARIES_ALL ${petsc_resolved_external_libs})
    list(REMOVE_DUPLICATES PETSC_LIBRARIES_ALL) # Remove any duplicate libraries
    message(STATUS "Complete PETSc libraries for linking: ${PETSC_LIBRARIES_ALL}")

    # --- Step 8: Test if a minimal PETSc program can compile and run ---
    # This is crucial for validating the found installation.
    # PETSc heavily relies on MPI, so we must find MPI first.
    find_package(MPI REQUIRED) # This will set MPI_INCLUDE_DIRS, MPI_LIBRARIES, etc.

    # Define a custom function to check if a source code snippet compiles and runs.
    function(petsc_check_source_runs _includes _libraries _source_code _result_var _language)
        # Set required CMake variables for the check
        set(CMAKE_REQUIRED_INCLUDES "${_includes}" PARENT_SCOPE)
        set(CMAKE_REQUIRED_LIBRARIES "${_libraries}" PARENT_SCOPE)

        if ("${_language}" STREQUAL "C")
            include(CheckCSourceRuns)
            check_c_source_runs("${_source_code}" ${_result_var} RUN_OUTPUT_VARIABLE _run_output)
        elseif ("${_language}" STREQUAL "CXX")
            include(CheckCXXSourceRuns)
            check_cxx_source_runs("${_source_code}" ${_result_var} RUN_OUTPUT_VARIABLE _run_output)
        else()
            message(FATAL_ERROR "Unsupported language for PETSc check: ${_language}. Must be C or CXX.")
        endif()

        # Clean up the required CMake variables
        set(CMAKE_REQUIRED_INCLUDES "" PARENT_SCOPE)
        set(CMAKE_REQUIRED_LIBRARIES "" PARENT_SCOPE)
    endfunction()

    set(_PETSC_ERR_FUNC "CHKERRQ(ierr)") # Modern PETSc error handling macro
    set(_PETSC_TEST_SOURCE "
static const char help[] = \"PETSc test program.\";
#include <petscts.h>   // Required for TS (Time Stepping) objects
#include <petscsys.h>  // Required for PetscInitialize, PetscFinalize

int main(int argc,char *argv[]) {
    PetscErrorCode ierr;
    TS ts;

    // Initialize PETSc. MPI_COMM_WORLD is typically used here.
    ierr = PetscInitialize(&argc,&argv,0,help);${_PETSC_ERR_FUNC};
    // Create a simple TS object
    ierr = TSCreate(PETSC_COMM_WORLD,&ts);${_PETSC_ERR_FUNC};
    // Set options for TS (e.g., from command line)
    ierr = TSSetFromOptions(ts);${_PETSC_ERR_FUNC};
    // Destroy the TS object
    ierr = TSDestroy(&ts);${_PETSC_ERR_FUNC};
    // Finalize PETSc
    ierr = PetscFinalize();${_PETSC_ERR_FUNC};
    return 0;
}
")

    # Combine all necessary include directories for the test program (PETSc + MPI)
    set(petsc_includes_for_test "${petsc_includes_all}" "${MPI_INCLUDE_DIRS}")
    list(REMOVE_DUPLICATES petsc_includes_for_test)

    # Combine all necessary libraries for the test program (PETSc + MPI)
    set(petsc_libraries_for_test "${PETSC_LIBRARIES_ALL}" "${MPI_LIBRARIES}")
    list(REMOVE_DUPLICATES petsc_libraries_for_test)

    petsc_check_source_runs(
        "${petsc_includes_for_test}"
        "${petsc_libraries_for_test}"
        "${_PETSC_TEST_SOURCE}"
        petsc_executable_runs_check # Output variable for the check result
        "${PETSC_LANGUAGE_BINDINGS}" # Language (C or CXX) for the test
    )

    # Set the cache variable indicating if PETSc executable can run
    set(PETSC_EXECUTABLE_RUNS "${petsc_executable_runs_check}" CACHE BOOL
        "Can the system successfully run a PETSc executable? This variable can be manually set to \"YES\" to force CMake to accept a given PETSc configuration, but this will almost always result in a broken build. If you change PETSC_DIR, PETSC_ARCH, or PETSC_CURRENT you would have to reset this variable." FORCE)

    # --- Step 9: Set output variables for this module and create imported target ---
    set(PETSc_INCLUDE_DIRS ${petsc_includes_all} CACHE PATH "PETSc include directories")
    set(PETSc_LIBRARIES ${PETSC_LIBRARIES_ALL} CACHE STRING "PETSc libraries to link against")
    set(PETSc_COMPILER ${petsc_cc} CACHE FILEPATH "Compiler used by PETSc")
    set(PETSc_MPIEXEC ${petsc_mpiexec} CACHE FILEPATH "Executable for running MPI programs (from PETSc config)")
    # Set standard PETSc definitions; assuming __INSDIR__ for PETSc >= 3.1
    set(PETSc_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc compiler definitions")

    # Create an imported INTERFACE library target for PETSc.
    # This target will propagate include directories, link libraries, and compile definitions
    # to any target that links against it.
    if (NOT TARGET PETSc::PETSc)
        add_library(PETSc::PETSc INTERFACE IMPORTED)
    endif()
    target_include_directories(PETSc::PETSc INTERFACE ${PETSc_INCLUDE_DIRS})
    target_link_libraries(PETSc::PETSc INTERFACE ${PETSc_LIBRARIES})
    target_compile_definitions(PETSc::PETSc INTERFACE ${PETSc_DEFINITIONS})
    set(PETSc_TARGET PETSc::PETSc) # Set a variable for convenience in root CMakeLists.txt

endif() # End of PETSC_DIR and PETSC_ARCH check

# --- Step 10: Handle standard arguments for find_package ---
# This ensures REQUIRED, QUIET, COMPONENTS, and VERSION handling are consistent.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSC_EXECUTABLE_RUNS
    VERSION_VAR PETSc_VERSION
    FAIL_MESSAGE "PETSc could not be found. Please ensure PETSC_DIR and PETSC_ARCH are set correctly, or that PETScConfig.cmake is discoverable."
)

# Mark internal variables as advanced so they don't clutter the GUI by default
mark_as_advanced(PETSC_DIR PETSC_ARCH PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_COMPILER PETSc_DEFINITIONS PETSc_MPIEXEC PETSC_EXECUTABLE_RUNS)

