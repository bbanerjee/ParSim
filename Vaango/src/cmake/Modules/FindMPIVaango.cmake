#
# The MIT License
#
# Copyright (c) 2014-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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
# - FindMPIVaango
# This module leverages CMake's built-in FindMPI module for robust MPI discovery.
# It provides project-specific variables for MPI configuration within Vaango.
#
# The following variables are set:
#
# TXMPI_FOUND: TRUE if MPI is found and configured.
# HAVE_MPI: TRUE if MPI is found and configured (same as TXMPI_FOUND).
#
# If MPI is found (i.e., TXMPI_FOUND is TRUE), the following standard MPI variables
# provided by FindMPI will be available and used by Vaango:
#
# MPI_C_FOUND: TRUE if C MPI is found
# MPI_CXX_FOUND: TRUE if C++ MPI is found
# MPI_Fortran_FOUND: TRUE if Fortran MPI is found
#
# MPI_INCLUDE_DIRS: the directories containing the C/C++ MPI header files.
# MPI_MODULE_DIRS: the directories containing the Fortran module files
#                  and the Fortran include files.
# MPI_LIBRARY_DIRS: the directories containing the MPI libraries.
# MPI_LIBRARIES: the list of MPI libraries to link against.
# MPI_LINK_FLAGS: MPI-specific linker flags.
# MPI_C_COMPILER: The C MPI compiler wrapper (e.g., mpicc).
# MPI_CXX_COMPILER: The C++ MPI compiler wrapper (e.g., mpicxx).
# MPI_Fortran_COMPILER: The Fortran MPI compiler wrapper (e.g., mpifort).
# MPI_EXECUTABLE: Path to mpiexec or mpirun.

# Make sure we don't try to find MPI multiple times
if (TXMPI_FOUND)
    return()
endif()

# Standard find_package(MPI) handles most cases including compiler wrappers.
# It will check for C, C++, and Fortran MPI.
# FindMPI will automatically try compiling test programs with compiler wrappers
# (like mpicxx, mpifort) before falling back to system-wide searches.

# The REQUIRED argument is passed directly from the parent find_package call.
# This variable is set by the standard FindPackageHandleStandardArgs module
# which FindMPI uses.
# IF MPIVaango_FIND_REQUIRED is true then we set the REQUIRED argument.
# This handles the case where a higher-level CMakeLists.txt calls
# find_package(MPIVaango REQUIRED)
set(MPI_FIND_QUIETLY ${MPIVaango_FIND_QUIETLY})
set(MPI_FIND_REQUIRED ${MPIVaango_FIND_REQUIRED})
set(MPI_FIND_COMPONENTS ${MPIVaango_FIND_COMPONENTS})
set(MPI_FIND_VERSION ${MPIVaango_FIND_VERSION})

# Use CMake's standard FindMPI module.
# FindMPI will set MPI_FOUND, MPI_C_FOUND, MPI_CXX_FOUND, MPI_Fortran_FOUND,
# and all the other MPI_... variables.
message(STATUS "Searching for MPI using CMake's FindMPI module...")
find_package(MPI)

# After find_package(MPI) completes, MPI_FOUND will be set.
# We map MPI_FOUND to TXMPI_FOUND and HAVE_MPI for consistency with existing Vaango code.
if (MPI_FOUND)
    set(TXMPI_FOUND TRUE CACHE INTERNAL "MPI found by FindMPIVaango.")
    set(HAVE_MPI 1 CACHE BOOL "Whether have MPI.") # For consistency with original code
    message(STATUS "MPI found and configured.")

    # In modern CMake, you should typically use the imported targets provided by FindMPI.
    # e.g., target_link_libraries(my_app PRIVATE MPI::MPI_CXX)
    #
    # However, if existing Vaango code explicitly uses MPI_INCLUDE_DIRS, MPI_LIBRARIES, etc.,
    # then these variables are already set by FindMPI.
    # No need to re-construct MPI_INCLUDE_DIRS, MPI_LIBRARY_DIRS, or MPI_EXECUTABLES.
    # FindMPI already provides these.

    # Example: Print some key variables set by FindMPI for verification
    if (MPI_C_FOUND)
        message(STATUS "MPI C compiler: ${MPI_C_COMPILER}")
        message(STATUS "MPI C includes: ${MPI_C_INCLUDE_PATH}")
        message(STATUS "MPI C libraries: ${MPI_C_LIBRARIES}")
    endif()
    if (MPI_CXX_FOUND)
        message(STATUS "MPI CXX compiler: ${MPI_CXX_COMPILER}")
        message(STATUS "MPI CXX includes: ${MPI_CXX_INCLUDE_PATH}")
        message(STATUS "MPI CXX libraries: ${MPI_CXX_LIBRARIES}")
    endif()
    if (MPI_Fortran_FOUND)
        message(STATUS "MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
        message(STATUS "MPI Fortran includes: ${MPI_Fortran_INCLUDE_PATH}")
        message(STATUS "MPI Fortran module dirs: ${MPI_Fortran_MODULE_DIR}")
        message(STATUS "MPI Fortran libraries: ${MPI_Fortran_LIBRARIES}")
    endif()

    message(STATUS "MPI include directories: ${MPI_INCLUDE_DIRS}")
    message(STATUS "MPI library directories: ${MPI_LIBRARY_DIRS}")
    message(STATUS "MPI libraries to link: ${MPI_LIBRARIES}")
    message(STATUS "MPI executables: ${MPI_EXECUTABLE}")
    message(STATUS "MPI link flags: ${MPI_LINK_FLAGS}")

else ()
    # MPI was not found by find_package(MPI)
    set(TXMPI_FOUND FALSE CACHE INTERNAL "MPI not found by FindMPIVaango.")
    set(HAVE_MPI 0 CACHE BOOL "Whether have MPI.") # For consistency with original code

    if (MPIVaango_FIND_REQUIRED)
        message(FATAL_ERROR "MPI required but not found by FindMPIVaango.")
    else ()
        message(STATUS "MPI not found. Vaango will be built without MPI support.")
    endif ()
endif ()

# Handle the `TX_HAVE_MPICXX_COMPILER_WRAPPER` and `TX_HAVE_MPIFC_COMPILER_WRAPPER`
# if they are truly needed *outside* the scope of standard MPI variables.
# Modern FindMPI (since CMake 3.0) does internal checks with try_compile.
# If these flags are merely indicators that a compiler wrapper was used,
# MPI_C_COMPILER, MPI_CXX_COMPILER, and MPI_Fortran_COMPILER
# would already reflect the wrapper paths (e.g., /usr/bin/mpicxx)
# if CMake found MPI via those wrappers.

# However, if for some reason Vaango's build system specifically needs to know
# if `CC` or `ftn` were the *exact* compilers that supplied MPI,
# we can add a lightweight check *after* FindMPI has run:

if (TXMPI_FOUND)
    # Check if the primary CXX compiler used for MPI is a Cray-like wrapper
    if (MPI_CXX_COMPILER AND ${MPI_CXX_COMPILER} MATCHES "CC$")
        set(TX_HAVE_MPICXX_COMPILER_WRAPPER TRUE CACHE INTERNAL "MPI C++ automatically included/linked by CC compiler wrapper.")
        message(STATUS "Detected MPI C++ compiler is 'CC' wrapper.")
    endif()

    # Check if the primary Fortran compiler used for MPI is a Cray-like wrapper
    if (MPI_Fortran_COMPILER AND ${MPI_Fortran_COMPILER} MATCHES "ftn$")
        set(TX_HAVE_MPIFC_COMPILER_WRAPPER TRUE CACHE INTERNAL "MPI Fortran automatically included/linked by ftn compiler wrapper.")
        message(STATUS "Detected MPI Fortran compiler is 'ftn' wrapper.")
    endif()
endif()

# Cleanup variables used internally by this module if not needed externally
# unset(MPI_FIND_QUIETLY)
# unset(MPI_FIND_REQUIRED)
# unset(MPI_FIND_COMPONENTS)
# unset(MPI_FIND_VERSION)

include(FindPackageHandleStandardArgs)
# Use FPHSA to set TXMPI_FOUND (final result) and issue messages.
# This assumes that MPIVaango_FIND_REQUIRED is passed as a REQUIRED argument
# to the parent find_package call, or derived from it as done above.
find_package_handle_standard_args(MPIVaango
    REQUIRED_VARS MPI_C_FOUND MPI_CXX_FOUND MPI_Fortran_FOUND # Or just MPI_FOUND
    FAIL_MESSAGE "Could not find a complete MPI installation."
    VERSION_VAR MPI_VERSION
)