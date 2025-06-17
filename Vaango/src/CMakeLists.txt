#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

# This is the root Vaango CMakeLists file.

# Flags for cmake:
#  -DCMAKE_BUILD_TYPE=Debug -- for the debug build
#  -DVISIT_DIR=aa/bb/cc -- to build Visit
#  -DBUILD_UNIT_TESTS=1   -- to build unit tests
#  -DUSE_CLANG=1  -- use clang compilers

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

#----------------------------------------------------------------------------
# Prevent in-source builds 
#----------------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR 
        "Vaango cannot be configured in-source. "
        "Please create a separate build directory!")
endif()

#----------------------------------------------------------------------------
# Set the project with languages
#----------------------------------------------------------------------------
project(Vaango 
    VERSION 1.1.0
    LANGUAGES C CXX Fortran
    DESCRIPTION "Vaango project"
)

# Optional: Create legacy variables for compatibility
set(MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(PATCH_VERSION ${PROJECT_VERSION_PATCH})

# Optional: Enable verbose makefile (uncomment if needed)
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Create interface library for common compile options
add_library(vaango_options INTERFACE)
add_library(vaango_warnings INTERFACE)

# Create aliases for cleaner usage
add_library(Vaango::Options ALIAS vaango_options)
add_library(Vaango::Warnings ALIAS vaango_warnings)

# Backward compatibility
include (CMakeBackwardCompatibilityC)

#----------------------------------------------------------------------------
# Default to release builds
#----------------------------------------------------------------------------
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )

#--------------------------------------------------
# Set the local module path and load macros
#--------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Macros)
message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

#----------------------------------------------------------------------------------
# Configure options
#----------------------------------------------------------------------------------
option(MPM             "Build MPM and MPMICE" ON)
option(ICE             "Build ICE and MPMICE" ON)
option(FORTRAN         "Build with fortran" ON)
option(PETSC           "Build with PETSc" OFF)
option(HYPRE           "Build with HYPRE" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(VISIT             "Build with Visit" OFF) 

# Explicitly setting VISIT ON if VISIT_DIR is defined,
# but consider making VISIT an option if it's user-configurable.
if (DEFINED VISIT_DIR AND NOT VISIT)
    set(VISIT ON CACHE BOOL "Build with Visit" FORCE)
endif ()

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.") # Use bin for executables
mark_as_advanced(CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_RUNTIME_OUTPUT_DIRECTORY)

#----------------------------------------------------------------------------
# C/C++ compiler configuration
#----------------------------------------------------------------------------
# Set C++ standard requirements
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable compile commands export
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler selection (if needed)
if(USE_CLANG)
    find_program(CLANG_C_COMPILER clang REQUIRED)
    find_program(CLANG_CXX_COMPILER clang++ REQUIRED)
    find_program(GFORTRAN_COMPILER gfortran REQUIRED)
    
    set(CMAKE_C_COMPILER ${CLANG_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${CLANG_CXX_COMPILER})
    set(CMAKE_Fortran_COMPILER ${GFORTRAN_COMPILER})
endif()

#----------------------------------------------------------------------------
# Compile options
#----------------------------------------------------------------------------
# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(vaango_options INTERFACE
        $<$<CONFIG:Debug>:-g -O0>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
    )
    target_compile_options(vaango_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(vaango_options INTERFACE
        -rdynamic
        $<$<CONFIG:Debug>:-g -O0>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
    )
    target_compile_options(vaango_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        # Uncomment these for stricter checking:
        # -Wsuggest-override
        # -Wconversion
        # -Wsign-conversion
    )
    target_link_options(vaango_options INTERFACE
        -rdynamic
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_compile_options(vaango_options INTERFACE
        $<$<CONFIG:Debug>:-g -O0>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
    )
    # Intel compiler specific options can be added here
endif()

#----------------------------------------------------------------------------
# Architecture-specific options
#----------------------------------------------------------------------------
if(NOT DEFINED NBITS)
    # Auto-detect architecture if not specified
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(NBITS 64)
    else()
        set(NBITS 32)
    endif()
endif()

if(NBITS EQUAL 64)
    message(STATUS "Building for 64-bit architecture")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(vaango_options INTERFACE
            -m64
            -msse2
        )
    endif()
else()
    message(STATUS "Building for 32-bit architecture")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(vaango_options INTERFACE
            -m32
            -msse2
        )
    endif()
endif()

#----------------------------------------------------------------------------
# Platform-specific configuration
#----------------------------------------------------------------------------
if(WIN32)
    target_compile_definitions(vaango_options INTERFACE
        _USE_MATH_DEFINES          # Defines M_PI for Windows
        _WIN32_WINNT=0x0500       # Windows 2000 features
        _CRT_SECURE_NO_DEPRECATE  # Turn off VS warnings about strcpy etc.
        WIN32_LEAN_AND_MEAN       # Prevent rad1/grp1 definitions
        NOMINMAX                  # Prevent min/max macro definitions
    )
    
    if(NOT BUILD_SHARED_LIBS)
        target_compile_definitions(vaango_options INTERFACE
            BUILD_CORE_STATIC
            GLEW_STATIC
        )
    endif()
    
    # Set warning level for MSVC
    if(MSVC)
        target_compile_options(vaango_warnings INTERFACE
            /W1  # Warning level 1 (instead of default W3/W4)
        )
        # Remove any existing /W flags and replace with /W1
        string(REGEX REPLACE "/W[0-4]" "/W1" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    endif()
endif()

#----------------------------------------------------------------------------
# Platform detection and system information
#----------------------------------------------------------------------------
include(ProcessorCount)
include(CheckTypeSize)

# Detect number of processors
# Set build parallel level if not set
ProcessorCount(NUM_CPUS)
if(NUM_CPUS EQUAL 0)
    set(NUM_CPUS 1)
    set(CMAKE_BUILD_PARALLEL_LEVEL 1)
else()
    set(CMAKE_BUILD_PARALLEL_LEVEL 4)
endif()

# Platform-specific settings with modern CMake approach
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")  # macOS
    set(IS_OSX TRUE)
    target_compile_definitions(vaango_options INTERFACE IS_OSX=1)
    message(STATUS "Platform: macOS")
    
elseif(CMAKE_SYSTEM_NAME STREQUAL "AIX")
    set(IS_AIX TRUE)
    target_compile_definitions(vaango_options INTERFACE IS_AIX=1)
    message(STATUS "Platform: AIX")
    
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(IS_LINUX TRUE)
    target_compile_definitions(vaango_options INTERFACE IS_LINUX=1)
    message(STATUS "Platform: Linux")
    
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    if(MINGW)
        set(IS_WIN TRUE)
        target_compile_definitions(vaango_options INTERFACE IS_WIN=1)
        message(STATUS "Platform: Windows (MinGW)")
    else()
        message(FATAL_ERROR "Currently Windows only configures with MinGW. See install instructions")
    endif()
    
else()
    message(WARNING "Untested system: ${CMAKE_SYSTEM_NAME} - You are on your own!")
endif()

message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Architecture: ${NBITS}-bit")
message(STATUS "CPU Count: ${NUM_CPUS}")
if(CMAKE_BUILD_TYPE)
    message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
endif()

#----------------------------------------------------------------------------
# Fortran settings
#----------------------------------------------------------------------------
include(fortranDependency)

# NO_FORTRAN and DEF_NO_FORTRAN settings
# NO_FORTRAN is used in configVars.mk[.in] and 
# DEF_NO_FORTRAN is used in uintah_testdefs.h.in
if (CMAKE_Fortran_COMPILER)
    # Fortran compiler was found and enabled by the project() command
    set(NO_FORTRAN "no")
    set(DEF_NO_FORTRAN "") # No #define NO_FORTRAN needed if Fortran is enabled
    message(STATUS "Fortran compiler found: ${CMAKE_Fortran_COMPILER}")
else ()
    # Fortran compiler was NOT found, even though requested in project()
    # This scenario is less common if Fortran is a strict requirement,
    # as project() would typically error out if a language isn't found.
    # However, it's good to be robust.
    set(NO_FORTRAN "yes")
    set(DEF_NO_FORTRAN "#define NO_FORTRAN")
    message(FATAL_ERROR "Fortran compiler not found, even though specified in project() command. Fortran-dependent features will be disabled.")
endif ()

message(STATUS "Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")

# Set default values for Fortran underscore style and library.
# These defaults can be overridden based on specific compiler IDs.
set(FORTRAN_UNDERSCORE "end" CACHE STRING "What type of Fortran underscore style - linux,end,none")
set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
set(F_LIBRARY "" CACHE STRING "Fortran runtime library to link")

# Modern approach using CMAKE_Fortran_COMPILER_ID
# CMake's internal modules often handle linking the correct runtime library
# automatically when you enable the Fortran language for targets.
# However, if explicit linking is still required, use find_library.
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # GFortran typically uses an "end" underscore style.
    # CMake's Fortran language support usually links libgfortran automatically.
    # If not, you could explicitly find it:
    # find_library(F_LIBRARY NAMES gfortran)
    set(FORTRAN_UNDERSCORE "end" CACHE STRING "What type of Fortran underscore style - linux,end,none")
    set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
    message(STATUS "Configured for GFortran.")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # Intel Fortran (ifort) also typically uses an "end" underscore style.
    # CMake usually handles ifort runtime libraries automatically.
    # If not, you might need to find ifcore or other Intel libraries.
    # find_library(F_LIBRARY NAMES ifcore)
    set(FORTRAN_UNDERSCORE "end" CACHE STRING "What type of Fortran underscore style - linux,end,none")
    set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
    message(STATUS "Configured for Intel Fortran.")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC") # PGI is now NVHPC SDK
    # PGI/NVHPC Fortran typically uses an "end" underscore style.
    # find_library(F_LIBRARY NAMES pgftnrtl)
    set(FORTRAN_UNDERSCORE "end" CACHE STRING "What type of Fortran underscore style - linux,end,none")
    set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
    message(STATUS "Configured for PGI/NVHPC Fortran.")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL") # IBM XL Fortran
    # XL Fortran typically uses a "none" underscore style.
    # find_library(F_LIBRARY NAMES xlf90)
    set(FORTRAN_UNDERSCORE "none" CACHE STRING "What type of Fortran underscore style - linux,end,none")
    set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_NONE")
    message(STATUS "Configured for IBM XL Fortran.")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    # Cray Fortran (ftn) might handle things differently, often no explicit linking needed.
    set(FORTRAN_UNDERSCORE "end" CACHE STRING "What type of Fortran underscore style - linux,end,none") # Common for Cray
    set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
    message(STATUS "Configured for Cray Fortran.")
else ()
    # Fallback for unknown compilers.
    message(WARNING "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}. "
                    "Assuming default underscore style 'end' and automatic runtime library linking.")
    # If Fortran is essential and you can't proceed without a known compiler,
    # consider uncommenting the FATAL_ERROR below:
    # message(FATAL_ERROR "Unsupported Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

message(STATUS "Fortran underscore style: ${FORTRAN_UNDERSCORE}")
message(STATUS "Fortran underscore define: ${DEF_FORTRAN_UNDERSCORE}")
message(STATUS "Fortran library: ${F_LIBRARY}")

#----------------------------------------------------------------------------
# Type size checks and endianness detection
#----------------------------------------------------------------------------
check_type_size("long long" SIZEOF_LONG_LONG)
message(STATUS "sizeof(long long): ${SIZEOF_LONG_LONG}")

# Long long size validation and swap determination
if(SIZEOF_LONG_LONG EQUAL 8)
    set(LONG_LONG_SWAP "SWAP_8")
elseif(SIZEOF_LONG_LONG EQUAL 4)
    set(LONG_LONG_SWAP "SWAP_4")
else()
    message(FATAL_ERROR 
        "Unexpected sizeof(long long): ${SIZEOF_LONG_LONG}\n"
        "This usually indicates an incorrect LD_LIBRARY_PATH for your compiler.\n"
        "To debug: check the compile command and C++ program from config.log,\n"
        "compile and run it manually to see the actual error.\n"
        "Fix your LD_LIBRARY_PATH or compiler setup.")
endif()

# Endianness detection (modern CMake way)
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    target_compile_definitions(vaango_options INTERFACE IS_BIG_ENDIAN=1)
    message(STATUS "Endianness: Big Endian")
else()
    target_compile_definitions(vaango_options INTERFACE IS_LITTLE_ENDIAN=1)
    message(STATUS "Endianness: Little Endian")
endif()

#----------------------------------------------------------------------------
# Configure generated header files
#----------------------------------------------------------------------------
# Generate Endian.h
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/Util/testEndian.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Core/Util/Endian.h
    @ONLY
)

# Generate mem_init.h
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/Malloc/testmem_init.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Core/Malloc/mem_init.h
    @ONLY
)

# Add binary directory to include path for generated headers
target_include_directories(vaango_options INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Call the function to perform all the feature checks.
# This populates CMake cache variables (e.g., HAVE_SEMAPHORE_H)
# and sets global compile definitions (e.g., HAVE_C11_HASHMAP).
include(CheckCppFeatures)
vaango_check_features()

# Optional: You can access the results of the checks directly here.
# These messages provide immediate feedback during the CMake configure step.
message(STATUS "Selected RefCounted implementation: ${VAANGO_REFCOUNT_IMPL}")
message(STATUS "Selected AtomicCounter implementation: ${VAANGO_ATOMIC_IMPL}")
message(STATUS "Selected Hashmap Macro: ${VAANGO_HASHMAP_IMPL_MACRO}")
message(STATUS "Time implementation: ${VAANGO_TIME_IMPL}")

#----------------------------------------------------------------------------
# Optional variable intialization  
#----------------------------------------------------------------------------
set(HAVE_HDF5 "no")
set(HAVE_TABPROPS "no")
set(HAVE_SPATIALOPS "no")
set(HAVE_EXPRLIB "no")

#----------------------------------------------------------------------------------
# Turn off ICE or MPM (for sci_defs/uintah_defs.h)
#----------------------------------------------------------------------------------
set(DEF_MPM "")
set(DEF_ICE "")
set(DEF_NO_FORTRAN "")

if (MPM)
  set(BUILD_MPM "yes")
  set(USE_MPM   "yes")
else ()
  set(DEF_MPM "#define NO_MPM")
  set(BUILD_MPM "no")
  set(USE_MPM   "no")
endif()

if (ICE)
  set(BUILD_ICE "yes")
  set(USE_ICE   "yes")
else()
  set(DEF_ICE "#define NO_ICE")
  set(BUILD_ICE "no")
  set(USE_ICE   "no")
endif()

#----------------------------------------------------------------------------------
# Find the Perl package to use
#----------------------------------------------------------------------------------
find_program(PERL perl)
set(FSPEC ${CMAKE_CURRENT_BINARY_DIR}/genfspec.pl CACHE STRING "Path to fspec")
configure_file(${Vaango_SOURCE_DIR}/tools/fspec.pl.in ${FSPEC})

#----------------------------------------------------------------------------------
# Find the MPI package to use
#----------------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
if (MPI_CXX_FOUND)
  message(STATUS "MPI found. Version: ${MPI_CXX_VERSION}")
  #message(STATUS "MPI CXX Libraries: ${MPI_CXX_LIBRARIES}")
  #message(STATUS "MPI Fortran Libraries: ${MPI_Fortran_LIBRARIES}")
  #message(STATUS "MPI Include: ${MPI_CXX_INCLUDE_PATH}")

  set(DEF_MPI_CONST_WORKS "#define MPICONST const")
  set(DEF_MPI_MAX_THREADS "#define MPI_MAX_THREADS 64")
  set(DEF_MPI3_ENABLED "#define UINTAH_ENABLE_MPI3 true")
endif()

#----------------------------------------------------------------------------------
# Find the LibXML2 package to use
#----------------------------------------------------------------------------------
find_package(LibXml2 REQUIRED)
if (LIBXML2_FOUND)
  #include_directories(${LIBXML2_INCLUDE_DIR})
  message(STATUS "LibXml2 found. Version: ${LIBXML2_VERSION_STRING}")
endif()

#-----------------------------------------------------------------------------
# Configure time cross-platform header file generation
#-----------------------------------------------------------------------------
# SCI_MALLOC on?
IF (CMAKE_BUILD_TYPE MATCHES DEBUG)
  set(DEF_MALLOC "")
  set(DEF_MALLOC_LINENUMBERS "")
else()
  set(DEF_MALLOC "#define DISABLE_SCI_MALLOC 1")
  set(DEF_MALLOC_LINENUMBERS "")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

# SCI Assertions
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(THROW_LEVEL  0)
  set(ASSERTION_LEVEL  0)
else()
  set(THROW_LEVEL  0)
  set(ASSERTION_LEVEL  3)
endif()

# SCI Template definitions
set(TEMPLATE_TAG "")
set(TEMPLATE_BOX "<>")

#-----------------------------------------------------------------------------
# Threads library (for populating sci_defs/threads_def.h)
#-----------------------------------------------------------------------------
find_package(Threads REQUIRED)

if(TARGET Threads::Threads)
    set(DEF_THREAD "#define SCI_PTHREAD 1")
    
    # Modern way to get thread library properties
    get_target_property(THREAD_LIBRARY Threads::Threads INTERFACE_LINK_LIBRARIES)
    if(NOT THREAD_LIBRARY)
        set(THREAD_LIBRARY "")
    endif()
    set(THREAD_LIB_FLAG "${THREAD_LIBRARY}")
    
    # Determine thread type based on platform
    if(CMAKE_USE_PTHREADS_INIT OR UNIX)
        set(THREADS_TYPE "pthreads")
    elseif(WIN32)
        set(THREADS_TYPE "win32")
    else()
        set(THREADS_TYPE "unix")
    endif()
    
    message(STATUS "Threads library ${THREADS_TYPE} found")
else()
    set(DEF_THREAD "#define SCI_NOTHREAD 1")
    set(THREAD_LIB_FLAG "")
    set(THREAD_LIBRARY "")
    message(WARNING "No threads library found. Software may work incorrectly.")
endif()

set(THREAD_IMPL "Thread_${THREADS_TYPE}.cc")

#-----------------------------------------------------------------------------
# Find dl library (for dlopen())
#-----------------------------------------------------------------------------
find_library(DL_LIBRARY NAMES dl)
if(DL_LIBRARY)
    message(STATUS "DL library found: ${DL_LIBRARY}")
else()
    set(DL_LIBRARY "")
    message(STATUS "DL library not found")
endif()

#------------------------------------------------------------------------------
# Find the Eigen3 libraries 
#------------------------------------------------------------------------------
find_package(Eigen3 REQUIRED)

# Modern target-based approach
if(TARGET Eigen3::Eigen)
    message(STATUS "Eigen3 found as imported target")
    # No need for manual include_directories with modern targets
    
    if(USE_CLANG)
        # Modern way to add system includes to a target
        target_compile_options(your_target_name PRIVATE 
            $<$<COMPILE_LANGUAGE:CXX>:-isystem>
            $<$<COMPILE_LANGUAGE:CXX>:$<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>>
        )
    endif()
else()
    # Fallback for older Eigen3 installations
    include_directories(${EIGEN3_INCLUDE_DIR})
    message(STATUS "Eigen3 include directory is ${EIGEN3_INCLUDE_DIR}")
    
    if(USE_CLANG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN3_INCLUDE_DIR}")
    endif()
endif()

#------------------------------------------------------------------------------
# Find the HDF5 libraries 
#------------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX)

if(HDF5_FOUND)
    if(TARGET hdf5::hdf5_cpp OR TARGET HDF5::HDF5)
        # Modern target-based approach
        message(STATUS "HDF5 found as imported target")
        if(TARGET hdf5::hdf5_cpp)
            set(HDF5_TARGET "hdf5::hdf5_cpp;hdf5::hdf5")
        else()
            set(HDF5_TARGET HDF5::HDF5)
        endif()
        # Target will handle includes and linking automatically
    else()
        # Fallback for older HDF5 installations
        include_directories(${HDF5_INCLUDE_DIRS})
        message(STATUS "HDF5 include directory is ${HDF5_INCLUDE_DIRS}")
        set(HDF5_TARGET ${HDF5_CXX_LIBRARIES})
    endif()
else()
    message(STATUS "No HDF5 found")
    set(HDF5_TARGET "")
endif()

#-----------------------------------------------------------------------------
# Build PETSc
#-----------------------------------------------------------------------------

# Configure PETSc with custom options
if(PETSC)
  include(ConfigurePetsc)
  configure_petsc(
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/petsc
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/petsc-install
      DEBUGGING OFF                    # Enable/disable debugging
      SHARED_LIBS ON                   # Build shared libraries
      WITH_MPI ON                      # Use OpenMPI
      DOWNLOAD_DEPS ON                 # Download BLAS/LAPACK
      EXTRA_CONFIGURE_ARGS             # Additional PETSc configure options
          --with-scalar-type=real
          --with-precision=double
          --download-superlu_dist=1
  )
  set(HAVE_PETSC TRUE)
  set(DEF_PETSC "#define HAVE_PETSC 1")
endif()

#-----------------------------------------------------------------------------
# Build HYPRE (for multigrid implicit linear solvers)
#-----------------------------------------------------------------------------

if(HYPRE)
    include(ConfigureHypre)
    configure_hypre()

    message("-- HYPRE: include_dir: ${HYPRE_INCLUDE_DIR}")
    message("-- HYPRE: library: ${HYPRE_LIBRARY}")
    message("-- HYPRE_FOUND : ${HYPRE_FOUND}")

    if (HYPRE_FOUND) 
      set(HAVE_HYPRE TRUE)
      set(DEF_HYPRE "#define HAVE_HYPRE 1")
      
      message(STATUS "HYPRE_VERSION_MAJOR: ${HYPRE_VERSION_MAJOR}")
      message(STATUS "HYPRE_VERSION_MINOR: ${HYPRE_VERSION_MINOR}")
      message(STATUS "HYPRE_VERSION_PATCH: ${HYPRE_VERSION_PATCH}")

      # Version information
      if(HYPRE_VERSION_MAJOR)
          set(DEF_HYPRE_VERSION_MAJOR "#define HYPRE_VERSION_MAJOR ${HYPRE_VERSION_MAJOR}")
      endif()
      if(HYPRE_VERSION_MINOR)
          set(DEF_HYPRE_VERSION_MINOR "#define HYPRE_VERSION_MINOR ${HYPRE_VERSION_MINOR}")
      endif()
      if(HYPRE_VERSION_PATCH)
          set(DEF_HYPRE_VERSION_PATCH "#define HYPRE_VERSION_PATCH ${HYPRE_VERSION_PATCH}")
      endif()

      message(STATUS "HYPRE found: ${HYPRE_VERSION_MAJOR}.${HYPRE_VERSION_MINOR}; Target: ${HYPRE_TARGET}")
    else()
      message(STATUS "HYPRE not found")
    endif()
    
endif()

#-----------------------------------------------------------------------------
# Find VisIt and build plugin
#-----------------------------------------------------------------------------
if(VISIT)
    if(NOT VISIT_DIR)
        message(FATAL_ERROR "VISIT_DIR must be specified when VISIT=ON")
    endif()
    
    if(NOT EXISTS "${VISIT_DIR}")
        message(FATAL_ERROR "VisIt directory does not exist: ${VISIT_DIR}")
    endif()
    
    message(STATUS "Checking for VisIt in ${VISIT_DIR}")
    
    set(VISIT_SOURCE_DIR "${VISIT_DIR}")
    message(STATUS "VISIT_SOURCE_DIR = ${VISIT_SOURCE_DIR}")
    
    # Use generator expressions for more robust path handling
    set(UINTAH_UDA_TO_VIS_LIB "$<TARGET_FILE_DIR:StandAlone_tools_uda2vis>/libStandAlone_tools_uda2vis.so")
    message(STATUS "Setting UDA_TO_VISIT_LIB_PATH = ${UINTAH_UDA_TO_VIS_LIB}")
    
    # Modern configure_file usage
    configure_file(
        "${CMAKE_SOURCE_DIR}/VisIt/udaReaderMTMD/avtudaReaderMTMDFileFormat.C.in"
        "${CMAKE_BINARY_DIR}/VisIt/udaReaderMTMD/avtudaReaderMTMDFileFormat.C"
        @ONLY
    )
    
    message(STATUS "Configuring avtudaReaderMTMDFileFormat.C")
    add_subdirectory(VisIt)
endif()

# Example of how to link these libraries to your targets (modern approach)
#-----------------------------------------------------------------------------
# When you create your executable or library targets, use:
#
# target_link_libraries(your_target_name PRIVATE
#     Threads::Threads
#     $<$<BOOL:${DL_LIBRARY}>:${DL_LIBRARY}>
#     Eigen3::Eigen  # if available as target
#     ${HDF5_TARGET}  # if available as target, otherwise use ${HDF5_LIBRARY}
#     $<$<BOOL:${PETSC}>:PETSc::PETSc>
#     $<$<BOOL:${HYPRE}>:HYPRE::HYPRE>
# )
#-----------------------------------------------------------------------
# Create the cross-platform header files
#-----------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs)

include(ConfigureHeaders)
configure_standard_sci_headers()

#-----------------------------------------------------------------------------
# For reading JSON input files
#-----------------------------------------------------------------------------
#set(JSON_INCLUDES ${Vaango_SOURCE_DIR}/submodules/json/src)
#include_directories(${JSON_INCLUDES})

#-----------------------------------------------------------------------------
# Packages to enable (must enable packages when building core for now)
set(LOAD_PACKAGE "Vaango" CACHE STRING "Building Vaango package by default.")
set(HANDLED_PACKAGE_Vaango "YES")

include_directories(
  ${PETSC_INCLUDES}
  ${HYPRE_INCLUDE}
  ${XML2_INCLUDE}
)

include_directories(
  ${Vaango_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

#=========================================================================
# For compliation, installation, and regression testing
#=========================================================================
MESSAGE(STATUS "")
MESSAGE(STATUS "Configuring done. Adding subdirectories.")

#------------------------------------------------------------------------
# Add main sources
#------------------------------------------------------------------------
add_subdirectory(Core)
add_subdirectory(CCA)
add_subdirectory(StandAlone)

#----------------------------------------------------------------------------------
# Regression testing
#----------------------------------------------------------------------------------
if (EXISTS ${CMAKE_SOURCE_DIR}/R_Tester/CMakeLists.txt)
  set(HAVE_REGRESSION_TESTS TRUE)
endif()
if (HAVE_REGRESSION_TESTS)
  add_subdirectory(R_Tester)
endif()

#----------------------------------------------------------------------------------
# Unit testing (CTest)
#----------------------------------------------------------------------------------
#enable_testing()
#add_test(NAME create_gold_standards 
#  COMMAND ${CMAKE_COMMAND} 
#   -D Vaango_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
#   -D Vaango_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
#   -P ${Vaango_SOURCE_DIR}/cmake/Macros/goldstandards.cmake)
#set_tests_properties(create_gold_standards PROPERTIES ENVIRONMENT
#  "PYTHONPATH=${Vaango_SOURCE_DIR}/R_Tester/toplevel;${Vaango_SOURCE_DIR}/R_Tester")

#----------------------------------------------------------------------------------
# Unit testing (googletests)
#----------------------------------------------------------------------------------
if (USE_CLANG)
  set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
  set(CMAKE_CXX_FLAGS "-std=c++11")
endif ()

#include(FetchContent)
#FetchContent_Declare(
#  googletest
#  GIT_REPOSITORY https://github.com/google/googletest.git
#  GIT_TAG release-1.12.1
#)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#FetchContent_MakeAvailable(googletest)

add_subdirectory(submodules)

