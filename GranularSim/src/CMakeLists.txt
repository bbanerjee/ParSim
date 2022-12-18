#
# The MIT License
#
# Copyright (c) 2017-2022 Parresia Research Limited, New Zealand
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
#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

option(USE_CLANG "Use Clang compiler" OFF)
set(CLANG_ROOT_DIR "/usr/local")

message("USE_CLANG = ${USE_CLANG}")
if (USE_CLANG)
  message("CLANG_ROOT_DIR = ${CLANG_ROOT_DIR}")
  set(CMAKE_C_COMPILER "${CLANG_ROOT_DIR}/bin/clang")
  set(CMAKE_CXX_COMPILER "${CLANG_ROOT_DIR}/bin/clang++")
else()
  option(USE_CLANG "Use Clang compiler" OFF)
endif ()

project(PARAELLIP3D_DEM_PD_SPH CXX)

set(MAJOR_VERSION "0")
set(MINOR_VERSION "1")
set(PATCH_VERSION "0")

#---
# Set module path
#---
set(ParSim_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})

#---
# Build type/compiler etc.
#---
#set(CMAKE_VERBOSE_MAKEFILE 1)
#set(CMAKE_BUILD_TYPE "Release")
#set(BUILD_SHARED_LIBS ON)
#set(CMAKE_BUILD_TYPE "Debug")

if (USE_CLANG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-inconsistent-missing-override -fopenmp")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 -Wall -fno-omit-frame-pointer -fsanitize=address -pedantic -Werror")
  #set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -fopenmp")
endif ()
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

#-------------------------------------------------
# Don't allow in-source builds
#-------------------------------------------------
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message( FATAL_ERROR "In-source builds not allowed. Please create a build directory and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

set(BASEPATH ${CMAKE_SOURCE_DIR})
include_directories(${BASEPATH})

#----------------------------------------------------------------------------
# Set module path to search for local Find<XXX>.cmake files
#----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#----------------------------------------------------------------------------
# Set thirdparty package directory
#----------------------------------------------------------------------------
set(THIRDPARTY_DIR /home/banerjee/Packages)

#----------------------------------------------------------------------------
# Find the MPI package to use
#----------------------------------------------------------------------------
find_package(MPI REQUIRED)
#find_package(MPIParSim REQUIRED)
message(STATUS "MPI_C_FOUND = ${MPI_C_FOUND}")
message(STATUS "MPI_CXX_FOUND = ${MPI_CXX_FOUND}")
message(STATUS "MPI_FOUND = ${MPI_FOUND}")
if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  # Take the MPI include path and split it into -I's
  #SET(MPI_INCLUDE_PATH_CONV "")
  #FOREACH(I ${MPI_INCLUDE_PATH})
  #  SET(MPI_INCLUDE_PATH_CONV "${MPI_INCLUDE_PATH_CONV} -I${I}")
  #ENDFOREACH(I)
  set(MPI_DIR "" CACHE PATH "Path to MPI")
  set(MPI_LIBRARY ${MPI_LIBRARIES})
else()
  message(STATUS "No mpi found")
  set(MPI_DIR "")
  set(MPI_LIBRARY "")
endif()

#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.53
#----------------------------------------------------------------------------
set(BOOST_LIBRARYDIR /usr/lib/x86_64-linux-gnu)
find_package(Boost 1.71 COMPONENTS mpi serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Add requirements for Output in VTK format
#----------------------------------------------------------------------------
#find_package(VTK 5.10 REQUIRED)
#find_package(VTK 6 REQUIRED)
#find_package(VTK 7 REQUIRED)
set(VTK_DIR ${THIRDPARTY_DIR}/VTK-9.1.0/build CACHE PATH "vtk build directory")
find_package(VTK 9.1 REQUIRED 
             COMPONENTS
             CommonCore
             CommonDataModel
             IOXML
             QUIET)
message(STATUS "VTK header directory = ${VTK_INCLUDE_DIRS}")
message(STATUS "VTK libraries directory = ${VTK_LIBRARIES}")
include_directories(${VTK_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Add requirements for zen xml parser
#----------------------------------------------------------------------------
set(ZEN_LINUX 1)

#----------------------------------------------------------------------------
# Add requirements for Zlib compression library 
#----------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  message(STATUS "Zlib compression library found")
  include_directories(${ZLIB_INCLUDE_DIRS})
else()
  message(STATUS "Zlib compression library not found")
  set(ZLIB_DIR "")
  set(ZLIB_LIBRARIES "")
  set(ZLIB_INCLUDE_DIRS "")
endif()

#------------------------------------------------------------------------------
# Find the Eigen3 libraries 
#------------------------------------------------------------------------------
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "Eigen3 include directory is ${EIGEN3_INCLUDE_DIR}")
if (USE_CLANG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN3_INCLUDE_DIR}")
endif()

#------------------------------------------------------------------------------
# Find the Qhull libraries 
# (Make sure the qhull libraries and includes are in qhull/build)
#------------------------------------------------------------------------------
set(QHULL_USE_STATIC ON)
if (USE_CLANG)
  set(QHULL_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/qhull/build_clang)
else (USE_CLANG)
  set(QHULL_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/qhull/build)
endif (USE_CLANG)
find_package(Qhull REQUIRED)


#----------------------------------------------------------------------------
# Set up library sources
#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(ELLIP3D_SRCS
    ${ELLIP3D_SRCS}
)

add_subdirectory(Core)
add_subdirectory(InputOutput)
add_subdirectory(Boundary)
add_subdirectory(DiscreteElements)
add_subdirectory(FluidDynamics)
add_subdirectory(Peridynamics)
add_subdirectory(SmoothParticleHydro)
add_subdirectory(Simulations)

add_library(ellip3d_lib 
  ${ELLIP3D_SRCS})
target_link_libraries(ellip3d_lib
  ${VTK_LIBRARIES}
)

#----------------------------------------------------------------------------
# Create executable
#----------------------------------------------------------------------------
add_executable(granularsim main.cpp)
target_link_libraries(granularsim 
  ellip3d_lib
  ${MPI_LIBRARY}
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ZLIB_LIBRARIES}
)

#----------------------------------------------------------------------------
# Including extra cmake rules for clang devtools
#----------------------------------------------------------------------------
include(cmake/clang-dev-tools.cmake)

#----------------------------------------------------------------------------
# For googletest
#----------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# For populating VTK includes and libraries
vtk_module_autoinit(
  TARGETS granularsim 
          ellip3d_lib
  MODULES ${VTK_LIBRARIES}
)
