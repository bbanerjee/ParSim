#
# The MIT License
#
# Copyright (c) 2017- Callaghan Innovation, New Zealand
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
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(PARAELLIP3D_DEM_PD CXX)

set(MAJOR_VERSION "0")
set(MINOR_VERSION "1")
set(PATCH_VERSION "0")

if (USE_CLANG)
  set(CMAKE_CXX_COMPILER "clang++")
endif ()

#set(CMAKE_VERBOSE_MAKEFILE 1)
#set(CMAKE_BUILD_TYPE "Release")
set(BUILD_SHARED_LIBS ON)
set(CMAKE_BUILD_TYPE "Debug")

if (USE_CLANG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wno-deprecated -fopenmp")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -Wall -Wno-deprecated -fopenmp")
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
# Find the MPI package to use
#----------------------------------------------------------------------------
find_package(MPI REQUIRED)
if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  # Take the MPI include path and split it into -I's
  SET(MPI_INCLUDE_PATH_CONV "")
  FOREACH(I ${MPI_INCLUDE_PATH})
    SET(MPI_INCLUDE_PATH_CONV "${MPI_INCLUDE_PATH_CONV} -I${I}")
  ENDFOREACH(I)
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
find_package(Boost 1.53 COMPONENTS mpi serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Add requirements for Output in VTK format
#----------------------------------------------------------------------------
find_package(VTK 5.10 REQUIRED)
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
add_subdirectory(Simulations)

add_library(ELLIP3D_LIB 
  ${ELLIP3D_SRCS})

#----------------------------------------------------------------------------
# Create executable
#----------------------------------------------------------------------------
add_executable(paraEllip3D main.cpp)
target_link_libraries(paraEllip3D 
  ELLIP3D_LIB
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
if (USE_CLANG)
  set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
endif ()
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory(googletest)

