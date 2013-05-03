# This is the root EMUNE CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(EMU2DC CXX)

set(MAJOR_VERSION "1")
set(MINOR_VERSION "0")
set(PATCH_VERSION "0")

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_BUILD_TYPE "Release")
set(BUILD_SHARED_LIBS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall")
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

set(BASEPATH ${CMAKE_SOURCE_DIR})
include_directories(${BASEPATH})

#----------------------------------------------------------------------------
# Set module path to search for local Find<XXX>.cmake files
#----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

#----------------------------------------------------------------------------------
# Find the LibXML2 package to use
#----------------------------------------------------------------------------------
find_package(LibXml2 REQUIRED)
if (LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(LIBXML2_DIR "" CACHE PATH "Path to MPI")
  set(XML2_INCLUDE ${LIBXML2_INCLUDE_DIR})
  set(XML2_LIBRARY ${LIBXML2_LIBRARIES})
else()
  message(FATAL_ERROR "No libxml2 found")
  set(XML2_INCLUDE "")
  set(XML2_LIBRARY "")
endif()

#----------------------------------------------------------------------------
# Find Uintah::ProblemSpec 
#----------------------------------------------------------------------------
find_package(UintahProblemSpec REQUIRED)
include_directories(${PROBLEMSPEC_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.49
#----------------------------------------------------------------------------
find_package(Boost 1.49 COMPONENTS REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Find Triangle library
#----------------------------------------------------------------------------
find_package(Triangle REQUIRED)
include_directories(${TRIANGLE_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Just test
#----------------------------------------------------------------------------
set(EMU2DC_SRC
  Domain.cc
)
add_library(EMU2DC_MODULES ${EMU2DC_SRC})

add_executable(emu2dc emu2dc.cc)
target_link_libraries(emu2dc 
  EMU2DC_MODULES 
  ${PROBLEMSPEC_LIBRARY}
  ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Test map
#----------------------------------------------------------------------------
set(testmap_src
  Node.cc
)

add_library(testmap_lib ${testmap_src})
add_executable(test_map test_map.cc)
target_link_libraries(test_map 
  testmap_lib)

#----------------------------------------------------------------------------
# Test point search
#----------------------------------------------------------------------------
set(test_point_search_src
  Domain.cc
  Node.cc
  FamilyComputer.cc
)

add_library(test_point_search_lib ${test_point_search_src})
add_executable(test_point_search test_point_search.cc)
target_link_libraries(test_point_search 
  test_point_search_lib
  ${PROBLEMSPEC_LIBRARY})

#----------------------------------------------------------------------------
# Test input file reader
#----------------------------------------------------------------------------
#foreach(dir ${PROBLEMSPEC_INCLUDE_DIR})
# message(STATUS "Dir = ${dir}")
#endforeach()

set(GEOMETRY_SRC
  Geometry/Point3D.cc
  Geometry/Vector3D.cc
  Geometry/Polygon3D.cc
)

set(test_reader_src
  Body.cc
  Crack.cc
  DamageModel.cc
  Domain.cc
  Element.cc
  FamilyComputer.cc
  HorizonComputer.cc
  Material.cc
  Node.cc
  ProblemSpecReader.cc
  SimulationState.cc
)

add_library(test_reader_lib 
  ${GEOMETRY_SRC}
  ${test_reader_src})
add_executable(test_reader test_reader.cc)
target_link_libraries(test_reader 
  test_reader_lib
  ${PROBLEMSPEC_LIBRARY}
  ${Boost_LIBRARIES}
  ${TRIANGLE_LIBRARY})

#----------------------------------------------------------------------------
# Test surface triangulation
#----------------------------------------------------------------------------
# Add requirements for PCL
#find_package(PCL 1.6 REQUIRED)

#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})
#add_definitions(${PCL_DEFINITIONS})

# Create an executable for the PLY converter
#add_executable(test_surf_tri test_surf_tri.cc)
#target_link_libraries (test_surf_tri 
#  ${PROBLEMSPEC_LIBRARY}
#  ${PCL_LIBRARIES})
