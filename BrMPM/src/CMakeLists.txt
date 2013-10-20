# This is the root BrMPM CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(BrMPM CXX)

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
# Add requirements for Output in VTK format
#----------------------------------------------------------------------------
find_package(VTK 5.8 REQUIRED)
# message(STATUS "VTK header directory = ${VTK_INCLUDE_DIRS}")
# message(STATUS "VTK libraries directory = ${VTK_LIBRARIES}")
include_directories(${VTK_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Test variant
#----------------------------------------------------------------------------
add_executable(test_variant test_variant.cc)
target_link_libraries(test_variant 
  ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Set up library sources
#----------------------------------------------------------------------------
set(BrMPM_SRC
  Matrix3D.cc
  MPMUtils.cc
  MPMDatawarehouse.cc
  MPMContact.cc
  MPMVelocityContact.cc
#  MPM2D.cc
  MPMFreeContact.cc
  MPMFrictionContact.cc
  MPMFrictionlessContact.cc
  MPMMaterial.cc
  MPMMatrix.cc
  MPMPatch.cc
  MPMShapeFunction.cc
  MPMTime.cc
  Body.cc
  Domain.cc
  ForceBC.cc
  InitialConditions.cc
  Material.cc
  Node.cc
  Output.cc
  OutputVTK.cc
  ProblemSpecReader.cc
  ProblemSpecUtil.cc
  SimulationState.cc
  Time.cc
  MPMTime.cc
  VelocityBC.cc
)
set(GEOMETRY_SRC
  Geometry/Box3D.cc
  Geometry/Point3D.cc
  Geometry/Polygon3D.cc
  Geometry/Vector3D.cc
)
set(GEOMETRY_PIECE_SRC
  GeometryPiece/BoxGeometryPiece.cc
  GeometryPiece/GeometryPiece.cc
  GeometryPiece/GeometryPieceFactory.cc
  GeometryPiece/GeometryReader.cc
)

add_library(BrMPM_LIB 
  ${BrMPM_SRC}
  ${GEOMETRY_SRC}
  ${GEOMETRY_PIECE_SRC})

#----------------------------------------------------------------------------
# Test input file reader
#----------------------------------------------------------------------------
#foreach(dir ${PROBLEMSPEC_INCLUDE_DIR})
# message(STATUS "Dir = ${dir}")
#endforeach()

#add_executable(test_reader test_reader.cc)
#target_link_libraries(test_reader 
#  BrMPM_LIB
#  ${PROBLEMSPEC_LIBRARY}
#  ${Boost_LIBRARIES}
#  ${TRIANGLE_LIBRARY}
#  ${VTK_LIBRARIES}
#  ${XML2_LIBRARY})

#----------------------------------------------------------------------------
# Test mpm
#----------------------------------------------------------------------------
# add_executable(test_two_contact test_two_contact.cc)
# target_link_libraries(test_two_contact 
#   BrMPM_LIB
#   ${PROBLEMSPEC_LIBRARY}
#   ${Boost_LIBRARIES}
#   ${VTK_LIBRARIES}
#   ${XML2_LIBRARY})

