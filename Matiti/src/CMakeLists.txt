# This is the root MATITI CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(MATITI CXX)

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
# Add requirements for Output in VTK format
#----------------------------------------------------------------------------
find_package(VTK 5.8 REQUIRED)
# message(STATUS "VTK header directory = ${VTK_INCLUDE_DIRS}")
# message(STATUS "VTK libraries directory = ${VTK_LIBRARIES}")
include_directories(${VTK_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Set up library sources
#----------------------------------------------------------------------------
set(MATITI_SRC
  Peridynamics.cc
  Core/Bond.cc
  Core/Body.cc
  Core/Crack.cc
  Core/Domain.cc
  Core/Element.cc
  Core/FamilyComputer.cc
  Core/HorizonComputer.cc
  Core/Node.cc
  Core/SimulationState.cc
  Core/Time.cc
  BoundaryConditions/LoadBC.cc
  BoundaryConditions/LoadBCFactory.cc
  BoundaryConditions/ForceBC.cc
  BoundaryConditions/TractionBC.cc
  BoundaryConditions/InitialConditions.cc
  BoundaryConditions/VelocityBC.cc
  InputOutput/Output.cc
  InputOutput/OutputVTK.cc
  InputOutput/ProblemSpecReader.cc
  InputOutput/ProblemSpecUtil.cc
  MaterialModels/DamageModel.cc
  MaterialModels/Material.cc
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
  GeometryPiece/PlaneGeometryReader.cc
)

add_library(MATITI_LIB 
  ${MATITI_SRC}
  ${GEOMETRY_SRC}
  ${GEOMETRY_PIECE_SRC})

#----------------------------------------------------------------------------
# Just test
#----------------------------------------------------------------------------
#add_executable(emu2dc emu2dc.cc)
#target_link_libraries(emu2dc 
#  MATITI_LIB 
#  ${PROBLEMSPEC_LIBRARY}
#  ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Test map
#----------------------------------------------------------------------------
add_executable(test_map StandAlone/test_map.cc)
target_link_libraries(test_map 
  MATITI_LIB
  ${PROBLEMSPEC_LIBRARY}
  ${TRIANGLE_LIBRARY}
  ${VTK_LIBRARIES}
  ${XML2_LIBRARY})

#----------------------------------------------------------------------------
# Test point search
#----------------------------------------------------------------------------
add_executable(test_point_search StandAlone/test_point_search.cc)
target_link_libraries(test_point_search 
  MATITI_LIB 
  ${PROBLEMSPEC_LIBRARY}
  ${TRIANGLE_LIBRARY}
  ${VTK_LIBRARIES}
  ${XML2_LIBRARY})

#----------------------------------------------------------------------------
# Test input file reader
#----------------------------------------------------------------------------
#foreach(dir ${PROBLEMSPEC_INCLUDE_DIR})
# message(STATUS "Dir = ${dir}")
#endforeach()

add_executable(test_reader StandAlone/test_reader.cc)
target_link_libraries(test_reader 
  MATITI_LIB
  ${PROBLEMSPEC_LIBRARY}
  ${Boost_LIBRARIES}
  ${TRIANGLE_LIBRARY}
  ${VTK_LIBRARIES}
  ${XML2_LIBRARY})

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

#----------------------------------------------------------------------------
# Test peridynamics
#----------------------------------------------------------------------------
add_executable(test_peri StandAlone/test_peri.cc)
target_link_libraries(test_peri 
  MATITI_LIB
  ${PROBLEMSPEC_LIBRARY}
  ${Boost_LIBRARIES}
  ${TRIANGLE_LIBRARY}
  ${VTK_LIBRARIES}
  ${XML2_LIBRARY})

