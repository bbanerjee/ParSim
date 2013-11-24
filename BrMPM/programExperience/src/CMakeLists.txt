# This is the root BrMPM CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(programExperience CXX)

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


#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.49
#----------------------------------------------------------------------------
find_package(Boost 1.49 COMPONENTS REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Test variant
#----------------------------------------------------------------------------
add_executable(test_data test_data.cc)
target_link_libraries(test_data
  programExperience_LIB  
  ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Set up library sources
#----------------------------------------------------------------------------
set(programExperience_SRC
  Vector3D.cc
  Point3D.cc
  IntVector3D.cc
  Matrix3D.cc
)
#set(GEOMETRY_MATH_SRC
 # GeometryMath/Box3D.cc
 # GeometryMath/IntVector3D.cc
#  GeometryMath/Point3D.cc
 # GeometryMath/Polygon3D.cc
#  GeometryMath/Vector3D.cc
#  GeometryMath/Matrix3D.cc
#  GeometryMath/MPMMatrix.cc
#)
#set(GEOMETRY_PIECE_SRC
#  GeometryPiece/BoxGeometryPiece.cc
#  GeometryPiece/GeometryPiece.cc
#  GeometryPiece/GeometryPieceFactory.cc
#  GeometryPiece/GeometryReader.cc
#)
#set(SHAPE_FUNCTION_SRC
#  ShapeFunctions/MPMShapeFunctionFactory.cc
 # ShapeFunctions/LinearShapeFunction.cc
#)
#set(CONTACT_SRC
 # Contact/MPMContactFactory.cc
#  Contact/MPMContact.cc
#  Contact/MPMFreeContact.cc
#  Contact/MPMFrictionContact.cc
#  Contact/MPMFrictionlessContact.cc
#  Contact/MPMVelocityContact.cc
#)

#add_library(BrMPM_LIB 
add_library(programExperience_LIB 
   ${programExperience_SRC})

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

