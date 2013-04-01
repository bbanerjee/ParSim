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
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall")

set(BASEPATH ${CMAKE_SOURCE_DIR})
include_directories(${BASEPATH})

#----------------------------------------------------------------------------
# Needs boost libraries: at least 1.49
#----------------------------------------------------------------------------
find_package(Boost 1.49 COMPONENTS REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

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
