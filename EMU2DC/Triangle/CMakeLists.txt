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
set(BASE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/include")
set(BASE_SOURCE_PATH "${CMAKE_SOURCE_DIR}/src")

#----------------------------------------------------------------------------
# Build triangle++
#----------------------------------------------------------------------------
set(TFLAGS "-DREDUCED -DANSI_DECLARATORS -DTRILIBRARY -DCDT_ONLY -DLINUX")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TFLAGS}")

set(TPP_SRC
  ${BASE_SOURCE_PATH}/assert.cpp
  ${BASE_SOURCE_PATH}/del_impl.cpp
)
include_directories(${BASE_INCLUDE_PATH})

add_library(Triangle++ ${TPP_SRC})

add_executable(tpp ${BASE_SOURCE_PATH}/main.cpp)
target_link_libraries(tpp 
  Triangle++)

