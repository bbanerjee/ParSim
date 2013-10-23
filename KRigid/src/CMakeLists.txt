# This is the root KRigid CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(KRigid Fortran)

set(MAJOR_VERSION "1")
set(MINOR_VERSION "1")
set(PATCH_VERSION "0")

set(CMAKE_VERBOSE_MAKEFILE 1)
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_Fortran_COMPILER_INIT = ${CMAKE_Fortran_COMPILER_INIT}")
message("CMAKE_Fortran_COMPILER_FULLPATH = ${CMAKE_Fortran_COMPILER_FULLPATH}")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffree-line-length-none -ffixed-line-length-none -Wall")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffree-line-length-none -Wall")

# Needs lapack and blas
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(KRigid_SRC
  inertial_properties.f	
  initial_conditions.f		
  external_forces.f		
  body_forces.f		
  time_parameters.f		
  translational_motion.f		
  rotational_motion.f		
  atimesb.f		
  output.f		
)
add_library(KRigidLib ${KRigid_SRC})

add_executable(KRigid main.f)
target_link_libraries(KRigid 
  KRigidLib
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES})

