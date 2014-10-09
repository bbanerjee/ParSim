# This is the root EMUNE CMakeLists file.

#----------------------------------------------------------------------------
# Create the project and check cmake version
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(EMUNE Fortran)

set(MAJOR_VERSION "1")
set(MINOR_VERSION "3")
set(PATCH_VERSION "0")

#set(CMAKE_VERBOSE_MAKEFILE 1)
message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_Fortran_COMPILER_INIT = ${CMAKE_Fortran_COMPILER_INIT}")
message("CMAKE_Fortran_COMPILER_FULLPATH = ${CMAKE_Fortran_COMPILER_FULLPATH}")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffree-line-length-none -ffixed-line-length-none -Wall")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffree-line-length-none -Wall")

set(EMUNE_SRC
  Damage.f90
  Dynamics.f90
  Get_family.f90
  Global_variables.f90
  Horizon.f90
  Input_subroutines.f90
  Module_precision.f90
  Objects.f90
  Peri_motion.f90
  Solver.f90
  Volume_partition.f90
)
add_library(EMUNE_MODULES ${EMUNE_SRC})

add_executable(emune EMUNE_new.f90)
target_link_libraries(emune EMUNE_MODULES)

