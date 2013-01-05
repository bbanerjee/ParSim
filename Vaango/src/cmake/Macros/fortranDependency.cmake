# Adds fortran file to list of srcs and creates the .h file for it
macro(FORTRAN_DEPENDENCY fortran_file lib)

  # should pass in a fortran_file and package it belongs to.  I.e.,
  # FORTRAN_DEPENDENCY(fortran/bcscalar_fort.h CCA_Components_Arches)
  string(REPLACE .F _fort.h header ${CMAKE_CURRENT_BINARY_DIR}/${fortran_file})
  string(REPLACE .F .fspec fspec_file ${fortran_file})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran)
  set(${lib}_FORTRAN_SRCS ${${lib}_FORTRAN_SRCS} ${fortran_file})
  # if everybody has cmake 2.4.5, we can do IS_NEWER_THAN, and not create these files on every cmake
  if(${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} IS_NEWER_THAN ${header})
    if(NOT EXISTS ${header})
      execute_process(COMMAND perl ${FSPEC} ${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} ${header})
    endif(NOT EXISTS ${header})
  endif(${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} IS_NEWER_THAN ${header})

endmacro(FORTRAN_DEPENDENCY)


