#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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


