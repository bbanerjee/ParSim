#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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


# CMakeLists.txt for Vaango/src/CCA/Components/MPM/ConstitutiveModel/fortran

set(FORTRAN_CM_FLAGS "${CMAKE_Fortran_FLAGS}") # Only general compiler flags here
message(STATUS "Fortran general flags: ${FORTRAN_CM_FLAGS}")

set(FORTRAN_CM_DEFINITIONS I1SOFTLIM NEWXSV DILAT_LIM)
message(STATUS "Fortran definitions: ${FORTRAN_CM_DEFINITIONS}")

set(FORTRAN_CM_SRC 
  ${CMAKE_CURRENT_SOURCE_DIR}/Hooke.F90  
  ${CMAKE_CURRENT_SOURCE_DIR}/Kayenta_calcs.F  
  ${CMAKE_CURRENT_SOURCE_DIR}/MIGUtilsF.F  
  ${CMAKE_CURRENT_SOURCE_DIR}/diamm.F
  ${CMAKE_CURRENT_SOURCE_DIR}/visco.F90
)

set(FORTRAN_CM_CPP_SRC 
  ${CMAKE_CURRENT_SOURCE_DIR}/MIGUtils.cc   
)

set_source_files_properties(
  ${FORTRAN_CM_SRC}
  PROPERTIES
    COMPILE_OPTIONS "${FORTRAN_CM_FLAGS};-Wno-unused-dummy-argument;-Wno-unused-parameter"
)

add_library(MPM_ConstitutiveModel_FORTRAN ${FORTRAN_CM_SRC})

target_compile_definitions(MPM_ConstitutiveModel_FORTRAN PRIVATE
    ${FORTRAN_CM_DEFINITIONS}
)

target_link_libraries(MPM_ConstitutiveModel_FORTRAN PRIVATE
  ${F_LIBRARY}
  Vaango::Options
  Vaango::Warnings
)

SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${FORTRAN_CM_CPP_SRC}
  PARENT_SCOPE
)

