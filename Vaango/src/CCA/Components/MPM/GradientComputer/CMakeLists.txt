#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2014-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

# CMakeLists.txt for Vaango/src/CCA/Components/MPM/GradientComputer

set(Vaango_CCA_Components_MPM_GradientComputer_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/DeformationGradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/DisplacementGradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/GradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/VelocityGradientComputer.cc
)

add_library(MPM_GradientComputer ${Vaango_CCA_Components_MPM_GradientComputer_SRCS})

target_link_libraries(MPM_GradientComputer PRIVATE
  ${LIBS}
  MPI::MPI_CXX
  Vaango::Options
  Vaango::Warnings
)
if (PETSC)
   target_link_petsc(MPM_GradientComputer PRIVATE)
endif()

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_MPM_GradientComputer)
ENDIF(BUILD_SHARED_LIBS)

install(TARGETS
  MPM_GradientComputer
)