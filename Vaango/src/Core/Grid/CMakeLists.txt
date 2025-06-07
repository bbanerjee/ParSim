#  The MIT License
# 
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#  Copyright (c) 2015-2024 Biswajit Banerjee
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

# CMakeLists.txt for Vaango/src/Core/Grid

SET(Vaango_Core_Grid_MPMInterpolator_SRCS
    MPMInterpolators/AxiCPDIInterpolator.cc
    MPMInterpolators/AxiGIMPInterpolator.cc
    MPMInterpolators/BSplineInterpolator.cc
    MPMInterpolators/CPDIInterpolator.cc
    MPMInterpolators/CPTIInterpolator.cc
    MPMInterpolators/FastAxiCPDIInterpolator.cc
    MPMInterpolators/FastCPDIInterpolator.cc
    MPMInterpolators/GIMPInterpolator.cc
    MPMInterpolators/LinearInterpolator.cc
    MPMInterpolators/AxiLinearInterpolator.cc
    MPMInterpolators/ParticleInterpolator.cc
    MPMInterpolators/TOBSplineInterpolator.cc
)

SET(Vaango_Core_Grid_SRCS
  Box.h
  AMR.cc
  AMR_CoarsenRefine.cc
  Box.cc
  DbgOutput.cc
  Ghost.cc
  Grid.cc
  Level.cc
  Material.cc
  MaterialManager.cc
  Patch.cc
  PatchRangeTree.cc
  Region.cc
  EmptyMaterial.cc
  Task.cc
  UnknownVariable.cc
  ${Vaango_Core_Grid_MPMInterpolator_SRCS}
)

add_subdirectory(BoundaryConditions)
add_subdirectory(Variables)
add_subdirectory(PatchBVH)
add_subdirectory(UnitTests)

#message(STATUS "Core/Grid sources are ${Vaango_Core_Grid_SRCS}")

add_library(Vaango_Core_Grid ${Vaango_Core_Grid_SRCS})

target_link_libraries(Vaango_Core_Grid
  Vaango_Core_Exceptions
  Vaango_Core_Util
  Vaango_Core_Containers
  Vaango_Core_Parallel
  Vaango_Core_ProblemSpec
  Vaango_Core_Math
  ${MPI_LIBRARY}
  ${Z_LIBRARY}
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_Core_Grid)
ENDIF(BUILD_SHARED_LIBS)

