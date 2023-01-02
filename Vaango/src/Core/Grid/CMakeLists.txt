#  The MIT License
# 
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#  Copyright (c) 2015-2023 Biswajit Banerjee
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

SET(Vaango_Core_Grid_SRCS
  Box.h
  AMR.cc
  AMR_CoarsenRefine.cc
  axiCpdiInterpolator.cc
  AxiGIMPInterpolator.cc
  Box.cc
  BSplineInterpolator.cc
  cpdiInterpolator.cc
  cptiInterpolator.cc
  DbgOutput.cc
  fastAxiCpdiInterpolator.cc
  fastCpdiInterpolator.cc
  Ghost.cc
  GIMPInterpolator.cc
  Grid.cc
  Level.cc
  LinearInterpolator.cc
  AxiLinearInterpolator.cc
  Material.cc
  MaterialManager.cc
  ParticleInterpolator.cc
  Patch.cc
  PatchRangeTree.cc
  Region.cc
  EmptyMaterial.cc
  SimulationState.cc
  SimulationTime.cc
  Task.cc
  TOBSplineInterpolator.cc
  UnknownVariable.cc
)

add_subdirectory(BoundaryConditions)
add_subdirectory(Variables)
add_subdirectory(PatchBVH)

#message(STATUS "Core/Grid sources are ${Vaango_Core_Grid_SRCS}")

ADD_LIBRARY(Vaango_Core_Grid ${Vaango_Core_Grid_SRCS})

TARGET_LINK_LIBRARIES(Vaango_Core_Grid
  Vaango_Core_Exceptions
  Vaango_Core_Thread
  Vaango_Core_Util
  Vaango_Core_Containers
  Vaango_Core_Parallel
  Vaango_Core_ProblemSpec
  Vaango_Core_Math
  ${MPI_LIBRARY}
  ${Z_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_Core_Grid)
ENDIF(BUILD_SHARED_LIBS)

