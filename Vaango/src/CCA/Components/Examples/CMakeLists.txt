#
# The MIT License
#
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

SET(Vaango_CCA_Components_Examples_SRCS
  AMRWave.cc           
  Wave.cc              
  Poisson1.cc          
  Poisson2.cc          
  Poisson3.cc          
  Poisson4.cc          
  Burger.cc            
  HeatEquation.cc      
  ParticleTest1.cc     
  Interpolator.cc      
  ExamplesLabel.cc     
  RegridderTest.cc     
  SolverTest1.cc       
  RegionDB.cc          
  Benchmark.cc         
)

set(Examples_LIBS
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_ProblemSpec
  Vaango_Core_GeometryPiece
  Vaango_Core_Exceptions
  Vaango_Core_Math
  Vaango_Core_Parallel
)

if (HAVE_HYPRE)
  SET(Vaango_CCA_Components_Examples_SRCS
    ${Vaango_CCA_Components_Examples_SRCS}
    SolverTest2.cc       
  )
endif (HAVE_HYPRE)


if (CUDA)
  SET(Vaango_CCA_Components_Examples_SRCS
    ${Vaango_CCA_Components_Examples_SRCS}
    PoissonGPU1.cc            
    PoissonGPU1Kernel.cu      
    GPUSchedulerTest.cc       
    UnifiedSchedulerTest.cc   
    GPUSchedulerTestKernel.cu          
  )
endif (CUDA)

add_library(Vaango_CCA_Components_Examples ${Vaango_CCA_Components_Examples_SRCS})

target_link_libraries(Vaango_CCA_Components_Examples PRIVATE
  ${Examples_LIBS}
  ${CUDA_LIBRARY}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
  $<$<BOOL:${HAVE_HYPRE}>:${HYPRE_TARGET}>
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_CCA_Components_Examples)
ENDIF(BUILD_SHARED_LIBS)

install(TARGETS
  Vaango_CCA_Components_Examples
)