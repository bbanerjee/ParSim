#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Examples

SET(Packages_Uintah_CCA_Components_Examples_SRCS
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
  Packages_Uintah_CCA_Ports
  Packages_Uintah_Core_Grid
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_ProblemSpec
  Packages_Uintah_Core_GeometryPiece
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_Disclosure
  Packages_Uintah_Core_Math
  Packages_Uintah_Core_Parallel
)

if (RADIATION)
  SET(Packages_Uintah_CCA_Components_Examples_SRCS
    ${Packages_Uintah_CCA_Components_Examples_SRCS}
    RMCRT_Test.cc       
  )
  set(Examples_LIBS
    ${Examples_LIBS}
    Packages_Uintah_CCA_Components_Models_Radiation
    Packages_Uintah_CCA_Components_Arches_Radiation
  )
endif (RADIATION)

if (CUDA)
  SET(Packages_Uintah_CCA_Components_Examples_SRCS
    ${Packages_Uintah_CCA_Components_Examples_SRCS}
    PoissonGPU1.cc            
    PoissonGPU1Kernel.cu      
    GPUSchedulerTest.cc       
    UnifiedSchedulerTest.cc   
    GPUSchedulerTestKernel.cu          
  )
endif (CUDA)

ADD_LIBRARY(Packages_Uintah_CCA_Components_Examples ${Packages_Uintah_CCA_Components_Examples_SRCS})

TARGET_LINK_LIBRARIES(Packages_Uintah_CCA_Components_Examples
  ${Examples_LIBS}
  ${XML2_LIBRARY}
  ${M_LIBRARY}
  ${CUDA_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_CCA_Components_Examples)
ENDIF(BUILD_SHARED_LIBS)
