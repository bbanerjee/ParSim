#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2023 Biswajit Banerjee
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
# CMakeLists.txt for Vaango/src/CCA/Components/Schedulers

SET(Vaango_CCA_Components_Schedulers_SRCS
  DependencyBatch.cc 
  DependencyException.cc 
  DetailedDependency.cc 
  DetailedTask.cc 
  DetailedTasks.cc 
  DynamicMPIScheduler.cc 
  MemoryLog.cc 
  MPIScheduler.cc 
  OnDemandDataWarehouse.cc 
  Relocate.cc 
  RelocateScatterRecord.cc 
  RelocateMPIScatterRecord.cc 
  RuntimeStats.cc 
  SchedulerCommon.cc 
  SchedulerFactory.cc 
  SendState.cc 
  TaskGraph.cc 
  TaskGraphCompTable.cc 
  UnifiedScheduler.cc 
)

if (CUDA)
  set(Vaango_CCA_Components_Schedulers_SRCS
    ${Vaango_CCA_Components_Schedulers_SRCS}
    GPUDataWarehouse.cu
    GPUGridVariableInfo.cc
    GPUGridVariableGhosts.cc
    GPUMemoryPool.cc
  )
  set(CUDA_OBJECT
      GPUDataWarehouse.o
     )
endif (CUDA)

ADD_LIBRARY(Vaango_CCA_Components_Schedulers ${Vaango_CCA_Components_Schedulers_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_Schedulers
  Vaango_CCA_Components_ProblemSpecification
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_ProblemSpec
  Vaango_CCA_Ports
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  
  ${XML2_LIBRARY}
  ${TAU_LIBRARY}
  ${MPI_LIBS}
  ${VAMPIR_LIBRARY}
  ${PERFEX_LIBRARY}
  ${CUDA_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_CCA_Components_Schedulers)
ENDIF(BUILD_SHARED_LIBS)
