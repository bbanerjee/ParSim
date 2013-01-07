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

# CMakeLists.txt for Vaango/src/CCA/Components/Schedulers

SET(Vaango_CCA_Components_Schedulers_SRCS
  DetailedTasks.cc 
  MemoryLog.cc 
  CommRecMPI.cc 
  MPIScheduler.cc 
  DynamicMPIScheduler.cc 
  MessageLog.cc 
  MixedScheduler.cc 
  OnDemandDataWarehouse.cc 
  Relocate.cc 
  SchedulerCommon.cc 
  SchedulerFactory.cc 
  SendState.cc 
  SingleProcessorScheduler.cc 
  TaskGraph.cc 
  ThreadedMPIScheduler.cc 
  ThreadedMPIScheduler2.cc 
  UnifiedScheduler.cc 
  ThreadPool.cc 
  DependencyException.cc 
  IncorrectAllocation.cc 
  Util.cc 
  templates.cc
)

if (CUDA)
  set(Vaango_CCA_Components_Schedulers_SRCS
    ${Vaango_CCA_Components_Schedulers_SRCS}
    GPUThreadedMPIScheduler.cc
  )
endif (CUDA)

ADD_LIBRARY(Vaango_CCA_Components_Schedulers ${Vaango_CCA_Components_Schedulers_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_Schedulers
  Vaango_CCA_Components_ProblemSpecification
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_Disclosure
  Vaango_Core_ProblemSpec
  Vaango_CCA_Ports
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_Core_Thread
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
