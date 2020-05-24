#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#  Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

# CMakeLists.txt for Vaango/src/StandAlone

message(STATUS "Binary dir = ${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

ADD_EXECUTABLE(vaango vaango.cc)
ADD_EXECUTABLE(compare_uda compare_uda.cc)
ADD_EXECUTABLE(slb slb.cc)
ADD_EXECUTABLE(restart_merger restart_merger.cc)

SET(VAANGO_LIBS 
        Vaango_Core_Thread       
        Vaango_Core_DataArchive 
        Vaango_Core_Grid        
        Vaango_Core_Parallel    
        Vaango_Core_Labels      
        Vaango_Core_Util        
        Vaango_Core_Math        
        Vaango_Core_Disclosure  
        Vaango_Core_Exceptions  
        Vaango_Core_OS  
        Vaango_CCA_Ports        
        Vaango_CCA_Components_Parent 
        Vaango_CCA_Components_DataArchiver  
        Vaango_CCA_Components_LoadBalancers 
        Vaango_CCA_Components_Regridder     
        Vaango_Core_ProblemSpec             
        Vaango_CCA_Components_SimulationController 
        Vaango_CCA_Components_Schedulers           
        Vaango_CCA_Components_ProblemSpecification 
        Vaango_CCA_Components_Solvers
)

IF(ICE)
  SET(VAANGO_LIBS ${VAANGO_LIBS} Vaango_CCA_Components_Models)
ENDIF(ICE)
IF(HYPRE)
  SET(VAANGO_LIBS ${VAANGO_LIBS} ${HYPRE_LIBRARY})
ENDIF()
IF(PETSC)
  SET(VAANGO_LIBS ${VAANGO_LIBS} ${PETSC_LIBRARY})
ENDIF()

TARGET_LINK_LIBRARIES(vaango ${VAANGO_LIBS})

TARGET_LINK_LIBRARIES(compare_uda
  Vaango_Core_DataArchive
  Vaango_Core_Grid
  Vaango_Core_Malloc
  Vaango_Core_Thread
  Vaango_Core_Math
  Vaango_Core_Util
)

set(RESTART_LIBS
  Vaango_Core_DataArchive
#  Vaango_Core_Grid
#  Vaango_Core_OS
  Vaango_CCA_Components_Parent 
#  Vaango_Core_ProblemSpec
  Vaango_CCA_Components_DataArchiver
)
IF(HYPRE)
  SET(RESTART_LIBS ${RESTART_LIBS} ${HYPRE_LIBRARY})
ENDIF()
TARGET_LINK_LIBRARIES(restart_merger ${RESTART_LIBS})
  
TARGET_LINK_LIBRARIES(slb
  Vaango_Core_GeometryPiece
  Vaango_Core_Grid
  Vaango_Core_ProblemSpec
  Vaango_CCA_Components_ProblemSpecification
)  

# Go to subdirectories
add_subdirectory(tools)
add_subdirectory(Benchmarks)
