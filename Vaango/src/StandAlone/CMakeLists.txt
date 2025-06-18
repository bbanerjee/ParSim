#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#  Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

message(STATUS "Binary dir = ${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

set(VAANGO_UTILS_SRC
    Utils/vaango_utils.cc
    Utils/vaango_options.cc
)

set(COMPARE_UDA_UTILS_SRC
    Utils/compare_uda_utils.cc
    Utils/compare_uda_options.cc
    Utils/MaterialParticleVarData.cc
    Utils/MaterialParticleData.cc
    Utils/FieldComparator.cc
)

add_executable(vaango vaango.cc ${VAANGO_UTILS_SRC})
add_executable(compare_uda compare_uda.cc ${COMPARE_UDA_UTILS_SRC})
add_executable(slb slb.cc)
add_executable(restart_merger restart_merger.cc)

SET(VAANGO_LIBS 
        Vaango_Core_DataArchive 
        Vaango_Core_Grid        
        Vaango_Core_Parallel    
        Vaango_Core_Util        
        Vaango_Core_Math        
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
        $<$<BOOL:${ICE}>:Vaango_CCA_Components_Models>
        $<$<BOOL:${HAVE_HYPRE}>:${HYPRE_TARGET}>
        $<$<BOOL:${HAVE_PETSC}>:${PETSC_TARGET}>
        MPI::MPI_CXX
        Vaango::Options
        Vaango::Warnings
)

target_link_libraries(vaango PRIVATE
  ${VAANGO_LIBS}
)
target_compile_definitions(vaango PRIVATE
  COMPILE_COMMANDS_PATH="${CMAKE_BINARY_DIR}/compile_commands.json"
  CMAKE_SRC_DIR="${CMAKE_SOURCE_DIR}"
)
target_link_libraries(compare_uda PRIVATE
  ${VAANGO_LIBS}
)


set(RESTART_LIBS
  Vaango_Core_DataArchive
  Vaango_Core_OS  
  Vaango_CCA_Ports
  Vaango_CCA_Components_Parent 
  Vaango_CCA_Components_DataArchiver
  Vaango_CCA_Components_ProblemSpecification 
  MPM_ParticleCreator
  $<$<BOOL:${HAVE_HYPRE}>:${HYPRE_TARGET}>
  MPI::MPI_CXX
  Vaango::Options
  Vaango::Warnings
)
target_link_libraries(restart_merger PRIVATE
  ${RESTART_LIBS}
)
  
target_link_libraries(slb PRIVATE
  Vaango_Core_GeometryPiece
  Vaango_Core_Grid
  Vaango_Core_ProblemSpec
  Vaango_CCA_Components_ProblemSpecification
  Vaango::Options
  Vaango::Warnings
)  

# Go to subdirectories
add_subdirectory(tools)
add_subdirectory(Benchmarks)
