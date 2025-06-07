#
# The MIT License
#
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

SET(Vaango_MPM_Solvers_SRCS
  ImpMPMSolvers/SimpleSolver.cc 
  ImpMPMSolvers/Solver.cc       
)

IF (PETSC)
  SET(Vaango_MPM_Solvers_SRCS 
      ${Vaango_MPM_Solvers_SRCS} 
      ImpMPMSolvers/PetscSolver.cc)
  SET(LIBS ${LIBS} ${PETSC_LIBRARIES})
ELSE ()
  SET(Vaango_MPM_Solvers_SRCS 
      ${Vaango_MPM_Solvers_SRCS} 
      ImpMPMSolvers/FakePetscSolver.cc)
ENDIF ()

SET(Vaango_MPM_DiffusionTasks_SRCS
  ReactionDiffusion/ScalarDiffusionTasks.cc
)

SET(Vaango_CCA_Components_MPM_SRCS
  ${Vaango_MPM_Solvers_SRCS}
  ${Vaango_MPM_DiffusionTasks_SRCS}
  SerialMPM.cc    
  MPM_UpdateStressLast.cc    
  RigidMPM.cc     
  FractureMPM.cc  
  ImpMPM.cc       
  ShellMPM.cc     
  AMRMPM.cc       
  UofU_MPM.cc    
)

add_subdirectory(Core)
add_subdirectory(Contact)
add_subdirectory(ConstitutiveModel)
add_subdirectory(ThermalContact)
add_subdirectory(PhysicalBC)
add_subdirectory(ParticleCreator)
add_subdirectory(CohesiveZone)
add_subdirectory(Crack)
add_subdirectory(HeatConduction)
add_subdirectory(MMS)
add_subdirectory(GradientComputer)
add_subdirectory(ReactionDiffusion)

SET(LIBS
  ${LIBS}
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_ProblemSpec
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_Core_GeometryPiece
  Vaango_Core_Math
  MPM_ParticleCreator
  MPM_GradientComputer
  MPM_ConstitutiveModel
  MPM_PhysicalBC
  MPM_MMS
  MPM_Diffusion
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  MPI::MPI_CXX
)

add_library(Vaango_CCA_Components_MPM ${Vaango_CCA_Components_MPM_SRCS})

target_link_libraries(Vaango_CCA_Components_MPM 
  ${LIBS}
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_CCA_Components_MPM)
ENDIF(BUILD_SHARED_LIBS)
