#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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
SET(Vaango_CCA_Components_MPM_SRCS
  SerialMPM.cc    
  MPM_UpdateStressLast.cc    
  RigidMPM.cc     
  MPMCommon.cc    
  FractureMPM.cc  
  ImpMPM.cc       
  ShellMPM.cc     
  AMRMPM.cc       
  SimpleSolver.cc 
  Solver.cc       
  MPMBoundCond.cc 
  MPMFlags.cc	  
  ImpMPMFlags.cc
  UofU_MPM.cc    
)

add_subdirectory(ConstitutiveModel)
add_subdirectory(Contact)
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
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_Labels
  Vaango_Core_ProblemSpec
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_Core_GeometryPiece
  Vaango_Core_Math
  MPM_ConstitutiveModel
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

IF (PETSC)
  SET(Vaango_CCA_Components_MPM_SRCS ${Vaango_CCA_Components_MPM_SRCS} PetscSolver.cc)
  SET(LIBS ${LIBS} ${PETSC_LIBRARIES})
ELSE ()
  SET(Vaango_CCA_Components_MPM_SRCS ${Vaango_CCA_Components_MPM_SRCS} FakePetscSolver.cc)
ENDIF ()


ADD_LIBRARY(Vaango_CCA_Components_MPM ${Vaango_CCA_Components_MPM_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_MPM ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_CCA_Components_MPM)
ENDIF(BUILD_SHARED_LIBS)
