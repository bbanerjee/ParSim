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

# CMakeLists.txt for Packages/Uintah/CCA/Components/MPM

SET(Packages_Uintah_CCA_Components_MPM_SRCS
  SerialMPM.cc    
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

SET(LIBS
  Packages_Uintah_CCA_Ports
  Packages_Uintah_Core_Grid
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_Disclosure
  Packages_Uintah_Core_Labels
  Packages_Uintah_Core_ProblemSpec
  Packages_Uintah_Core_Parallel
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_GeometryPiece
  Packages_Uintah_Core_Math
  MPM_ConstitutiveModel
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

IF (PETSC)
  SET(Packages_Uintah_CCA_Components_MPM_SRCS ${Packages_Uintah_CCA_Components_MPM_SRCS} PetscSolver.cc)
  SET(LIBS ${LIBS} ${PETSC_LIBRARIES})
ELSE ()
  SET(Packages_Uintah_CCA_Components_MPM_SRCS ${Packages_Uintah_CCA_Components_MPM_SRCS} FakePetscSolver.cc)
ENDIF ()


ADD_LIBRARY(Packages_Uintah_CCA_Components_MPM ${Packages_Uintah_CCA_Components_MPM_SRCS})

TARGET_LINK_LIBRARIES(Packages_Uintah_CCA_Components_MPM ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_CCA_Components_MPM)
ENDIF(BUILD_SHARED_LIBS)
