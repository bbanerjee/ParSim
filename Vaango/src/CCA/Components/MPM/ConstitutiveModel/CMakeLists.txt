#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2020 Parresia Reaserch Limited, New Zealand
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
# CMakeLists.txt for Vaango/src/CCA/Components/MPM/ConstitutiveModel

SET(MPM_ConstitutiveModel_SRCS
  ConstitutiveModelFactory.cc   
  ConstitutiveModel.cc          
  ImplicitCM.cc                 
  MPMMaterial.cc                
)

#message(STATUS "Compiling FORTRAN code? ${FORTRAN}")
if (FORTRAN)
  add_subdirectory(fortran)
endif(FORTRAN)

add_subdirectory(DamageModels)
add_subdirectory(ElasticModels)
add_subdirectory(ExplosiveModels)
add_subdirectory(J2PlasticModels)
add_subdirectory(ManufacturedSolutions)
add_subdirectory(PorousModels)
add_subdirectory(RockSoilModels)
add_subdirectory(SpecialPurposeModels)
add_subdirectory(TabularModels)
add_subdirectory(ViscoElasticModels)

add_subdirectory(InternalVarModels)
add_subdirectory(KinHardeningModels)
add_subdirectory(PressureModels)
add_subdirectory(ShearModulusModels)
add_subdirectory(YieldCondModels)

add_subdirectory(ModelState)
add_subdirectory(Utilities)

add_subdirectory(Models)
add_subdirectory(J2PlasticSubmodels)
#add_subdirectory(ModelsUnderDevelopment)

SET(LIBS
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_Disclosure
  Vaango_Core_Labels
  Vaango_Core_ProblemSpec
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_Core_GeometryPiece
  Vaango_Core_Math
  MPM_ConstitutiveModel_PlasticityModels
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
  ${HDF5_LIBRARY}
)

#message(STATUS "Sources: ${MPM_ConstitutiveModel_SRCS}")

ADD_LIBRARY(MPM_ConstitutiveModel ${MPM_ConstitutiveModel_SRCS})

TARGET_LINK_LIBRARIES(MPM_ConstitutiveModel ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_MPM_ConstitutiveModel)
ENDIF(BUILD_SHARED_LIBS)

IF(NOT CMAKE_COMPILER_IS_GNUCXX) 
  add_subdirectory(UnitTests)
ENDIF()


