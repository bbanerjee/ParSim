#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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
# CMakeLists.txt for Vaango/src/CCA/Components/MPM/ConstitutiveModel/Biswajit/Models

#-----------------------------------------------------------------------------
# For reading JSON input files
#-----------------------------------------------------------------------------
#set(JSON_INCLUDES ${Vaango_SOURCE_DIR}/submodules/json/src)
#message(STATUS "JSON include dir = ${JSON_INCLUDES}")
#include_directories(${JSON_INCLUDES})
#set(CMAKE_VERBOSE_MAKEFILE 1)

#-----------------------------------------------------------------------------
# The models
#-----------------------------------------------------------------------------
SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelStateBase.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_SoilModelBrannon.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Arenisca3.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_CamClay.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Default.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Tabular.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_TabularCap.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Arenisca.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Arenisca.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Constant.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Arena.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_ArenaMixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Tabular.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_NeuralNet.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuliModel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuliModelFactory.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/PressureModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/PressureModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Borja.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Hypoelastic.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Hyperelastic.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_MieGruneisen.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Air.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Water.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/Pressure_Granite.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCondition.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldConditionFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCondUtils.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Arenisca3.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_CamClay.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Gurson.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_ArenaMixture.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Tabular.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_TabularCap.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_vonMises.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ShearModulusModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ShearModulusModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ShearModulus_Borja.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ShearModulus_Constant.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ShearModulus_Nadal.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardeningModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardeningModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardening_None.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardening_Prager.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardening_Armstrong.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardening_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_SoilModelBrannonKappa.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_BorjaPressure.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_TabularCap.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TensorUtils.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TableUtils.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TabularData.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombState.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombTypes.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombBase.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombSheng.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombClassic.cc 
  PARENT_SCOPE
)

message(STATUS "Compiler is ${CMAKE_CXX_COMPILER_ID}")
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_COMPILER_IS_INTEL 1)
endif()
IF(NOT CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_COMPILER_IS_INTEL) 
  add_subdirectory(UnitTests)
ENDIF()
