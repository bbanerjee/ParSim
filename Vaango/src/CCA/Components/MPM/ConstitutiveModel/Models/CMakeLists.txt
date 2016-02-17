#
# The MIT License
#
# Copyright (c) 2004 Scientific Computing and Imaging Institute,
# University of Utah.
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015 Parresia Research Limited, New Zealand
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

SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelStateBase.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Arenisca3.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_CamClay.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_Default.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ModelState_MasonSand.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Arenisca.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Arenisca.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_Constant.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ElasticModuli_MasonSand.cc
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
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Arenisca3.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_CamClay.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_Gurson.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/YieldCond_MasonSand.cc 
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
  ${CMAKE_CURRENT_SOURCE_DIR}/KinematicHardening_MasonSand.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_ArenaKappa.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_BorjaPressure.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_MasonSand.cc 
  PARENT_SCOPE
)

add_subdirectory(UnitTests)
