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
SET(MPM_ConstitutiveModel_PlasticityModels_SRCS
  PlasticityState.cc 
  DamageModel.cc 
  DamageModelFactory.cc 
  NullDamage.cc 
  JohnsonCookDamage.cc 
  HancockMacKenzieDamage.cc 
  MPMEquationOfState.cc 
  MPMEquationOfStateFactory.cc 
  DefaultHypoElasticEOS.cc 
  HyperElasticEOS.cc 
  MieGruneisenEOS.cc 
  MieGruneisenEOSEnergy.cc 
  FlowModel.cc 
  FlowStressModelFactory.cc 
  IsoHardeningFlow.cc 
  JohnsonCookFlow.cc 
  ZAFlow.cc
  ZAPolymerFlow.cc
  MTSFlow.cc 
  SCGFlow.cc 
  PTWFlow.cc 
  YieldCondition.cc 
  YieldConditionFactory.cc 
  GursonYield.cc 
  VonMisesYield.cc 
  StabilityCheck.cc 
  StabilityCheckFactory.cc 
  AcousticTensorCheck.cc 
  BeckerCheck.cc 
  DruckerCheck.cc 
  NoneCheck.cc 
  DruckerBeckerCheck.cc 
  ShearModulusModel.cc 
  ShearModulusModelFactory.cc 
  ConstantShear.cc 
  MTSShear.cc 
  NPShear.cc 
  PTWShear.cc 
  SCGShear.cc 
  SpecificHeatModel.cc 
  SpecificHeatModelFactory.cc 
  ConstantCp.cc 
  CubicCp.cc 
  CopperCp.cc 
  SteelCp.cc 
  MeltingTempModel.cc 
  MeltingTempModelFactory.cc 
  ConstantMeltTemp.cc 
  LinearMeltTemp.cc 
  SCGMeltTemp.cc 
  BPSMeltTemp.cc 
  ViscoPlasticityModel.cc 
  ViscoPlasticityModelFactory.cc 
  SuvicI.cc 
  KinematicHardeningModel.cc 
  KinematicHardeningModelFactory.cc 
  NoKinematicHardening.cc 
  PragerKinematicHardening.cc 
  ArmstrongFrederickKinematicHardening.cc 
  DevStressModel.cc 
  DevStressModelFactory.cc 
  HypoViscoElasticDevStress.cc 
  HypoElasticDevStress.cc 
)

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
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

ADD_LIBRARY(MPM_ConstitutiveModel_PlasticityModels ${MPM_ConstitutiveModel_PlasticityModels_SRCS})

TARGET_LINK_LIBRARIES(MPM_ConstitutiveModel_PlasticityModels ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_MPM_ConstitutiveModel_PlasticityModels)
ENDIF(BUILD_SHARED_LIBS)

