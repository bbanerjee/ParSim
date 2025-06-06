/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef UINTAH_HOMEBREW_MPMLABEL_H
#define UINTAH_HOMEBREW_MPMLABEL_H

#include <memory>
#include <vector>

namespace Uintah {

using std::vector;

class MPMDiffusionLabel;
class VarLabel;

class MPMLabel
{
public:
  MPMLabel();
  ~MPMLabel();

  const VarLabel* timeStepLabel;
  const VarLabel* simulationTimeLabel;
  const VarLabel* delTLabel;
  const VarLabel* doMechLabel;

  // Label containing subclasses.
  std::unique_ptr<MPMDiffusionLabel> diffusion;

  const VarLabel* partCountLabel;

  // Heat flux from fire
  const VarLabel* heatRate_CCLabel;

  // non PermanentParticleState
  const VarLabel* pTemperatureGradientLabel; // for heat conduction
  const VarLabel* pPressureLabel;
  const VarLabel* pScratchVecLabel;
  const VarLabel* pScratchLabel;
  const VarLabel* pLocalizedMPMLabel;
  const VarLabel* pLocalizedMPMLabel_preReloc;
  const VarLabel* pRemoveLabel;
  const VarLabel* pRemoveLabel_preReloc;
  const VarLabel* TotalVolumeDeformedLabel;
  const VarLabel* pXXLabel;
  const VarLabel* pPartitionUnityLabel;

  // Two more labels for velGrad and defGrad (with least disruption
  // of existing code in mind).  These labels are used in the MPM
  // task that compute the velocity gradient and the defromation
  // gradient.
  const VarLabel* pVelGradLabel;
  const VarLabel* pVelGradLabel_preReloc;
  const VarLabel* pDispGradLabel;
  const VarLabel* pDispGradLabel_preReloc;
  const VarLabel* pDeformRateMidLabel;
  const VarLabel* pDefGradLabel;
  const VarLabel* pDefGradMidLabel;
  const VarLabel* pDefGradLabel_preReloc;
  const VarLabel* pPolarDecompRLabel;
  const VarLabel* pPolarDecompRLabel_preReloc;
  const VarLabel* pPolarDecompRMidLabel;

  // PermanentParticleState
  const VarLabel* pStressLabel;
  const VarLabel* pStressLabel_preReloc;
  const VarLabel* pStressUnrotatedLabel;
  const VarLabel* pVolumeLabel;
  const VarLabel* pVolumeMidLabel;
  const VarLabel* pVolumeLabel_preReloc;
  const VarLabel* pMassLabel;
  const VarLabel* pMassLabel_preReloc;
  const VarLabel* pDispLabel;
  const VarLabel* pDispLabel_preReloc;
  const VarLabel* pVelocityLabel;
  const VarLabel* pVelocityLabel_preReloc;
  const VarLabel* pAccelerationLabel;
  const VarLabel* pAccelerationLabel_preReloc;
  const VarLabel* pVelocityXPICLabel;
  const VarLabel* pVelocityXPICLabel_preReloc;
  const VarLabel* pCoriolisImportanceLabel;
  const VarLabel* pCoriolisImportanceLabel_preReloc;
  const VarLabel* pBodyForceAccLabel;
  const VarLabel* pBodyForceAccLabel_preReloc;
  const VarLabel* pExternalForceLabel;
  const VarLabel* pExternalForceCorner1Label;
  const VarLabel* pExternalForceCorner2Label;
  const VarLabel* pExternalForceCorner3Label;
  const VarLabel* pExternalForceCorner4Label;
  const VarLabel* pExtForceLabel_preReloc;
  const VarLabel* pExtForceCorner1Label_preReloc;
  const VarLabel* pExtForceCorner2Label_preReloc;
  const VarLabel* pExtForceCorner3Label_preReloc;
  const VarLabel* pExtForceCorner4Label_preReloc;
  const VarLabel* pXLabel;
  const VarLabel* pXLabel_preReloc;
  const VarLabel* pSurfLabel;
  const VarLabel* pSurfLabel_preReloc;
  const VarLabel* pTemperatureLabel;               // for heat conduction
  const VarLabel* pTemperatureLabel_preReloc;      // for heat conduction
  const VarLabel* pTempCurrentLabel;               // for thermal stress
  const VarLabel* pTempPreviousLabel;              // for thermal stress
  const VarLabel* pTempPreviousLabel_preReloc;     // for thermal stress
  const VarLabel* pdTdtLabel;                      // for heat conduction
  const VarLabel* pdTdtLabel_preReloc;             // for heat conduction
  const VarLabel* pExternalHeatRateLabel;          // for heat conduction
  const VarLabel* pExternalHeatRateLabel_preReloc; // for heat conduction
  const VarLabel* pExternalHeatFluxLabel;          // for heat conduction
  const VarLabel* pExternalHeatFluxLabel_preReloc; // for heat conduction
  const VarLabel* pParticleIDLabel;
  const VarLabel* pParticleIDLabel_preReloc;
  const VarLabel* pSizeLabel;
  const VarLabel* pCurSizeLabel;
  const VarLabel* pSizeLabel_preReloc;

  const VarLabel* pFiberDirLabel;
  const VarLabel* pFiberDirLabel_preReloc;

  const VarLabel* pScaleFactorLabel;
  const VarLabel* pScaleFactorLabel_preReloc;

  const VarLabel* gColorLabel;
  const VarLabel* gLambdaDotLabel;
  const VarLabel* gMassLabel;
  const VarLabel* gMassAllLabel;
  const VarLabel* gPositionLabel;
  const VarLabel* gVelocityLabel;
  const VarLabel* gAccelerationLabel;
  const VarLabel* gVelocityXPICLabel;
  const VarLabel* gVelocityBCLabel;
  const VarLabel* gVelocityStarLabel;
  const VarLabel* gInternalForceLabel;
  const VarLabel* gBodyForceLabel;
  const VarLabel* gExternalForceLabel;
  const VarLabel* NC_CCweightLabel;
  const VarLabel* gContactLabel;
  const VarLabel* gTemperatureRateLabel; // for heat conduction
  const VarLabel* gTemperatureLabel;     // for heat conduction
  const VarLabel* gSpecificVolumeLabel;          // specific volume
  const VarLabel* gSp_vol_srcLabel;      // specific volume
  const VarLabel* gTemperatureNoBCLabel; // for heat conduction
  const VarLabel* gTemperatureStarLabel; // for heat conduction
  const VarLabel* gdTdtLabel;
  const VarLabel* gHeatFluxLabel;
  const VarLabel* gExternalHeatRateLabel;
  const VarLabel* gExternalHeatFluxLabel;
  const VarLabel* gThermalContactTemperatureRateLabel;
  const VarLabel* gNormTractionLabel;
  const VarLabel* gSurfNormLabel;
  const VarLabel* gStressLabel;
  const VarLabel* gStressForSavingLabel;
  const VarLabel* gVolumeLabel;
  const VarLabel* gZOILabel;
  const VarLabel* cVolumeLabel;
  const VarLabel* numLocInCellLabel;
  const VarLabel* numInCellLabel;
  const VarLabel* gradPAccNCLabel;
  const VarLabel* dTdt_NCLabel;          // for heat conduction
  const VarLabel* massBurnFractionLabel; // for burn modeling
  const VarLabel* frictionalWorkLabel;
  const VarLabel* gNumNearParticlesLabel;

  const VarLabel* StrainEnergyLabel;
  const VarLabel* AccStrainEnergyLabel;
  const VarLabel* KineticEnergyLabel;
  const VarLabel* ThermalEnergyLabel;
  const VarLabel* TotalMassLabel;
  const VarLabel* NeedAddMPMMaterialLabel;
  const VarLabel* BndyForceLabel[6];
  const VarLabel* BndyTractionLabel[6];
  const VarLabel* BndyContactAreaLabel[6];
  const VarLabel* BndyContactCellAreaLabel[6];
  const VarLabel* CenterOfMassPositionLabel;
  const VarLabel* TotalMomentumLabel;
  const VarLabel* RigidReactionForceLabel;
  const VarLabel* RigidReactionTorqueLabel{ nullptr };
  const VarLabel* TotalLocalizedParticleLabel;

  // Needs to be modified (flagged as a reduction task in runtime)
  VarLabel* SumTransmittedForceLabel{ nullptr };

  const VarLabel* pCellNAPIDLabel;

  // Implicit MPM labels
  const VarLabel* gVelocityOldLabel;
  const VarLabel* dispNewLabel;
  const VarLabel* dispIncLabel;
  const VarLabel* dispIncQNorm0;
  const VarLabel* dispIncNormMax;
  const VarLabel* dispIncQNorm;
  const VarLabel* dispIncNorm;

  // Labels for particle erosion
  const VarLabel* pErosionLabel;
  const VarLabel* pErosionLabel_preReloc;

  // MPM Physical BC labels (permanent particle state)
  const VarLabel* materialPointsPerLoadCurveLabel;
  const VarLabel* pLoadCurveIDLabel;
  const VarLabel* pLoadCurveIDLabel_preReloc;

  const VarLabel* p_qLabel;
  const VarLabel* p_qLabel_preReloc;

  // for Fracture ----------
  const VarLabel* pDispGradsLabel;
  const VarLabel* pDispGradsLabel_preReloc;
  const VarLabel* pStrainEnergyDensityLabel;
  const VarLabel* pStrainEnergyDensityLabel_preReloc;

  const VarLabel* pgCodeLabel;
  const VarLabel* pKineticEnergyDensityLabel;
  const VarLabel* pVelGradsLabel;

  const VarLabel* gNumPatlsLabel;
  const VarLabel* GNumPatlsLabel;
  const VarLabel* gDisplacementLabel;
  const VarLabel* GDisplacementLabel;
  const VarLabel* gGridStressLabel;
  const VarLabel* GGridStressLabel;
  const VarLabel* gDispGradsLabel;
  const VarLabel* GDispGradsLabel;
  const VarLabel* gVelGradsLabel;
  const VarLabel* GVelGradsLabel;
  const VarLabel* gStrainEnergyDensityLabel;
  const VarLabel* GStrainEnergyDensityLabel;
  const VarLabel* gKineticEnergyDensityLabel;
  const VarLabel* GKineticEnergyDensityLabel;

  const VarLabel* GCrackNormLabel;
  const VarLabel* GMassLabel;
  const VarLabel* GVolumeLabel;
  const VarLabel* GVelocityLabel;
  const VarLabel* GTemperatureLabel;
  const VarLabel* GTemperatureNoBCLabel;
  const VarLabel* GExternalForceLabel;
  const VarLabel* GExternalHeatRateLabel;
  const VarLabel* GThermalContactTemperatureRateLabel;
  const VarLabel* GInternalForceLabel;
  const VarLabel* GdTdtLabel;
  const VarLabel* GTemperatureRateLabel;
  const VarLabel* GTemperatureStarLabel;
  const VarLabel* GVelocityStarLabel;
  const VarLabel* GAccelerationLabel;
  const VarLabel* GSp_volLabel;
  const VarLabel* GSp_vol_srcLabel;
  // ------------------------------

  // Labels for shell materials
  const VarLabel* pThickTopLabel;
  const VarLabel* pInitialThickTopLabel;
  const VarLabel* pThickBotLabel;
  const VarLabel* pInitialThickBotLabel;
  const VarLabel* pNormalLabel;
  const VarLabel* pInitialNormalLabel;
  const VarLabel* pThickTopLabel_preReloc;
  const VarLabel* pInitialThickTopLabel_preReloc;
  const VarLabel* pThickBotLabel_preReloc;
  const VarLabel* pInitialThickBotLabel_preReloc;
  const VarLabel* pNormalLabel_preReloc;
  const VarLabel* pInitialNormalLabel_preReloc;
  const VarLabel* pTypeLabel;
  const VarLabel* pTypeLabel_preReloc;

  const VarLabel* gNormalRotRateLabel;
  const VarLabel* gNormalRotMomentLabel;
  const VarLabel* gNormalRotMassLabel;
  const VarLabel* gNormalRotAccLabel;

  // Debugging Labels
  const VarLabel* pColorLabel;
  const VarLabel* pColorLabel_preReloc;

  // For Cohesive Zones
  const VarLabel* czAreaLabel;
  const VarLabel* czAreaLabel_preReloc;
  const VarLabel* czNormLabel;
  const VarLabel* czNormLabel_preReloc;
  const VarLabel* czTangLabel;
  const VarLabel* czTangLabel_preReloc;
  const VarLabel* czDispTopLabel;
  const VarLabel* czDispTopLabel_preReloc;
  const VarLabel* czDispBottomLabel;
  const VarLabel* czDispBottomLabel_preReloc;
  const VarLabel* pCZSeparationLabel;
  const VarLabel* pCZSeparationLabel_preReloc;
  const VarLabel* pCZForceLabel;
  const VarLabel* pCZForceLabel_preReloc;
  const VarLabel* pCZTopMatLabel;
  const VarLabel* pCZTopMatLabel_preReloc;
  const VarLabel* pCZBotMatLabel;
  const VarLabel* pCZBotMatLabel_preReloc;
  const VarLabel* pCZFailedLabel;
  const VarLabel* pCZFailedLabel_preReloc;
  const VarLabel* pCZIDLabel;
  const VarLabel* pCZIDLabel_preReloc;
  const VarLabel* pCellNACZIDLabel;

  // For adaptive mesh refinement
  const VarLabel* pRefinedLabel;
  const VarLabel* pRefinedLabel_preReloc;
  const VarLabel* pLastLevelLabel;
  const VarLabel* pLastLevelLabel_preReloc;

  const VarLabel* MPMRefineCellLabel;

  // For friction contact
  const VarLabel* gMatlProminenceLabel;
  const VarLabel* gAlphaMaterialLabel;
  const VarLabel* gNormAlphaToBetaLabel;

  // Hydro-mechanical coupling
  const VarLabel* ccPorosity;
  const VarLabel* ccPorePressure;
  const VarLabel* ccPorePressureOld;
  const VarLabel* ccRHS_FlowEquation;
  const VarLabel* ccTransmissivityMatrix;
  const VarLabel* pFluidMassLabel;
  const VarLabel* pFluidVelocityLabel;
  const VarLabel* pFluidAccelerationLabel;
  const VarLabel* pSolidMassLabel;
  const VarLabel* pPorosityLabel;
  const VarLabel* pPorosityLabel_preReloc;
  const VarLabel* pPrescribedPorePressureLabel;
  const VarLabel* pPorePressureLabel;
  const VarLabel* pPorePressureFilterLabel;

  const VarLabel* pStressRateLabel;
  const VarLabel* pStressRateLabel_preReloc;
  const VarLabel* gFluidMassBarLabel;
  const VarLabel* gFluidMassLabel;
  const VarLabel* gFluidVelocityLabel;
  const VarLabel* FluidVelInc;
  const VarLabel* gFluidVelocityStarLabel;
  const VarLabel* gFluidAccelerationLabel;
  const VarLabel* gInternalFluidForceLabel;
  const VarLabel* gExternalFluidForceLabel;
  const VarLabel* gInternalDragForceLabel;
  const VarLabel* gFlowInertiaForceLabel;
  const VarLabel* gPorePressureLabel;
  const VarLabel* gPorePressureFilterLabel;

  const VarLabel* pFluidMassLabel_preReloc;
  const VarLabel* pFluidVelocityLabel_preReloc;
  const VarLabel* pFluidAccelerationLabel_preReloc;
  const VarLabel* pSolidMassLabel_preReloc;
  const VarLabel* pPorePressureLabel_preReloc;
  const VarLabel* pPorePressureFilterLabel_preReloc;
  const VarLabel* gFluidMassBarLabel_preReloc;
  const VarLabel* gFluidMassLabel_preReloc;
  const VarLabel* gFluidVelocityLabel_preReloc;
  const VarLabel* gFluidVelocityStarLabel_preReloc;
  const VarLabel* gFluidAccelerationLabel_preReloc;

  // MPM Hydrostatic BC label
  const VarLabel* boundaryPointsPerCellLabel;
};
} // End namespace Uintah

#endif
