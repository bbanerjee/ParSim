#ifndef __VAANGO_PERIDYNAMICSLABEL_H__
#define __VAANGO_PERIDYNAMICSLABEL_H__


#include <vector>

namespace Uintah 
{
  class VarLabel;
}

namespace Vaango 
{
  class PeridynamicsLabel {

    public:

      PeridynamicsLabel();
      ~PeridynamicsLabel();

      const Uintah::VarLabel* delTLabel;
      
      // non PermanentParticleState
      const Uintah::VarLabel* pPressureLabel;
      const Uintah::VarLabel* pVolumeDeformedLabel;

      // PermanentParticleState
      // Velocity, displacement, and deformation gradient
      const Uintah::VarLabel* pStressLabel;
      const Uintah::VarLabel* pStressLabel_preReloc;
      const Uintah::VarLabel* pVolumeLabel;
      const Uintah::VarLabel* pVolumeLabel_preReloc;
      const Uintah::VarLabel* pMassLabel;
      const Uintah::VarLabel* pMassLabel_preReloc;
      const Uintah::VarLabel* pVelocityLabel;
      const Uintah::VarLabel* pVelocityLabel_preReloc;
      const Uintah::VarLabel* pVelocityStarLabel;
      const Uintah::VarLabel* pAccelerationLabel;
      const Uintah::VarLabel* pAccelerationLabel_preReloc;
      const Uintah::VarLabel* pInternalForceLabel;
      const Uintah::VarLabel* pInternalForceLabel_preReloc;
      const Uintah::VarLabel* pExternalForceLabel;
      const Uintah::VarLabel* pExternalForceLabel_preReloc;

      const Uintah::VarLabel* pParticleIDLabel;
      const Uintah::VarLabel* pParticleIDLabel_preReloc;
      const Uintah::VarLabel* pSizeLabel;
      const Uintah::VarLabel* pSizeLabel_preReloc;
      const Uintah::VarLabel* pSurfLabel;
      const Uintah::VarLabel* pSurfLabel_preReloc;

      const Uintah::VarLabel* gMassLabel;
      const Uintah::VarLabel* gAccelerationLabel;
      const Uintah::VarLabel* gVelocityLabel;
      const Uintah::VarLabel* gVelocityStarLabel;
      const Uintah::VarLabel* gExternalForceLabel;
      const Uintah::VarLabel* gInternalForceLabel;
      const Uintah::VarLabel* gContactLabel;
      const Uintah::VarLabel* gNormTractionLabel;
      const Uintah::VarLabel* gSurfNormLabel;
      const Uintah::VarLabel* gStressLabel;
      const Uintah::VarLabel* gVolumeLabel;
      
      const Uintah::VarLabel* StrainEnergyLabel;
      const Uintah::VarLabel* AccStrainEnergyLabel;
      const Uintah::VarLabel* KineticEnergyLabel;
      const Uintah::VarLabel* CenterOfMassPositionLabel;
      const Uintah::VarLabel* TotalMomentumLabel;

      const Uintah::VarLabel* partCountLabel;
      const Uintah::VarLabel* pCellNAPIDLabel;

      // MPM Physical BC labels (permanent particle state)
      const Uintah::VarLabel* particlesPerLoadCurveLabel;
      const Uintah::VarLabel* pLoadCurveIDLabel;
      const Uintah::VarLabel* pLoadCurveIDLabel_preReloc;

      // These labels are for Peridynamics
      const Uintah::VarLabel* pPositionLabel;               // Particle position
      const Uintah::VarLabel* pPositionStarLabel;
      const Uintah::VarLabel* pPositionLabel_preReloc;
      const Uintah::VarLabel* pHorizonLabel;                // Store the horizon size for each particle
      const Uintah::VarLabel* pHorizonLabel_preReloc;    
      const Uintah::VarLabel* pDamageLabel;                 // Store the horizon size for each particle
      const Uintah::VarLabel* pDamageLabel_preReloc;    
      const Uintah::VarLabel* pDisplacementLabel;           // Displacement
      const Uintah::VarLabel* pDisplacementStarLabel;           
      const Uintah::VarLabel* pDisplacementLabel_preReloc;
      const Uintah::VarLabel* pDefGradLabel;                // State-based peridynamics deformation gradient
      const Uintah::VarLabel* pDefGradLabel_preReloc;
      const Uintah::VarLabel* pShapeTensorInvLabel;         // Inverse of the shape tensor
      const Uintah::VarLabel* pShapeTensorInvLabel_preReloc;
      const Uintah::VarLabel* pPK1StressLabel;              // First P-K stress
      const Uintah::VarLabel* pPK1StressLabel_preReloc;
      const Uintah::VarLabel* pNeighborListLabel;           // Store the neighbor list for each particle
      const Uintah::VarLabel* pNeighborListLabel_preReloc;  
      const Uintah::VarLabel* pNeighborConnLabel;           // Store the neighbor connectivity for each particle
      const Uintah::VarLabel* pNeighborConnLabel_preReloc; 
      const Uintah::VarLabel* pNeighborCountLabel;          // Store the neighbor count for each particle
      const Uintah::VarLabel* pNeighborCountLabel_preReloc; 
      const Uintah::VarLabel* pNeighborBondEnergyLabel;     // Store the neighbor strain energy for each particle
      const Uintah::VarLabel* pNeighborBondEnergyLabel_preReloc; 
      const Uintah::VarLabel* pNeighborBondForceLabel;     // Store the neighbor internal force for each particle
      const Uintah::VarLabel* pNeighborBondForceLabel_preReloc; 

    };

} // End namespace Vaango

#endif
