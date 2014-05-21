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
      const Uintah::VarLabel* pVelGradLabel;
      const Uintah::VarLabel* pVelGradLabel_preReloc;
      const Uintah::VarLabel* pDispGradLabel;
      const Uintah::VarLabel* pDispGradLabel_preReloc;
      const Uintah::VarLabel* pDefGradLabel;
      const Uintah::VarLabel* pDefGradLabel_preReloc;
      const Uintah::VarLabel* pShapeTensInvLabel;
      const Uintah::VarLabel* pShapeTensInvLabel_preReloc;

      const Uintah::VarLabel* pStressLabel;
      const Uintah::VarLabel* pStressLabel_preReloc;
      const Uintah::VarLabel* pVolumeLabel;
      const Uintah::VarLabel* pVolumeLabel_preReloc;
      const Uintah::VarLabel* pMassLabel;
      const Uintah::VarLabel* pMassLabel_preReloc;
      const Uintah::VarLabel* pVelocityLabel;
      const Uintah::VarLabel* pVelocityLabel_preReloc;
      const Uintah::VarLabel* pDispLabel;
      const Uintah::VarLabel* pDispLabel_preReloc;
      const Uintah::VarLabel* pVelocityStarLabel;
      const Uintah::VarLabel* pAccelerationLabel;
      const Uintah::VarLabel* pAccelerationLabel_preReloc;
      const Uintah::VarLabel* pInternalForceLabel;
      const Uintah::VarLabel* pInternalForceLabel_preReloc;
      const Uintah::VarLabel* pExternalForceLabel;
      const Uintah::VarLabel* pExternalForceLabel_preReloc;

      // These labels are for storing neighbor information
      const Uintah::VarLabel* pHorizonLabel;                // Store the horizon size for each particle
      const Uintah::VarLabel* pHorizonLabel_preReloc;    
      const Uintah::VarLabel* pDamageLabel;                 // Store the horizon size for each particle
      const Uintah::VarLabel* pDamageLabel_preReloc;    
      const Uintah::VarLabel* pNeighborListLabel;           // Store the neighbor list for each particle
      const Uintah::VarLabel* pNeighborListLabel_preReloc;  
      const Uintah::VarLabel* pNeighborConnLabel;           // Store the neighbor connectivity for each particle
      const Uintah::VarLabel* pNeighborConnLabel_preReloc; 
      const Uintah::VarLabel* pNeighborCountLabel;          // Store the neighbor count for each particle
      const Uintah::VarLabel* pNeighborCountLabel_preReloc; 
      const Uintah::VarLabel* pNeighborBondEnergyLabel;     // Store the neighbor strain energy for each particle
      const Uintah::VarLabel* pNeighborBondEnergyLabel_preReloc; 

      const Uintah::VarLabel* pParticleIDLabel;
      const Uintah::VarLabel* pParticleIDLabel_preReloc;
      const Uintah::VarLabel* pSizeLabel;
      const Uintah::VarLabel* pSizeLabel_preReloc;
      const Uintah::VarLabel* pPositionLabel;
      const Uintah::VarLabel* pPositionLabel_preReloc;
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

    };

} // End namespace Vaango

#endif
