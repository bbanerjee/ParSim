#ifndef __VAANGO_PERIDYNAMICSLABEL_H__
#define __VAANGO_PERIDYNAMICSLABEL_H__

#include <vector>

namespace Vaango 
{

  class Uintah::VarLabel;

  class PeridynamicsLabel {

    public:

      PeridynamicsLabel();
      ~PeridynamicsLabel();

      const VarLabel* delTLabel;
      
      // non PermanentParticleState
      const VarLabel* pPressureLabel;
      const VarLabel* pVolumeDeformedLabel;

      // PermanentParticleState
      // Velocity, displacement, and deformation gradient
      const VarLabel* pVelGradLabel;
      const VarLabel* pVelGradLabel_preReloc;
      const VarLabel* pDispGradLabel;
      const VarLabel* pDispGradLabel_preReloc;
      const VarLabel* pDefGradLabel;
      const VarLabel* pDefGradLabel_preReloc;

      const VarLabel* pStressLabel;
      const VarLabel* pStressLabel_preReloc;
      const VarLabel* pVolumeLabel;
      const VarLabel* pVolumeLabel_preReloc;
      const VarLabel* pMassLabel;
      const VarLabel* pMassLabel_preReloc;
      const VarLabel* pVelocityLabel;
      const VarLabel* pVelocityLabel_preReloc;
      const VarLabel* pExternalForceLabel;
      const VarLabel* pExternalForceLabel_preReloc;

      const VarLabel* pXLabel;
      const VarLabel* pXLabel_preReloc;
      const VarLabel* pSurfLabel;
      const VarLabel* pSurfLabel_preReloc;

      const VarLabel* pParticleIDLabel;
      const VarLabel* pParticleIDLabel_preReloc;

      const VarLabel* gMassLabel;
      const VarLabel* gAccelerationLabel;
      const VarLabel* gVelocityLabel;
      const VarLabel* gExternalForceLabel;
      const VarLabel* gInternalForceLabel;
      const VarLabel* gContactLabel;
      const VarLabel* gNormTractionLabel;
      const VarLabel* gSurfNormLabel;
      const VarLabel* gStressLabel;
      const VarLabel* gVolumeLabel;
      
      const VarLabel* StrainEnergyLabel;
      const VarLabel* AccStrainEnergyLabel;
      const VarLabel* KineticEnergyLabel;
      const VarLabel* CenterOfMassPositionLabel;
      const VarLabel* TotalMomentumLabel;
    };

} // End namespace Vaango

#endif
