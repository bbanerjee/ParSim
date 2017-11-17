#ifndef __ELLIP3D_DEM_BOUNDARY_CONDITIONS_H__
#define __ELLIP3D_DEM_BOUNDARY_CONDITIONS_H__

#include <DiscreteElements/DEMPeriodicParticleBC.h>
#include <DiscreteElements/DEMBoundaryConditionUtils.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>
#include <string>

namespace dem {

  class DEMBoundaryConditions
  {
  private:

    BCUtils::DEM_DomainBCType d_domainBCType;
    BCUtils::DEM_ParticleBCType d_particleBCType;

    DEMPeriodicParticleBC<Displacement, 3> d_particleDispBC;
    DEMPeriodicParticleBC<DeformationGradient, 9> d_particleDefGradBC;
    DEMPeriodicParticleBC<AxisymmetricStrain, 4> d_particleAxiStrainBC;

  public:

    DEMBoundaryConditions() 
      : d_domainBCType(BCUtils::DEM_DomainBCType::FIXED)
      , d_particleBCType(BCUtils::DEM_ParticleBCType::NONE)
    {
    }

    bool read(const std::string& fileName);

    Displacement getDisplacement(double time) {
      return d_particleDispBC.getBCValue(time);
    }

    DeformationGradient getDeformationGradient(double time) {
      return d_particleDefGradBC.getBCValue(time);
    }

    AxisymmetricStrain getAxisymmetricStrain(double time) {
      return d_particleAxiStrainBC.getBCValue(time);
    }
  };

} // end namespace dem

#endif //ELLIP3D_DEM_BOUNDARY_CONDITIONS_H