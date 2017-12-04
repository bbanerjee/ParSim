#ifndef __ELLIP3D_DEM_BOUNDARY_CONDITIONS_H__
#define __ELLIP3D_DEM_BOUNDARY_CONDITIONS_H__

#include <DiscreteElements/DEMPeriodicParticleBC.h>
#include <DiscreteElements/DEMBoundaryConditionUtils.h>
#include <DiscreteElements/DEMContainers.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
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

    Displacement getDisplacement(double time) const {
      return d_particleDispBC.getBCValue(time);
    }

    DeformationGradient getDeformationGradient(double time) const {
      return d_particleDefGradBC.getBCValue(time);
    }

    AxisymmetricStrain getAxisymmetricStrain(double time) const {
      return d_particleAxiStrainBC.getBCValue(time);
    }

    DisplacementRate getDisplacementRate(double time) const {
      return d_particleDispBC.getBCRate(time);
    }

    DeformationGradientRate getDeformationGradientRate(double time) const {
      return d_particleDefGradBC.getBCRate(time);
    }

    AxisymmetricStrainRate getAxisymmetricStrainRate(double time) const {
      return d_particleAxiStrainBC.getBCRate(time);
    }

    std::pair<Displacement, DisplacementRate>
    getDisplacementBC(double time) const
    {
      return std::make_pair(getDisplacement(time),
                            getDisplacementRate(time));
    }

    std::pair<DeformationGradient, DeformationGradientRate>
    getDeformationGradientBC(double time) const
    {
      return std::make_pair(getDeformationGradient(time),
                            getDeformationGradientRate(time));
    }

    std::pair<AxisymmetricStrain, AxisymmetricStrainRate>
    getAxisymmetricStrainBC(double time) const
    {
      return std::make_pair(getAxisymmetricStrain(time),
                            getAxisymmetricStrainRate(time));
    }

    void applyParticleBC(double time, 
                         const Box& spatialDomain,
                         OrientedBox& modifiableDomain,
                         DEMParticlePArray& particles);
    void applyDisplacementBC(double time, 
                             const Box& spatialDomain,
                             OrientedBox& modifiableDomain,
                             DEMParticlePArray& particles);
    void applyDeformationGradientBC(double time, 
                                    const Box& spatialDomain,
                                    OrientedBox& modifiableDomain,
                                    DEMParticlePArray& particles);
    void applyAxisymmetricStrainBC(double time, 
                                   const Box& spatialDomain,
                                   OrientedBox& modifiableDomain,
                                   DEMParticlePArray& particles);
  };

} // end namespace dem

#endif //ELLIP3D_DEM_BOUNDARY_CONDITIONS_H