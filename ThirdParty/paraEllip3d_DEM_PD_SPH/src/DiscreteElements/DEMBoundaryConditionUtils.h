#ifndef __ELLIP3D_DEM_BOUNDARY_CONDITION_UTILS_H__
#define __ELLIP3D_DEM_BOUNDARY_CONDITION_UTILS_H__

#include <Core/Math/Vec.h>

#include <vector>
#include <string>

namespace dem {

  namespace BCUtils {

  enum class DEM_DomainBCType {
    FIXED,
    DISPLACEMENT,
    VELOCITY,
    TRACTION,
    PERIODIC,
    PERIODIC_PARTICLE
  };

  enum class DEM_ParticleBCType {
    NONE,
    FIXED,
    DISPLACEMENT,
    TRACTION,
    DEFORMATION_GRADIENT,
    AXISYMMETRIC_STRAIN
  };

  inline
  DEM_DomainBCType getDEMDomainBCType(const std::string& type)
  {
    if (type == "displacement") return DEM_DomainBCType::DISPLACEMENT;
    else if (type == "velocity") return DEM_DomainBCType::VELOCITY;
    else if (type == "traction") return DEM_DomainBCType::TRACTION;
    else if (type == "periodic") return DEM_DomainBCType::PERIODIC;
    else if (type == "periodic_particle") {
      return DEM_DomainBCType::PERIODIC_PARTICLE;
    }
    return DEM_DomainBCType::FIXED;
  }

  inline
  DEM_ParticleBCType getDEMParticleBCType(const std::string& type)
  {
    if (type == "fixed") return DEM_ParticleBCType::FIXED;
    else if (type == "displacement") return DEM_ParticleBCType::DISPLACEMENT;
    else if (type == "traction") return DEM_ParticleBCType::TRACTION;
    else if (type == "deformation_gradient") {
      return DEM_ParticleBCType::DEFORMATION_GRADIENT;
    } else if (type == "axisymmetric_strain") {
      return DEM_ParticleBCType::AXISYMMETRIC_STRAIN;
    }
    return DEM_ParticleBCType::NONE;
  }

  } // end namespace BCUtils
}
#endif