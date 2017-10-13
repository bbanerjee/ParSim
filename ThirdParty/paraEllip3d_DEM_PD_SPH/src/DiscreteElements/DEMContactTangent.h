#ifndef ELLIP3D_DEM_CONTACT_TANGENT_H
#define ELLIP3D_DEM_CONTACT_TANGENT_H

#include <Core/Math/Vec.h>

namespace dem {

class DEMContactTangent
{

public:
  std::size_t ct_particle1;
  std::size_t ct_particle2;
  Vec ct_tangentForce;
  Vec ct_tangentDisplacement;
  bool ct_tangentLoadingActive;
  Vec ct_tangentDispStart;
  REAL ct_tangentForcePeak;
  bool ct_tangentSlidingActive;

  DEMContactTangent()
    : ct_particle1(0)
    , ct_particle2(0)
    , ct_tangentForce(0)
    , ct_tangentDisplacement(0)
    , ct_tangentLoadingActive(0)
    , ct_tangentDispStart(0)
    , ct_tangentForcePeak(0)
    , ct_tangentSlidingActive(false)
  {
  }

  DEMContactTangent(std::size_t _ptcl1, std::size_t _ptcl2, Vec _tf, Vec _td, bool _tl,
             Vec _tds, REAL _tp, bool _ts)
    : ct_particle1(_ptcl1)
    , ct_particle2(_ptcl2)
    , ct_tangentForce(_tf)
    , ct_tangentDisplacement(_td)
    , ct_tangentLoadingActive(_tl)
    , ct_tangentDispStart(_tds)
    , ct_tangentForcePeak(_tp)
    , ct_tangentSlidingActive(_ts)
  {
  }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& ct_particle1;
    ar& ct_particle2;
    ar& ct_tangentForce;
    ar& ct_tangentDisplacement;
    ar& ct_tangentLoadingActive;
    ar& ct_tangentDispStart;
    ar& ct_tangentForcePeak;
    ar& ct_tangentSlidingActive;
  }
};

} // end namespace dem
#endif