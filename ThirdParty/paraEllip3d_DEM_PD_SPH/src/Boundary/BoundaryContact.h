#ifndef BOUNDARY_CONTACT_H
#define BOUNDARY_CONTACT_H

#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <boost/serialization/base_object.hpp>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

/////////////////////////////////////
class BoundaryContact
{
public:
  DEMParticle* bc_particle;
  Vec bc_point;
  Vec bc_normalForce;
  Vec bc_tangentForce;
  REAL bc_penetration;

public:
  BoundaryContact()
    : bc_particle(nullptr)
    , bc_point(0)
    , bc_normalForce(0)
    , bc_tangentForce(0)
    , bc_penetration(0)
  {
  }

  BoundaryContact(DEMParticle* p, Vec pt, Vec nm, Vec tg, REAL pntr)
    : bc_particle(p)
    , bc_point(std::move(pt))
    , bc_normalForce(std::move(nm))
    , bc_tangentForce(std::move(tg))
    , bc_penetration(pntr)
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << bc_point.x() << std::setw(OWID) << bc_point.y()
       << std::setw(OWID) << bc_point.z() << std::setw(OWID) << bc_normalForce.x()
       << std::setw(OWID) << bc_normalForce.y() << std::setw(OWID) << bc_normalForce.z()
       << std::setw(OWID) << bc_tangentForce.x() << std::setw(OWID) << bc_tangentForce.y()
       << std::setw(OWID) << bc_tangentForce.z() << std::setw(OWID) << bc_penetration
       << std::endl;
  }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& bc_particle;
    ar& bc_point;
    ar& bc_normalForce;
    ar& bc_tangentForce;
    ar& bc_penetration;
  }
};

} // end namespace dem

#endif
