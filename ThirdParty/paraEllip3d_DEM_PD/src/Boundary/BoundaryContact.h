#ifndef BOUNDARY_CONTACT_H
#define BOUNDARY_CONTACT_H

#include <Core/Const/const.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <boost/serialization/base_object.hpp>
#include <iostream>

namespace dem {

class Particle; // forward declaration, only use pointer to class Particle

/////////////////////////////////////
class BoundaryContact
{
public:
  Particle* ptcl;
  Vec point;
  Vec normal;
  Vec tangt;
  REAL penetr;

public:
  BoundaryContact()
    : ptcl(NULL)
    , point(0)
    , normal(0)
    , tangt(0)
    , penetr(0)
  {
  }

  BoundaryContact(Particle* p, Vec pt, Vec nm, Vec tg, REAL pntr)
    : ptcl(p)
    , point(pt)
    , normal(nm)
    , tangt(tg)
    , penetr(pntr)
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << point.getX() << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ() << std::setw(OWID) << normal.getX()
       << std::setw(OWID) << normal.getY() << std::setw(OWID) << normal.getZ()
       << std::setw(OWID) << tangt.getX() << std::setw(OWID) << tangt.getY()
       << std::setw(OWID) << tangt.getZ() << std::setw(OWID) << penetr
       << std::endl;
  }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& ptcl;
    ar& point;
    ar& normal;
    ar& tangt;
    ar& penetr;
  }
};

} // end namespace dem

#endif
