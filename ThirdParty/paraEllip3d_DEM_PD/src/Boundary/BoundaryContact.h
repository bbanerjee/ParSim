#ifndef BOUNDARY_CONTACT_H
#define BOUNDARY_CONTACT_H

#include <Core/Const/const.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <boost/serialization/base_object.hpp>
#include <iostream>
#include <utility>

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
    : ptcl(nullptr)
    , point(0)
    , normal(0)
    , tangt(0)
    , penetr(0)
  {
  }

  BoundaryContact(Particle* p, Vec pt, Vec nm, Vec tg, REAL pntr)
    : ptcl(p)
    , point(std::move(pt))
    , normal(std::move(nm))
    , tangt(std::move(tg))
    , penetr(pntr)
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << point.x() << std::setw(OWID) << point.y()
       << std::setw(OWID) << point.z() << std::setw(OWID) << normal.x()
       << std::setw(OWID) << normal.y() << std::setw(OWID) << normal.z()
       << std::setw(OWID) << tangt.x() << std::setw(OWID) << tangt.y()
       << std::setw(OWID) << tangt.z() << std::setw(OWID) << penetr
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
