#ifndef CORE_GEOMETRY_PLANE_H
#define CORE_GEOMETRY_PLANE_H

#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>

namespace dem {

/////////////////////////////////////
class Plane
{
public:
  Vec direction;
  Vec point;

public:
  Plane()
    : direction(0)
    , point(0)
  {
  }

  Plane(Vec dir, Vec pt)
    : direction(std::move(dir))
    , point(std::move(pt))
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << direction.x() << std::setw(OWID) << direction.y()
       << std::setw(OWID) << direction.z() << std::setw(OWID) << point.x()
       << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
       << std::endl;
  }

  Vec getDirection() const { return direction; }
  Vec getPosition() const { return point; }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& direction;
    ar& point;
  }
};

} // end namespace dem

#endif
