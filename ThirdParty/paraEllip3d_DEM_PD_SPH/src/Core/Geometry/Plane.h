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
  Vec normal;
  Vec point;

public:
  Plane()
    : normal(0)
    , point(0)
  {
  }

  Plane(Vec dir, Vec pt)
    : normal(std::move(dir))
    , point(std::move(pt))
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << normal.x() << std::setw(OWID) << normal.y()
       << std::setw(OWID) << normal.z() << std::setw(OWID) << point.x()
       << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
       << std::endl;
  }

  Vec getDirection() const { return normal; }
  Vec getPosition() const { return point; }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& normal;
    ar& point;
  }
};

} // end namespace dem

#endif
