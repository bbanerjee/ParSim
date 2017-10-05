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
  Vec direc;
  Vec point;

public:
  Plane()
    : direc(0)
    , point(0)
  {
  }

  Plane(Vec dir, Vec pt)
    : direc(std::move(dir))
    , point(std::move(pt))
  {
  }

  void print(std::ostream& os)
  {
    os << std::setw(OWID) << direc.x() << std::setw(OWID) << direc.y()
       << std::setw(OWID) << direc.z() << std::setw(OWID) << point.x()
       << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
       << std::endl;
  }

  Vec getDirec() const { return direc; }
  Vec getPoint() const { return point; }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& direc;
    ar& point;
  }
};

} // end namespace dem

#endif
