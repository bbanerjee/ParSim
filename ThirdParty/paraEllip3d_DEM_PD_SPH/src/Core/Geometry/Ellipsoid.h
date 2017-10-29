#ifndef __ELLIP3D_CORE_GEOM_ELLIPSOID_H__
#define __ELLIP3D_CORE_GEOM_ELLIPSOID_H__

#include <Core/Math/Vec.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class Matrix3;

class Ellipsoid
{

private:

  Vec d_center;
  Vec d_axis_a;
  Vec d_axis_b;
  Vec d_axis_c;
  REAL d_radius_a;
  REAL d_radius_b;
  REAL d_radius_c;

  void normalize_axes();

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_center;
    ar& d_axis_a;
    ar& d_axis_b;
    ar& d_axis_c;
    ar& d_radius_a;
    ar& d_radius_b;
    ar& d_radius_c;
  }

public:

  Ellipsoid(const Vec& center, const Vec& ax_a, const Vec& ax_b, const Vec& ax_c,
            REAL rad_a, REAL rad_b, REAL rad_c)
    : d_center(center)
    , d_axis_a(ax_a)
    , d_axis_b(ax_b)
    , d_axis_c(ax_c)
    , d_radius_a(rad_a)
    , d_radius_b(rad_b)
    , d_radius_c(rad_c)
  {
    normalize_axes();
  }

  OrientedBox getOrientedBoundingBox() const;

  Matrix3 toSphereTransformationMatrix() const;

  bool containsPoint(const Vec& pt) const;

  friend std::ostream& operator<<(std::ostream& os, const Ellipsoid& box);

};

} // namespace dem

#endif
