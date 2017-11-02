#ifndef __ELLIP3D_CORE_GEOM_ELLIPSOID_H__
#define __ELLIP3D_CORE_GEOM_ELLIPSOID_H__

#include <Core/Math/Vec.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Geometry/FaceEdge.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class Matrix3;

class Ellipsoid
{

private:

  Vec d_center;
  std::array<Vec, 3> d_axes;
  std::array<REAL, 3> d_radii;

  void normalize_axes();

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_center;
    ar& d_axes;
    ar& d_radii;
  }

public:

  Ellipsoid(const Vec& center, const Vec& ax_a, const Vec& ax_b, const Vec& ax_c,
            REAL rad_a, REAL rad_b, REAL rad_c)
    : d_center(center)
    , d_axes({{ax_a, ax_b, ax_c}})
    , d_radii({{rad_a, rad_b, rad_c}})
  {
    normalize_axes();
  }

  const Vec& center() const {return d_center;}
  const Vec& axis(int idx) const {return d_axes.at(idx);}
  REAL radii(int idx) const {return d_radii.at(idx);}

  OrientedBox getOrientedBoundingBox() const;

  Matrix3 toUnitSphereTransformationMatrix() const;

  bool containsPoint(const Vec& pt) const;

  bool intersects(const OrientedBox& box) const;
  bool intersects(const Face& face) const;

  void rotate(REAL angle, const Vec& axis);
  void translate(const Vec& dist);

  friend std::ostream& operator<<(std::ostream& os, const Ellipsoid& ellipsoid);

};

} // namespace dem

#endif
