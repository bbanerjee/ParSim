#include <Core/Const/Constants.h>
#include <Core/Geometry/Ellipsoid.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Math/Matrix3.h>
#include <iostream>

namespace dem {

void 
Ellipsoid::normalize_axes() {
  d_axis_a.normalizeInPlace();
  d_axis_b.normalizeInPlace();
  d_axis_c.normalizeInPlace();
}

OrientedBox
Ellipsoid::getOrientedBoundingBox() const
{
  OrientedBox box(d_center, d_axis_a, d_axis_b, d_axis_c,
                  d_radius_a, d_radius_b, r_radius_c);
  return box;
}

Matrix3 
Ellipsoid::toUnitSphereTransformationMatrix() const
{
  Matrix3 mat;
  mat(0,0) = d_axis_a[0]/d_radius_a;
  mat(0,1) = d_axis_a[1]/d_radius_a;
  mat(0,2) = d_axis_a[2]/d_radius_a;
  mat(1,0) = d_axis_b[0]/d_radius_b;
  mat(1,1) = d_axis_b[1]/d_radius_b;
  mat(1,2) = d_axis_b[2]/d_radius_b;
  mat(2,0) = d_axis_c[0]/d_radius_c;
  mat(2,1) = d_axis_c[1]/d_radius_c;
  mat(2,2) = d_axis_c[2]/d_radius_c;
}

bool 
Ellipsoid::containsPoint(const Vec& point) const
{
  Matrix3 N = toUnitSphereTransformationMatrix();
  Vec pp = N*(point - d_center) ;
  if (pp.length() > 1.0) {
    return false;
  }
  return true;
}

std::ostream&
operator<<(std::ostream& os, const Ellipsoid& b)
{
  os << " center = " << b.d_center 
     << " axis_a = " << b.d_axis_a 
     << " axis_b = " << b.d_axis_b 
     << " axis_c = " << b.d_axis_c 
     << " radius_a = " << b.d_radius_a 
     << " radius_b = " << b.d_radius_b 
     << " radius_c = " << b.d_radius_c 
     << "\n";
  return os;
}

} // namespace dem
