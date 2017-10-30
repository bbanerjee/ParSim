#include <Core/Const/Constants.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/ran.h>
#include <iostream>

namespace dem {

void
Box::print() const
{
  std::cout << d_v1.x() << ' ' << d_v1.y() << ' ' << d_v1.z() << ' ' << d_v2.x() << ' '
            << d_v2.y() << ' ' << d_v2.z() << std::endl;
}

Vec
Box::randomPoint() const
{
  REAL rand1 = ran(&idum);
  REAL rand2 = ran(&idum);
  REAL rand3 = ran(&idum);
  REAL x =
    rand1 * (d_center.x() - d_dimx / 2) + (1 - rand1) * (d_center.x() + d_dimx / 2);
  REAL y =
    rand2 * (d_center.y() - d_dimy / 2) + (1 - rand2) * (d_center.y() + d_dimy / 2);
  REAL z =
    rand3 * (d_center.z() - d_dimz / 2) + (1 - rand3) * (d_center.z() + d_dimz / 2);
  return Vec(x, y, z);
}

// Special test:  Point is "inside" only if (x-xmin) >= -eps
// and (x-xmax) < -eps.
// Needed because of restrictive assumptions in the DEM code
bool 
Box::inside(const Vec& pt, REAL tol) const
{
  if (pt.x() - d_v1.x() >= -tol && pt.x() - d_v2.x() < -tol &&
      pt.y() - d_v1.y() >= -tol && pt.y() - d_v2.y() < -tol &&
      pt.z() - d_v1.z() >= -tol && pt.z() - d_v2.z() < -tol) {
    return true;
  }
  return false;
}

bool
Box::outside(const Vec& pt) const
{
  if ((pt.x() < d_v1.x() || pt.x() > d_v2.x()) ||
      (pt.y() < d_v1.y() || pt.y() > d_v2.y()) ||
      (pt.z() < d_v1.z() || pt.z() > d_v2.z())) {
    return true;
  }
  return false;
}

std::ostream&
operator<<(std::ostream& os, const Box& b)
{
  os << "d_dimx = " << b.d_dimx << " d_dimy = " << b.d_dimy << " d_dimz = " << b.d_dimz
     << " d_center = " << b.d_center << " d_v1 = " << b.d_v1 << " d_v2 = " << b.d_v2
     << "\n";
  return os;
}

} // namespace dem
