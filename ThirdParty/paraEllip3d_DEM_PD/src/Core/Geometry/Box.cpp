#include <Core/Const/const.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/ran.h>
#include <iostream>

namespace dem {

void
Box::print() const
{
  std::cout << v1.x() << ' ' << v1.y() << ' ' << v1.z() << ' ' << v2.x() << ' '
            << v2.y() << ' ' << v2.z() << std::endl;
}

Vec
Box::randomPoint() const
{
  REAL rand1 = ran(&idum);
  REAL rand2 = ran(&idum);
  REAL rand3 = ran(&idum);
  REAL x =
    rand1 * (center.x() - dimx / 2) + (1 - rand1) * (center.x() + dimx / 2);
  REAL y =
    rand2 * (center.y() - dimy / 2) + (1 - rand2) * (center.y() + dimy / 2);
  REAL z =
    rand3 * (center.z() - dimz / 2) + (1 - rand3) * (center.z() + dimz / 2);
  return Vec(x, y, z);
}

std::ostream&
operator<<(std::ostream& os, const Box& b)
{
  os << "dimx = " << b.dimx << " dimy = " << b.dimy << " dimz = " << b.dimz
     << " center = " << b.center << " v1 = " << b.v1 << " v2 = " << b.v2
     << "\n";
  return os;
}

} // namespace dem
