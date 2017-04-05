#include <Core/Const/const.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/ran.h>
#include <iostream>

namespace dem {

void
Box::print() const
{
  std::cout << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '
            << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
}

Vec
Box::randomPoint() const
{
  REAL rand1 = ran(&idum);
  REAL rand2 = ran(&idum);
  REAL rand3 = ran(&idum);
  REAL x = rand1 * (center.getX() - dimx / 2) +
           (1 - rand1) * (center.getX() + dimx / 2);
  REAL y = rand2 * (center.getY() - dimy / 2) +
           (1 - rand2) * (center.getY() + dimy / 2);
  REAL z = rand3 * (center.getZ() - dimz / 2) +
           (1 - rand3) * (center.getZ() + dimz / 2);
  return Vec(x, y, z);
}

std::ostream& operator<<(std::ostream& os, const Box& b) {
  os << "dimx = " << b.dimx << " dimy = " << b.dimy << " dimz = " << b.dimz
     << " center = " << b.center << " v1 = " << b.v1 << " v2 = " << b.v2 << "\n";
  return os;
}

} // namespace dem
