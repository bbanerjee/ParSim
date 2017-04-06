#include <Core/Geometry/Cylinder.h>
#include <Core/Math/ran.h>
#include <iostream>

namespace dem {

void
Cylinder::print() const
{
  std::cout << "radius=" << radius << std::endl;
  std::cout << "height=" << height << std::endl;
  std::cout << "center=";
  center.print(std::cout);
}

Vec
Cylinder::randomPoint() const
{
  REAL rand1 = ran(&idum);
  REAL rand2 = ran(&idum);
  REAL rand3 = ran(&idum);
  REAL z = (center.z() + height / 2) * rand1 +
           (center.z() - height / 2) * (1 - rand1);
  REAL theta = 2 * Pi * rand2;
  REAL r = radius * rand3;
  REAL x = center.x() + r * cos(theta);
  REAL y = center.y() + r * sin(theta);
  return Vec(x, y, z);
}
}
