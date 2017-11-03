#include <Core/Geometry/Cylinder.h>
#include <Core/Math/ran.h>
#include <iostream>

namespace dem {

void
Cylinder::print() const
{
  std::cout << "radius=" << d_radius << std::endl;
  std::cout << "height=" << d_height << std::endl;
  std::cout << "center=";
  d_center.print(std::cout);
}

Vec
Cylinder::randomPoint() const
{
  REAL rand1 = ran(&idum);
  REAL rand2 = ran(&idum);
  REAL rand3 = ran(&idum);
  REAL z = (d_center.z() + d_height / 2) * rand1 + 
    (d_center.z() - d_height / 2) * (1 - rand1);
  REAL theta = 2 * Pi * rand2;
  REAL r = d_radius * rand3;
  REAL x = d_center.x() + r * cos(theta);
  REAL y = d_center.y() + r * sin(theta);
  return Vec(x, y, z);
}
}
