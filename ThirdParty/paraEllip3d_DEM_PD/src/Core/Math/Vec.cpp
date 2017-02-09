#include <Core/Const/const.h>
#include <Core/Math/Vec.h>
#include <iostream>

namespace dem {

bool
Vec::operator==(const Vec v)
{
  return x == v.x && y == v.y && z == v.z;
}

bool
Vec::operator==(const REAL d)
{
  return x == d && y == d && z == d;
}

bool
Vec::operator!=(const Vec v)
{
  return x != v.x || y != v.y || z != v.z;
}

void
Vec::operator+=(Vec v)
{
  x += v.x;
  y += v.y;
  z += v.z;
}

void
Vec::operator-=(Vec v)
{
  x -= v.x;
  y -= v.y;
  z -= v.z;
}

void
Vec::operator*=(REAL d)
{
  x *= d;
  y *= d;
  z *= d;
}

void
Vec::operator/=(REAL d)
{
  x /= d;
  y /= d;
  z /= d;
}

Vec
Vec::operator+(Vec v) const
{
  return Vec(x + v.x, y + v.y, z + v.z);
}

Vec
Vec::operator-(Vec v) const
{
  return Vec(x - v.x, y - v.y, z - v.z);
}

Vec
Vec::operator%(Vec p) const
{
  return Vec(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
}

Vec Vec::operator*(REAL d) const
{
  return Vec(x * d, y * d, z * d);
}

REAL Vec::operator*(Vec p) const
{
  return (x * p.x + y * p.y + z * p.z);
}

void
Vec::print(std::ostream& ofs) const
{
  ofs << std::setw(OWID) << x << std::setw(OWID) << y << std::setw(OWID) << z;
}

Vec
Vec::fromString(const std::string& str)
{

  // Parse out the [num,num,num]
  std::string::size_type i1 = str.find("[");
  std::string::size_type i2 = str.find_first_of(",");
  std::string::size_type i3 = str.find_last_of(",");
  std::string::size_type i4 = str.find("]");

  std::string x_val(str, i1 + 1, i2 - i1 - 1);
  std::string y_val(str, i2 + 1, i3 - i2 - 1);
  std::string z_val(str, i3 + 1, i4 - i3 - 1);

  // Validate the values (**TODO** Add check for repeated ".")
  std::string validChars(" -+.0123456789eE");
  std::string::size_type xpos = x_val.find_first_not_of(validChars);
  std::string::size_type ypos = y_val.find_first_not_of(validChars);
  std::string::size_type zpos = z_val.find_first_not_of(validChars);

  if (xpos != std::string::npos || ypos != std::string::npos ||
      zpos != std::string::npos) {
    std::cout << "Unexpected non-real value in Vector."
              << " Converting to double if possible\n";
  }

  // Create Vec from string
  Vec result;
  result.set(std::stod(x_val), std::stod(y_val), std::stod(z_val));
  return result;
}

Vec operator*(REAL d, Vec v)
{
  return Vec(v.getX() * d, v.getY() * d, v.getZ() * d);
}

Vec
operator/(Vec v, REAL d)
{
  return Vec(v.getX() / d, v.getY() / d, v.getZ() / d);
}

REAL
vfabs(Vec v)
{
  REAL x = v.getX();
  REAL y = v.getY();
  REAL z = v.getZ();
  return sqrt(x * x + y * y + z * z);
}

Vec
vcos(Vec v)
{
  return Vec(cos(v.getX()), cos(v.getY()), cos(v.getZ()));
}

Vec
vacos(Vec v)
{
  return Vec(acos(v.getX()), acos(v.getY()), acos(v.getZ()));
}

Vec
operator-(Vec v)
{
  return -1.0 * v;
}

Vec
normalize(Vec v)
{
  REAL alf = vfabs(v);
  if (alf < EPS) // important, otherwise may cause numerical instability
    return v;
  return v / (vfabs(v));
}

// return what vector vec is rotated to by vector rot.
// note: that vec rotates along x, y, z axis by rot.x, rot.y, rot.z
// is equivalent to that vec rotates along vector rot by vfabs(rot)
Vec
rotateVec(Vec vec, Vec rot)
{
  REAL alf = vfabs(rot);
  if (alf < EPS) // important, otherwise may cause numerical instability
    return vec;

  Vec nx = rot / alf;
  Vec vx = (vec * nx) * nx;
  Vec vy = vec - vx;

  REAL theta = atan(vfabs(vy) / vfabs(vx));
#ifndef NDEBUG
  debugInf << "Vec.cpp: iter=" << iteration << " alf=" << alf
           << " theta=" << theta << std::endl;
#endif
  if (theta < EPS) // important, otherwise my cause numerical instability
    return vec;

  Vec ny = normalize(vy);
  Vec nz = normalize(nx % ny); // normalize for higher precision
  REAL radius = vfabs(vy);
  return radius * sin(alf) * nz + radius * cos(alf) * ny + vx;
}

} // namespace dem
