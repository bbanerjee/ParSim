#include <Core/Const/const.h>
#include <Core/Math/Vec.h>
#include <Core/Math/IntVec.h>
#include <iostream>

namespace dem {

bool
Vec::operator==(const Vec v)
{
  return d_data[0] == v.d_data[0] && d_data[1] == v.d_data[1] && d_data[2] == v.d_data[2];
}

bool
Vec::operator==(const REAL val)
{
  return d_data[0] == val && d_data[1] == val && d_data[2] == val;
}

bool
Vec::operator!=(const Vec v)
{
  return d_data[0] != v.d_data[0] || d_data[1] != v.d_data[1] || d_data[2] != v.d_data[2];
}

void
Vec::operator+=(Vec v)
{
  d_data[0] += v.d_data[0];
  d_data[1] += v.d_data[1];
  d_data[2] += v.d_data[2];
}

void
Vec::operator-=(Vec v)
{
  d_data[0] -= v.d_data[0];
  d_data[1] -= v.d_data[1];
  d_data[2] -= v.d_data[2];
}

void
Vec::operator*=(REAL val)
{
  d_data[0] *= val;
  d_data[1] *= val;
  d_data[2] *= val;
}

void
Vec::operator/=(REAL val)
{
  d_data[0] /= val;
  d_data[1] /= val;
  d_data[2] /= val;
}

Vec
Vec::operator+(Vec v) const
{
  return Vec(d_data[0] + v.d_data[0], d_data[1] + v.d_data[1], d_data[2] + v.d_data[2]);
}

Vec
Vec::operator-(Vec v) const
{
  return Vec(d_data[0] - v.d_data[0], d_data[1] - v.d_data[1], d_data[2] - v.d_data[2]);
}

Vec
Vec::operator%(Vec p) const
{
  return Vec(d_data[1] * p.d_data[2] - d_data[2] * p.d_data[1], d_data[2] * p.d_data[0] - d_data[0] * p.d_data[2],
             d_data[0] * p.d_data[1] - d_data[1] * p.d_data[0]);
}

Vec 
Vec::operator*(REAL d) const
{
  return Vec(d_data[0] * d, d_data[1] * d, d_data[2] * d);
}

REAL 
Vec::operator*(Vec p) const
{
  return (d_data[0] * p.d_data[0] + d_data[1] * p.d_data[1] + d_data[2] * p.d_data[2]);
}

REAL
Vec::lengthSq() const
{
  return d_data[0] * d_data[0] + d_data[1] * d_data[1] + d_data[2] * d_data[2];
}

REAL
Vec::length() const
{
  return std::sqrt(lengthSq());
}

// Divides a vector by an int vector and produces a new vector
Vec 
Vec::operator/(const IntVec& intVec) const
{
  REAL x = d_data[0]/static_cast<REAL>(intVec.x());
  REAL y = d_data[1]/static_cast<REAL>(intVec.y());
  REAL z = d_data[2]/static_cast<REAL>(intVec.z());
  return Vec(x,y,z);
}

// Multiplies a vector by an int vector and produces a new vector
Vec 
Vec::operator*(const IntVec& intVec) const
{
  REAL x = d_data[0]*static_cast<REAL>(intVec.x());
  REAL y = d_data[1]*static_cast<REAL>(intVec.y());
  REAL z = d_data[2]*static_cast<REAL>(intVec.z());
  return Vec(x,y,z);
}

void
Vec::print(std::ostream& ofs) const
{
  ofs << std::setw(OWID) << d_data[0] << std::setw(OWID) << d_data[1] << std::setw(OWID)
      << d_data[2];
}

Vec
Vec::fromString(const std::string& str)
{

  // Parse out the [num,num,num] formatted string
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
    std::cerr << "Unexpected non-real value in Vector."
              << " Converting to double if possible\n";
  }

  // Create Vec from string
  Vec result;
  result.set(std::stod(x_val), std::stod(y_val), std::stod(z_val));
  return result;
}

Vec operator*(REAL d, Vec v)
{
  return Vec(v.x() * d, v.y() * d, v.z() * d);
}

Vec
operator/(Vec v, REAL d)
{
  return Vec(v.x() / d, v.y() / d, v.z() / d);
}

REAL
vfabs(Vec v)
{
  REAL x = v.x();
  REAL y = v.y();
  REAL z = v.z();
  return sqrt(x * x + y * y + z * z);
}

Vec
vcos(Vec v)
{
  return Vec(cos(v.x()), cos(v.y()), cos(v.z()));
}

Vec
vacos(Vec v)
{
  return Vec(acos(v.x()), acos(v.y()), acos(v.z()));
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

Vec& 
operator/=(Vec& realVec, const IntVec& intVec) {
  realVec.d_data[0] /= static_cast<REAL>(intVec.x());
  realVec.d_data[1] /= static_cast<REAL>(intVec.y());
  realVec.d_data[2] /= static_cast<REAL>(intVec.z());
  return realVec;
}

Vec& 
operator*=(Vec& realVec, const IntVec& intVec) {
  realVec.d_data[0] *= static_cast<REAL>(intVec.x());
  realVec.d_data[1] *= static_cast<REAL>(intVec.y());
  realVec.d_data[2] *= static_cast<REAL>(intVec.z());
  return realVec;
}

std::ostream&
operator<<(std::ostream& os, const Vec& v)
{
  os << ' ' << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
  return os;
}

} // namespace dem
