/*
 * IntIntVector3D.cc
 *
 *  Created on: 22/10/2013
 *      Author: banerjee
 */

#include <Geometry/IntVector3D.h>
#include <cmath>

using namespace BrMPM;

IntVector3D::IntVector3D(const IntVector3D& vec)
{
  d_vec[0] = vec[0];
  d_vec[1] = vec[1];
  d_vec[2] = vec[2];
}

bool
IntVector3D::operator==(const IntVector3D& vec) const
{
  return (vec.x() == d_vec[0] && vec.y() == d_vec[1] && vec.z() == d_vec[2]);
}

bool
IntVector3D::operator!=(const IntVector3D& vec) const
{
  return (vec.x() != d_vec[0] && vec.y() != d_vec[1] && vec.z() != d_vec[2]);
}

void
IntVector3D::operator=(const IntVector3D& vec)
{
  d_vec[0] = vec.x();
  d_vec[1] = vec.y();
  d_vec[2] = vec.z();
}

IntVector3D
IntVector3D::operator*(const int val) const
{
  return IntVector3D(d_vec[0]*val, d_vec[1]*val, d_vec[2]*val);
}

IntVector3D
IntVector3D::operator*(const IntVector3D& vec) const
{
  return IntVector3D(d_vec[0]*vec.x(), d_vec[1]*vec.y(), d_vec[2]*vec.z());
}

IntVector3D
IntVector3D::operator+(const IntVector3D& vec) const
{
  return IntVector3D(d_vec[0]+vec.x(), d_vec[1]+vec.y(), d_vec[2]+vec.z());
}

IntVector3D
IntVector3D::operator-(const IntVector3D& vec) const
{
  return IntVector3D(d_vec[0]-vec.x(), d_vec[1]-vec.y(), d_vec[2]-vec.z());
}

IntVector3D&
IntVector3D::operator*=(const int val)
{
  d_vec[0] *= val;
  d_vec[1] *= val;
  d_vec[2] *= val;
  return *this;
}

IntVector3D&
IntVector3D::operator+=(const IntVector3D& vec)
{
  d_vec[0] += vec.x();
  d_vec[1] += vec.y();
  d_vec[2] += vec.z();
  return *this;
}

IntVector3D&
IntVector3D::operator-=(const IntVector3D& vec)
{
  d_vec[0] -= vec.x();
  d_vec[1] -= vec.y();
  d_vec[2] -= vec.z();
  return *this;
}

int
IntVector3D::max() const
{
  return std::max(std::max(d_vec[0],d_vec[1]), d_vec[2]);
}

int
IntVector3D::min() const
{
  return std::min(std::min(d_vec[0],d_vec[1]), d_vec[2]);
}

namespace BrMPM {

  IntVector3D
  min(const IntVector3D& v1, const IntVector3D& v2)
  {
    return IntVector3D(std::min(v1.x(),v2.x()), std::min(v1.y(),v2.y()), std::min(v1.z(), v2.z()));
  }

  IntVector3D
  max(const IntVector3D& v1, const IntVector3D& v2)
  {
    return IntVector3D(std::max(v1.x(),v2.x()), std::max(v1.y(),v2.y()), std::max(v1.z(), v2.z()));
  }

  std::ostream& operator<<(std::ostream& out, const IntVector3D& vec)
  {
    out << "[" << vec.d_vec[0] << " " << vec.d_vec[1] << " " << vec.d_vec[2] << "]";
    return out;
  }
}

