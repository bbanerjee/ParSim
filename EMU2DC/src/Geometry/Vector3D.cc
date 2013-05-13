#include <Geometry/Vector3D.h>
#include <cmath>

using namespace Emu2DC;

Vector3D::Vector3D(const Vector3D& vec)
{
  d_vec[0] = vec[0];
  d_vec[1] = vec[1];
  d_vec[2] = vec[2];
}

Vector3D::Vector3D(const Point3D& start, const Point3D& end)
{
 d_vec[0] = end.x() - start.x();
 d_vec[1] = end.y() - start.y();
 d_vec[2] = end.z() - start.z();
}

double 
Vector3D::length() const
{
  return std::sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
}

double 
Vector3D::lengthSq() const
{
  return (d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
}

double 
Vector3D::dot(const Vector3D& vec) const
{
  return d_vec[0]*vec[0] + d_vec[1]*vec[1] + d_vec[2]*vec[2];
}

Vector3D 
Vector3D::cross(const Vector3D& vec) const
{
  return Vector3D(
        d_vec[1]*vec[2]-d_vec[2]*vec[1],
        d_vec[2]*vec[0]-d_vec[0]*vec[2],
        d_vec[0]*vec[1]-d_vec[1]*vec[0]);
}

bool 
Vector3D::operator==(const Vector3D& vec) const
{
  return (vec.x() == d_vec[0] && vec.y() == d_vec[1] && vec.z() == d_vec[2]);
}

bool 
Vector3D::operator!=(const Vector3D& vec) const
{
  return (vec.x() != d_vec[0] && vec.y() != d_vec[1] && vec.z() != d_vec[2]);
}

Vector3D& 
Vector3D::operator=(const Vector3D& vec)
{
  d_vec[0] = vec.x();
  d_vec[1] = vec.y();
  d_vec[2] = vec.z();
  return *this;
}

Vector3D 
Vector3D::operator*(const double val) const
{
  return Vector3D(d_vec[0]*val, d_vec[1]*val, d_vec[2]*val);
}

Vector3D 
Vector3D::operator/(const double val) const
{
  return Vector3D(d_vec[0]/val, d_vec[1]/val, d_vec[2]/val);
}

Vector3D 
Vector3D::operator+(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]+vec.x(), d_vec[1]+vec.y(), d_vec[2]+vec.z());
}

Vector3D 
Vector3D::operator-(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]-vec.x(), d_vec[1]-vec.y(), d_vec[2]-vec.z());
}

Vector3D& 
Vector3D::operator*=(const double val)
{
  d_vec[0] *= val;
  d_vec[1] *= val;
  d_vec[2] *= val;
  return *this;
}

Vector3D& 
Vector3D::operator/=(const double val)
{
  d_vec[0] /= val;
  d_vec[1] /= val;
  d_vec[2] /= val;
  return *this;
}

Vector3D& 
Vector3D::operator+=(const Vector3D& vec)
{
  d_vec[0] += vec.x();
  d_vec[1] += vec.y();
  d_vec[2] += vec.z();
  return *this;
}

Vector3D& 
Vector3D::operator-=(const Vector3D& vec)
{
  d_vec[0] -= vec.x();
  d_vec[1] -= vec.y();
  d_vec[2] -= vec.z();
  return *this;
}

namespace Emu2DC {
  std::ostream& operator<<(std::ostream& out, const Vector3D& vec) 
  {
    out << "[" << vec.d_vec[0] << " " << vec.d_vec[1] << " " << vec.d_vec[2] << "]";
    return out;
  }
}


