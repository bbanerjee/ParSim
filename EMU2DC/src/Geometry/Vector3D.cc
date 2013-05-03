#include <Geometry/Vector3D.h>
#include <cmath>

using namespace Emu2DC;

Vector3D::Vector3D(const Vector3D& vec)
{
  d_x = vec.d_x;
  d_y = vec.d_y;
  d_z = vec.d_z;
}

Vector3D::Vector3D(const Point3D& start, const Point3D& end)
{
 d_x = end.x() - start.x();
 d_y = end.y() - start.y();
 d_z = end.z() - start.z();
}

double 
Vector3D::length() const
{
  return std::sqrt(d_x*d_x+d_y*d_y+d_z*d_z);
}

double 
Vector3D::lengthSq() const
{
  return (d_x*d_x+d_y*d_y+d_z*d_z);
}

double 
Vector3D::dot(const Vector3D& vec) const
{
  return d_x*vec.d_x + d_y*vec.d_y + d_z*vec.d_z;
}

Vector3D 
Vector3D::cross(const Vector3D& vec) const
{
  return Vector3D(
        d_y*vec.d_z-d_z*vec.d_y,
        d_z*vec.d_x-d_x*vec.d_z,
        d_x*vec.d_y-d_y*vec.d_x);
}

bool 
Vector3D::operator==(const Vector3D& vec) const
{
  return (vec.d_x == d_x && vec.d_y == d_y && vec.d_z == d_z);
}

bool 
Vector3D::operator!=(const Vector3D& vec) const
{
  return (vec.d_x != d_x && vec.d_y != d_y && vec.d_z != d_z);
}

Vector3D& 
Vector3D::operator=(const Vector3D& vec)
{
  d_x = vec.d_x;
  d_y = vec.d_y;
  d_z = vec.d_z;
  return *this;
}

Vector3D 
Vector3D::operator*(const double val) const
{
  return Vector3D(d_x*val, d_y*val, d_z*val);
}

Vector3D 
Vector3D::operator+(const Vector3D& vec) const
{
  return Vector3D(d_x+vec.x(), d_y+vec.y(), d_z+vec.z());
}

Vector3D 
Vector3D::operator-(const Vector3D& vec) const
{
  return Vector3D(d_x-vec.x(), d_y-vec.y(), d_z-vec.z());
}

Vector3D& 
Vector3D::operator*=(const double val)
{
  d_x *= val;
  d_y *= val;
  d_z *= val;
  return *this;
}

Vector3D& 
Vector3D::operator+=(const Vector3D& vec)
{
  d_x += vec.x();
  d_y += vec.y();
  d_z += vec.z();
  return *this;
}

Vector3D& 
Vector3D::operator-=(const Vector3D& vec)
{
  d_x -= vec.x();
  d_y -= vec.y();
  d_z -= vec.z();
  return *this;
}

namespace Emu2DC {
  std::ostream& operator<<(std::ostream& out, const Vector3D& vec) 
  {
    out << "[" << vec.d_x << " " << vec.d_y << " " << vec.d_z << "]";
    return out;
  }
}


