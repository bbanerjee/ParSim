#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

using namespace Emu2DC;

Point3D::Point3D(const Point3D& pt)
{
  d_x = pt.d_x;
  d_y = pt.d_y;
  d_z = pt.d_z;
}

bool 
Point3D::operator==(const Point3D& pt) const
{
  return (pt.d_x == d_x && pt.d_y == d_y && pt.d_z == d_z);
}

bool 
Point3D::operator!=(const Point3D& pt) const
{
  return (pt.d_x != d_x && pt.d_y != d_y && pt.d_z != d_z);
}

Point3D& 
Point3D::operator=(const Point3D& pt)
{
  d_x = pt.d_x;
  d_y = pt.d_y;
  d_z = pt.d_z;
  return *this;
}

Point3D 
Point3D::operator+(const Vector3D& vec) const
{
  return Point3D(d_x+vec.x(), d_y+vec.y(), d_z+vec.z());
}

Point3D 
Point3D::operator-(const Vector3D& vec) const
{
  return Point3D(d_x-vec.x(), d_y-vec.y(), d_z-vec.z());
}

Point3D& Point3D::operator+=(const Vector3D& vec)
{
  d_x += vec.x();
  d_y += vec.y();
  d_z += vec.z();
  return *this;
}

Point3D& Point3D::operator-=(const Vector3D& vec)
{
  d_x -= vec.x();
  d_y -= vec.y();
  d_z -= vec.z();
  return *this;
}

Vector3D
Point3D::operator+(const Point3D& pt) const
{
  return Vector3D(d_x+pt.d_x, d_y+pt.d_y, d_z+pt.d_z);
}

Vector3D
Point3D::operator-(const Point3D& pt) const
{
  return Vector3D(d_x-pt.d_x, d_y-pt.d_y, d_z-pt.d_z);
}

namespace Emu2DC {
  std::ostream& operator<<(std::ostream& out, const Point3D& pt) 
  {
    out << "[" << pt.d_x << " " << pt.d_y << " " << pt.d_z << "]";
    return out;
  }
}


