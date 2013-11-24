#include <Point3D.h>
#include <Vector3D.h>

using namespace BrMPM;

Point3D::Point3D(const Point3D& pt)
{
  d_pos[0] = pt.d_pos[0];
  d_pos[1] = pt.d_pos[1];
  d_pos[2] = pt.d_pos[2];
}

bool 
Point3D::operator==(const Point3D& pt) const
{
  return (pt.d_pos[0] == d_pos[0] && pt.d_pos[1] == d_pos[1] && pt.d_pos[2] == d_pos[2]);
}

bool 
Point3D::operator!=(const Point3D& pt) const
{
  return (pt.d_pos[0] != d_pos[0] && pt.d_pos[1] != d_pos[1] && pt.d_pos[2] != d_pos[2]);
}

void
Point3D::operator=(const Point3D& pt)
{
  d_pos[0] = pt.d_pos[0];
  d_pos[1] = pt.d_pos[1];
  d_pos[2] = pt.d_pos[2];
}

Point3D 
Point3D::operator+(const double& shift) const
{
  return Point3D(d_pos[0]+shift, d_pos[1]+shift, d_pos[2]+shift);
}

Point3D 
Point3D::operator-(const double& shift) const
{
  return Point3D(d_pos[0]-shift, d_pos[1]-shift, d_pos[2]-shift);
}

Point3D 
Point3D::operator+(const Vector3D& vec) const
{
  return Point3D(d_pos[0]+vec.x(), d_pos[1]+vec.y(), d_pos[2]+vec.z());
}

Point3D 
Point3D::operator-(const Vector3D& vec) const
{
  return Point3D(d_pos[0]-vec.x(), d_pos[1]-vec.y(), d_pos[2]-vec.z());
}

Point3D& Point3D::operator+=(const Vector3D& vec)
{
  d_pos[0] += vec.x();
  d_pos[1] += vec.y();
  d_pos[2] += vec.z();
  return *this;
}

Point3D& Point3D::operator-=(const Vector3D& vec)
{
  d_pos[0] -= vec.x();
  d_pos[1] -= vec.y();
  d_pos[2] -= vec.z();
  return *this;
}

Vector3D
Point3D::operator+(const Point3D& pt) const
{
  return Vector3D(d_pos[0]+pt.d_pos[0], d_pos[1]+pt.d_pos[1], d_pos[2]+pt.d_pos[2]);
}

Vector3D
Point3D::operator-(const Point3D& pt) const
{
  return Vector3D(d_pos[0]-pt.d_pos[0], d_pos[1]-pt.d_pos[1], d_pos[2]-pt.d_pos[2]);
}

namespace BrMPM {
  std::ostream& operator<<(std::ostream& out, const Point3D& pt) 
  {
    out << "[" << pt.d_pos[0] << " " << pt.d_pos[1] << " " << pt.d_pos[2] << "]";
    return out;
  }

}


