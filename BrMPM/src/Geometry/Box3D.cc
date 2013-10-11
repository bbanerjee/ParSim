#include <Geometry/Box3D.h>

using namespace BrMPM;

Box3D::Box3D()
{
}

Box3D::Box3D(const Box3D& box)
  : d_lower(box.d_lower), d_upper(box.d_upper)
{
  // If wrong lower/upper fix.
  correctBoundingBox();
}

Box3D::Box3D(const Point3D& lower, const Point3D& upper)
  : d_lower(lower), d_upper(upper)
{
  // If wrong lower/upper fix.
  correctBoundingBox();
}

Box3D::~Box3D()
{
}

Box3D& 
Box3D::operator=(const Box3D& box)
{
  d_lower = box.d_lower;
  d_upper = box.d_upper;
  return *this;
}

bool 
Box3D::overlaps(const Box3D& box, double epsilon) const
{
  if(d_lower.x()+epsilon > box.d_upper.x() || d_upper.x() < box.d_lower.x()+epsilon) {
    return false;
  }
  if(d_lower.y()+epsilon > box.d_upper.y() || d_upper.y() < box.d_lower.y()+epsilon) {
    return false;
  }
  if(d_lower.z()+epsilon > box.d_upper.z() || d_upper.z() < box.d_lower.z()+epsilon) {
    return false;
  }
  return true;
}

bool 
Box3D::contains(const Point3D& pt) const
{
  return pt.x() >= d_lower.x() && pt.y() >= d_lower.y()
      && pt.z() >= d_lower.z() && pt.x() <= d_upper.x()
      && pt.y() <= d_upper.y() && pt.z() <= d_upper.z();
}

Point3D 
Box3D::lower() const
{
  return d_lower;
}

Point3D 
Box3D::upper() const
{
  return d_upper;
}

bool 
Box3D::isDegenerate()
{
  return d_lower.x() >= d_upper.x() || d_lower.y() >= d_upper.y() || d_lower.z() >= d_upper.z();
}

void 
Box3D::correctBoundingBox()
{
  for( int index = 0; index < 3; index++ ) {
    if(d_upper[index] < d_lower[index] ) {
      double temp = d_upper[index];
      d_upper[index] = d_lower[index];
      d_lower[index] = temp;
    }
  }
}


namespace BrMPM {
  std::ostream& 
  operator<<(std::ostream& os, const Box3D& box)
  {
    os << box.lower() << " - " << box.upper();
    return os;
  }
}
