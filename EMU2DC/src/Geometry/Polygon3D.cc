#include <Geometry/Polygon3D.h>

using namespace Emu2DC;


Polygon3D::Polygon3D()
{
}

Polygon3D::Polygon3D(const std::vector<Point3D>& poly)
{
  for (auto iter = poly.begin(); iter != poly.end(); ++iter) {
    d_vertices.emplace_back(Point3D(*iter));
  } 
}

Polygon3D::~Polygon3D()
{
}

bool 
Polygon3D::operator==(const Polygon3D& poly) const
{
  if (d_vertices.size() != poly.numVertices()) return false;
  for (auto this_iter = d_vertices.begin(), iter = poly.begin(); iter != poly.end(); ++this_iter, ++iter) {
    if (*iter != *this_iter) return false;
  } 
  return true;
}


unsigned int 
Polygon3D::numVertices() const
{
  return d_vertices.size();
}

void 
Polygon3D::addVertex(const Point3D& pt)
{
  d_vertices.emplace_back(Point3D(pt));
}

Polygon3D& 
Polygon3D::operator+=(const Point3D& pt)
{
  d_vertices.emplace_back(Point3D(pt));
  return *this;
}
    
const Point3D& 
Polygon3D::vertex(const int& index) const
{
  // TODO: check index first ?
  return d_vertices[index];
}


std::ostream& operator<<(std::ostream& out, const Polygon3D& poly)
{
  out.setf(std::ios::floatfield);
  out.precision(6);
  for (auto iter = poly.begin(); iter != poly.end(); ++iter) {
      out << *iter ;
  }
  out << std::endl;
  return out;
}
