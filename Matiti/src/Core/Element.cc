#include <Element.h>
#include <Node.h>
#include <Types.h>
#include <Exception.h>
#include <Geometry/Vector3D.h>

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace Matiti;

// The element constructor
Element::Element()
 : d_id(0), d_volume(0.0), d_nodes(0)
{
}
  
Element::Element(const int& id, const NodePArray& nodes)
 : d_id(id), d_nodes(nodes)
{
  computeVolume(); 
}

// The element destructor
Element::~Element()
{
}

void 
Element::initialize(const int id, const NodePArray& nodes)
{
  d_id = id;
  d_nodes = nodes;
  computeVolume();
}

void
Element::computeVolume()
{
  computeVolume3D();
}

void
Element::computeVolume2D()
{
  // Set up coordinate array
  std::vector<double> xcoords, ycoords;
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    NodeP node = *iter;
    const Point3D& pos = node->position();
    xcoords.push_back(pos.x());
    ycoords.push_back(pos.y());
  }
  xcoords.push_back(xcoords[0]);
  ycoords.push_back(ycoords[0]);

  // Compute area
  double area = 0.0;
  for (unsigned int ii = 0; ii < xcoords.size()-1; ii++) {
    area += (xcoords[ii]*ycoords[ii+1]-xcoords[ii+1]*ycoords[ii]);
  }
  area = std::abs(0.5*area);
  d_volume = area;
}

void
Element::computeVolume3D()
{
  std::vector<Point3D> points;
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    NodeP node = *iter;
    const Point3D& pos = node->position();
    points.push_back(Point3D(pos));
  }

  if (d_nodes.size() == 4) {           // Tets
    d_volume = computeVolumeTetrahedron(points[0], points[1], points[2], points[3]);
  } else if (d_nodes.size() == 6) {    // Prism
    d_volume = computeVolumeTetrahedron(points[0], points[1], points[2], points[3]);
    d_volume += computeVolumeTetrahedron(points[1], points[2], points[3], points[4]);
    d_volume += computeVolumeTetrahedron(points[2], points[3], points[4], points[5]);
  } else if (d_nodes.size() == 8) {    // Hex
    d_volume = computeVolumeTetrahedron(points[0], points[1], points[2], points[5]);
    d_volume += computeVolumeTetrahedron(points[0], points[2], points[7], points[5]);
    d_volume += computeVolumeTetrahedron(points[0], points[2], points[3], points[7]);
    d_volume += computeVolumeTetrahedron(points[0], points[4], points[5], points[7]);
    d_volume += computeVolumeTetrahedron(points[2], points[7], points[5], points[6]);
  } else {
    std::ostringstream out;
    out << "**ERROR** Elements can only be tetrahedra, prisms, or hexahedra." << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  if (d_volume <= 0.0) {
    std::ostringstream out;
    out << "**ERROR** Zero volume element found." << *this;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}

double
Element::computeVolumeTetrahedron(const Point3D& p0, const Point3D& p1, const Point3D& p2,
                                  const Point3D& p3) const
{
  Vector3D vec01(p0, p1);
  Vector3D vec02(p0, p2);
  Vector3D vec03(p0, p3);
  Vector3D vec1x2 = vec01.cross(vec02);
  double volume = 1.0/6.0*std::abs(vec1x2.dot(vec03));
  return volume;
}

void Element::computeGeometry2D(double& area, double& xlength, double& ylength) const
{
  // Set up coordinate array
  std::vector<double> xcoords, ycoords;
  for (constNodePIterator iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    NodeP node = *iter;
    const Point3D& pos = node->position();
    xcoords.push_back(pos.x());
    ycoords.push_back(pos.y());
  }
  xcoords.push_back(xcoords[0]);
  ycoords.push_back(ycoords[0]);

  // Compute area
  area = 0.0;
  for (unsigned int ii = 0; ii < xcoords.size()-1; ii++) {
    area += (xcoords[ii]*ycoords[ii+1]-xcoords[ii+1]*ycoords[ii]);
  }
  area = std::abs(0.5*area);

  // Compute xlength, ylength
  double xmax = *(std::max_element(xcoords.begin(), xcoords.end()));
  double xmin = *(std::min_element(xcoords.begin(), xcoords.end()));
  double ymax = *(std::max_element(ycoords.begin(), ycoords.end()));
  double ymin = *(std::min_element(ycoords.begin(), ycoords.end()));
  xlength = std::abs(xmax-xmin);
  ylength = std::abs(ymax-ymin);
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Element& elem)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Element " << elem.d_id << " = [";
    for (auto iter = (elem.d_nodes).begin(); iter != (elem.d_nodes).end(); ++iter) {
      out << (*iter)->getID() << " ";
    }
    out << "]" << std::endl;
    //for (auto iter = (elem.d_nodes).begin(); iter != (elem.d_nodes).end(); ++iter) {
    //  out << (*iter) << std::endl;
    //}
    return out;
  }
}
