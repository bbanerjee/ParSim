/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Element2D.h>
#include <Core/Node.h>
#include <Types/Types.h>
#include <Core/Exception.h>
#include <Geometry/Vector3D.h>

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace Matiti;

// The element constructor
Element2D::Element2D()
 :  d_area(0.0), d_surface_area(0.0), d_nodes(0)
{
}
  
Element2D::Element2D(const NodeP node1, const NodeP node2, const NodeP node3)
{
  d_nodes.push_back(node1);
  d_nodes.push_back(node2);
  d_nodes.push_back(node3);
  computeArea(); 
}

Element2D::Element2D(const NodeP node1, const NodeP node2, const NodeP node3, 
                     const NodeP node4)
{
  d_nodes.push_back(node1);
  d_nodes.push_back(node2);
  d_nodes.push_back(node3);
  d_nodes.push_back(node4);
  computeArea(); 
}

Element2D::Element2D(const NodePArray& nodes)
 : d_nodes(nodes)
{
  computeArea(); 
}

// The element destructor
Element2D::~Element2D()
{
}

void 
Element2D::initialize(const NodePArray& nodes)
{
  d_nodes = nodes;
  computeArea();
}

bool
Element2D::hasNode(const NodeP node) const
{
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    if (node == *iter) return true;
  }
  return false;
}

bool 
Element2D::isSubset(const NodePArray& nodeSet) const
{
  unsigned int count = 0;
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    for (auto set = nodeSet.begin(); set != nodeSet.end(); set++) {
      if (*iter == *set) {
        count++;
        break;
      }
    }
  }
  if (count == d_nodes.size()) return true;
  return false;
}

void
Element2D::computeArea()
{
  std::vector<Point3D> points;
  unsigned int count = 0;
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    NodeP node = *iter;
    if (node->onSurface()) {
      count++;
    }
    const Point3D& pos = node->position();
    points.push_back(Point3D(pos));
  }

  if (points.size() == 3) { // Triangle
    d_area = computeAreaTriangle(points[0], points[1], points[2]);
  } else if (points.size() == 4) { // Quadrilateral
    d_area = computeAreaQuadrilateral(points[0], points[1], points[2], points[3]);
  } else {
    std::ostringstream out;
    out << "**ERROR** Elements can only have 4, 6, or 8 nodes.";
    throw Exception(out.str(), __FILE__, __LINE__);
  }
  
  if (count != d_nodes.size()) {
    d_surface_area = 0.0;
  } else {
    d_surface_area = d_area;
  }
}

double
Element2D::computeAreaTriangle(const Point3D& p1,
                               const Point3D& p2,
                               const Point3D& p3) const
{
  Vector3D vec1(p1, p2);
  Vector3D vec2(p1, p3);
  Vector3D vec1x2 = vec1.cross(vec2);
  double area = 0.5*vec1x2.length();
  return area;
}

double
Element2D::computeAreaQuadrilateral(const Point3D& p1,
                                    const Point3D& p2,
                                    const Point3D& p3,
                                    const Point3D& p4) const
{
  double area = computeAreaTriangle(p1, p2, p3);
  area += computeAreaTriangle(p1, p3, p4);
  return area;
}


namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Element2D& elem)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Element2D " << " = [";
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
