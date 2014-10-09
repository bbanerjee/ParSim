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

#include <GeometryMath/Polygon3D.h>

using namespace BrMPM;


Polygon3D::Polygon3D()
{
}

Polygon3D::Polygon3D(const std::vector<Point3D>& poly)
{
  for (auto iter = poly.begin(); iter != poly.end(); ++iter) {
    d_vertices.emplace_back(Point3D(*iter));
  } 
}

Polygon3D::Polygon3D(const Polygon3D& poly)
{
  for (auto iter = (poly.d_vertices).begin(); iter != (poly.d_vertices).end(); ++iter) {
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

const Point3D& 
Polygon3D::operator[](const int& index) const
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
