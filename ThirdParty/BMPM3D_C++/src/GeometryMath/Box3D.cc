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

#include <GeometryMath/Box3D.h>

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
