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

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

using namespace Matiti;

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

Point3D& 
Point3D::operator=(const Point3D& pt)
{
  d_pos[0] = pt.d_pos[0];
  d_pos[1] = pt.d_pos[1];
  d_pos[2] = pt.d_pos[2];
  return *this;
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

namespace Matiti {
  std::ostream& operator<<(std::ostream& out, const Point3D& pt) 
  {
    out << "[" << pt.d_pos[0] << " " << pt.d_pos[1] << " " << pt.d_pos[2] << "]";
    return out;
  }
}


