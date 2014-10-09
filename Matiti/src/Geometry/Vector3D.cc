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

#include <Geometry/Vector3D.h>
#include <cmath>

using namespace Matiti;

Vector3D::Vector3D(const Vector3D& vec)
{
  d_vec[0] = vec[0];
  d_vec[1] = vec[1];
  d_vec[2] = vec[2];
}

Vector3D::Vector3D(const Point3D& start, const Point3D& end)
{
 d_vec[0] = end.x() - start.x();
 d_vec[1] = end.y() - start.y();
 d_vec[2] = end.z() - start.z();
}

double 
Vector3D::length() const
{
  return std::sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
}

double 
Vector3D::lengthSq() const
{
  return (d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
}

Vector3D
Vector3D::invDirection() const
{
  double v1 = (d_vec[0] == 0.0) ? std::numeric_limits<double>::infinity() : 1.0/d_vec[0] ;
  double v2 = (d_vec[1] == 0.0) ? std::numeric_limits<double>::infinity() : 1.0/d_vec[1] ;
  double v3 = (d_vec[2] == 0.0) ? std::numeric_limits<double>::infinity() : 1.0/d_vec[2] ;
  return Vector3D(v1, v2, v3);
}

double 
Vector3D::dot(const Vector3D& vec) const
{
  return d_vec[0]*vec[0] + d_vec[1]*vec[1] + d_vec[2]*vec[2];
}

Vector3D 
Vector3D::cross(const Vector3D& vec) const
{
  return Vector3D(
        d_vec[1]*vec[2]-d_vec[2]*vec[1],
        d_vec[2]*vec[0]-d_vec[0]*vec[2],
        d_vec[0]*vec[1]-d_vec[1]*vec[0]);
}

bool 
Vector3D::operator==(const Vector3D& vec) const
{
  return (vec.x() == d_vec[0] && vec.y() == d_vec[1] && vec.z() == d_vec[2]);
}

bool 
Vector3D::operator!=(const Vector3D& vec) const
{
  return (vec.x() != d_vec[0] && vec.y() != d_vec[1] && vec.z() != d_vec[2]);
}

Vector3D& 
Vector3D::operator=(const Vector3D& vec)
{
  d_vec[0] = vec.x();
  d_vec[1] = vec.y();
  d_vec[2] = vec.z();
  return *this;
}

Vector3D 
Vector3D::operator*(const double val) const
{
  return Vector3D(d_vec[0]*val, d_vec[1]*val, d_vec[2]*val);
}

Vector3D 
Vector3D::operator/(const double val) const
{
  return Vector3D(d_vec[0]/val, d_vec[1]/val, d_vec[2]/val);
}

Vector3D 
Vector3D::operator*(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]*vec.x(), d_vec[1]*vec.y(), d_vec[2]*vec.z());
}

Vector3D 
Vector3D::operator/(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]/vec.x(), d_vec[1]/vec.y(), d_vec[2]/vec.z());
}

Vector3D 
Vector3D::operator+(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]+vec.x(), d_vec[1]+vec.y(), d_vec[2]+vec.z());
}

Vector3D 
Vector3D::operator-(const Vector3D& vec) const
{
  return Vector3D(d_vec[0]-vec.x(), d_vec[1]-vec.y(), d_vec[2]-vec.z());
}

Vector3D& 
Vector3D::operator*=(const double val)
{
  d_vec[0] *= val;
  d_vec[1] *= val;
  d_vec[2] *= val;
  return *this;
}

Vector3D& 
Vector3D::operator/=(const double val)
{
  d_vec[0] /= val;
  d_vec[1] /= val;
  d_vec[2] /= val;
  return *this;
}

Vector3D& 
Vector3D::operator+=(const Vector3D& vec)
{
  d_vec[0] += vec.x();
  d_vec[1] += vec.y();
  d_vec[2] += vec.z();
  return *this;
}

Vector3D& 
Vector3D::operator-=(const Vector3D& vec)
{
  d_vec[0] -= vec.x();
  d_vec[1] -= vec.y();
  d_vec[2] -= vec.z();
  return *this;
}

double 
Vector3D::max() const
{
  return std::max(std::max(d_vec[0],d_vec[1]), d_vec[2]);
}

double 
Vector3D::min() const
{
  return std::min(std::min(d_vec[0],d_vec[1]), d_vec[2]);
}

namespace Matiti {

  Vector3D 
  min(const Vector3D& v1, const Vector3D& v2)
  {
    return Vector3D(std::min(v1.x(),v2.x()), std::min(v1.y(),v2.y()), std::min(v1.z(), v2.z()));
  }

  Vector3D 
  max(const Vector3D& v1, const Vector3D& v2)
  {
    return Vector3D(std::max(v1.x(),v2.x()), std::max(v1.y(),v2.y()), std::max(v1.z(), v2.z()));
  }

  std::ostream& operator<<(std::ostream& out, const Vector3D& vec) 
  {
    out << "[" << vec.d_vec[0] << " " << vec.d_vec[1] << " " << vec.d_vec[2] << "]";
    return out;
  }
}


