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

#ifndef __MATITI_POINT3D_H__
#define __MATITI_POINT3D_H__

#include <iostream>
#include <limits>
#include <cmath>

namespace BrMPM {

  class Vector3D;

  class Point3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const BrMPM::Point3D& p);

  public:
    Point3D(): d_pos{std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max()} {}
    Point3D(double x, double y, double z): d_pos{x, y, z} {}
    Point3D(const Point3D& pt);
    ~Point3D() {}

    bool operator==(const Point3D& pt) const;
    bool operator!=(const Point3D& pt) const;

    void operator=(const Point3D& pt);
    inline void set(const double& val) {d_pos[0] = val; d_pos[1] = val; d_pos[2] = val;}

    Point3D operator+(const double& shift) const;
    Point3D operator-(const double& shift) const;

    Point3D operator+(const Vector3D& vec) const;
    Point3D operator-(const Vector3D& vec) const;
    Point3D& operator+=(const Vector3D& vec);
    Point3D& operator-=(const Vector3D& vec);

    Vector3D operator+(const Point3D& pt) const;
    Vector3D operator-(const Point3D& pt) const;

    double operator[](int index) const {return d_pos[index];}
    double& operator[](int index) {return d_pos[index];}

    void x(const double xx) {d_pos[0] = xx;}
    double x() const {return d_pos[0];}
    void y(const double yy) {d_pos[1] = yy;}
    double y() const {return d_pos[1];}
    void z(const double zz) {d_pos[2] = zz;}
    double z() const {return d_pos[2];}
  
    inline double distance(const Point3D& pt) const
    {
      double dist_sq = 0.0, dx = 0.0;
      for (unsigned int ii=0; ii < 3; ++ii) {
        dx = d_pos[ii] - pt.d_pos[ii];
        dist_sq += (dx*dx); 
      }
      return std::sqrt(dist_sq);
    }

  private:
    double d_pos[3];

  }; // end class Point3D
 
} // End namespace 

#endif //ifndef __MATITI_POINT3D_H__
