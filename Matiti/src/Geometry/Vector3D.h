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

#ifndef __MATITI_VECTOR3D_H__
#define __MATITI_VECTOR3D_H__

#include <Types/Types.h>
#include <Geometry/Point3D.h>

namespace Matiti {

  class Vector3D {

  public:

    friend std::ostream& operator<<(std::ostream& os, const Matiti::Vector3D& p);

  public:

    Vector3D(): d_vec({{std::numeric_limits<double>::max(), 
                        std::numeric_limits<double>::max(), 
                        std::numeric_limits<double>::max()}}) {}
    Vector3D(double x, double y, double z): d_vec({{x, y, z}}) {}
    Vector3D(const Vector3D& vec);
    Vector3D(const Point3D& start, const Point3D& end);

    ~Vector3D() {}

    //  **WARNING** Not checking index.  Do checking outside.
    double operator[](int index) const {return d_vec[index];}
    double& operator[](int index) {return d_vec[index];}

    void x(const double xx) {d_vec[0] = xx;}
    double x() const {return d_vec[0];}
    void y(const double yy) {d_vec[1] = yy;}
    double y() const {return d_vec[1];}
    void z(const double zz) {d_vec[2] = zz;}
    double z() const {return d_vec[2];}
    
    double length() const;
    double lengthSq() const;

    Vector3D invDirection() const;

    double dot(const Vector3D& vec) const;
    Vector3D cross(const Vector3D& vec) const;

    bool operator==(const Vector3D& vec) const;
    bool operator!=(const Vector3D& vec) const;

    Vector3D& operator=(const Vector3D& vec);

    Vector3D operator*(const double fac) const;
    Vector3D operator/(const double fac) const;
    Vector3D operator*(const Vector3D& vec) const;
    Vector3D operator/(const Vector3D& vec) const;

    Vector3D operator+(const Vector3D& vec) const;
    Vector3D operator-(const Vector3D& vec) const;

    Vector3D& operator*=(const double);
    Vector3D& operator/=(const double);
    Vector3D& operator+=(const Vector3D& vec);
    Vector3D& operator-=(const Vector3D& vec);

    void reset() {d_vec[0] = 0.0; d_vec[1] = 0.0; d_vec[2] = 0.0;}
    bool isnan() const {return std::isnan(d_vec[0]) || std::isnan(d_vec[1]) || std::isnan(d_vec[2]);}

    double max() const;
    double min() const;

  private:

    //double d_x, d_y, d_z;
    Array3 d_vec;

  }; // end class

} // end namespace

namespace Matiti {

 Vector3D min(const Vector3D& v1, const Vector3D& v2);
 Vector3D max(const Vector3D& v1, const Vector3D& v2);

} // end namespace


#endif

