#ifndef __EMU2DC_VECTOR3D_H__
#define __EMU2DC_VECTOR3D_H__

#include <Geometry/Point3D.h>

namespace Emu2DC {

  class Vector3D {

  public:

    friend std::ostream& operator<<(std::ostream& os, const Emu2DC::Vector3D& p);

  public:

    Vector3D(): d_x(std::numeric_limits<double>::max()), 
                d_y(std::numeric_limits<double>::max()), 
                d_z(std::numeric_limits<double>::max()) {}
    Vector3D(double x, double y, double z): d_x(x), d_y(y), d_z(z) {}
    Vector3D(const Vector3D& vec);
    Vector3D(const Point3D& start, const Point3D& end);

    ~Vector3D() {}

    void x(const double xx) {d_x = xx;}
    double x() const {return d_x;}
    void y(const double yy) {d_y = yy;}
    double y() const {return d_y;}
    void z(const double zz) {d_z = zz;}
    double z() const {return d_z;}
    
    double length() const;
    double lengthSq() const;

    double dot(const Vector3D& vec) const;
    Vector3D cross(const Vector3D& vec) const;

    bool operator==(const Vector3D& vec) const;
    bool operator!=(const Vector3D& vec) const;

    Vector3D& operator=(const Vector3D& vec);

    Vector3D operator*(const double) const;
    Vector3D operator+(const Vector3D& vec) const;
    Vector3D operator-(const Vector3D& vec) const;

    Vector3D& operator*=(const double);
    Vector3D& operator+=(const Vector3D& vec);
    Vector3D& operator-=(const Vector3D& vec);

  private:

    double d_x, d_y, d_z;

  }; // end class

} // end namespace

#endif

