#ifndef __EMU2DC_POINT3D_H__
#define __EMU2DC_POINT3D_H__

#include <iostream>
#include <limits>

namespace Emu2DC {

  class Vector3D;

  class Point3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Emu2DC::Point3D& p);

  public:
    Point3D(): d_x(std::numeric_limits<double>::max()), 
               d_y(std::numeric_limits<double>::max()), 
               d_z(std::numeric_limits<double>::max()) {}
    Point3D(double x, double y, double z): d_x(x), d_y(y), d_z(z) {}
    Point3D(const Point3D& pt);
    ~Point3D() {}

    bool operator==(const Point3D& pt) const;
    bool operator!=(const Point3D& pt) const;
    Point3D& operator=(const Point3D& pt);

    Point3D operator+(const Vector3D& vec) const;
    Point3D operator-(const Vector3D& vec) const;
    Point3D& operator+=(const Vector3D& vec);
    Point3D& operator-=(const Vector3D& vec);

    void x(const double xx) {d_x = xx;}
    double x() const {return d_x;}
    void y(const double yy) {d_y = yy;}
    double y() const {return d_y;}
    void z(const double zz) {d_z = zz;}
    double z() const {return d_z;}
  
  private:
    double d_x, d_y, d_z;

  }; // end class Point3D
 
} // End namespace 

#endif //ifndef __EMU2DC_POINT3D_H__
