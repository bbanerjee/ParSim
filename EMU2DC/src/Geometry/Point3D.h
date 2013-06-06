#ifndef __EMU2DC_POINT3D_H__
#define __EMU2DC_POINT3D_H__

#include <iostream>
#include <limits>
#include <cmath>

namespace Emu2DC {

  class Vector3D;

  class Point3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Emu2DC::Point3D& p);

  public:
    Point3D(): d_pos{std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max()} {}
    Point3D(double x, double y, double z): d_pos{x, y, z} {}
    Point3D(const Point3D& pt);
    ~Point3D() {}

    bool operator==(const Point3D& pt) const;
    bool operator!=(const Point3D& pt) const;
    Point3D& operator=(const Point3D& pt);

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

#endif //ifndef __EMU2DC_POINT3D_H__
