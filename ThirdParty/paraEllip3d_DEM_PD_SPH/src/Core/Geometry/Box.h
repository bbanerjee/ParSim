#ifndef BOX_H
#define BOX_H

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class Box
{
private:
  REAL d_dimx;
  REAL d_dimy;
  REAL d_dimz;
  Vec d_center;
  Vec d_v1; // lower corner of the box, minimum x,y,z value
  Vec d_v2; // lower corner of the box, maximum x,y,z value

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_dimx;
    ar& d_dimy;
    ar& d_dimz;
    ar& d_center;
    ar& d_v1;
    ar& d_v2;
  }

public:
  Box()
    : d_dimx(0)
    , d_dimy(0)
    , d_dimz(0)
    , d_center(0)
    , d_v1(0)
    , d_v2(0)
  {
  }

  Box(REAL dx, REAL dy, REAL dz, Vec c)
    : d_dimx(dx)
    , d_dimy(dy)
    , d_dimz(dz)
    , d_center(c)
  {
    d_v1.setX(c.x() - dx / 2.0);
    d_v1.setY(c.y() - dy / 2.0);
    d_v1.setZ(c.z() - dz / 2.0);
    d_v2.setX(c.x() + dx / 2.0);
    d_v2.setY(c.y() + dy / 2.0);
    d_v2.setZ(c.z() + dz / 2.0);
  }

  Box(REAL dx, REAL dy, REAL dz, Vec ref, int i)
    : d_dimx(dx)
    , d_dimy(dy)
    , d_dimz(dz)
    , d_v1(std::move(ref))
  {
    d_center.setX(d_v1.x() + dx / 2.0);
    d_center.setY(d_v1.y() + dy / 2.0);
    d_center.setZ(d_v1.z() + dz / 2.0);
    d_v2.setX(d_v1.x() + dx);
    d_v2.setY(d_v1.y() + dy);
    d_v2.setZ(d_v1.z() + dz);
  }

  Box(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2)
    : d_dimx(x2 - x1)
    , d_dimy(y2 - y1)
    , d_dimz(z2 - z1)
    , d_v1(x1, y1, z1)
    , d_v2(x2, y2, z2)
  {
    d_center = (d_v1 + d_v2) / 2;
  }

  Box(Vec _d_v1, Vec _d_v2)
    : d_v1(std::move(_d_v1))
    , d_v2(std::move(_d_v2))
  {
    d_center = (d_v1 + d_v2) / 2;
    Vec vt = d_v2 - d_v1;
    d_dimx = vt.x();
    d_dimy = vt.y();
    d_dimz = vt.z();
  }

  // Create a new box using an existing box +- some buffer distance
  Box(const Box& box, const double& buffer)
  {
    d_v1 = box.d_v1 - buffer;
    d_v2 = box.d_v2 + buffer;
    d_center = (d_v1 + d_v2) / 2;
    Vec d_v21 = d_v2 - d_v1;
    d_dimx = d_v21.x();
    d_dimy = d_v21.y();
    d_dimz = d_v21.z();
  }


  REAL getDimx() const { return d_dimx; }
  REAL getDimy() const { return d_dimy; }
  REAL getDimz() const { return d_dimz; }
  Vec getCenter() const { return d_center; }
  Vec getMinCorner() const { return d_v1; }
  Vec getMaxCorner() const { return d_v2; }
  REAL getVolume() const { return d_dimx * d_dimy * d_dimz; }

  void setDimx(REAL dx) { d_dimx = dx; }
  void setDimy(REAL dy) { d_dimy = dy; }
  void setDimz(REAL dz) { d_dimz = dz; }
  void setCenter(Vec v) { d_center = v; }
  void setV1(Vec v) { d_v1 = v; }
  void setV2(Vec v) { d_v2 = v; }
  Vec randomPoint() const;
  void print() const;
  void set(REAL dx, REAL dy, REAL dz, Vec c)
  {
    d_dimx = dx;
    d_dimy = dy;
    d_dimz = dz;
    d_center = c;
    d_v1.setX(c.x() - dx / 2.0);
    d_v1.setY(c.y() - dy / 2.0);
    d_v1.setZ(c.z() - dz / 2.0);
    d_v2.setX(c.x() + dx / 2.0);
    d_v2.setY(c.y() + dy / 2.0);
    d_v2.setZ(c.z() + dz / 2.0);
  }

  void set(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2)
  {
    d_dimx = x2 - x1;
    d_dimy = y2 - y1;
    d_dimz = z2 - z1;
    d_v1.set(x1, y1, z1);
    d_v2.set(x2, y2, z2);
    d_center = (d_v1 + d_v2) / 2;
  }

  bool inside(const Vec& pt, REAL tolerance) const;

  bool outside(const Vec& pt) const;

  friend std::ostream& operator<<(std::ostream& os, const Box& box);

};

} // namespace dem

#endif
