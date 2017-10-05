#ifndef BOX_H
#define BOX_H

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class Box
{

public:
  Box()
    : dimx(0)
    , dimy(0)
    , dimz(0)
    , center(0)
    , v1(0)
    , v2(0)
  {
  }

  Box(REAL dx, REAL dy, REAL dz, Vec c)
    : dimx(dx)
    , dimy(dy)
    , dimz(dz)
    , center(c)
  {
    v1.setX(c.x() - dx / 2.0);
    v1.setY(c.y() - dy / 2.0);
    v1.setZ(c.z() - dz / 2.0);
    v2.setX(c.x() + dx / 2.0);
    v2.setY(c.y() + dy / 2.0);
    v2.setZ(c.z() + dz / 2.0);
  }

  Box(REAL dx, REAL dy, REAL dz, Vec ref, int i)
    : dimx(dx)
    , dimy(dy)
    , dimz(dz)
    , v1(std::move(ref))
  {
    center.setX(v1.x() + dx / 2.0);
    center.setY(v1.y() + dy / 2.0);
    center.setZ(v1.z() + dz / 2.0);
    v2.setX(v1.x() + dx);
    v2.setY(v1.y() + dy);
    v2.setZ(v1.z() + dz);
  }

  Box(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2)
    : dimx(x2 - x1)
    , dimy(y2 - y1)
    , dimz(z2 - z1)
    , v1(x1, y1, z1)
    , v2(x2, y2, z2)
  {
    center = (v1 + v2) / 2;
  }

  Box(Vec _v1, Vec _v2)
    : v1(std::move(_v1))
    , v2(std::move(_v2))
  {
    center = (v1 + v2) / 2;
    Vec vt = v2 - v1;
    dimx = vt.x();
    dimy = vt.y();
    dimz = vt.z();
  }

  // Create a new box using an existing box +- some buffer distance
  Box(const Box& box, const double& buffer)
  {
    v1 = box.v1 - buffer;
    v2 = box.v2 + buffer;
    center = (v1 + v2) / 2;
    Vec v21 = v2 - v1;
    dimx = v21.x();
    dimy = v21.y();
    dimz = v21.z();
  }


  REAL getDimx() const { return dimx; }
  REAL getDimy() const { return dimy; }
  REAL getDimz() const { return dimz; }
  Vec getCenter() const { return center; }
  Vec getMinCorner() const { return v1; }
  Vec getMaxCorner() const { return v2; }
  REAL getVolume() const { return dimx * dimy * dimz; }

  void setDimx(REAL dx) { dimx = dx; }
  void setDimy(REAL dy) { dimy = dy; }
  void setDimz(REAL dz) { dimz = dz; }
  void setCenter(Vec v) { center = v; }
  void setV1(Vec v) { v1 = v; }
  void setV2(Vec v) { v2 = v; }
  Vec randomPoint() const;
  void print() const;
  void set(REAL dx, REAL dy, REAL dz, Vec c)
  {
    dimx = dx;
    dimy = dy;
    dimz = dz;
    center = c;
    v1.setX(c.x() - dx / 2.0);
    v1.setY(c.y() - dy / 2.0);
    v1.setZ(c.z() - dz / 2.0);
    v2.setX(c.x() + dx / 2.0);
    v2.setY(c.y() + dy / 2.0);
    v2.setZ(c.z() + dz / 2.0);
  }

  void set(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2)
  {
    dimx = x2 - x1;
    dimy = y2 - y1;
    dimz = z2 - z1;
    v1.set(x1, y1, z1);
    v2.set(x2, y2, z2);
    center = (v1 + v2) / 2;
  }

  bool inside(const Vec& pt, REAL tolerance) const;

  bool outside(const Vec& pt) const;

  friend std::ostream& operator<<(std::ostream& os, const Box& box);

private:
  REAL dimx;
  REAL dimy;
  REAL dimz;
  Vec center;
  Vec v1; // lower corner of the box, minimum x,y,z value
  Vec v2; // lower corner of the box, maximum x,y,z value

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& dimx;
    ar& dimy;
    ar& dimz;
    ar& center;
    ar& v1;
    ar& v2;
  }
};

} // namespace dem

#endif
