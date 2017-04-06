#ifndef VEC_H
#define VEC_H

#include <Core/Types/realtypes.h>
#include <boost/mpi.hpp>
#include <iomanip>
#include <iostream>

namespace dem {

class Vec
{

public:
  Vec()
    : d_x(0)
    , d_y(0)
    , d_z(0)
  {
  }
  Vec(REAL val)
    : d_x(val)
    , d_y(val)
    , d_z(val)
  {
  }
  Vec(REAL _x, REAL _y, REAL _z)
    : d_x(_x)
    , d_y(_y)
    , d_z(_z)
  {
  }
  REAL x() const { return d_x; }
  REAL y() const { return d_y; }
  REAL z() const { return d_z; }
  void setX(REAL _x) { d_x = _x; }
  void setY(REAL _y) { d_y = _y; }
  void setZ(REAL _z) { d_z = _z; }
  void set(REAL _x, REAL _y, REAL _z)
  {
    d_x = _x;
    d_y = _y;
    d_z = _z;
  }
  void set(Vec v)
  {
    d_x = v.x();
    d_y = v.y();
    d_z = v.z();
  }

  bool operator==(const Vec v);
  bool operator==(const REAL d);
  bool operator!=(const Vec v);
  void operator+=(const Vec v);
  void operator-=(const Vec v);
  void operator*=(REAL d);
  void operator/=(REAL d);
  Vec operator+(Vec v) const;
  Vec operator-(Vec v) const;
  Vec operator%(Vec p) const;  // cross product of this vector and p
  REAL operator*(Vec p) const; // dot product of this vector and p
  Vec operator*(REAL d) const;

  REAL lengthSq() const;
  REAL length() const;

  void print(std::ostream& ofs) const;

  // Creates a Vec from a string that looks like "[1.0, -2.1e-3, 3]".
  static Vec fromString(const std::string& str);

  friend std::ostream& operator<<(std::ostream& os, const Vec& vec);

private:
  REAL d_x;
  REAL d_y;
  REAL d_z;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_x;
    ar& d_y;
    ar& d_z;
  }
};

// non-member functions
Vec operator*(REAL d, Vec v);
Vec operator/(Vec v, REAL d);
Vec operator-(Vec v);
REAL vfabs(Vec v);
Vec vcos(Vec v);
Vec vacos(Vec v);
Vec rotateVec(Vec v, Vec alf);
Vec normalize(Vec v);

} // namespace dem

#endif
