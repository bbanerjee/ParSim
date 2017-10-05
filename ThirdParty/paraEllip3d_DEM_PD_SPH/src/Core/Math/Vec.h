#ifndef VEC_H
#define VEC_H

#include <Core/Types/RealTypes.h>
#include <boost/mpi.hpp>
#include <iomanip>
#include <iostream>
#include <array>
#include <vector>

namespace dem {

class IntVec;

class Vec
{

public:
  Vec()
    : d_data({{0, 0, 0}})
  {
  }
  Vec(REAL val)
    : d_data({{val, val, val}})
  {
  }
  Vec(REAL _x, REAL _y, REAL _z)
    : d_data({{_x, _y, _z}})
  {
  }
  Vec(std::vector<REAL>::const_iterator begin,
      std::vector<REAL>::const_iterator end) 
  {
    if (!(end - begin < 2)) {
      d_data[0] = *begin;
      d_data[1] = *(begin+1);
      d_data[2] = *(begin+2);
    } else {
      d_data[0] = 0.0;
      d_data[1] = 0.0;
      d_data[2] = 0.0;
    }
  }

  REAL x() const { return d_data[0]; }
  REAL y() const { return d_data[1]; }
  REAL z() const { return d_data[2]; }

  const REAL* data() const { return d_data.data(); }
  REAL* data() { return d_data.data(); }

  void setX(REAL _x) { d_data[0] = _x; }
  void setY(REAL _y) { d_data[1] = _y; }
  void setZ(REAL _z) { d_data[2] = _z; }
  void set(REAL _x, REAL _y, REAL _z)
  {
    d_data[0] = _x;
    d_data[1] = _y;
    d_data[2] = _z;
  }
  void set(Vec v)
  {
    d_data = v.d_data;
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

  friend double dot(const Vec& p1, const Vec& p2);
  friend Vec cross(const Vec& p1, const Vec& p2);

  //Vec operator%(Vec p) const;  // cross product of this vector and p
  //REAL operator*(Vec p) const; // dot product of this vector and p
  Vec operator*(REAL d) const;

  // Divides a vector by an int vector and produces a new vector
  Vec operator/(const IntVec& intVec) const;

  // Multiplies a vector by an int vector and produces a new vector
  Vec operator*(const IntVec& intVec) const;

  // Divides a vector by a real vector component-wise and produces a new vector
  Vec operator/(const Vec& intVec) const;

  // Multiplies a vector by a real vector component-wise and produces a new vector
  Vec operator*(const Vec& vec) const;


  REAL lengthSq() const;
  REAL length() const;

  void print(std::ostream& ofs) const;

  // Creates a Vec from a string that looks like "[1.0, -2.1e-3, 3]".
  static Vec fromString(const std::string& str);

  // Divides a vector by an int vector in place
  friend Vec& operator/=(Vec& realVec, const IntVec& intVec);

  // Multiplies a vector by an int vector in place
  friend Vec& operator*=(Vec& realVec, const IntVec& intVec);

  friend std::ostream& operator<<(std::ostream& os, const Vec& vec);

private:

  std::array<REAL, 3> d_data;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_data;
  }
};

// non-member functions
Vec operator*(REAL d, Vec v);
Vec operator/(Vec v, REAL d);
Vec operator-(Vec v);
REAL vnormL2(Vec v);
Vec vcos(Vec v);
Vec vacos(Vec v);
Vec rotateVec(Vec v, Vec alf);
Vec normalize(Vec v);

} // namespace dem

#endif
