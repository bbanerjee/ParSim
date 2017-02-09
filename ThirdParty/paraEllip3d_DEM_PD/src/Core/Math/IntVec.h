#ifndef ELLIP3D_INTVEC_H
#define ELLIP3D_INTVEC_H

#include <boost/mpi.hpp>
#include <iomanip>
#include <iostream>

namespace dem {

class IntVec
{

public:
  IntVec()
    : x(0)
    , y(0)
    , z(0)
  {
  }
  IntVec(int d)
    : x(d)
    , y(d)
    , z(d)
  {
  }
  IntVec(int _x, int _y, int _z)
    : x(_x)
    , y(_y)
    , z(_z)
  {
  }
  int getX() const { return x; }
  int getY() const { return y; }
  int getZ() const { return z; }
  void setX(int _x) { x = _x; }
  void setY(int _y) { y = _y; }
  void setZ(int _z) { z = _z; }
  void set(int _x, int _y, int _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }
  void set(IntVec v)
  {
    x = v.getX();
    y = v.getY();
    z = v.getZ();
  }

  bool operator==(const IntVec v);
  bool operator==(const int d);
  bool operator!=(const IntVec v);
  void operator+=(const IntVec v);
  void operator-=(const IntVec v);
  void operator*=(int d);
  IntVec operator+(IntVec v) const;
  IntVec operator-(IntVec v) const;
  int operator*(IntVec p) const; // dot product of this vector and p
  IntVec operator*(int d) const;
  void print(std::ostream& ofs) const;

  // Creates an IntVector from a string that looks like "[1, 2, 3]".
  static IntVec fromString(const std::string& str);

private:
  int x;
  int y;
  int z;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& x;
    ar& y;
    ar& z;
  }
};

// non-member functions
IntVec operator*(int d, IntVec v);
IntVec operator-(IntVec v);

} // namespace dem

#endif
