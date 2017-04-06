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
    : d_x(0)
    , d_y(0)
    , d_z(0)
  {
  }
  IntVec(int d)
    : d_x(d)
    , d_y(d)
    , d_z(d)
  {
  }
  IntVec(int _x, int _y, int _z)
    : d_x(_x)
    , d_y(_y)
    , d_z(_z)
  {
  }
  int x() const { return d_x; }
  int y() const { return d_y; }
  int z() const { return d_z; }
  void setX(int _x) { d_x = _x; }
  void setY(int _y) { d_y = _y; }
  void setZ(int _z) { d_z = _z; }
  void set(int _x, int _y, int _z)
  {
    d_x = _x;
    d_y = _y;
    d_z = _z;
  }
  void set(IntVec v)
  {
    d_x = v.x();
    d_y = v.y();
    d_z = v.z();
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
  int d_x;
  int d_y;
  int d_z;

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
IntVec operator*(int d, IntVec v);
IntVec operator-(IntVec v);

} // namespace dem

#endif
