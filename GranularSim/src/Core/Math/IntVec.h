#ifndef ELLIP3D_INTVEC_H
#define ELLIP3D_INTVEC_H

#include <boost/mpi.hpp>
#include <iomanip>
#include <iostream>
#include <array>

namespace dem {

class IntVec
{

public:
  IntVec()
    : d_data({{0,0,0}})
  {
  }
  IntVec(int d)
    : d_data({{d, d, d}})
  {
  }
  IntVec(int _x, int _y, int _z)
    : d_data({{_x, _y, _z}})
  {
  }

  IntVec(const IntVec& vec)
    : d_data(vec.d_data)
  {
  }

  IntVec& operator=(const IntVec& vec)
  {
    this->d_data = vec.d_data;
    return *this;
  }

  IntVec& operator=(const std::array<int,3>& vec)
  {
    this->d_data = vec;
    return *this;
  }

  int x() const { return d_data[0]; }
  int y() const { return d_data[1]; }
  int z() const { return d_data[2]; }

  int& x() { return d_data[0]; }
  int& y() { return d_data[1]; }
  int& z() { return d_data[2]; }

  const int* data() const { return d_data.data(); }
  int* data() { return d_data.data(); }

  void setX(int _x) { d_data[0] = _x; }
  void setY(int _y) { d_data[1] = _y; }
  void setZ(int _z) { d_data[2] = _z; }

  void set(int _x, int _y, int _z)
  {
    d_data = {{_x, _y, _z}};
  }
  void set(IntVec v)
  {
    d_data = v.d_data;
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

  friend std::ostream& operator<<(std::ostream& os, const IntVec& vec);

private:
  std::array<int, 3> d_data;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_data;
  }
};

// non-member functions
IntVec operator*(int d, IntVec v);
IntVec operator-(IntVec v);

} // namespace dem

#endif
