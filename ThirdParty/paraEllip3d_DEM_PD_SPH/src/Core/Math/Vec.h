#ifndef VEC_H
#define VEC_H

#include <Core/Types/RealTypes.h>
#include <Core/Const/Constants.h>
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
  {
    d_data[0] = 0; d_data[1] = 0; d_data[2] = 0;
  }

  Vec(REAL val) 
  {
    d_data[0] = val; d_data[1] = val; d_data[2] = val;
  }

  Vec(REAL _x, REAL _y, REAL _z) 
  {
    d_data[0] = _x; d_data[1] = _y; d_data[2] = _z;
  }

  Vec(std::vector<REAL>::const_iterator begin,
      std::vector<REAL>::const_iterator end) 
  {
    if (!(end - begin < 2)) {
      d_data[0] = *begin;
      d_data[1] = *(begin+1);
      d_data[2] = *(begin+2);
    } else {
      d_data[0] = 0.0; d_data[1] = 0.0; d_data[2] = 0.0;
    }
  }

  inline REAL x() const { return d_data[0]; }
  inline REAL y() const { return d_data[1]; }
  inline REAL z() const { return d_data[2]; }

  const REAL* data() const { return d_data; }
  REAL* data() { return d_data; }

  inline void setX(REAL _x) { d_data[0] = _x; }
  inline void setY(REAL _y) { d_data[1] = _y; }
  inline void setZ(REAL _z) { d_data[2] = _z; }
  inline void set(REAL _x, REAL _y, REAL _z) 
  {
    d_data[0] = _x; d_data[1] = _y; d_data[2] = _z;
  }

  inline
  void set(Vec v) 
  {
    d_data[0] = v.d_data[0]; d_data[1] = v.d_data[1]; d_data[2] = v.d_data[2];
  }

  inline double& operator[](int idx) {
    return d_data[idx];
  }

  inline double operator[](int idx) const {
    return d_data[idx];
  }

  inline
  bool operator==(const Vec v) const
  {
    return d_data[0] == v.d_data[0] && 
           d_data[1] == v.d_data[1] && d_data[2] == v.d_data[2];
  }

  inline
  bool operator==(const REAL val)
  {
    return d_data[0] == val && d_data[1] == val && d_data[2] == val;
  }

  inline
  bool operator!=(const Vec v)
  {
    return d_data[0] != v.d_data[0] || 
           d_data[1] != v.d_data[1] || d_data[2] != v.d_data[2];
  }

  inline
  void operator+=(const Vec v)
  {
    d_data[0] += v.d_data[0]; d_data[1] += v.d_data[1]; d_data[2] += v.d_data[2];
  }

  inline
  void operator-=(const Vec v)
  {
    d_data[0] -= v.d_data[0]; d_data[1] -= v.d_data[1]; d_data[2] -= v.d_data[2];
  }

  inline
  void operator*=(REAL val)
  {
    d_data[0] *= val; d_data[1] *= val; d_data[2] *= val;
  }

  inline
  void operator/=(REAL val)
  {
    d_data[0] /= val; d_data[1] /= val; d_data[2] /= val;
  }

  inline
  Vec operator+(Vec v) const
  {
    return Vec(d_data[0] + v.d_data[0], d_data[1] + v.d_data[1], 
               d_data[2] + v.d_data[2]);
  }

  inline
  Vec operator-(Vec v) const
  {
    return Vec(d_data[0] - v.d_data[0], d_data[1] - v.d_data[1], 
               d_data[2] - v.d_data[2]);
  }

  //Vec operator%(Vec p) const;  // cross product of this vector and p
  //REAL operator*(Vec p) const; // dot product of this vector and p
  inline
  Vec operator*(REAL d) const
  {
    return Vec(d_data[0] * d, d_data[1] * d, d_data[2] * d);
  }


  // Divides a vector by an int vector and produces a new vector
  Vec operator/(const IntVec& intVec) const;

  // Multiplies a vector by an int vector and produces a new vector
  Vec operator*(const IntVec& intVec) const;

  // Divides a vector by a real vector component-wise and produces a new vector
  Vec operator/(const Vec& intVec) const;

  // Multiplies a vector by a real vector component-wise and produces a new vector
  Vec operator*(const Vec& vec) const;

  inline
  REAL lengthSq() const
  {
    return d_data[0] * d_data[0] + d_data[1] * d_data[1] + d_data[2] * d_data[2];
  }

  inline
  REAL length() const
  {
    return std::sqrt(lengthSq());
  }

  // Normalize in place
  inline
  void normalizeInPlace()
  {
    REAL length = this->length();
    if (length < EPS) return;
    d_data[0] /= length;
    d_data[1] /= length;
    d_data[2] /= length;
  }

  void print(std::ostream& ofs) const;

  // Creates a Vec from a string that looks like "[1.0, -2.1e-3, 3]".
  static Vec fromString(const std::string& str);

  // Divides a vector by an int vector in place
  friend Vec& operator/=(Vec& realVec, const IntVec& intVec);

  // Multiplies a vector by an int vector in place
  friend Vec& operator*=(Vec& realVec, const IntVec& intVec);

  friend std::ostream& operator<<(std::ostream& os, const Vec& vec);

  friend inline
  REAL 
  dot(const Vec& v1, const Vec& v2)
  {
    return (v1.d_data[0] * v2.d_data[0] + v1.d_data[1] * v2.d_data[1] + 
            v1.d_data[2] * v2.d_data[2]);
  }

  friend inline 
  Vec
  cross(const Vec& v1, const Vec& v2)
  {
    return Vec(v1.d_data[1] * v2.d_data[2] - v1.d_data[2] * v2.d_data[1], 
              v1.d_data[2] * v2.d_data[0] - v1.d_data[0] * v2.d_data[2],
              v1.d_data[0] * v2.d_data[1] - v1.d_data[1] * v2.d_data[0]);
  }

  friend inline 
  Vec vcos(Vec v)
  {
    return Vec(cos(v.d_data[0]), cos(v.d_data[1]), cos(v.d_data[2]));
  }

  friend inline 
  Vec vacos(Vec v)
  {
    return Vec(acos(v.d_data[0]), acos(v.d_data[1]), acos(v.d_data[2]));
  }

  friend inline 
  Vec operator*(REAL d, Vec v)
  {
    return Vec(v.d_data[0] * d, v.d_data[1] * d, v.d_data[2] * d);
  }

  friend inline 
  Vec operator/(Vec v, REAL d)
  {
    return Vec(v.d_data[0] / d, v.d_data[1] / d, v.d_data[2] / d);
  }

  friend inline 
  REAL vnormL2(Vec v)
  {
    REAL x = v.d_data[0]; REAL y = v.d_data[1]; REAL z = v.d_data[2];
    return sqrt(x * x + y * y + z * z);
  }

  friend inline 
  Vec operator-(Vec v)
  {
    return -1.0 * v;
  }

  friend inline 
  Vec normalize(Vec v)
  {
    REAL alf = vnormL2(v);
    if (alf < EPS) // important, otherwise may cause numerical instability
      return v;
    return v / (vnormL2(v));
  }

private:

  //std::array<REAL, 3> d_data;
  REAL d_data[3];

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_data;
  }
};

Vec rotateVec(Vec v, Vec alf);

} // namespace dem

#endif
