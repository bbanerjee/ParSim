#ifndef VEC_H
#define VEC_H

#include <Core/Types/realtypes.h>
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>

namespace dem {
  
  class Vec {
    
  public:
    Vec() : x(0), y(0), z(0) {}
    Vec(REAL d) : x(d), y(d), z(d) {}
    Vec(REAL _x, REAL _y, REAL _z) : x(_x),y(_y),z(_z) {}
    REAL getX() const {return x;}
    REAL getY() const {return y;}
    REAL getZ() const {return z;}
    void setX(REAL _x) {x = _x;}
    void setY(REAL _y) {y = _y;}
    void setZ(REAL _z) {z = _z;}
    void set(REAL _x, REAL _y, REAL _z) {x = _x; y = _y; z = _z;}
    void set(Vec v) {x = v.getX(); y = v.getY(); z = v.getZ();}
    
    bool operator==(const Vec v);
    bool operator==(const REAL d);   
    bool operator!=(const Vec v); 
    void operator+=(const Vec v);
    void operator-=(const Vec v);
    void operator*=(REAL d);
    void operator/=(REAL d);
    Vec  operator+(Vec v) const;
    Vec  operator-(Vec v) const;
    Vec  operator%(Vec p) const;   // cross product of this vector and p
    REAL operator*(Vec p) const;   // dot product of this vector and p
    Vec  operator*(REAL d) const;
    void print(std::ostream &ofs) const;
    
  private:
    REAL x;
    REAL y;
    REAL z;
    
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & x;
      ar & y;
      ar & z;
    }
    
  };
  
  // non-member functions
  Vec  operator*(REAL d, Vec v);
  Vec  operator/(Vec v, REAL d);
  Vec  operator-(Vec v);
  REAL vfabs(Vec v);
  Vec  vcos(Vec v);
  Vec  vacos(Vec v);
  Vec  rotateVec(Vec v, Vec alf);
  Vec  normalize(Vec v);
  
} // namespace dem

#endif
