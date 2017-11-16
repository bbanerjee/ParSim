#ifndef CORE_MECHANICS_CONCEPTS_DEFORMATIONS_H
#define CORE_MECHANICS_CONCEPTS_DEFORMATIONS_H

#include <Core/Math/Matrix3.h>
#include <Core/Math/Vec.h>
#include <iostream>
#include <vector>

namespace dem {

struct DeformationGradient
{
  std::array<REAL, 9> _value;

  DeformationGradient() = default;

  DeformationGradient(const Matrix3& F)
  {
    _value[0] = F(1,1); _value[1] = F(1,2); _value[2] = F(1,3);
    _value[3] = F(2,1); _value[4] = F(2,2); _value[5] = F(2,3);
    _value[6] = F(3,1); _value[7] = F(3,2); _value[8] = F(3,3);
  }

  DeformationGradient(const std::array<REAL, 9>& F)
  {
    _value = F;
  }

  DeformationGradient(REAL F11, REAL F12, REAL F13,
                      REAL F21, REAL F22, REAL F23,
                      REAL F31, REAL F32, REAL F33)
  {
    _value[0] = F11; _value[1] = F12; _value[2] = F13;
    _value[3] = F21; _value[4] = F22; _value[5] = F23;
    _value[6] = F31; _value[7] = F32; _value[8] = F33;
  }

  inline
  DeformationGradient operator*(REAL scalar) const 
  {
    std::array<REAL, 9> val = _value;
    for (auto& v : val) v *= scalar;
    return DeformationGradient(val);
  }

  inline
  std::array<REAL, 9> data() const
  {
    return _value;
  }

  friend 
  DeformationGradient operator*(REAL scalar, const DeformationGradient& F) 
  {
    return F*scalar;
  }

  friend 
  DeformationGradient operator+(const DeformationGradient& F1, 
                                const DeformationGradient& F2) 
  {
    std::array<REAL, 9> val;
    for (auto ii = 0u; ii < 9; ++ii) val[ii] = F1._value[ii] + F2._value[ii];
    return DeformationGradient(val);
  }

  friend std::ostream& operator<<(std::ostream& os, const DeformationGradient& F)
  {
    os << "F11: " << F._value[0] << " F12: " << F._value[1] 
       << " F13: " << F._value[2] << " F21: " << F._value[3] 
       << " F22: " << F._value[4] << " F23: " << F._value[5]
       << " F31: " << F._value[6] << " F32: " << F._value[7]
       << " F33: " << F._value[8];
    return os;
  }
};

struct Displacement
{
private:
  REAL _u, _v, _w;

public:
  Displacement() = default;

  Displacement(const Vec& disp)
   : _u(disp.x()), _v(disp.y()), _w(disp.z())
  {
  }

  Displacement(REAL u1, REAL u2, REAL u3) 
   : _u(u1), _v(u2), _w(u3)
  {
  }

  inline
  Displacement operator*(REAL scalar) const 
  {
    return Displacement(_u*scalar, _v*scalar, _w*scalar);
  }

  inline
  Vec data() const
  {
    return Vec(_u, _v, _w);
  }

  friend 
  Displacement operator*(REAL scalar, const Displacement& disp) 
  {
    return disp*scalar;
  }

  friend 
  Displacement operator+(const Displacement& u1, const Displacement& u2) 
  {
    return Displacement(u1._u + u2._u, u1._v + u2._v, u1._w + u2._w);
  }

  friend std::ostream& operator<<(std::ostream& os, const Displacement& disp)
  {
    os << "u: (" << disp._u << " " << disp._v << " " << disp._w << ")";
    return os;
  }
};

struct Velocity
{
  std::array<REAL, 3> _value;

  Velocity(const Vec& vel)
  {
    _value[0] = vel.x();
    _value[1] = vel.y();
    _value[2] = vel.z();
  }

  Velocity(REAL v1, REAL v2, REAL v3)
  {
    _value[0] = v1;
    _value[1] = v2;
    _value[2] = v3;
  }

  friend std::ostream& operator<<(std::ostream& os, const Velocity& vel)
  {
    os << "v: (" << vel._value[0] << " " << vel._value[1] 
       << " " << vel._value[2] << ")";
    return os;
  }
};

struct Acceleration
{
  std::array<REAL, 3> _value;

  Acceleration(const Vec& acc)
  {
    _value[0] = acc.x();
    _value[1] = acc.y();
    _value[2] = acc.z();
  }

  Acceleration(REAL a1, REAL a2, REAL a3)
  {
    _value[0] = a1;
    _value[1] = a2;
    _value[2] = a3;
  }

  friend std::ostream& operator<<(std::ostream& os, const Acceleration& acc)
  {
    os << "a: (" << acc._value[0] << " " << acc._value[1] 
       << " " << acc._value[2] << ")";
    return os;
  }
};

} // end namespace dem

#endif
