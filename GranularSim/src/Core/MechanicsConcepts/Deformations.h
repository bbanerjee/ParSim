#ifndef CORE_MECHANICS_CONCEPTS_DEFORMATIONS_H
#define CORE_MECHANICS_CONCEPTS_DEFORMATIONS_H

#include <Core/Math/Matrix3.h>
#include <Core/Math/Vec.h>
#include <iostream>
#include <vector>
#include <array>

namespace dem {

class DeformationGradient
{
private:
  Matrix3 _value;

public:

  DeformationGradient()  = default;
  DeformationGradient(REAL val) 
  {
    _value = val;
  }

  DeformationGradient(const Matrix3& F)
  {
    _value = F;
  }

  DeformationGradient(const std::array<REAL, 9>& F)
  {
    _value(0,0) = F[0]; _value(0,1) = F[1]; _value(0,2) = F[2];
    _value(1,0) = F[3]; _value(1,1) = F[4]; _value(1,2) = F[5];
    _value(2,0) = F[6]; _value(2,1) = F[7]; _value(2,2) = F[8];
  }

  DeformationGradient(REAL F11, REAL F12, REAL F13,
                      REAL F21, REAL F22, REAL F23,
                      REAL F31, REAL F32, REAL F33)
  {
    _value(0,0) = F11; _value(0,1) = F12; _value(0,2) = F13;
    _value(1,0) = F21; _value(1,1) = F22; _value(1,2) = F23;
    _value(2,0) = F31; _value(2,1) = F32; _value(2,2) = F33;
  }

  inline
  DeformationGradient operator-(const DeformationGradient& F) const 
  {
    return DeformationGradient(_value - F._value); 
  }

  inline
  DeformationGradient operator*(REAL scalar) const 
  {
    return DeformationGradient(_value * scalar);
  }

  inline
  const Matrix3& data() const
  {
    return _value;
  }

  DeformationGradient inverse() const
  {
    return _value.Inverse();
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
    return F1._value + F2._value;
  }

  inline 
  double operator()(int i, int j) const
  {
    return _value(i,j);
  }

  //inline 
  //double& operator()(int i, int j)
  //{
  //  return _value(i,j);
  //}

  friend std::ostream& operator<<(std::ostream& os, const DeformationGradient& F)
  {
    os << "(F11, F12, F13, F21, F22, F23, F31, F32, F33): " 
       << " " << F._value(0,0) << " " << F._value(0,1) 
       << " " << F._value(0,2) << " " << F._value(1,0) 
       << " " << F._value(1,1) << " " << F._value(1,2)
       << " " << F._value(2,0) << " " << F._value(2,1)
       << " " << F._value(2,2);
    return os;
  }
};

using DeformationGradientRate = DeformationGradient;

class Displacement
{
private:
  std::array<REAL, 3> _value;

public:
  Displacement() = default;
  Displacement(REAL val) 
  {
    _value[0] = val; _value[1] = val; _value[2] = val;
  }

  Displacement(const Vec& disp)
  {
   _value[0] = disp.x(); _value[1] = disp.y(); _value[2] = disp.z();
  }

  Displacement(REAL u1, REAL u2, REAL u3) 
  {
   _value[0] = u1; _value[1] = u2; _value[2] = u3;
  }

  inline
  Displacement operator-(const Displacement& disp) const 
  {
    return Displacement(_value[0] - disp._value[0], 
                        _value[1] - disp._value[1], 
                        _value[2] - disp._value[2]);
  }

  inline
  Displacement operator*(REAL scalar) const 
  {
    return Displacement(_value[0]*scalar, _value[1]*scalar, _value[2]*scalar);
  }

  inline
  const std::array<REAL, 3>& data() const
  {
    return _value;
  }

  inline 
  double operator[](int i) const
  {
    return _value[i];
  }

  //inline 
  //double& operator[](int i) 
  //{
  //  return _value[i];
  //}

  friend 
  Displacement operator*(REAL scalar, const Displacement& disp) 
  {
    return disp*scalar;
  }

  friend 
  Displacement operator+(const Displacement& u1, const Displacement& u2) 
  {
    return Displacement(u1._value[0] + u2._value[0], 
                        u1._value[1] + u2._value[1], 
                        u1._value[2] + u2._value[2]);
  }

  friend std::ostream& operator<<(std::ostream& os, const Displacement& disp)
  {
    os << "u: (" << disp._value[0] << " " << disp._value[1] << " " << disp._value[2] << ")";
    return os;
  }
};

using DisplacementRate = Displacement;

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
