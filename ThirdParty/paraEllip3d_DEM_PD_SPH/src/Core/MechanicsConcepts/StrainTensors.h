#ifndef CORE_MECHANICS_CONCEPTS_STRAIN_TENSORS_H
#define CORE_MECHANICS_CONCEPTS_STRAIN_TENSORS_H

#include <Core/Math/Matrix3.h>
#include <iostream>
#include <vector>

namespace dem {

struct SmallStrain
{
private:
  std::array<REAL, 6> _value;

public:
  SmallStrain(const Matrix3& strain)
  {
    _value[0] = strain(1,1);
    _value[1] = strain(2,2);
    _value[2] = strain(3,3);
    _value[3] = strain(2,3);
    _value[4] = strain(1,3);
    _value[5] = strain(1,2);
  }

  SmallStrain(REAL epsilon11, REAL epsilon22, REAL epsilon33,
              REAL epsilon23, REAL epsilon13, REAL epsilon12)
  {
    _value[0] = epsilon11;
    _value[1] = epsilon22;
    _value[2] = epsilon33;
    _value[3] = epsilon23;
    _value[4] = epsilon13;
    _value[5] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const SmallStrain& strain)
  {
    os << "epsilon11: " << strain._value[0] << "epsilon22: " << strain._value[1] 
       << "epsilon33: " << strain._value[2] << "epsilon23: " << strain._value[3] 
       << "epsilon13: " << strain._value[4] << "epsilon12: " << strain._value[5];
    return os;
  }
};

struct PlaneStressStrain
{
private:
  std::array<REAL, 4> _value;

public:
  PlaneStressStrain(const Matrix3& epsilon)
  {
    _value[0] = epsilon(1,1);
    _value[1] = epsilon(2,2);
    _value[2] = epsilon(3,3);
    _value[3] = epsilon(1,2);
  }

  PlaneStressStrain(REAL epsilon11, REAL epsilon22, REAL epsilon33, REAL epsilon12)
  {
    _value[0] = epsilon11;
    _value[1] = epsilon22;
    _value[2] = epsilon33;
    _value[3] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStressStrain& strain)
  {
    os << "epsilon11: " << strain._value[0] << "epsilon22: " << strain._value[1] 
       << "epsilon33: " << strain._value[2] << "epsilon12: " << strain._value[3];
    return os;
  }
};

struct PlaneStrain
{
private:
  std::array<REAL, 3> _value;

public:
  PlaneStrain(const Matrix3& epsilon)
  {
    _value[0] = epsilon(1,1);
    _value[1] = epsilon(2,2);
    _value[2] = epsilon(1,2);
  }

  PlaneStrain(REAL epsilon11, REAL epsilon22, REAL epsilon12)
  {
    _value[0] = epsilon11;
    _value[1] = epsilon22;
    _value[2] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStrain& strain)
  {
    os << "epsilon11: " << strain._value[0] << "epsilon22: " << strain._value[1] 
       << "epsilon12: " << strain._value[2];
    return os;
  }
};

struct AxisymmetricStrain
{
private:
  std::array<REAL, 4> _value;

public:
  AxisymmetricStrain() = default;

  AxisymmetricStrain(const Matrix3& epsilon)
  {
    _value[0] = epsilon(1,1);
    _value[1] = epsilon(2,2);
    _value[2] = epsilon(3,3);
    _value[3] = epsilon(1,2);
  }

  AxisymmetricStrain(REAL epsilonrr, REAL epsilonzz, REAL epsilontt, REAL epsilonrz)
  {
    _value[0] = epsilonrr;
    _value[1] = epsilonzz;
    _value[2] = epsilontt;
    _value[3] = epsilonrz;
  }

  inline
  AxisymmetricStrain operator*(REAL scalar) const 
  {
    return AxisymmetricStrain(_value[0]*scalar, _value[1]*scalar, 
                              _value[2]*scalar, _value[3]*scalar);
  }

  inline
  const std::array<REAL, 4>& data() const
  {
    return _value;
  }

  friend 
  AxisymmetricStrain operator*(REAL scalar, const AxisymmetricStrain& strain) 
  {
    return strain*scalar;
  }

  friend 
  AxisymmetricStrain operator+(const AxisymmetricStrain& eps1, 
                               const AxisymmetricStrain& eps2) 
  {
    return AxisymmetricStrain(eps1._value[0] + eps2._value[0],
                              eps1._value[1] + eps2._value[1],
                              eps1._value[2] + eps2._value[2],
                              eps1._value[3] + eps2._value[3]);
  }

  friend std::ostream& operator<<(std::ostream& os, const AxisymmetricStrain& strain)
  {
    os << "epsilon(rr, zz, tt, rz): " << strain._value[0] 
       << " " << strain._value[1] 
       << " " << strain._value[2] << " " << strain._value[3];
    return os;
  }
};
} // end namespace dem

#endif
