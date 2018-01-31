#ifndef CORE_MECHANICS_CONCEPTS_STRESS_TENSORS_H
#define CORE_MECHANICS_CONCEPTS_STRESS_TENSORS_H

#include <Core/Math/Matrix3.h>
#include <iostream>
#include <vector>

namespace dem {

struct CauchyStress
{
  std::array<REAL, 6> component;

  CauchyStress(const Matrix3& stress)
  {
    component[0] = stress(1,1);
    component[1] = stress(2,2);
    component[2] = stress(3,3);
    component[3] = stress(2,3);
    component[4] = stress(1,3);
    component[5] = stress(1,2);
  }

  CauchyStress(REAL sigma11, REAL sigma22, REAL sigma33,
               REAL sigma23, REAL sigma13, REAL sigma12)
  {
    component[0] = sigma11;
    component[1] = sigma22;
    component[2] = sigma33;
    component[3] = sigma23;
    component[4] = sigma13;
    component[5] = sigma12;
  }

  friend std::ostream& operator<<(std::ostream& os, const CauchyStress& stress)
  {
    os << "stress11: " << stress.component[0] << "stress22: " << stress.component[1] 
       << "stress33: " << stress.component[2] << "stress23: " << stress.component[3] 
       << "stress13: " << stress.component[4] << "stress12: " << stress.component[5];
    return os;
  }
};

struct PlaneStress
{
  std::array<REAL, 3> component;

  PlaneStress(const Matrix3& stress)
  {
    component[0] = stress(1,1);
    component[1] = stress(2,2);
    component[2] = stress(1,2);
  }

  PlaneStress(REAL sigma11, REAL sigma22, REAL sigma12)
  {
    component[0] = sigma11;
    component[1] = sigma22;
    component[2] = sigma12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStress& stress)
  {
    os << "stress11: " << stress.component[0] << "stress22: " << stress.component[1] 
       << "stress12: " << stress.component[2];
    return os;
  }
};

struct PlaneStrainStress
{
  std::array<REAL, 4> component;

  PlaneStrainStress(const Matrix3& stress)
  {
    component[0] = stress(1,1);
    component[1] = stress(2,2);
    component[2] = stress(3,3);
    component[3] = stress(1,2);
  }

  PlaneStrainStress(REAL sigma11, REAL sigma22, REAL sigma33, REAL sigma12)
  {
    component[0] = sigma11;
    component[1] = sigma22;
    component[2] = sigma33;
    component[3] = sigma12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStrainStress& stress)
  {
    os << "stress11: " << stress.component[0] << "stress22: " << stress.component[1] 
       << "stress33: " << stress.component[2] << "stress12: " << stress.component[3];
    return os;
  }
};

struct AxisymmetricStress
{
  std::array<REAL, 4> component;

  AxisymmetricStress(const Matrix3& stress)
  {
    component[0] = stress(1,1);
    component[1] = stress(2,2);
    component[2] = stress(3,3);
    component[3] = stress(1,2);
  }

  AxisymmetricStress(REAL sigmarr, REAL sigmazz, REAL sigmatt, REAL sigmarz)
  {
    component[0] = sigmarr;
    component[1] = sigmazz;
    component[2] = sigmatt;
    component[3] = sigmarz;
  }

  friend std::ostream& operator<<(std::ostream& os, const AxisymmetricStress& stress)
  {
    os << "stressrr: " << stress.component[0] << "stresszz: " << stress.component[1] 
       << "stresstt: " << stress.component[2] << "stressrz: " << stress.component[3];
    return os;
  }
};
} // end namespace dem

#endif
