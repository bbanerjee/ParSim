#ifndef CORE_MECHANICS_CONCEPTS_STRAIN_TENSORS_H
#define CORE_MECHANICS_CONCEPTS_STRAIN_TENSORS_H

#include <Core/Math/Matrix3.h>
#include <iostream>
#include <vector>

namespace dem {

struct SmallStrain
{
  std::array<REAL, 6> component;

  SmallStrain(const Matrix3& strain)
  {
    component[0] = strain(1,1);
    component[1] = strain(2,2);
    component[2] = strain(3,3);
    component[3] = strain(2,3);
    component[4] = strain(1,3);
    component[5] = strain(1,2);
  }

  SmallStrain(REAL epsilon11, REAL epsilon22, REAL epsilon33,
              REAL epsilon23, REAL epsilon13, REAL epsilon12)
  {
    component[0] = epsilon11;
    component[1] = epsilon22;
    component[2] = epsilon33;
    component[3] = epsilon23;
    component[4] = epsilon13;
    component[5] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const SmallStrain& strain)
  {
    os << "epsilon11: " << strain.component[0] << "epsilon22: " << strain.component[1] 
       << "epsilon33: " << strain.component[2] << "epsilon23: " << strain.component[3] 
       << "epsilon13: " << strain.component[4] << "epsilon12: " << strain.component[5];
    return os;
  }
};

struct PlaneStressStrain
{
  std::array<REAL, 4> component;

  PlaneStressStrain(const Matrix3& epsilon)
  {
    component[0] = epsilon(1,1);
    component[1] = epsilon(2,2);
    component[2] = epsilon(3,3);
    component[3] = epsilon(1,2);
  }

  PlaneStressStrain(REAL epsilon11, REAL epsilon22, REAL epsilon33, REAL epsilon12)
  {
    component[0] = epsilon11;
    component[1] = epsilon22;
    component[2] = epsilon33;
    component[3] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStressStrain& strain)
  {
    os << "epsilon11: " << strain.component[0] << "epsilon22: " << strain.component[1] 
       << "epsilon33: " << strain.component[2] << "epsilon12: " << strain.component[3];
    return os;
  }
};

struct PlaneStrain
{
  std::array<REAL, 3> component;

  PlaneStrain(const Matrix3& epsilon)
  {
    component[0] = epsilon(1,1);
    component[1] = epsilon(2,2);
    component[2] = epsilon(1,2);
  }

  PlaneStrain(REAL epsilon11, REAL epsilon22, REAL epsilon12)
  {
    component[0] = epsilon11;
    component[1] = epsilon22;
    component[2] = epsilon12;
  }

  friend std::ostream& operator<<(std::ostream& os, const PlaneStrain& strain)
  {
    os << "epsilon11: " << strain.component[0] << "epsilon22: " << strain.component[1] 
       << "epsilon12: " << strain.component[2];
    return os;
  }
};

struct AxisymmetricStrain
{
  std::array<REAL, 4> component;

  AxisymmetricStrain(const Matrix3& epsilon)
  {
    component[0] = epsilon(1,1);
    component[1] = epsilon(2,2);
    component[2] = epsilon(3,3);
    component[3] = epsilon(1,2);
  }

  AxisymmetricStrain(REAL epsilonrr, REAL epsilonzz, REAL epsilontt, REAL epsilonrz)
  {
    component[0] = epsilonrr;
    component[1] = epsilonzz;
    component[2] = epsilontt;
    component[3] = epsilonrz;
  }

  friend std::ostream& operator<<(std::ostream& os, const AxisymmetricStrain& strain)
  {
    os << "epsilonrr: " << strain.component[0] << "epsilonzz: " << strain.component[1] 
       << "epsilontt: " << strain.component[2] << "epsilonrz: " << strain.component[3];
    return os;
  }
};
} // end namespace dem

#endif
