#ifndef CORE_MECHANICS_CONCEPTS_FORCES_H
#define CORE_MECHANICS_CONCEPTS_FORCES_H

#include <Core/Math/Vec.h>
#include <iostream>
#include <vector>

namespace dem {

struct Traction
{
  std::array<REAL, 3> component;

  Traction(const Vec& traction)
  {
    component[0] = traction.x();
    component[1] = traction.y();
    component[2] = traction.z();
  }

  Traction(REAL t1, REAL t2, REAL t3)
  {
    component[0] = t1;
    component[1] = t2;
    component[2] = t3;
  }

  friend std::ostream& operator<<(std::ostream& os, const Traction& traction)
  {
    os << "t1: " << traction.component[0] << "t2: " << traction.component[1] 
       << "t3: " << traction.component[2] ;
    return os;
  }
};

struct Force
{
  std::array<REAL, 3> component;

  Force(const Vec& force)
  {
    component[0] = force.x();
    component[1] = force.y();
    component[2] = force.z();
  }

  Force(REAL f1, REAL f2, REAL f3)
  {
    component[0] = f1;
    component[1] = f2;
    component[2] = f3;
  }

  friend std::ostream& operator<<(std::ostream& os, const Force& force)
  {
    os << "f1: " << force.component[0] << "f2: " << force.component[1] 
       << "f3: " << force.component[2];
    return os;
  }
};

struct Moment
{
  std::array<REAL, 3> component;

  Moment(const Vec& moment)
  {
    component[0] = moment.x();
    component[1] = moment.y();
    component[2] = moment.z();
  }

  Moment(REAL m1, REAL m2, REAL m3)
  {
    component[0] = m1;
    component[1] = m2;
    component[2] = m3;
  }

  friend std::ostream& operator<<(std::ostream& os, const Moment& moment)
  {
    os << "m1: " << moment.component[0] << "m2: " << moment.component[1] 
       << "m3: " << moment.component[2];
    return os;
  }
};

} // end namespace dem

#endif
