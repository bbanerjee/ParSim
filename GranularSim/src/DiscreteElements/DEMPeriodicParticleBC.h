#ifndef ELLIP3D_PERIODIC_PARTICLE_BC_H
#define ELLIP3D_PERIODIC_PARTICLE_BC_H

#include <Boundary/BoundaryConditionCurve.h>

namespace dem {

template <typename T, int N>
class DEMPeriodicParticleBC
{

private:

  BoundaryConditionCurve<T, N> d_bcCurve;

public:

  DEMPeriodicParticleBC() = default;

  DEMPeriodicParticleBC(const XMLProblemSpec& ps)
  {
    d_bcCurve.read(ps);
  }

  void read(const XMLProblemSpec& ps) 
  {
    d_bcCurve.read(ps);
  }

  T getBCValue(double time) const
  {
    T value = d_bcCurve.getBCValue(time);
    return value;
  }

  T getBCRate(double time) const
  {
    T value = d_bcCurve.getBCRate(time);
    return value;
  }

  friend std::ostream& 
  operator<<(std::ostream& os, const DEMPeriodicParticleBC<T, N>& bc)
  {
    os << " BCs: \n" << bc.d_bcCurve;
    return os;
  }

};

} // namespace dem ends

#endif
