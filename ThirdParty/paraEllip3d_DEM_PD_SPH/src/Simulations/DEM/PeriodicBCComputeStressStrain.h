#ifndef ELLIP3D_PERIODIC_BC_COMPUTE_STRESS_STRAIN_H
#define ELLIP3D_PERIODIC_BC_COMPUTE_STRESS_STRAIN_H

#include <Simulations/Command.h>
#include <DiscreteElements/DEMContainers.h>
#include <vector>
#include <ostream>

namespace dem {
class PeriodicBCComputeStressStrain : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in PeriodicBCComputeStressStrain.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* sph)
  {
    std::cout << "**ERROR** Execute with DEM + SPH. "
              << "Should not be called in PeriodicBCComputeStressStrain.\n";
  }

private:

  DEMTetrahedronPArray createTessellation(const DEMParticlePArray& particles);

  void computeBoundaryTractions();

};
}

#endif
