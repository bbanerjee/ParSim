#ifndef SPH_INTERACTIONS_H
#define SPH_INTERACTIONS_H

#include <SmoothParticleHydro/SPHParticle.h>
#include <SmoothParticleHydro/SPHContainers.h>

#include <Core/Geometry/Box.h>
#include <Core/Math/Vec.h>
#include <Core/Types/integertypes.h>
#include <Core/Types/realtypes.h>

#include <iostream>
#include <math.h>
#include <vector>

namespace sph {

struct SPHParticleData {
  bool     isSet;
  REAL     rho;
  REAL     mu;
  REAL     press;
  REAL     R;
  REAL     p_over_rhoSq;
  REAL     mass;
  REAL     coeff;
  dem::Vec pos;
  dem::Vec vel;

  void setData(const SPHParticleP& particle);
};

/**
 *  Particle b is always free.  Particle a can be free, ghost, or boundary
 */
struct SPHInteractionData {
  REAL     dist_a_b_free;
  REAL     dist_a_b_free_sq;
  REAL     vel_dot_gradW;
  REAL     vel_radial;
  REAL     rho_av;
  REAL     p_rhoSq_av;
  dem::Vec pos_a_b_free;
  dem::Vec vel_a_b_free;
  dem::Vec vel_b_free_a;
};

struct SPHInteractionKernel {
  REAL     Wqmin;
  REAL     Wa_b_free;
  REAL     Wb_free_a;
  dem::Vec gradWa_b_free_a;
  dem::Vec gradWb_free_a_b;
  REAL     dWa_b_free_da;
  REAL     dWb_free_a_db;
  REAL     phi;
};

class SPHInteractions
{
  public:

  SPHInteractions();
  SPHInteractions (const SPHInteractions&) = delete;
  SPHInteractions& operator= (const SPHInteractions&) = delete;

  void setDataParticleA(const SPHParticleP& particle);
  void setDataParticleB(const SPHParticleP& particle);
  void swapParticles();

  void initializeInteractionData();

  template <int dim>
  void computeInteractionKernel(const REAL& smoothLength);

  void updateParticleCoeffs();
  void updateMomentumExchangeCoeffs(const REAL& smoothLength,
                                    const REAL& alpha,
                                    const REAL& soundSpeed);

  void doGhostBoundaryVelocityCorrection(const SPHParticleP& ghostParticle);
  void doDomainBoundaryVelocityCorrection(const dem::Box& domain);

  bool areFarApart(const REAL& kernelSize) const;

  void updateInteractionFreeFree(SPHParticleP& free_part_a,
                                 SPHParticleP& free_part_b,
                                 const REAL& epsilon) const;

  void updateInteractionGhostFree(SPHParticleP& ghost_part_a,
                                  SPHParticleP& free_part_b,
                                  const REAL& epsilon) const;

  void updateInteractionBoundaryFree(SPHParticleP& boundary_part_a,
                                     SPHParticleP& free_part_b,
                                     const REAL& epsilon) const;

  private:

  SPHParticleData d_part_a;
  SPHParticleData d_part_b_free;
  SPHInteractionData d_interact;
  SPHInteractionKernel d_kernel;

  void doBoundaryVelocityCorrection(const REAL& da, const REAL& dB);
}; 

} // end namespace sph

#endif
