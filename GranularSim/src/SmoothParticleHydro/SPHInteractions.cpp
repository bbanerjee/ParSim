#include <SmoothParticleHydro/SPHInteractions.h>
#include <SmoothParticleHydro/SPHKernels.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>

using namespace sph;

using dem::Vec;
using dem::Box;

void
SPHParticleData::setData(const SPHParticleP& particle) {
  pos = particle->currentPosition();
  rho = particle->density();
  mu = particle->getViscosity();
  press = particle->getPressure();
  R = (press >= 0) ? 0.006 : 0.6;
  p_over_rhoSq = press/(rho*rho);
  coeff = 1 + R;

  mass = particle->mass();
  vel = particle->getVelocity();

  isSet = true;
}

SPHInteractions::SPHInteractions()
{
  d_part_a.isSet = false;
  d_part_b_free.isSet = false;
  d_interact.dist_a_b_free = -1.0;
}

void 
SPHInteractions::setDataParticleA(const SPHParticleP& particle)
{
  d_part_a.setData(particle);
}

void 
SPHInteractions::setDataParticleB(const SPHParticleP& particle)
{
  d_part_b_free.setData(particle);
}

// If particle b is not free, we have to swap the two particles
// for the logic to work (assuming particle a is now free)
void 
SPHInteractions::swapParticles()
{
  SPHParticleData tmp(d_part_a);
  d_part_a = d_part_b_free;
  d_part_b_free = tmp;

  d_interact.pos_a_b_free = -d_interact.pos_a_b_free;
  d_interact.vel_b_free_a = d_interact.vel_a_b_free;
  d_interact.vel_a_b_free = -d_interact.vel_a_b_free;

  d_kernel.gradWb_free_a_b = d_kernel.gradWa_b_free_a;
  d_kernel.gradWa_b_free_a = -d_kernel.gradWa_b_free_a;
}

void
SPHInteractions::initializeInteractionData() {
  if (d_part_a.isSet && d_part_b_free.isSet) {
    d_interact.pos_a_b_free = d_part_a.pos - d_part_b_free.pos;
    d_interact.dist_a_b_free_sq = d_interact.pos_a_b_free.lengthSq();
    d_interact.dist_a_b_free = std::sqrt(d_interact.dist_a_b_free_sq);
    d_interact.vel_a_b_free = d_part_a.vel - d_part_a.vel;
    d_interact.vel_b_free_a = -d_interact.vel_a_b_free;
  }
}

bool
SPHInteractions::areFarApart(const REAL& kernelSize) const
{
  if (d_interact.dist_a_b_free > kernelSize) return true;
  return false;
}

template <int dim>
void
SPHInteractions::computeInteractionKernel(const REAL& smoothLength)
{
  SPHKernels kernelFactory;
  d_kernel.Wqmin = kernelFactory.minQuinticSplineKernel<dim>(smoothLength);
  d_kernel.Wa_b_free = kernelFactory.quinticSplineKernel<dim>(
    d_interact.dist_a_b_free, smoothLength);
  d_kernel.Wb_free_a = d_kernel.Wa_b_free;
  d_kernel.gradWa_b_free_a = kernelFactory.gradientQuinticSplineKernel<dim>(
    d_interact.pos_a_b_free, d_interact.dist_a_b_free, smoothLength); 
  d_kernel.gradWb_free_a_b = -d_kernel.gradWa_b_free_a;
  d_kernel.dWa_b_free_da = d_kernel.gradWa_b_free_a.length();
  d_kernel.dWb_free_a_db = d_kernel.dWa_b_free_da;
  auto Wratio = d_kernel.Wa_b_free/d_kernel.Wqmin;
  d_kernel.phi = Wratio * Wratio * Wratio * Wratio;
}

void
SPHInteractions::updateParticleCoeffs() 
{
  d_part_a.coeff = 1 + d_part_a.R * d_kernel.phi;
  d_part_b_free.coeff = 1 + d_part_b_free.R * d_kernel.phi;
}

void 
SPHInteractions::updateMomentumExchangeCoeffs(const REAL& smoothLength,
                                              const REAL& alpha,
                                              const REAL& soundSpeed)
{
  d_interact.vel_dot_gradW = dot(d_interact.vel_a_b_free, d_kernel.gradWa_b_free_a);
  d_interact.vel_radial = dot(d_interact.vel_a_b_free, d_interact.pos_a_b_free);

  d_interact.rho_av = (d_part_a.rho + d_part_a.rho)/2;
  REAL Gamma_ab = 0;
  if (d_interact.vel_radial < 0) {
    REAL mu_ab = smoothLength * d_interact.vel_radial /
      (d_interact.dist_a_b_free_sq + 0.01 * smoothLength * smoothLength);
    Gamma_ab = (-alpha*soundSpeed*mu_ab)/d_interact.rho_av;
  }
  d_interact.p_rhoSq_av = d_part_a.p_over_rhoSq*d_part_a.coeff + 
                          d_part_b_free.p_over_rhoSq*d_part_b_free.coeff + Gamma_ab;
}


// Calculate Vab using the method shown in Morris's paper, 1996
void 
SPHInteractions::doGhostBoundaryVelocityCorrection(const SPHParticleP& particle)
{
  // Note that dB is an approximation to the distance
  // from the ghost point to the tangent plane
  // used in Morris
  auto dem_part = particle->getDEMParticle();
  auto da = dem_part->shortestDistToBoundary(d_part_b_free.pos);
  auto dB = dem_part->shortestDistToBoundary(d_part_a.pos);
  doBoundaryVelocityCorrection(da, dB);
}

// Calculate Vab using the method shown in Morris's paper, 1996
void 
SPHInteractions::doDomainBoundaryVelocityCorrection(const Box& domain)
{
  Vec domain_min = domain.minCorner();
  Vec domain_max = domain.maxCorner();

  auto xmin = domain_min.x();
  auto xmax = domain_max.x();
  auto ymin = domain_min.y();
  auto ymax = domain_max.y();
  auto zmin = domain_min.z();
  auto zmax = domain_max.z();

  auto x_a = d_part_a.pos.x();
  auto y_a = d_part_a.pos.y();
  auto z_a = d_part_a.pos.z();
  auto x_b = d_part_b_free.pos.x();
  auto y_b = d_part_b_free.pos.y();
  auto z_b = d_part_b_free.pos.z();

  auto da = 0.0;
  auto dB = 0.0;
  if (x_a < xmin) { 
    da = x_b - xmin;
    dB = xmin - x_a;
  } else if (x_a > xmax) {
    da = xmax - x_b;
    dB = x_a - xmax;
  } else if (y_a < ymin) {
    da = y_b - ymin;
    dB = ymin - y_a;
  } else if (y_a > ymax) {
    da = ymax - y_b;
    dB = y_a - ymax;
  } else if (z_a < zmin) {
    da = z_b - zmin;
    dB = zmin - z_a;
  } else if (z_a > zmax) {
    da = zmax - z_b;
    dB = z_a - zmax;
  }
  doBoundaryVelocityCorrection(da, dB);
}

void 
SPHInteractions::doBoundaryVelocityCorrection(const REAL& da, const REAL& dB)
{
  auto beta = 1 + dB / da;
  if (beta > 2.0 || beta < 0 || std::isnan(beta)) {
    beta = 2.0;
  }
  d_interact.vel_b_free_a = beta * (d_part_b_free.vel - d_part_a.vel);
  d_interact.vel_a_b_free = -d_interact.vel_b_free_a;
}

void 
SPHInteractions::updateInteractionFreeFree(SPHParticleP& free_part_a,
                                           SPHParticleP& free_part_b,
                                           const REAL& epsilon) const
{
  free_part_a->incDensityRate(d_part_b_free.mass*d_interact.vel_dot_gradW);
  free_part_b->incDensityRate(d_part_a.mass*d_interact.vel_dot_gradW);

  auto dva_dt = -d_part_b_free.mass * d_interact.p_rhoSq_av * d_kernel.gradWa_b_free_a;
  auto dvb_dt = -d_part_a.mass * d_interact.p_rhoSq_av * d_kernel.gradWb_free_a_b;
  free_part_a->incAcceleration(dva_dt);
  free_part_b->incAcceleration(dvb_dt);

  auto delta_a = epsilon * d_part_b_free.mass * 
    (-d_interact.vel_a_b_free) * d_kernel.Wa_b_free / d_interact.rho_av;
  auto delta_b = epsilon * d_part_a.mass * 
    (d_interact.vel_a_b_free) * d_kernel.Wb_free_a / d_interact.rho_av;
  free_part_a->incVelocityCorrection(delta_a);
  free_part_b->incVelocityCorrection(delta_b);
}

void 
SPHInteractions::updateInteractionGhostFree(SPHParticleP& ghost_part_a,
                                            SPHParticleP& free_part_b,
                                            const REAL& epsilon) const
{
  free_part_b->incDensityRate(d_part_a.mass*d_interact.vel_dot_gradW);

  auto dva_dt = -d_part_b_free.mass * d_interact.p_rhoSq_av * d_kernel.gradWa_b_free_a;
  auto dvb_dt = -d_part_a.mass * d_interact.p_rhoSq_av * d_kernel.gradWb_free_a_b;
  free_part_b->incAcceleration(dvb_dt);

  auto dem_part = ghost_part_a->getDEMParticle();
  dem_part->addForce(d_part_a.mass * dva_dt);
  dem_part->addMoment(cross(d_part_a.pos - dem_part->currentPosition(),
                            d_part_a.mass * dva_dt));

  auto delta_b = epsilon * d_part_a.mass * 
    (d_interact.vel_a_b_free) * d_kernel.Wb_free_a / d_interact.rho_av;
  free_part_b->incVelocityCorrection(delta_b);
}

void 
SPHInteractions::updateInteractionBoundaryFree(SPHParticleP& boundary_part_a,
                                               SPHParticleP& free_part_b,
                                               const REAL& epsilon) const
{
  boundary_part_a->incDensityRate(d_part_b_free.mass*d_interact.vel_dot_gradW);
  free_part_b->incDensityRate(d_part_a.mass*d_interact.vel_dot_gradW);

  auto dvb_dt = -d_part_a.mass * d_interact.p_rhoSq_av * d_kernel.gradWb_free_a_b;
  free_part_b->incAcceleration(dvb_dt); 

  auto delta_b = epsilon * d_part_a.mass * 
    (d_interact.vel_a_b_free) * d_kernel.Wb_free_a / d_interact.rho_av;
  free_part_b->incVelocityCorrection(delta_b);
}

namespace sph {
  template void 
  SPHInteractions::computeInteractionKernel<1>(const REAL& smoothLength);
  template void 
  SPHInteractions::computeInteractionKernel<2>(const REAL& smoothLength);
  template void 
  SPHInteractions::computeInteractionKernel<3>(const REAL& smoothLength);
}