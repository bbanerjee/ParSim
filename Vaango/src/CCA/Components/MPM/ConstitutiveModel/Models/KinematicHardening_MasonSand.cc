/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <Core/Exceptions/InternalError.h>
#include <cmath>

#define NUMERICALLY_INTEGRATE_BACKSTRESS

using namespace Uintah;
using namespace Vaango;

const Matrix3 KinematicHardening_MasonSand::Identity(1.0, 0.0, 0.0,
                                                     0.0, 1.0, 0.0,
                                                     0.0, 0.0, 1.0);

KinematicHardening_MasonSand::KinematicHardening_MasonSand(ProblemSpecP& ps,
                                                           InternalVariableModel* intvar)
{
  // Copy of internal variables object
  d_intvar = intvar;

  // Initial fluid pressure
  ps->require("fluid_pressure_initial", d_cm.fluid_pressure_initial);

  // Local labels
  initializeLocalMPMLabels();
}
         
KinematicHardening_MasonSand::KinematicHardening_MasonSand(const KinematicHardening_MasonSand* cm)
{
  // Copy of internal variables object
  d_intvar = cm->d_intvar;

  // Initial fluid pressure
  d_cm.fluid_pressure_initial = cm->d_cm.fluid_pressure_initial;

  // Local labels
  initializeLocalMPMLabels();
}
         
KinematicHardening_MasonSand::~KinematicHardening_MasonSand()
{
  //VarLabel::destroy(pZetaLabel);
  //VarLabel::destroy(pZetaLabel_preReloc);
}

void KinematicHardening_MasonSand::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP plastic_ps = ps->appendChild("kinematic_hardening_model");
  plastic_ps->setAttribute("type","mason_sand_pore_pressure");
  plastic_ps->appendElement("fluid_pressure_initial",d_cm.fluid_pressure_initial);
}

void 
KinematicHardening_MasonSand::computeBackStress(const ModelStateBase* state_input,
                                                Matrix3& backStress_new)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // If the state is tensile the back stress does not change. Return old backtress.
  if (state->I1 > 0.0) {
    backStress_new = Identity*(state->zeta/3.0);
    return;
  }

  // Get the variables of interest
  double p0            = d_cm.fluid_pressure_initial;
  double phi0          = state->phi0;
  double Sw0           = state->Sw0;
  double phat_old      = -state->zeta/3.0;
  double ep_v_bar_old  = -state->ep_v;
  double dep_v_bar     = -state->dep_v;

#ifndef NUMERICALLY_INTEGRATE_BACKSTRESS

  double phat_new = computePressureUnloaded(ev_p_bar_old, Sw0, phi0, p0);

#else
  // Compute volumetric strains in air, water, and matrix material at p = zeta
  double exp_ev = 0.0;
  double dexp_ev_air    = d_air.computeDerivExpElasticVolumetricStrain(phat_old, 0.0, exp_ev);
  double dexp_ev_water  = d_water.computeDerivExpElasticVolumetricStrain(phat_old, p0, exp_ev);
  double dexp_ev_matrix = d_granite.computeDerivExpElasticVolumetricStrain(phat_old, 0.0, exp_ev);

  // Compute denominator of rate equation
  double BB = phi0*((1.0 - Sw0)*dexp_ev_air + Sw0*dexp_ev_water) + (1-phi0)*dexp_ev_matrix;

  // Compute derivative of zeta wrt plastic vol strain
  double dzeta_dep_v = -std::exp(-ep_v_bar_old)/BB;

  // Compute new zeta using forward Euler
  double phat_new = phat_old + dzeta_dep_v*dep_v_bar;
  std::cout << "\t Backstress: " 
            << " air = " << dexp_ev_air
            << " water = " << dexp_ev_water
            << " matrix = " << dexp_ev_matrix << std::endl
            << "\t\t BB = " << BB
            << " dzeta/dep_v = " << dzeta_dep_v
            << " delta ep_v = " << dep_v_bar
            << " ep_v = " << ep_v_bar_old
            << " phat_old = " << phat_old
            << " delta zeta = " << dzeta_dep_v*dep_v_bar
            << " phat_new = " << phat_new << std::endl;
#endif

  // Compute backstress tensor
  backStress_new = -phat_new*Identity;  

  return;
}

//------------------------------------------------------
// Newton solve for pressure
//------------------------------------------------------
double
KinematicHardening_MasonSand::computePressureUnloaded(const double& ev_p, 
                                                      const double& S0, 
                                                      const double& phi0, 
                                                      const double& p0) 
{
  double maxIter = 100;
  double tol = 1.0e-6;
  //double p = p0;
  double p = 0.0;
  double k = 1;
  while (k < maxIter) {
    double g_Dg = computeGpByDgp(p, ev_p, S0, phi0, p0);
    p -= g_Dg;
    k++;
    if (abs(g_Dg) < tol) {
      break;
    }
  }

  if (k >= maxIter) {
    std::ostringstream out;
    out << "**ERROR** Too many iterations needed to find backstress.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  return p;
}

double
KinematicHardening_MasonSand::computeGpByDgp(const double& p, 
                                             const double& ev_p, 
                                             const double& S0, 
                                             const double& phi0, 
                                             const double& p0) 
{
  std::cout << "p = " << p << " ev_p = " << ev_p << std::endl;

  // Compute volumetric strains in air, water, and matrix material at p = zeta
  double exp_ev_air = 0.0;
  double exp_ev_water = 0.0;
  double exp_ev_matrix = 0.0;
  double dexp_ev_air    = d_air.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_air);
  double dexp_ev_water  = d_water.computeDerivExpElasticVolumetricStrain(p, p0, exp_ev_water);
  double dexp_ev_matrix = d_granite.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_matrix);

  double air_term = (1.0-S0)*phi0;
  double water_term = S0*phi0;
  double solid_term = (1.0-phi0);

  double Gp = air_term*exp_ev_air + water_term*exp_ev_water + solid_term*exp_ev_matrix 
              - std::exp(-ev_p);

  double dGp = air_term*dexp_ev_air + water_term*dexp_ev_water + solid_term*dexp_ev_matrix;

  double Gp_by_dGp = Gp/dGp;
  if (!std::isfinite(Gp_by_dGp)) {
    std::ostringstream out;
    out << "**ERROR** Increment in Newton iteration is not finite."
        << " p = " << p << " Gp = " << Gp << " dGp = " << dGp << std::endl;
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  return Gp_by_dGp;
}

