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

  // Get the variables of interest
  double p0 = d_cm.fluid_pressure_initial;
  double I1            = state->I1;
  double phi0          = state->phi0;
  double Sw0           = state->Sw0;
  double zeta_bar_old  = -state->zeta;
  double ep_v_bar_old  = -state->ep_v;
  double dep_v_bar     = -state->dep_v;

  // If the state is tensile the back stress does not change. Return old backtress.
  if (I1 > 0.0) {
    backStress_new = -zeta_bar_old*Identity;
    return;
  }

  // Compute volumetric strains in air, water, and matrix material at p = zeta
  double dexp_ev_air    = d_air.computeDerivExpElasticVolumetricStrain(zeta_bar_old, 0.0);
  double dexp_ev_water  = d_water.computeDerivExpElasticVolumetricStrain(zeta_bar_old, p0);
  double dexp_ev_matrix = d_granite.computeDerivExpElasticVolumetricStrain(zeta_bar_old, 0.0);

  // Compute denominator of rate equation
  double BB = phi0*((1.0 - Sw0)*dexp_ev_air + Sw0*dexp_ev_water) + (1-phi0)*dexp_ev_matrix;

  // Compute derivative of zeta wrt plastic vol strain
  double dzeta_dep_v = std::exp(-ep_v_bar_old)/BB;

  // Compute new zeta using forward Euler
  double zeta_bar_new = zeta_bar_old + dzeta_dep_v*dep_v_bar;
 
  // Compute backstress tensor
  backStress_new = -zeta_bar_new*Identity;  

  return;
}

