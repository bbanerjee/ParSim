/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_Prager.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

KinematicHardening_Prager::KinematicHardening_Prager(ProblemSpecP& ps)
{
  d_cm.beta = 1.0;
  ps->get("beta", d_cm.beta);
  ps->require("hardening_modulus", d_cm.hardening_modulus);
}

KinematicHardening_Prager::KinematicHardening_Prager(
  const KinematicHardening_Prager* cm)
{
  d_cm.beta = cm->d_cm.beta;
  d_cm.hardening_modulus = cm->d_cm.hardening_modulus;
}

KinematicHardening_Prager::~KinematicHardening_Prager() = default;

void
KinematicHardening_Prager::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP plastic_ps = ps->appendChild("kinematic_hardening_model");
  plastic_ps->setAttribute("type", "prager_hardening");

  plastic_ps->appendElement("beta", d_cm.beta);
  plastic_ps->appendElement("hardening_modulus", d_cm.hardening_modulus);
}

/* Assumes von Mises plasticity and an associated flow rule.  The back stress
is given by the rate equation D/Dt(beta) = 2/3~gammadot~Hprime~df/dsigma */
void
KinematicHardening_Prager::computeBackStress(
  const ModelStateBase*, const double& delT, const particleIndex idx,
  const double& delLambda, const Matrix3& df_dsigma_normal_new,
  const Matrix3& backStress_old, Matrix3& backStress_new)
{
  // Get the hardening modulus (constant for Prager kinematic hardening)
  double H_prime = d_cm.beta * d_cm.hardening_modulus;
  double stt = sqrt(2.0 / 3.0);

  // Compute updated backstress
  backStress_new =
    backStress_old + df_dsigma_normal_new * (delLambda * H_prime * stt);

  return;
}

void
KinematicHardening_Prager::eval_h_beta(const Matrix3& df_dsigma,
                                       const ModelStateBase*, Matrix3& h_beta)
{
  double H_prime = d_cm.beta * d_cm.hardening_modulus;
  h_beta = df_dsigma * (2.0 / 3.0 * H_prime);
  return;
}
