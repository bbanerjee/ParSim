/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardening_Armstrong.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/InternalError.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

KinematicHardening_Armstrong::KinematicHardening_Armstrong(ProblemSpecP& ps)
{
  d_cm.beta = 1.0;
  ps->get("beta", d_cm.beta);
  ps->require("hardening_modulus_1", d_cm.hardening_modulus_1);
  ps->require("hardening_modulus_2", d_cm.hardening_modulus_2);
}

KinematicHardening_Armstrong::KinematicHardening_Armstrong(
  const KinematicHardening_Armstrong* cm)
{
  d_cm.beta = cm->d_cm.beta;
  d_cm.hardening_modulus_1 = cm->d_cm.hardening_modulus_1;
  d_cm.hardening_modulus_2 = cm->d_cm.hardening_modulus_2;
}

KinematicHardening_Armstrong::~KinematicHardening_Armstrong() = default;

void
KinematicHardening_Armstrong::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP plastic_ps = ps->appendChild("kinematic_hardening_model");
  plastic_ps->setAttribute("type", "armstrong_frederick_hardening");

  plastic_ps->appendElement("beta", d_cm.beta);
  plastic_ps->appendElement("hardening_modulus_1", d_cm.hardening_modulus_1);
  plastic_ps->appendElement("hardening_modulus_2", d_cm.hardening_modulus_2);
}

void
KinematicHardening_Armstrong::computeBackStress(
  const ModelStateBase*, [[maybe_unused]] const double& delT, [[maybe_unused]] const particleIndex idx,
  const double& delLambda, const Matrix3& df_dsigma_normal_new,
  const Matrix3& backStress_old, Matrix3& backStress_new)
{
  // Get the hardening modulus
  double H_1 = d_cm.beta * d_cm.hardening_modulus_1;
  double H_2 = d_cm.beta * d_cm.hardening_modulus_2;
  double stt = sqrt(3.0 / 2.0);
  double o_stt = 1.0 / stt;
  double denom = 1.0 / (1.0 + stt * H_2 * delLambda);

  // Compute updated backstress
  backStress_new =
    backStress_old + df_dsigma_normal_new * (delLambda * H_1 * o_stt);
  backStress_new = backStress_new * denom;

  return;
}

void
KinematicHardening_Armstrong::eval_h_beta(const Matrix3& df_dsigma,
                                          const ModelStateBase* state,
                                          Matrix3& h_beta)
{
  double H_1 = d_cm.beta * d_cm.hardening_modulus_1;
  double H_2 = d_cm.beta * d_cm.hardening_modulus_2;
  Matrix3 beta = state->backStress;
  double norm_r = df_dsigma.Norm();
  h_beta = df_dsigma * (2.0 / 3.0 * H_1) - beta * (H_2 * norm_r);
  return;
}
