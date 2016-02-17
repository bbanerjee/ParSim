/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

KinematicHardening_MasonSand::KinematicHardening_MasonSand(ProblemSpecP& ps,
                                                           InternalVariableModel* intvar)
{
  // Copy of internal variables object
  d_intvar = intvar;

  // Initial fluid pressure
  ps->require("fluid_pressure_initial",d_cm.fluid_pressure_initial);

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
  VarLabel::destroy(pPorePressureLabel);
  VarLabel::destroy(pPorePressureLabel_preReloc);
  VarLabel::destroy(pZetaLabel);
  VarLabel::destroy(pZetaLabel_preReloc);
}

void KinematicHardening_MasonSand::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP plastic_ps = ps->appendChild("kinematic_hardening_model");
  plastic_ps->setAttribute("type","mason_sand_pore_pressure");
  plastic_ps->appendElement("fluid_pressure_initial",d_cm.fluid_pressure_initial);
}

void 
KinematicHardening_MasonSand::computeBackStress(const ModelStateBase* ,
                                                const double& delT,
                                                const particleIndex idx,
                                                const double& delLambda,
                                                const Matrix3& df_dsigma_normal_new,
                                                const Matrix3& backStress_old,
                                                Matrix3& backStress_new)
{
  // **TODO** Compute updated backstress

  return;
}

void 
KinematicHardening_MasonSand::eval_h_beta(const Matrix3& df_dsigma,
                                          const ModelStateBase* state_input,
                                          Matrix3& h_beta)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // **TODO** Compute derivative of backstress
  return;
}

