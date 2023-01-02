/*
 * The MIT License
 *
 * Copyright (c) 1997-2014 The University of Utah
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/MPMBoundCond.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Components/MPM/ReactionDiffusion/JGConcentrationDiffusion.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>

#include <iostream>

namespace Uintah {

JGConcentrationDiffusion::JGConcentrationDiffusion(ProblemSpecP& ps,
                                                   MaterialManagerP& matManager,
                                                   MPMFlags* Mflag,
                                                   string diff_type)
  : ScalarDiffusionModel(ps, matManager, Mflag, diff_type)
{
}

void
JGConcentrationDiffusion::scheduleComputeFlux(Task* task,
                                              const MPMMaterial* matl,
                                              const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  Ghost::GhostType gnone        = Ghost::None;

  task->requires(Task::OldDW, d_lb->pConcGradientLabel, matlset, gnone);
  task->computes(d_lb->pFluxLabel, matlset);
}

void
JGConcentrationDiffusion::computeFlux(const Patch* patch,
                                      const MPMMaterial* matl,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  auto interpolator = d_Mflag->d_interpolator->clone(patch);
  vector<IntVector> ni(interpolator->size());
  vector<Vector> d_S(interpolator->size());

  Vector dx = patch->dCell();
  int dwi   = matl->getDWIndex();
  constParticleVariable<Vector> pConcGradient;
  ParticleVariable<Vector> pFlux;

  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

  old_dw->get(pConcGradient, d_lb->pConcGradientLabel, pset);
  new_dw->allocateAndPut(pFlux, d_lb->pFluxLabel, pset);

  double timestep = 10000.0;
  for (const auto& idx : *pset) {
    pFlux[idx] = diffusivity * pConcGradient[idx];

    timestep = std::min(timestep, computeStableTimeStep(diffusivity, dx));
  } // End of Particle Loop

  // cout << "Time Step: " << timestep << endl;

  // new_dw->put(delt_vartype(timestep), d_lb->delTLabel, patch->getLevel());
  // delete interpolator;
}

} // end namespace Uintah