/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ConstantRate.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Util/DebugStream.h>

#include <iostream>
using namespace Uintah;

static DebugStream cout_doing("AMRMPM", false);

ConstantRate::ConstantRate(ProblemSpecP& ps,
                           MaterialManagerP& sS,
                           MPMFlags* Mflag,
                           std::string diff_type)
  : ScalarDiffusionModel(ps, sS, Mflag, diff_type)
{
  ps->require("constant_rate", d_constant_rate);
  std::cout << "rate: " << d_constant_rate << std::endl;
}

ConstantRate::~ConstantRate() {}

void
ConstantRate::addInitialComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            [[maybe_unused]] const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(d_lb->diffusion->pFlux, matlset);
}

void
ConstantRate::addParticleState(std::vector<const VarLabel*>& from,
                               std::vector<const VarLabel*>& to) const
{
  from.push_back(d_lb->diffusion->pFlux);
  to.push_back(d_lb->diffusion->pFlux_preReloc);
}

void
ConstantRate::computeFlux(const Patch* patch,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  // Get the current simulation time
  // double simTime = d_materialManager->getElapsedSimTime();

  // simTime_vartype simTime;
  // old_dw->get(simTime, d_lb->simulationTimeLabel);

  auto interpolator = d_Mflag->d_interpolator->clone(patch);
  std::vector<IntVector> ni(interpolator->size());
  std::vector<double> S(interpolator->size());

  int dwi = matl->getDWIndex();
  //  Vector dx = patch->dCell();
  //  double comp_diffusivity;

  ParticleVariable<Vector> pFlux;

  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

  new_dw->allocateAndPut(pFlux, d_lb->diffusion->pFlux_preReloc, pset);

  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
       iter++) {
    particleIndex idx = *iter;

    pFlux[idx] = Vector(0.0, 0.0, 0.0);
  } // End of Particle Loop
}

void
ConstantRate::initializeSDMData(const Patch* patch,
                                const MPMMaterial* matl,
                                DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<Vector> pFlux;

  new_dw->allocateAndPut(pFlux, d_lb->diffusion->pFlux, pset);

  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
       iter++) {
    pFlux[*iter] = Vector(0, 0, 0);
  }
}

void
ConstantRate::scheduleComputeFlux(Task* task,
                                  const MPMMaterial* matl,
                                  [[maybe_unused]] const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // task->needs(Task::OldDW, d_lb->simulationTimeLabel,);

  task->computes(d_lb->diffusion->pFlux_preReloc, matlset);
}

void
ConstantRate::addSplitParticlesComputesAndRequires(
  [[maybe_unused]] Task* task,
  [[maybe_unused]] const MPMMaterial* matl,
  [[maybe_unused]] const PatchSet* patches) const
{
  // Do nothing for now
}

void
ConstantRate::splitSDMSpecificParticleData([[maybe_unused]] const Patch* patch,
                                           [[maybe_unused]] const int dwi,
                                           [[maybe_unused]] const int nDims,
                                           [[maybe_unused]] ParticleVariable<int>& prefOld,
                                           [[maybe_unused]] ParticleVariable<int>& pref,
                                           [[maybe_unused]] const unsigned int oldNumParts,
                                           [[maybe_unused]] const int numNewPartNeeded,
                                           [[maybe_unused]] DataWarehouse* old_dw,
                                           [[maybe_unused]] DataWarehouse* new_dw)
{
  // Do nothing for now
}

void
ConstantRate::scheduleComputeDivergence(Task* task,
                                        const MPMMaterial* matl,
                                        [[maybe_unused]] const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->needs(Task::NewDW, d_lb->gMassLabel, Ghost::None);

  task->computes(d_lb->diffusion->gConcentrationRate, matlset);
}

void
ConstantRate::computeDivergence(const Patch* patch,
                                const MPMMaterial* matl,
                                [[maybe_unused]] DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  int dwi             = matl->getDWIndex();
  Ghost::GhostType gn = Ghost::None;

  //*********Start - Used for testing purposes - CG *******
  // int timestep = d_materialManager->getCurrentTopLevelTimestep();
  //*********End   - Used for testing purposes - CG *******

  constNCVariable<double> gMass;
  NCVariable<double> gConcRate;

  new_dw->get(gMass, d_lb->gMassLabel, dwi, patch, gn, 0);
  new_dw->allocateAndPut(gConcRate,
                         d_lb->diffusion->gConcentrationRate,
                         dwi,
                         patch);

  for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
       iter++) {
    IntVector n  = *iter;
    gConcRate[n] = d_constant_rate * gMass[n];
  }
}

void
ConstantRate::scheduleComputeDivergence_CFI([[maybe_unused]] Task* t,
                                            [[maybe_unused]] const MPMMaterial* matl,
                                            [[maybe_unused]] const PatchSet* patch) const
{
}

void
ConstantRate::computeDivergence_CFI([[maybe_unused]] const PatchSubset* finePatches,
                                    [[maybe_unused]] const MPMMaterial* matl,
                                    [[maybe_unused]] DataWarehouse* old_dw,
                                    [[maybe_unused]] DataWarehouse* new_dw)
{
}

void
ConstantRate::outputProblemSpec(ProblemSpecP& ps, bool output_rdm_tag) const
{
  ProblemSpecP rdm_ps = ps;

  if (output_rdm_tag) {
    rdm_ps = ps->appendChild("diffusion_model");
    rdm_ps->setAttribute("type", "constant_rate");
  }
  ScalarDiffusionModel::baseOutputSDMProbSpec(rdm_ps);
  rdm_ps->appendElement("constant_rate", d_constant_rate);
}
