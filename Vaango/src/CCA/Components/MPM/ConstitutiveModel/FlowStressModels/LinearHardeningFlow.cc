/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include "LinearHardeningFlow.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Math/FastMatrix.h>
#include <cmath>

using namespace Uintah;
using Vaango::ModelStateBase;

LinearHardeningFlow::LinearHardeningFlow(ProblemSpecP& ps)
{
  ps->require("K", d_CM.K);
  ps->require("sigma_0", d_CM.sigma_0);

  // Initialize internal variable labels for evolution
  pAlphaLabel =
    VarLabel::create("p.alpha", ParticleVariable<double>::getTypeDescription());
  pAlphaLabel_preReloc = VarLabel::create(
    "p.alpha+", ParticleVariable<double>::getTypeDescription());
}

LinearHardeningFlow::LinearHardeningFlow(const LinearHardeningFlow* cm)
{
  d_CM.K = cm->d_CM.K;
  d_CM.sigma_0 = cm->d_CM.sigma_0;

  // Initialize internal variable labels for evolution
  pAlphaLabel =
    VarLabel::create("p.alpha", ParticleVariable<double>::getTypeDescription());
  pAlphaLabel_preReloc = VarLabel::create(
    "p.alpha+", ParticleVariable<double>::getTypeDescription());
}

LinearHardeningFlow::~LinearHardeningFlow()
{
  VarLabel::destroy(pAlphaLabel);
  VarLabel::destroy(pAlphaLabel_preReloc);
}

void
LinearHardeningFlow::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP flow_ps = ps->appendChild("flow_model");
  flow_ps->setAttribute("type", "linear");
  flow_ps->appendElement("K", d_CM.K);
  flow_ps->appendElement("sigma_0", d_CM.sigma_0);
}

void
LinearHardeningFlow::addInitialComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pAlphaLabel, matlset);
}

void
LinearHardeningFlow::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                         const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pAlphaLabel, matlset, Ghost::None);
  task->computes(pAlphaLabel_preReloc, matlset);
}

void
LinearHardeningFlow::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                         const PatchSet*, bool /*recurse*/,
                                         bool SchedParent)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  if (SchedParent) {
    task->requires(Task::ParentOldDW, pAlphaLabel, matlset, Ghost::None);
  } else {
    task->requires(Task::OldDW, pAlphaLabel, matlset, Ghost::None);
  }
}

void
LinearHardeningFlow::addParticleState(std::vector<const VarLabel*>& from,
                                   std::vector<const VarLabel*>& to)
{
  from.push_back(pAlphaLabel);
  to.push_back(pAlphaLabel_preReloc);
}

void
LinearHardeningFlow::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                            const PatchSet*, MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pAlphaLabel_preReloc, matlset, Ghost::None);
}

void
LinearHardeningFlow::allocateCMDataAdd(DataWarehouse* new_dw,
                                    ParticleSubset* addset,
                                    ParticleLabelVariableMap* newState,
                                    ParticleSubset* delset, DataWarehouse*)
{
  constParticleVariable<double> o_Alpha;
  new_dw->get(o_Alpha, pAlphaLabel_preReloc, delset);

  ParticleVariable<double> pAlpha;
  new_dw->allocateTemporary(pAlpha, addset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pAlpha[*n] = o_Alpha[*o];
  }

  (*newState)[pAlphaLabel] = pAlpha.clone();
}

void
LinearHardeningFlow::initializeInternalVars(ParticleSubset* pset,
                                         DataWarehouse* new_dw)
{
  new_dw->allocateAndPut(pAlpha_new, pAlphaLabel, pset);
  for (auto idx : *pset) {
    pAlpha_new[idx] = 0.0;
  }
}

void
LinearHardeningFlow::getInternalVars(ParticleSubset* pset, DataWarehouse* old_dw)
{
  old_dw->get(pAlpha, pAlphaLabel, pset);
}

void
LinearHardeningFlow::allocateAndPutInternalVars(ParticleSubset* pset,
                                             DataWarehouse* new_dw)
{
  new_dw->allocateAndPut(pAlpha_new, pAlphaLabel_preReloc, pset);
}

void
LinearHardeningFlow::allocateAndPutRigid(ParticleSubset* pset,
                                      DataWarehouse* new_dw)
{
  new_dw->allocateAndPut(pAlpha_new, pAlphaLabel_preReloc, pset);
  for (auto idx : *pset) {
    pAlpha_new[idx] = 0.0;
  }
}

void
LinearHardeningFlow::updateElastic(const particleIndex idx)
{
  pAlpha_new[idx] = pAlpha[idx];
}

void
LinearHardeningFlow::updatePlastic(const particleIndex idx, const double& delGamma)
{
  pAlpha_new[idx] = pAlpha[idx] + sqrt(2.0 / 3.0) * delGamma;
}

double
LinearHardeningFlow::computeFlowStress(const ModelStateBase* state, const double&,
                                    const double&, const MPMMaterial*,
                                    const particleIndex idx)
{
  //  double flowStress = d_CM.sigma_0 + d_CM.K*pAlpha[idx];
  double flowStress = d_CM.sigma_0 + d_CM.K * state->eqPlasticStrain;
  return flowStress;
}

double
LinearHardeningFlow::computeEpdot(const ModelStateBase* state, const double&,
                               const double&, const MPMMaterial*,
                               const particleIndex)
{
  return state->eqPlasticStrainRate;
}

void
LinearHardeningFlow::computeTangentModulus(const Matrix3& stress,
                                        const ModelStateBase*, const double&,
                                        const MPMMaterial*, const particleIndex,
                                        TangentModulusTensor&,
                                        TangentModulusTensor&)
{
  throw InternalError("Empty Function: LinearHardeningFlow::computeTangentModulus",
                      __FILE__, __LINE__);
}

void
LinearHardeningFlow::evalDerivativeWRTScalarVars(const ModelStateBase* state,
                                              const particleIndex idx,
                                              Vector& derivs) const
{
  derivs[0] = evalDerivativeWRTStrainRate(state, idx);
  derivs[1] = evalDerivativeWRTTemperature(state, idx);
  derivs[2] = evalDerivativeWRTPlasticStrain(state, idx);
}

double
LinearHardeningFlow::evalDerivativeWRTPlasticStrain(const ModelStateBase*,
                                                 const particleIndex) const
{
  return d_CM.K;
}

///////////////////////////////////////////////////////////////////////////
/*  Compute the shear modulus. */
///////////////////////////////////////////////////////////////////////////
double
LinearHardeningFlow::computeShearModulus(const ModelStateBase* state)
{
  return state->shearModulus;
}

///////////////////////////////////////////////////////////////////////////
/* Compute the melting temperature */
///////////////////////////////////////////////////////////////////////////
double
LinearHardeningFlow::computeMeltingTemp(const ModelStateBase* state)
{
  return state->meltingTemp;
}

double
LinearHardeningFlow::evalDerivativeWRTTemperature(const ModelStateBase*,
                                               const particleIndex) const
{
  return 0.0;
}

double
LinearHardeningFlow::evalDerivativeWRTStrainRate(const ModelStateBase*,
                                              const particleIndex) const
{
  return 0.0;
}

double
LinearHardeningFlow::evalDerivativeWRTAlpha(const ModelStateBase*,
                                         const particleIndex) const
{
  return d_CM.K;
}
