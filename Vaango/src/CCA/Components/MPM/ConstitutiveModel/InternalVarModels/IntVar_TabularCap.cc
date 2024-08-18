/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <errno.h>
#include <fenv.h>

//#define USE_TOTAL_STRAIN

using namespace Vaango;

using Uintah::VarLabel;
using Uintah::ProblemSpecP;
using Uintah::MPMMaterial;
using Uintah::PatchSet;
using Uintah::MaterialSubset;
using Uintah::Patch;
using Uintah::ParticleSubset;
using Uintah::DataWarehouse;
using Uintah::MPMLabel;
using Uintah::Task;
using Uintah::Ghost;
using Uintah::ParticleLabelVariableMap;
using Uintah::constParticleLabelVariableMap;
using Uintah::TabularCapIntVar;
using Uintah::constParticleVariable;
using Uintah::particleIndex;
using Uintah::ParticleVariable;

/*!-----------------------------------------------------*/
IntVar_TabularCap::IntVar_TabularCap(Uintah::ProblemSpecP& ps)
  : d_capX_fn(ps)
{
  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}

/*!-----------------------------------------------------*/
IntVar_TabularCap::IntVar_TabularCap(const IntVar_TabularCap* cm)
{
  d_capX_fn = cm->d_capX_fn;

  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}

/*!-----------------------------------------------------*/
IntVar_TabularCap::~IntVar_TabularCap()
{
  VarLabel::destroy(pCapXLabel);
  VarLabel::destroy(pCapXLabel_preReloc);
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "tabular_cap");
  d_capX_fn.table.outputProblemSpec(int_var_ps);
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::addInitialComputesAndRequires(Task* task,
                                                 const MPMMaterial* matl,
                                                 const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pCapXLabel, matlset);
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::initializeInternalVariable(Uintah::ParticleSubset* pset,
                                              Uintah::DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<double> pCapX;
  new_dw->allocateAndPut(pCapX, pCapXLabel, pset);

#ifdef USE_TOTAL_STRAIN
  // The table is of the form x -> ev_bar, y -> X_p_bar = X_bar/3
  double capX = -1000;
#else
  // The table is of the form x -> ev_p_bar, y -> X_p_bar = X_bar/3
  DoubleVec1D gg = d_capX_fn.table.interpolate<1>({ { 0.0 } });
  double capX    = -gg[0] * 3.0;
#endif

  for (auto particle : *pset) {
    pCapX[particle] = capX;
  }
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::addComputesAndRequires(Task* task,
                                          const MPMMaterial* matl,
                                          const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pCapXLabel, matlset, Ghost::None);
  task->computes(pCapXLabel_preReloc, matlset);
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::addParticleState(std::vector<const VarLabel*>& from,
                                    std::vector<const VarLabel*>& to)
{
  from.push_back(pCapXLabel);
  to.push_back(pCapXLabel_preReloc);
}

/* Get one (possibly composite) internal variable */
template <>
void
IntVar_TabularCap::getInternalVariable<double>(
  ParticleSubset* pset,
  DataWarehouse* old_dw,
  constParticleVariable<double>& pCapX)
{
  old_dw->get(pCapX, pCapXLabel, pset);
}

/* Allocate one (possibly composite) internal variable */
template <>
void
IntVar_TabularCap::allocateAndPutInternalVariable(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* new_dw,
  Uintah::ParticleVariable<double>& pCapX_new)
{
  new_dw->allocateAndPut(pCapX_new, pCapXLabel_preReloc, pset);
}

template <>
void
IntVar_TabularCap::evolveInternalVariable<TabularCapIntVar>(
  particleIndex pidx,
  const ModelStateBase* state,
  constParticleVariable<TabularCapIntVar>& var_old,
  ParticleVariable<TabularCapIntVar>& var)
{
}

//--------------------------------------------------------------------------------------
// Compute hydrostatic strength
//--------------------------------------------------------------------------------------
double
IntVar_TabularCap::computeInternalVariable(
  const std::string& label,
  const ModelStateBase* /*state_old*/,
  const ModelStateBase* state_cur) const
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_cur);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double ep_v_bar  = -state->ep_v;
  double X_bar_new = -state->capX;
  if (ep_v_bar < 0.0) { // tension
    X_bar_new = computeDrainedHydrostaticStrength(0.0);
  } else {
#ifdef USE_TOTAL_STRAIN
    // The table is of the form x -> ev_bar, y -> X_p_bar = X_bar/3
    double ev_bar = ep_v_bar - state->elasticStrainTensor.Trace();
    X_bar_new     = computeDrainedHydrostaticStrength(ev_bar);
#else
    // The table is of the form x -> ep_v_bar, y -> X_p_bar = X_bar/3
    X_bar_new       = computeDrainedHydrostaticStrength(ep_v_bar);
#endif
  }

  double X_new = -X_bar_new;

  return X_new;
}

/**
 *  Compute drained hydrostatic strength
 */
double
IntVar_TabularCap::computeDrainedHydrostaticStrength(
  const double& ep_v_bar) const
{
  DoubleVec1D gg       = d_capX_fn.table.interpolate<1>({ { ep_v_bar } });
  double X_bar_drained = gg[0] * 3.0;
  //std::cout << "CapX : ep_v = " << ep_v_bar << " p = " << gg[0] << "\n";

  return X_bar_drained;
}

/**
 *  Compute derivative of internal variable with respect to volumetric
 *  plastic strain
 */
double
IntVar_TabularCap::computeVolStrainDerivOfInternalVariable(
  const std::string& label,
  const ModelStateBase* state_input) const
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double epsilon = 1.0e-6;
#ifdef USE_TOTAL_STRAIN
  // The table is of the form x -> ev_bar, y -> X_p_bar = X_bar/3
  double ev_bar     = -state->elasticStrainTensor.Trace() - state->ep_v;
  double ev_bar_lo  = ev_bar - epsilon;
  double ev_bar_hi  = ev_bar + epsilon;
  double X_p_bar_lo = computeDrainedHydrostaticStrength(ev_bar_lo);
  double X_p_bar_hi = computeDrainedHydrostaticStrength(ev_bar_hi);
#else
  // The table is of the form x -> ev_p_bar, y -> X_p_bar = X_bar/3
  double ep_v_bar    = -state->ep_v;
  double ep_v_bar_lo = ep_v_bar - epsilon;
  double ep_v_bar_hi = ep_v_bar + epsilon;
  double X_p_bar_lo  = computeDrainedHydrostaticStrength(ep_v_bar_lo);
  double X_p_bar_hi  = computeDrainedHydrostaticStrength(ep_v_bar_hi);
#endif

  double dX_p_bar_dep_v_bar = (X_p_bar_hi - X_p_bar_lo) / (2 * epsilon);
  return dX_p_bar_dep_v_bar;
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::allocateCMDataAddRequires(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet*,
                                             MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pCapXLabel_preReloc, matlset, Ghost::None);
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::allocateCMDataAdd(DataWarehouse* old_dw,
                                     ParticleSubset* addset,
                                     ParticleLabelVariableMap* newState,
                                     ParticleSubset* delset,
                                     DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<double> pCapX;
  Uintah::constParticleVariable<double> o_capX;

  new_dw->allocateTemporary(pCapX, addset);

  new_dw->get(o_capX, pCapXLabel_preReloc, delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pCapX[*n] = o_capX[*o];
  }

  (*newState)[pCapXLabel] = pCapX.clone();
}

/*!-----------------------------------------------------*/
void
IntVar_TabularCap::allocateAndPutRigid(ParticleSubset* pset,
                                       DataWarehouse* new_dw,
                                       constParticleLabelVariableMap& var)
{
  Uintah::ParticleVariable<double> pCapX_new;
  new_dw->allocateAndPut(pCapX_new, pCapXLabel_preReloc, pset);
  for (int& iter : *pset) {
    pCapX_new[iter] = dynamic_cast<Uintah::constParticleVariable<double>&>(
      *var[pCapXLabel])[iter];
  }
}

namespace Vaango {

template void
IntVar_TabularCap::getInternalVariable<double>(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* old_dw,
  Uintah::constParticleVariable<double>& pCapX);
template void
IntVar_TabularCap::evolveInternalVariable<TabularCapIntVar>(
  Uintah::particleIndex pidx,
  const ModelStateBase* state,
  Uintah::constParticleVariable<TabularCapIntVar>& var_old,
  Uintah::ParticleVariable<TabularCapIntVar>& var);
}
