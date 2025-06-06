/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_BorjaPressure.h>

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModel.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DOUT.hpp>

#include <cmath>
#include <iostream>

namespace {

Uintah::Dout g_dbg("IVBorjaP",
                   "IntVar_BorjaPressure",
                   "general info on borja pressure internal var",
                   false);
}

using namespace Vaango;
using namespace Uintah;

IntVar_BorjaPressure::IntVar_BorjaPressure(ProblemSpecP& ps,
                                           ShearModulusModel* shear)
{
  d_elastic = nullptr;
  d_shear   = shear;

  ParameterDict eosParams = (d_shear->getPressureModel())->getParameters();
  d_kappatilde            = eosParams["kappatilde"];
  d_kappahat              = eosParams["kappahat"];

  ps->require("pc0", d_pc0);
  ps->require("lambdatilde", d_lambdatilde);

  // Compute lambda_hat for large deformation
  d_lambdahat = d_lambdatilde / (1.0 - d_lambdatilde);

  // Initialize internal variable labels for evolution
  pPcLabel =
    VarLabel::create("p.p_c", ParticleVariable<double>::getTypeDescription());
  pPcLabel_preReloc =
    VarLabel::create("p.p_c+", ParticleVariable<double>::getTypeDescription());
}

IntVar_BorjaPressure::IntVar_BorjaPressure(const IntVar_BorjaPressure* cm)
{
  d_elastic = cm->d_elastic;
  d_shear   = cm->d_shear;

  d_pc0         = cm->d_pc0;
  d_lambdatilde = cm->d_lambdatilde;
  d_lambdahat   = cm->d_lambdahat;
  d_kappatilde  = cm->d_kappatilde;
  d_kappahat    = cm->d_kappahat;

  // Initialize internal variable labels for evolution
  pPcLabel =
    VarLabel::create("p.p_c", ParticleVariable<double>::getTypeDescription());
  pPcLabel_preReloc =
    VarLabel::create("p.p_c+", ParticleVariable<double>::getTypeDescription());
}

IntVar_BorjaPressure::~IntVar_BorjaPressure()
{
  VarLabel::destroy(pPcLabel);
  VarLabel::destroy(pPcLabel_preReloc);
}

void
IntVar_BorjaPressure::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "borja_consolidation_pressure");

  int_var_ps->appendElement("pc0", d_pc0);
  int_var_ps->appendElement("lambdatilde", d_lambdatilde);
}

void
IntVar_BorjaPressure::addParticleState(std::vector<const VarLabel*>& from,
                                       std::vector<const VarLabel*>& to)
{
  from.push_back(pPcLabel);
  to.push_back(pPcLabel_preReloc);
}

void
IntVar_BorjaPressure::addInitialComputesAndRequires(Task* task,
                                                    const MPMMaterial* matl,
                                                    const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pPcLabel, matlset);
}

void
IntVar_BorjaPressure::initializeInternalVariable(ParticleSubset* pset,
                                                 DataWarehouse* new_dw)
{
  ParticleVariable<double> pPc;
  new_dw->allocateAndPut(pPc, pPcLabel, pset);
  for (auto idx : *pset) {
    pPc[idx] = d_pc0;
  }
}

void
IntVar_BorjaPressure::addComputesAndRequires(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->needs(Task::OldDW, pPcLabel, matlset, Ghost::None);
  task->computes(pPcLabel_preReloc, matlset);
}

/* Get one (possibly composite) internal variable */
template<>
void
IntVar_BorjaPressure::getInternalVariable(ParticleSubset* pset,
                                          DataWarehouse* old_dw,
                                          constParticleVariable<double>& pPc)
{
  old_dw->get(pPc, pPcLabel, pset);
}

/* Get multiple local <int/double/Vector/Matrix3> internal variables */
template<>
std::vector<constParticleVariable<double>>
IntVar_BorjaPressure::getInternalVariables(ParticleSubset* pset,
                                           DataWarehouse* old_dw)
{
  constParticleDouble pPc;
  old_dw->get(pPc, pPcLabel, pset);

  constParticleDoubleVec intVarVec;
  intVarVec.push_back(pPc);

  return intVarVec;
}

template<>
void
IntVar_BorjaPressure::allocateAndPutInternalVariable(
  ParticleSubset* pset,
  DataWarehouse* new_dw,
  ParticleVariable<double>& pPc_new)
{
  new_dw->allocateAndPut(pPc_new, pPcLabel_preReloc, pset);
}

template<>
void
IntVar_BorjaPressure::evolveInternalVariable(
  Uintah::particleIndex pidx,
  const ModelStateBase* state,
  Uintah::constParticleVariable<BorjaIntVar>& var_old,
  Uintah::ParticleVariable<BorjaIntVar>& var)
{
}

////////////////////////////////////////////////////////////////////////////////////////
//  Compute the internal variable
double
IntVar_BorjaPressure::computeInternalVariable(
  const std::string& label,
  const ModelStateBase* state_old,
  const ModelStateBase* state_cur) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_cur);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get old p_c
  double pc_n = state->p_c0; // Old Pc

  // Get the trial elastic strain and the updated elastic strain
  // (volumetric part)
  double strain_elast_v_tr = state->epse_v_tr;
  double strain_elast_v    = state->epse_v;

  // Calculate new p_c
  double pc = pc_n * std::exp((strain_elast_v_tr - strain_elast_v) /
                              (d_kappahat - d_lambdahat));

#ifdef CATCH_NOT_FINITE
  if (!std::isfinite(pc)) {
    std::ostringstream desc;
    desc << "pc = " << pc << " pc_n = " << pc_n
         << " strain_elast_tr = " << strain_elast_v_tr
         << " strain_elast_v = " << strain_elast_v
         << " kappahat = " << d_kappahat << " lambdahat = " << d_lambdahat
         << "\n";
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }
#endif

  return pc;
}

////////////////////////////////////////////////////////////////////////////////////////
// Compute derivative of internal variable with respect to volumetric
// elastic strain
double
IntVar_BorjaPressure::computeVolStrainDerivOfInternalVariable(
  const std::string& label,
  const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  // if (!state) {
  //   std::ostringstream out;
  //   out << "**ERROR** The correct ModelState object has not been passed."
  //       << " Need ModelState_CamClay.";
  //   throw InternalError(out.str(), __FILE__, __LINE__);
  // }

  // Get old p_c
  double pc_n = state->p_c0;

  // Get the trial elastic strain and the updated elastic strain
  // (volumetric part)
  double strain_elast_v_tr = state->epse_v_tr;
  double strain_elast_v    = state->epse_v;

  // Calculate  dp_c/depse_v
  double pc = pc_n * std::exp((strain_elast_v_tr - strain_elast_v) /
                              (d_kappahat - d_lambdahat));

  return pc / (d_lambdahat - d_kappahat);
}

void
IntVar_BorjaPressure::allocateCMDataAddRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet*,
                                                MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->needs(Task::NewDW, pPcLabel_preReloc, matlset, Ghost::None);
}

void
IntVar_BorjaPressure::allocateCMDataAdd(DataWarehouse* old_dw,
                                        ParticleSubset* addset,
                                        ParticleLabelVariableMap* newState,
                                        ParticleSubset* delset,
                                        DataWarehouse* new_dw)
{
  ParticleVariable<double> pPc;
  constParticleVariable<double> o_Pc;

  new_dw->allocateTemporary(pPc, addset);

  new_dw->get(o_Pc, pPcLabel_preReloc, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pPc[*n] = o_Pc[*o];
  }

  (*newState)[pPcLabel] = pPc.clone();
}

void
IntVar_BorjaPressure::allocateAndPutRigid(ParticleSubset* pset,
                                          DataWarehouse* new_dw,
                                          constParticleVariableBase& pPc)
{
  ParticleVariable<double> pPc_new;
  new_dw->allocateAndPut(pPc_new, pPcLabel_preReloc, pset);
  for (auto idx : *pset) {
    pPc_new[idx] = dynamic_cast<constParticleVariable<double>&>(pPc)[idx];
  }
}
