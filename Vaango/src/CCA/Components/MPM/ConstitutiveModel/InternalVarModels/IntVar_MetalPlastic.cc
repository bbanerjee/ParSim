/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_MetalPlastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Labels/MPMLabel.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;

using MetalIntVar = Uintah::MetalIntVar;

IntVar_MetalPlastic::IntVar_MetalPlastic(ProblemSpecP& ps,
                                         ShearModulusModel* shear,
                                         MPMEquationOfState* eos)
{
  d_elastic = nullptr;
  d_shear   = shear;
  d_eos     = eos;
  initializeLocalMPMLabels();
}

IntVar_MetalPlastic::IntVar_MetalPlastic(const IntVar_MetalPlastic* cm)
{
  d_elastic = cm->d_elastic;
  d_shear   = cm->d_shear;
  initializeLocalMPMLabels();
}

void
IntVar_MetalPlastic::initializeLocalMPMLabels()
{
  auto type_d = Uintah::ParticleVariable<Uintah::MetalIntVar>::getTypeDescription();
  pIntVarLabel = Uintah::VarLabel::create("p.intVarMetalPlastic", type_d);
  pIntVarLabel_preReloc =
    Uintah::VarLabel::create("p.intVarMetalPlastic+", type_d);
}

IntVar_MetalPlastic::~IntVar_MetalPlastic()
{
  VarLabel::destroy(pIntVarLabel);
  VarLabel::destroy(pIntVarLabel_preReloc);
}

void
IntVar_MetalPlastic::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "elastic_plastic_hp");
}

/* get internal variable labels */
std::vector<const Uintah::VarLabel*>
IntVar_MetalPlastic::getLabels() const
{
  std::vector<const Uintah::VarLabel*> labels;
  labels.push_back(pIntVarLabel);
  labels.push_back(pIntVarLabel_preReloc);
  return labels;
}

void
IntVar_MetalPlastic::addParticleState(std::vector<const VarLabel*>& from,
                                      std::vector<const VarLabel*>& to)
{
  from.push_back(pIntVarLabel);
  to.push_back(pIntVarLabel_preReloc);
}

void
IntVar_MetalPlastic::addInitialComputesAndRequires(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pIntVarLabel, matlset);
}

void
IntVar_MetalPlastic::initializeInternalVariable(Uintah::ParticleSubset* pset,
                                                Uintah::DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<Uintah::MetalIntVar> pIntVar;
  new_dw->allocateAndPut(pIntVar, pIntVarLabel, pset);

  for (auto pidx : *pset) {
    pIntVar[pidx] = {0.0, 0.0};
  }
}

void
IntVar_MetalPlastic::initializeInternalVariable(const Patch* patch,
                                                const MPMMaterial* matl,
                                                ParticleSubset* pset,
                                                DataWarehouse* new_dw,
                                                MPMLabel* lb,
                                                ParameterDict& params)
{
  initializeInternalVariable(pset, new_dw);
}

void
IntVar_MetalPlastic::addComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pIntVarLabel, matlset, Ghost::None);
  task->computes(pIntVarLabel_preReloc, matlset);
}

std::vector<Uintah::constParticleVariable<double>>
IntVar_MetalPlastic::getInternalVariables(Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* old_dw,
                                          const double& dummy)
{
  Uintah::constParticleVariable<double> pIntVar;
  old_dw->get(pIntVar, pIntVarLabel, pset);

  std::vector<Uintah::constParticleVariable<double>> pIntVars;
  pIntVars.emplace_back(pIntVar);

  return pIntVars;
}

std::vector<Uintah::constParticleVariable<Uintah::Matrix3>>
IntVar_MetalPlastic::getInternalVariables(Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* old_dw,
                                          const Uintah::Matrix3& dummy)
{
  std::vector<Uintah::constParticleVariable<Uintah::Matrix3>> pIntVars;
  return pIntVars;
}

// Allocate and put the local particle internal variables
void
IntVar_MetalPlastic::allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                                    Uintah::DataWarehouse* new_dw,
                                                    Uintah::ParticleVariableBase& var)
{
  new_dw->allocateAndPut(var, pIntVarLabel, pset);
}

void
IntVar_MetalPlastic::allocateAndPutInternalVariable(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* new_dw,
  ParticleDoublePVec& pVars)
{
}

// Allocate and put the local <Matrix3> particle variables
void
IntVar_MetalPlastic::allocateAndPutInternalVariable(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* new_dw,
  ParticleMatrix3PVec& pVars)
{
}

template <>
void
IntVar_MetalPlastic::evolveInternalVariable(const Uintah::VarLabel* label,
                                            Uintah::particleIndex idx,
                                            const ModelStateBase* state,
                                            Uintah::constParticleVariable<MetalIntVar>& var_old,
                                            Uintah::ParticleVariable<MetalIntVar>& var_new)
{
  var_new[idx].eqPlasticStrain = 
    computeEqPlasticStrain(var_old[idx].eqPlasticStrain, state);
  var_new[idx].plasticPorosity = 
    computePlasticPorosity(var_old[idx].plasticPorosity, state);
}

double
IntVar_MetalPlastic::computeInternalVariable(const Uintah::VarLabel* label,
                                             const ModelStateBase* state_input) const 
{
  return 0.0;
}

double
IntVar_MetalPlastic::computeEqPlasticStrain(double eqPlasticStrain_old,
                                            const ModelStateBase* state) const
{
  return 0.0;
}

double
IntVar_MetalPlastic::computePlasticPorosity(double plasticPorosity_old,
                                            const ModelStateBase* state) const
{
  return 0.0;
}

double
IntVar_MetalPlastic::computeVolStrainDerivOfInternalVariable(
    const Uintah::VarLabel* label,
    const ModelStateBase* state) const 
{
  return 0.0;
}

/*!-----------------------------------------------------*/
void
IntVar_MetalPlastic::allocateCMDataAddRequires(Task* task,
                                               const MPMMaterial* matl,
                                               const PatchSet*,
                                               MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(
    Task::NewDW, pIntVarLabel_preReloc, matlset, Ghost::None);
}

void
IntVar_MetalPlastic::allocateCMDataAdd(DataWarehouse* old_dw,
                                       ParticleSubset* addset,
                                       ParticleLabelVariableMap* newState,
                                       ParticleSubset* delset,
                                       DataWarehouse* new_dw)
{
  ParticleVariable<MetalIntVar> pIntVar;
  constParticleVariable<MetalIntVar> o_pIntVar;

  new_dw->allocateTemporary(pIntVar, addset);

  new_dw->get(o_pIntVar, pIntVarLabel_preReloc, delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pIntVar[*n] = o_pIntVar[*o];
  }

  (*newState)[pIntVarLabel] = pIntVar.clone();
}

/*!-----------------------------------------------------*/
void
IntVar_MetalPlastic::allocateAndPutRigid(ParticleSubset* pset,
                                         DataWarehouse* new_dw,
                                         constParticleLabelVariableMap& var)
{
  ParticleVariable<MetalIntVar> pIntVar_new;
  new_dw->allocateAndPut(
    pIntVar_new, pIntVarLabel_preReloc, pset);
  for (auto pidx : *pset) {
    pIntVar_new[pidx] = dynamic_cast<constParticleVariable<MetalIntVar>&>(
      *var[pIntVarLabel])[pidx];
  }
}
