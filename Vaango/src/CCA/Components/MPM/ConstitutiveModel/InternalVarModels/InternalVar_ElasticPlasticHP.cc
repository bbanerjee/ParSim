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

IntVar_MetalPlastic::IntVar_MetalPlastic(ProblemSpecP& ps,
                                         ElasticModuliModel* elastic)
{
  d_elastic = elastic;
  d_shear   = nullptr;
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
  auto type_d = Uintah::ParticleVariable<double>::getTypeDescription();
  pEqPlasticStrainLabel = Uintah::VarLabel::create("p.eqPlasticStrain", type_d);
  pEqPlasticStrainLabel_preReloc =
    Uintah::VarLabel::create("p.eqPlasticStrain+", type_d);

  pPlasticPorosityLabel = Uintah::VarLabel::create("p.plasticPorosity", type_d);
  pPlasticPorosityLabel_preReloc =
    Uintah::VarLabel::create("p.plasticPorosity+", type_d);
}

IntVar_MetalPlastic::~IntVar_MetalPlastic()
{
  VarLabel::destroy(pEqPlasticStrainLabel);
  VarLabel::destroy(pEqPlasticStrainLabel_preReloc);

  VarLabel::destroy(pPlasticPorosityLabel);
  VarLabel::destroy(pPlasticPorosityLabel_preReloc);
}

void
IntVar_MetalPlastic::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "elastic_plastic_hp");
}

/* get internal variable labels */
std::vector<const Uintah::VarLabel*>
IntVar_MetalPlastic::getLabels() const override
{
  std::vector<const Uintah::VarLabel*> labels;
  labels.push_back(pEqPlasticStrainLabel);
  labels.push_back(pEqPlasticStrainLabel_preReloc);

  labels.push_back(pPlasticPorosityLabel);
  labels.push_back(pPlasticPorosityLabel_preReloc);

  return labels;
}

void
IntVar_MetalPlastic::addParticleState(std::vector<const VarLabel*>& from,
                                      std::vector<const VarLabel*>& to)
{
  from.push_back(pEqPlasticStrainLabel);
  to.push_back(pEqPlasticStrainLabel_preReloc);

  from.push_back(pPlasticPorosityLabel);
  to.push_back(pPlasticPorosityLabel_preReloc);
}

void
IntVar_MetalPlastic::addInitialComputesAndRequires(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pEqPlasticStrainLabel, matlset);
  task->computes(pPlasticPorosityLabel, matlset);
}

void
IntVar_MetalPlastic::initializeInternalVariable(const Patch* patch,
                                                const MPMMaterial* matl,
                                                ParticleSubset* pset,
                                                DataWarehouse* new_dw,
                                                MPMLabel* lb,
                                                ParameterDict& params)
{
  Uintah::ParticleVariable<double> pEqPlasticStrain, pPlasticPorosity;
  new_dw->allocateAndPut(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticPorosity, pPlasticPorosityLabel, pset);

  for (auto pidx : *pset) {
    pEqPlasticStrain[pidx] = 0.0;
    pPlasticPorosity[pidx] = 0.0;
  }
}

void
IntVar_MetalPlastic::addComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pEqPlasticStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticPorosityLabel, matlset, Ghost::None);
  task->computes(pEqPlasticStrainLabel_preReloc, matlset);
  task->computes(pPlasticPorosityLabel_preReloc, matlset);
}

std::vector<Uintah::constParticleVariable<double>>
IntVar_MetalPlastic::getInternalVariables(Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* old_dw,
                                          const double& dummy) override
{
  Uintah::constParticleVariable<double> pEqPlasticStrain, pPlasticPorosity;
  old_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
  old_dw->get(pPlasticPorosity, pPlasticPorosityLabel, pset);

  std::vector<Uintah::constParticleVariable<double>> pIntVars;
  pIntVars.emplace_back(pEqPlasticStrain);
  pIntVars.emplace_back(pPlasticPorosity);

  return pIntVars;
}

std::vector<Uintah::constParticleVariable<Uintah::Matrix3>>
IntVar_MetalPlastic::getInternalVariables(Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* old_dw,
                                          const Uintah::Matrix3& dummy) override
{
  std::vector<Uintah::constParticleVariable<Uintah::Matrix3>> pIntVars;
  return pIntVars;
}

// Allocate and put the local particle internal variables
void
IntVar_MetalPlastic::allocateAndPutInternalVariable(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* new_dw,
  vectorParticleDoubleP& pVars) override
{
  new_dw->allocateAndPut(*pVars[0], pEqPlasticStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(*pVars[1], pPlasticPorosityLabel_preReloc, pset);
}

// Allocate and put the local <Matrix3> particle variables
void
IntVar_MetalPlastic::allocateAndPutInternalVariable(
  Uintah::ParticleSubset* pset,
  Uintah::DataWarehouse* new_dw,
  vectorParticleMatrix3P& pVars) override
{
}

double
IntVar_MetalPlastic::computeInternalVariable(const Uintah::VarLabel* label,
                                             const ModelStateBase* state) const
{
}

/**
 * Function: computeEqPlasticStrain
 */
double
IntVar_MetalPlastic::computeEqPlasticStrain(const ModelStateBase* state) const
{
  return 0.0;
}

/**
 * Function: computePorosity
 */
double
IntVar_MetalPlastic::computePorosity(const ModelStateBase* state) const
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
    Task::NewDW, pEqPlasticStrainLabel_preReloc, matlset, Ghost::None);
  task->requires(
    Task::NewDW, pPlasticPorosityLabel_preReloc, matlset, Ghost::None);
}

void
IntVar_MetalPlastic::allocateCMDataAdd(DataWarehouse* old_dw,
                                       ParticleSubset* addset,
                                       ParticleLabelVariableMap* newState,
                                       ParticleSubset* delset,
                                       DataWarehouse* new_dw)
{
  ParticleVariable<double> pEqPlasticStrain, pPlasticPorosity;
  constParticleVariable<double> o_pEqPlasticStrain, o_pPlasticPorosity;

  new_dw->allocateTemporary(pEqPlasticStrain, addset);
  new_dw->allocateTemporary(pPlasticPorosity, addset);

  new_dw->get(o_pEqPlasticStrain, pEqPlasticStrainLabel_preReloc, delset);
  new_dw->get(o_pPlasticPorosity, pPlasticPorosityLabel_preReloc, delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pEqPlasticStrain[*n] = o_pEqPlasticStrain[*o];
    pPlasticPorosity[*n] = o_pPlasticPorosity[*o];
  }

  (*newState)[pEqPlasticStrainLabel] = pEqPlasticStrain.clone();
  (*newState)[pPlasticPorosityLabel] = pPlasticPorosity.clone();
}

/*!-----------------------------------------------------*/
void
IntVar_MetalPlastic::allocateAndPutRigid(ParticleSubset* pset,
                                         DataWarehouse* new_dw,
                                         constParticleLabelVariableMap& var)
{
  ParticleVariable<double> pEqPlasticStrain_new, pPlasticPorosity_new;
  new_dw->allocateAndPut(
    pEqPlasticStrain_new, pEqPlasticStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(
    pPlasticPorosity_new, pPlasticPorosityLabel_preReloc, pset);
  for (auto pidx : *pset) {
    pEqPlasticStrain_new[pidx] = dynamic_cast<constParticleVariable<double>&>(
      *var[pEqPlasticStrainLabel])[pidx];
    pPlasticPorosity_new[pidx] = dynamic_cast<constParticleVariable<double>&>(
      *var[pPlasticPorosityLabel])[pidx];
  }
}
