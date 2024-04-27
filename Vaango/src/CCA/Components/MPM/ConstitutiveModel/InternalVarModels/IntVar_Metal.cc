/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_Metal.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Grid/Variables/MPMIntVarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;

using MetalIntVar = Uintah::MetalIntVar;

IntVar_Metal::IntVar_Metal(ProblemSpecP& ps)
{
  d_elastic = nullptr;
  d_shear   = nullptr;
  d_eos     = nullptr;
  initializeLocalMPMLabels();
  getPorosityModelParams(ps);
}

IntVar_Metal::IntVar_Metal(const IntVar_Metal* cm)
{
  d_elastic = cm->d_elastic;
  d_shear   = cm->d_shear;
  d_eos     = cm->d_eos;
  initializeLocalMPMLabels();
  copyPorosityModelParams(cm);
}

void
IntVar_Metal::initializeLocalMPMLabels()
{
  auto type_d =
    Uintah::ParticleVariable<Uintah::MetalIntVar>::getTypeDescription();
  pIntVarLabel          = Uintah::VarLabel::create("p.intVarMetal", type_d);
  pIntVarLabel_preReloc = Uintah::VarLabel::create("p.intVarMetal+", type_d);
}

void
IntVar_Metal::getPorosityModelParams(ProblemSpecP& ps)
{
  d_poreNucleation.phi_n = 0.1; // Volume fraction of void nucleating particles
  d_poreNucleation.eps_mean_n = 0.3; // Mean strain for nucleation
  d_poreNucleation.eps_std_n  = 0.1; // Standard deviation strain for nucleation
  ps->getWithDefault("vol_frac_nucleation", d_poreNucleation.phi_n, 0.1);
  ps->getWithDefault(
    "mean_strain_nucleation", d_poreNucleation.eps_mean_n, 0.3);
  ps->getWithDefault(
    "stddev_strain_nucleation", d_poreNucleation.eps_std_n, 0.1);
}

void
IntVar_Metal::copyPorosityModelParams(const IntVar_Metal* cm)
{
  d_poreNucleation.phi_n      = cm->d_poreNucleation.phi_n;
  d_poreNucleation.eps_mean_n = cm->d_poreNucleation.eps_mean_n;
  d_poreNucleation.eps_std_n  = cm->d_poreNucleation.eps_std_n;
}

IntVar_Metal::~IntVar_Metal()
{
  VarLabel::destroy(pIntVarLabel);
  VarLabel::destroy(pIntVarLabel_preReloc);
}

void
IntVar_Metal::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "metal_internal_var");
  int_var_ps->appendElement("vol_frac_nucleation", d_poreNucleation.phi_n);
  int_var_ps->appendElement("mean_strain_nucleation",
                            d_poreNucleation.eps_mean_n);
  int_var_ps->appendElement("stddev_strain_nucleation",
                            d_poreNucleation.eps_std_n);
}

/* get internal variable labels */
std::vector<const Uintah::VarLabel*>
IntVar_Metal::getLabels() const
{
  std::vector<const Uintah::VarLabel*> labels;
  labels.push_back(pIntVarLabel);
  labels.push_back(pIntVarLabel_preReloc);
  return labels;
}

void
IntVar_Metal::addParticleState(std::vector<const VarLabel*>& from,
                               std::vector<const VarLabel*>& to)
{
  from.push_back(pIntVarLabel);
  to.push_back(pIntVarLabel_preReloc);
}

void
IntVar_Metal::addInitialComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pIntVarLabel, matlset);
}

void
IntVar_Metal::initializeInternalVariable(Uintah::ParticleSubset* pset,
                                         Uintah::DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<Uintah::MetalIntVar> pIntVar;
  new_dw->allocateAndPut(pIntVar, pIntVarLabel, pset);

  for (auto pidx : *pset) {
    pIntVar[pidx] = { 0.0, 0.0 };
  }
}

void
IntVar_Metal::initializeInternalVariable(const Patch* /*patch*/,
                                         const MPMMaterial* /*matl*/,
                                         ParticleSubset* pset,
                                         DataWarehouse* new_dw,
                                         MPMLabel* /*lb*/,
                                         ParameterDict& /*params*/)
{
  initializeInternalVariable(pset, new_dw);
}

void
IntVar_Metal::addComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pIntVarLabel, matlset, Ghost::None);
  task->computes(pIntVarLabel_preReloc, matlset);
}

/* Get one (possibly composite) internal variable */
template <>
void
IntVar_Metal::getInternalVariable(ParticleSubset* pset,
                                  DataWarehouse* old_dw,
                                  constParticleVariable<MetalIntVar>& var)
{
  old_dw->get(var, pIntVarLabel, pset);
}

/* Allocate one (possibly composite) internal variable */
template <>
void
IntVar_Metal::allocateAndPutInternalVariable(ParticleSubset* pset,
                                             DataWarehouse* new_dw,
                                             ParticleVariable<MetalIntVar>& var)
{
  new_dw->allocateAndPut(var, pIntVarLabel, pset);
}

template <>
void
IntVar_Metal::evolveInternalVariable(
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
IntVar_Metal::computeInternalVariable(const std::string& label,
                                      const ModelStateBase* state_old,
                                      const ModelStateBase* state_cur) const
{
  if (label == "porosity") {
    return computePlasticPorosity(state_old->porosity, state_cur);
  } else if (label == "eqPlasticStrain") {
    return computeEqPlasticStrain(state_old->eqPlasticStrain, state_cur);
  }
  return 0.0;
}

double
IntVar_Metal::computeEqPlasticStrain(double eqPlasticStrain_old,
                                     const ModelStateBase* state) const
{
  double h_epdot = eqPlasticStrainHardeningModulus(state);
  return state->eqPlasticStrain + h_epdot * state->lambdaIncPlastic;
}

double
IntVar_Metal::computePlasticPorosity(double plasticPorosity_old,
                                     const ModelStateBase* state) const
{
  double h_phi = plasticPorosityHardeningModulus(state);
  return state->porosity + h_phi * state->lambdaIncPlastic;
}

template <>
void
IntVar_Metal::computeHardeningModulus(const ModelStateBase* state,
                                      MetalIntVar& hardeningModulus) const
{
  hardeningModulus.eqPlasticStrain = eqPlasticStrainHardeningModulus(state);
  hardeningModulus.plasticPorosity = plasticPorosityHardeningModulus(state);
}

double
IntVar_Metal::eqPlasticStrainHardeningModulus(const ModelStateBase* state) const
{
  return 1.0;
}

double
IntVar_Metal::plasticPorosityHardeningModulus(const ModelStateBase* state) const
{
  double h_grow = poreGrowthHardeningModulus(state);
  double h_nucl = poreNucleationHardeningModulus(state);
  return h_grow + h_nucl;
}

/* Calculate hardening modulus due to void growth */
double
IntVar_Metal::poreGrowthHardeningModulus(const ModelStateBase* state) const
{
  double h_grow = (1.0 - state->porosity) * state->plasticFlowDirection.Trace();
  return h_grow;
}

/* Calculate hardening modulus due to void nucleation */
double
IntVar_Metal::poreNucleationHardeningModulus(const ModelStateBase* state) const
{
  double A      = voidNucleationFactor(state->eqPlasticStrain);
  double h_nucl = Vaango::Util::sqrt_two_third * A;
  return h_nucl;
}

/* Calculate the void nucleation factor */
double
IntVar_Metal::voidNucleationFactor(double eqPlasticStrain) const
{
  double fac = (eqPlasticStrain - d_poreNucleation.eps_mean_n) /
               d_poreNucleation.eps_std_n;
  double A = d_poreNucleation.phi_n /
             (d_poreNucleation.eps_std_n * std::sqrt(2.0 * M_PI)) *
             std::exp(-0.5 * fac * fac);
  return A;
}

double
IntVar_Metal::computeVolStrainDerivOfInternalVariable(
  const std::string& label,
  const ModelStateBase* state) const
{
  return 0.0;
}

/*!-----------------------------------------------------*/
void
IntVar_Metal::allocateCMDataAddRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*,
                                        MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pIntVarLabel_preReloc, matlset, Ghost::None);
}

void
IntVar_Metal::allocateCMDataAdd(DataWarehouse* old_dw,
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
IntVar_Metal::allocateAndPutRigid(Uintah::ParticleSubset* pset,
                                  Uintah::DataWarehouse* new_dw,
                                  Uintah::constParticleVariableBase& var)
{
  ParticleVariable<MetalIntVar> pIntVar_new;
  new_dw->allocateAndPut(pIntVar_new, pIntVarLabel_preReloc, pset);
  for (auto pidx : *pset) {
    pIntVar_new[pidx] =
      static_cast<constParticleVariable<MetalIntVar>&>(var)[pidx];
  }
}

void
IntVar_Metal::allocateAndPutRigid(ParticleSubset* pset,
                                  DataWarehouse* new_dw,
                                  constParticleLabelVariableMap& var)
{
  ParticleVariable<MetalIntVar> pIntVar_new;
  new_dw->allocateAndPut(pIntVar_new, pIntVarLabel_preReloc, pset);
  for (auto pidx : *pset) {
    pIntVar_new[pidx] = dynamic_cast<constParticleVariable<MetalIntVar>&>(
      *var[pIntVarLabel])[pidx];
  }
}
