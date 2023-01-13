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

#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/DamageModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/FlowStressModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/ElasticPlasticHP.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondition/YieldConditionFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/MeltTempModels/MeltingTempModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/DeformationState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecHeatModels/SpecificHeatModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/StabilityModels/StabilityCheckFactory.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Gaussian.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/SymmMatrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <Core/Util/DebugStream.h>
#include <cmath>
#include <iostream>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;
using Vaango::ModelStateBase;

static DebugStream cout_EP("EP", false);
static DebugStream cout_EP1("EP1", false);
static DebugStream CSTi("EPi", false);
static DebugStream CSTir("EPir", false);

ElasticPlasticHP::ElasticPlasticHP(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("bulk_modulus", d_initialData.Bulk);
  ps->require("shear_modulus", d_initialData.Shear);

  d_initialData.alpha = 0.0; // default is per K.  Only used in implicit code
  ps->get("coeff_thermal_expansion", d_initialData.alpha);
  d_initialData.Chi = 0.9;
  ps->get("taylor_quinney_coeff", d_initialData.Chi);
  d_initialData.sigma_crit = 2.0e9; // default is Pa
  ps->get("critical_stress", d_initialData.sigma_crit);

  d_tol = 1.0e-10;
  ps->get("tolerance", d_tol);

  d_useModifiedEOS = false;
  ps->get("useModifiedEOS", d_useModifiedEOS);

  d_initialMaterialTemperature = 294.0;
  ps->get("initial_material_temperature", d_initialMaterialTemperature);

  d_checkTeplaFailureCriterion = true;
  ps->get("check_TEPLA_failure_criterion", d_checkTeplaFailureCriterion);

  d_doMelting = true;
  ps->get("do_melting", d_doMelting);

  d_checkStressTriax = true;
  ps->get("check_max_stress_failure", d_checkStressTriax);

  d_yield = YieldConditionFactory::create(ps);
  if (!d_yield) {
     std::ostringstream desc;
    desc << "An error occured in the YieldConditionFactory that has \n"
         << " slipped through the existing bullet proofing. Please tell \n"
         << " Biswajit.  "
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_stable = StabilityCheckFactory::create(ps);
  if (!d_stable)
    std::cerr << "Stability check disabled\n";

  d_flow = FlowStressModelFactory::create(ps);
  if (!d_flow) {
     std::ostringstream desc;
    desc << "An error occured in the FlowModelFactory that has \n"
         << " slipped through the existing bullet proofing. Please tell \n"
         << " Biswajit.  "
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_damage = DamageModelFactory::create(ps);
  if (!d_damage) {
     std::ostringstream desc;
    desc << "An error occured in the DamageModelFactory that has \n"
         << " slipped through the existing bullet proofing. Please tell \n"
         << " Biswajit.  "
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_eos = MPMEquationOfStateFactory::create(ps);
  d_eos->setBulkModulus(d_initialData.Bulk);
  if (!d_eos) {
     std::ostringstream desc;
    desc << "An error occured in the EquationOfStateFactory that has \n"
         << " slipped through the existing bullet proofing. Please tell \n"
         << " Jim.  "
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_shear = Vaango::ShearModulusModelFactory::create(ps);
  if (!d_shear) {
     std::ostringstream desc;
    desc << "ElasticPlasticHP::Error in shear modulus model factory"
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_melt = MeltingTempModelFactory::create(ps);
  if (!d_melt) {
     std::ostringstream desc;
    desc << "ElasticPlasticHP::Error in melting temp model factory"
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_computeSpecificHeat = false;
  ps->get("compute_specific_heat", d_computeSpecificHeat);
  d_Cp = SpecificHeatModelFactory::create(ps);

  setErosionAlgorithm();
  getInitialPorosityData(ps);
  getInitialDamageData(ps);
  initializeLocalMPMLabels();
}

ElasticPlasticHP::ElasticPlasticHP(const ElasticPlasticHP* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.Bulk       = cm->d_initialData.Bulk;
  d_initialData.Shear      = cm->d_initialData.Shear;
  d_initialData.alpha      = cm->d_initialData.alpha;
  d_initialData.Chi        = cm->d_initialData.Chi;
  d_initialData.sigma_crit = cm->d_initialData.sigma_crit;

  d_tol            = cm->d_tol;
  d_useModifiedEOS = cm->d_useModifiedEOS;

  d_initialMaterialTemperature = cm->d_initialMaterialTemperature;
  d_checkTeplaFailureCriterion = cm->d_checkTeplaFailureCriterion;
  d_doMelting                  = cm->d_doMelting;
  d_checkStressTriax           = cm->d_checkStressTriax;

  d_setStressToZero = cm->d_setStressToZero;
  d_allowNoTension  = cm->d_allowNoTension;
  d_allowNoShear    = cm->d_allowNoShear;

  d_evolvePorosity        = cm->d_evolvePorosity;
  d_porosity.f0           = cm->d_porosity.f0;
  d_porosity.f0_std       = cm->d_porosity.f0_std;
  d_porosity.fc           = cm->d_porosity.fc;
  d_porosity.fn           = cm->d_porosity.fn;
  d_porosity.en           = cm->d_porosity.en;
  d_porosity.sn           = cm->d_porosity.sn;
  d_porosity.porosityDist = cm->d_porosity.porosityDist;

  d_evolveDamage               = cm->d_evolveDamage;
  d_scalarDam.D0               = cm->d_scalarDam.D0;
  d_scalarDam.D0_std           = cm->d_scalarDam.D0_std;
  d_scalarDam.Dc               = cm->d_scalarDam.Dc;
  d_scalarDam.scalarDamageDist = cm->d_scalarDam.scalarDamageDist;

  d_computeSpecificHeat = cm->d_computeSpecificHeat;

  d_Cp     = SpecificHeatModelFactory::createCopy(cm->d_Cp);
  d_yield  = YieldConditionFactory::createCopy(cm->d_yield);
  d_stable = StabilityCheckFactory::createCopy(cm->d_stable);
  d_flow   = FlowStressModelFactory::createCopy(cm->d_flow);
  d_damage = DamageModelFactory::createCopy(cm->d_damage);
  d_eos    = MPMEquationOfStateFactory::createCopy(cm->d_eos);
  d_eos->setBulkModulus(d_initialData.Bulk);
  d_shear     = Vaango::ShearModulusModelFactory::createCopy(cm->d_shear);
  d_melt      = MeltingTempModelFactory::createCopy(cm->d_melt);

  initializeLocalMPMLabels();
}

ElasticPlasticHP::~ElasticPlasticHP()
{
  // Destructor
  VarLabel::destroy(pStrainRateLabel);
  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pPlasticStrainRateLabel);
  VarLabel::destroy(pDamageLabel);
  VarLabel::destroy(pPorosityLabel);
  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pEnergyLabel);

  VarLabel::destroy(pStrainRateLabel_preReloc);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);
  VarLabel::destroy(pPlasticStrainRateLabel_preReloc);
  VarLabel::destroy(pDamageLabel_preReloc);
  VarLabel::destroy(pPorosityLabel_preReloc);
  VarLabel::destroy(pLocalizedLabel_preReloc);
  VarLabel::destroy(pEnergyLabel_preReloc);

  delete d_flow;
  delete d_yield;
  delete d_stable;
  delete d_damage;
  delete d_eos;
  delete d_shear;
  delete d_melt;
  delete d_Cp;
}

//______________________________________________________________________
//
void
ElasticPlasticHP::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "elastic_plastic_hp");
  }

  cm_ps->appendElement("bulk_modulus", d_initialData.Bulk);
  cm_ps->appendElement("shear_modulus", d_initialData.Shear);
  cm_ps->appendElement("coeff_thermal_expansion", d_initialData.alpha);
  cm_ps->appendElement("taylor_quinney_coeff", d_initialData.Chi);
  cm_ps->appendElement("critical_stress", d_initialData.sigma_crit);
  cm_ps->appendElement("tolerance", d_tol);
  cm_ps->appendElement("useModifiedEOS", d_useModifiedEOS);
  cm_ps->appendElement("initial_material_temperature",
                       d_initialMaterialTemperature);
  cm_ps->appendElement("check_TEPLA_failure_criterion",
                       d_checkTeplaFailureCriterion);
  cm_ps->appendElement("do_melting", d_doMelting);
  cm_ps->appendElement("check_max_stress_failure", d_checkStressTriax);
  cm_ps->appendElement("compute_specific_heat", d_computeSpecificHeat);

  d_yield->outputProblemSpec(cm_ps);
  d_stable->outputProblemSpec(cm_ps);
  d_flow->outputProblemSpec(cm_ps);
  d_damage->outputProblemSpec(cm_ps);
  d_eos->outputProblemSpec(cm_ps);
  d_shear->outputProblemSpec(cm_ps);
  d_melt->outputProblemSpec(cm_ps);
  d_Cp->outputProblemSpec(cm_ps);

  cm_ps->appendElement("evolve_porosity", d_evolvePorosity);
  cm_ps->appendElement("initial_mean_porosity", d_porosity.f0);
  cm_ps->appendElement("initial_std_porosity", d_porosity.f0_std);
  cm_ps->appendElement("critical_porosity", d_porosity.fc);
  cm_ps->appendElement("frac_nucleation", d_porosity.fn);
  cm_ps->appendElement("meanstrain_nucleation", d_porosity.en);
  cm_ps->appendElement("stddevstrain_nucleation", d_porosity.sn);
  cm_ps->appendElement("initial_porosity_distrib", d_porosity.porosityDist);
  cm_ps->appendElement("porosity_seed", d_porosity.seed);

  cm_ps->appendElement("evolve_damage", d_evolveDamage);
  cm_ps->appendElement("initial_mean_scalar_damage", d_scalarDam.D0);
  cm_ps->appendElement("initial_std_scalar_damage", d_scalarDam.D0_std);
  cm_ps->appendElement("critical_scalar_damage", d_scalarDam.Dc);
  cm_ps->appendElement("initial_scalar_damage_distrib",
                       d_scalarDam.scalarDamageDist);
  cm_ps->appendElement("scalar_damage_seed", d_scalarDam.seed);
}

std::unique_ptr<ConstitutiveModel>
ElasticPlasticHP::clone()
{
  return std::make_unique<ElasticPlasticHP>(*this);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::initializeLocalMPMLabels()
{
  pRotationLabel = VarLabel::create(
    "p.rotation", ParticleVariable<Matrix3>::getTypeDescription());
  pStrainRateLabel = VarLabel::create(
    "p.strainRate", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel = VarLabel::create(
    "p.plasticStrain", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainRateLabel = VarLabel::create(
    "p.plasticStrainRate", ParticleVariable<double>::getTypeDescription());
  pDamageLabel = VarLabel::create(
    "p.damage", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel = VarLabel::create(
    "p.porosity", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pEnergyLabel = VarLabel::create(
    "p.energy", ParticleVariable<double>::getTypeDescription());

  pRotationLabel_preReloc = VarLabel::create(
    "p.rotation+", ParticleVariable<Matrix3>::getTypeDescription());
  pStrainRateLabel_preReloc = VarLabel::create(
    "p.strainRate+", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel_preReloc = VarLabel::create(
    "p.plasticStrain+", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainRateLabel_preReloc = VarLabel::create(
    "p.plasticStrainRate+", ParticleVariable<double>::getTypeDescription());
  pDamageLabel_preReloc = VarLabel::create(
    "p.damage+", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = VarLabel::create(
    "p.porosity+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
  pEnergyLabel_preReloc = VarLabel::create(
    "p.energy+", ParticleVariable<double>::getTypeDescription());
}
//______________________________________________________________________
//
void
ElasticPlasticHP::getInitialPorosityData(ProblemSpecP& ps)
{
  d_evolvePorosity = true;
  ps->get("evolve_porosity", d_evolvePorosity);
  d_porosity.f0           = 0.002; // Initial porosity
  d_porosity.f0_std       = 0.002; // Initial STD porosity
  d_porosity.fc           = 0.5;   // Critical porosity
  d_porosity.fn           = 0.1; // Volume fraction of void nucleating particles
  d_porosity.en           = 0.3; // Mean strain for nucleation
  d_porosity.sn           = 0.1; // Standard deviation strain for nucleation
  d_porosity.porosityDist = "constant";
  ps->get("initial_mean_porosity", d_porosity.f0);
  ps->get("initial_std_porosity", d_porosity.f0_std);
  ps->get("critical_porosity", d_porosity.fc);
  ps->get("frac_nucleation", d_porosity.fn);
  ps->get("meanstrain_nucleation", d_porosity.en);
  ps->get("stddevstrain_nucleation", d_porosity.sn);
  ps->get("initial_porosity_distrib", d_porosity.porosityDist);
  ps->get("porosity_seed", d_porosity.seed);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::getInitialDamageData(ProblemSpecP& ps)
{
  d_evolveDamage = true;
  ps->get("evolve_damage", d_evolveDamage);
  d_scalarDam.D0               = 0.0; // Initial scalar damage
  d_scalarDam.D0_std           = 0.0; // Initial STD scalar damage
  d_scalarDam.Dc               = 1.0; // Critical scalar damage
  d_scalarDam.scalarDamageDist = "constant";
  ps->get("initial_mean_scalar_damage", d_scalarDam.D0);
  ps->get("initial_std_scalar_damage", d_scalarDam.D0_std);
  ps->get("critical_scalar_damage", d_scalarDam.Dc);
  ps->get("initial_scalar_damage_distrib", d_scalarDam.scalarDamageDist);
  ps->get("scalar_damage_seed", d_scalarDam.seed);
}

/*! Compute specific heat

    double T = temperature;
    C_p = 1.0e3*(A + B*T + C/T^2)
    ** For steel **
    C_p = 1.0e3*(0.09278 + 7.454e-4*T + 12404.0/(T*T));
*/

//______________________________________________________________________
//
void
ElasticPlasticHP::setErosionAlgorithm()
{
  d_setStressToZero = false;
  d_allowNoTension  = false;
  d_allowNoShear    = false;
  if (flag->d_doErosion) {
    if (flag->d_erosionAlgorithm == "AllowNoTension")
      d_allowNoTension = true;
    else if (flag->d_erosionAlgorithm == "AllowNoShear")
      d_allowNoShear = true;
    else if (flag->d_erosionAlgorithm == "ZeroStress")
      d_setStressToZero = true;
  }
}
//______________________________________________________________________
//
void
ElasticPlasticHP::addParticleState(std::vector<const VarLabel*>& from,
                                   std::vector<const VarLabel*>& to)
{
  // This is an INCREMENTAL model. Needs polar decomp R to be saved.
  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);

  // Add the local particle state data for this constitutive model.
  from.push_back(pRotationLabel);
  from.push_back(pStrainRateLabel);
  from.push_back(pPlasticStrainLabel);
  from.push_back(pPlasticStrainRateLabel);
  from.push_back(pDamageLabel);
  from.push_back(pPorosityLabel);
  from.push_back(pLocalizedLabel);
  from.push_back(pEnergyLabel);

  to.push_back(pRotationLabel_preReloc);
  to.push_back(pStrainRateLabel_preReloc);
  to.push_back(pPlasticStrainLabel_preReloc);
  to.push_back(pPlasticStrainRateLabel_preReloc);
  to.push_back(pDamageLabel_preReloc);
  to.push_back(pPorosityLabel_preReloc);
  to.push_back(pLocalizedLabel_preReloc);
  to.push_back(pEnergyLabel_preReloc);

  // Add the particle state for the flow & deviatoric stress model
  d_flow->addParticleState(from, to);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::addInitialComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  task->computes(pStrainRateLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(pPlasticStrainRateLabel, matlset);
  task->computes(pDamageLabel, matlset);
  task->computes(pPorosityLabel, matlset);
  task->computes(pLocalizedLabel, matlset);
  task->computes(pEnergyLabel, matlset);

  // Add internal evolution variables computed by flow & deviatoric stress model
  d_flow->addInitialComputesAndRequires(task, matl, patch);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::initializeCMData(const Patch* patch,
                                   const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<double> pPlasticStrain, pDamage, pPorosity,
    pPlasticStrainRate, pStrainRate, pEnergy;
  ParticleVariable<int> pLocalized;

  new_dw->allocateAndPut(pStrainRate, pStrainRateLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticStrainRate, pPlasticStrainRateLabel, pset);
  new_dw->allocateAndPut(pDamage, pDamageLabel, pset);
  new_dw->allocateAndPut(pLocalized, pLocalizedLabel, pset);
  new_dw->allocateAndPut(pPorosity, pPorosityLabel, pset);
  new_dw->allocateAndPut(pEnergy, pEnergyLabel, pset);

  for (auto pidx : *pset) {
    pStrainRate[pidx]        = 0.0;
    pPlasticStrain[pidx]     = 0.0;
    pPlasticStrainRate[pidx] = 0.0;
    pDamage[pidx]            = d_damage->initialize();
    pPorosity[pidx]          = d_porosity.f0;
    pLocalized[pidx]         = 0;
    pEnergy[pidx]            = 0.;
  }

  // Do some extra things if the porosity or the damage distribution
  // is not uniform.
  if (d_porosity.porosityDist != "constant") {

    // Generate a Gaussian distributed random number given the mean
    // porosity and the std.
    auto patchID = patch->getID();
    auto patchID_mod = (patchID / 32) % 32;
    auto seed = ((d_porosity.seed + patchID_mod + 1) << patchID);
    Uintah::Gaussian gaussGen(d_porosity.f0, d_porosity.f0_std, seed, 1, DBL_MAX);
    for (auto pidx : *pset) {
      pPorosity[pidx] = std::abs(gaussGen.rand(1.0));
    }
  }

  if (d_scalarDam.scalarDamageDist != "constant") {

    // Generate a Gaussian distributed random number given the mean
    // damage and the std.
    auto patchID = patch->getID();
    auto patchID_mod = (patchID / 32) % 32;
    auto seed = ((d_porosity.seed + patchID_mod + 1) << patchID);
    Uintah::Gaussian gaussGen(
      d_scalarDam.D0, d_scalarDam.D0_std, seed, 1, DBL_MAX);
    for (auto pidx : *pset) {
      pDamage[pidx] = std::abs(gaussGen.rand(1.0));
    }
  }

  // Initialize the data for the flow model
  d_flow->initializeInternalVars(pset, new_dw);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::computeStableTimestep(const Patch* patch,
                                        const MPMMaterial* matl,
                                        DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;
  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  // Compute wave speed at each particle, store the maximum
  double modulus = initialData.Bulk + 4.0 / 3.0 * initialData.Shear;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (auto idx : *pset) {

    Vector pvelocity_idx = pVelocity[idx];
    if (pMass[idx] > 0) {
      double c_dil = std::sqrt(modulus * pVolume[idx] / pMass[idx]);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);
    } 
  }

  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

//______________________________________________________________________
//
void
ElasticPlasticHP::addComputesAndRequires(Task* task,
                                         const MPMMaterial* matl,
                                         const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  Ghost::GhostType gnone        = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);
  task->requires(Task::OldDW, lb->pPolarDecompRLabel, matlset, gnone);
  task->requires(Task::NewDW, lb->pPolarDecompRLabel_preReloc, matlset, gnone);

  // Other constitutive model and input dependent computes and requires
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);

  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);
  task->requires(Task::OldDW, pStrainRateLabel, matlset, gnone);
  task->requires(Task::OldDW, pPlasticStrainLabel, matlset, gnone);
  task->requires(Task::OldDW, pPlasticStrainRateLabel, matlset, gnone);
  task->requires(Task::OldDW, pDamageLabel, matlset, gnone);
  task->requires(Task::OldDW, pPorosityLabel, matlset, gnone);
  task->requires(Task::OldDW, pLocalizedLabel, matlset, gnone);
  task->requires(Task::OldDW, pEnergyLabel, matlset, gnone);

  task->computes(pStrainRateLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(pPlasticStrainRateLabel_preReloc, matlset);
  task->computes(pDamageLabel_preReloc, matlset);
  task->computes(pPorosityLabel_preReloc, matlset);
  task->computes(pLocalizedLabel_preReloc, matlset);
  task->computes(pEnergyLabel_preReloc, matlset);

  // Add internal evolution variables computed by flow model
  d_flow->addComputesAndRequires(task, matl, patches);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::computeStressTensor(const PatchSubset* patches,
                                      const MPMMaterial* matl,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  if (cout_EP.active()) {
    cout_EP << getpid() << " ElasticPlasticHP:ComputeStressTensor:Explicit"
            << " Matl = " << matl << " DWI = " << matl->getDWIndex()
            << " patch = " << (patches->get(0))->getID();
  }

  int matID = matl->getDWIndex();

  // Get the time increment (delT)
  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  // Viscous heating
  double include_AV_heating = 0.0;
  if (flag->d_artificialViscosityHeating) {
    include_AV_heating = 1.0;
  }

  // General stuff
  Matrix3 pDefRate_mid[idx](0.0); // Rate of deformation
  Matrix3 tensorW(0.0); // Spin
  Matrix3 tensorF;
  tensorF.Identity(); // Deformation gradient
  Matrix3 tensorU;
  tensorU.Identity(); // Right Cauchy-Green stretch
  Matrix3 tensorR;
  tensorR.Identity();     // Rotation
  Matrix3 sigma(0.0);     // The Cauchy stress
  Matrix3 devD(0.0); // Deviatoric part of tensor D
  Matrix3 tensorS(0.0);   // Devaitoric part of tensor Sig
  Matrix3 tensorF_new;
  tensorF_new.Identity(); // Deformation gradient


  double bulk         = d_initialData.Bulk;
  double shear        = d_initialData.Shear;
  double rho_0        = matl->getInitialDensity();
  double Tm           = matl->getMeltTemperature();
  double sqrtThreeTwo = sqrt(1.5);
  double Vaango::Util::sqrt_two_third = 1.0 / sqrtThreeTwo;

  double totalStrainEnergy = 0.0;

  // Loop thru patches
  for (int patchIndex = 0; patchIndex < patches->size(); patchIndex++) {
    const Patch* patch = patches->get(patchIndex);
    Vector dx = patch->dCell();
    double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;

    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Get particle data
    constParticleVariable<int> pLocalized;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass, pVolume_old, pTemperature_old;
                                  pPlasticStrain_old, pDamage_old, 
                                  pPorosity_old;
                                  pStrainRate_old, pPlasticStrainRate_old, 
                                  pEnergy_old;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pDefGrad_old, pRotation_old;

    old_dw->get(pLocalized, pLocalizedLabel, pset);
    old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume_old, lb->pVolumeLabel, pset);
    old_dw->get(pTemperature_old, lb->pTemperatureLabel, pset);
    old_dw->get(pPlasticStrain_old, pPlasticStrainLabel, pset);
    old_dw->get(pDamage_old, pDamageLabel, pset);
    old_dw->get(pStrainRate_old, pStrainRateLabel, pset);
    old_dw->get(pPlasticStrainRate_old, pPlasticStrainRateLabel, pset);
    old_dw->get(pPorosity_old, pPorosityLabel, pset);
    old_dw->get(pEnergy_old, pEnergyLabel, pset);
    old_dw->get(pRotation_old, lb->pPolarDecompRLabel, pset);
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

    constParticleVariable<double> pVolume_new;
    constParticleVariable<Matrix3> pRotation_new, pDefRate_mid, 
                                   pDefGrad_new, pStress_old;

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pRotation_new, lb->pPolarDecompRLabel_preReloc, pset);
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);

    ParticleVariable<int> pLocalized_new;
    ParticleVariable<double> pPlasticStrain_new, pDamage_new, pPorosity_new;
    ParticleVariable<double> pStrainRate_new, pPlasticStrainRate_new;
    ParticleVariable<double> pdTdt, p_q, pEnergy_new;
    ParticleVariable<Matrix3> pStress_new;

    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pDamage_new, pDamageLabel_preReloc, pset);
    new_dw->allocateAndPut(pPorosity_new, pPorosityLabel_preReloc, pset);
    new_dw->allocateAndPut(pStrainRate_new, pStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrainRate_new, pPlasticStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pEnergy_new, pEnergyLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    d_flow->getInternalVars(pset, old_dw);

    d_flow->allocateAndPutInternalVars(pset, new_dw);

    //______________________________________________________________________
    // Loop thru particles
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    for (auto idx: *pset) {

      /* Delete particle after setting everything to zero if 
         unphysical deformation gradient */
      double J = pDefGrad_new[idx].Determinant();
      if (!(J > 0.) || J > 1.e5) {
        std::cerr
          << "**ERROR** Negative (or huge) Jacobian of deformation gradient."
          << "  Deleting particle " << pParticleID[idx] << "\n";
        std::cerr << "matID = " << matl->getDWIndex() << "\n";
        std::cerr << "J = " << J << "\n";
        std::cerr << "F_old = " << pDefGrad_old[idx] << "\n";
        std::cerr << "F_new = " << pDefGrad_new[idx] << "\n";
        std::cerr << "Temp = " << pTemperature_old[idx] << "\n";
        std::cerr << "Tm = " << Tm << "\n";
        pLocalized_new[idx] = -999;
        pPlasticStrain_new[idx] = 0.0;
        pDamage_new[idx] = 0.0;
        pPorosity_new[idx] = 0.0;
        pStrainRate_new[idx] = 0.0;
        pPlasticStrainRate_new[idx] = 0.0;
        pdTdt[idx] = 0.0;
        p_q[idx] = 0.0;
        pEnergy_new[idx] = 0.0;
        pStress_new[idx] = Vaango::Util::Zero;
        auto defState = scinew DeformationState();
        defState->D    = Vaango::Util::Zero;
        defState->devD = Vaango::Util::Zero;
        continue;
      } 

      // Assign zero int. heating by default, modify with appropriate sources
      // This has units (in MKS) of K/s  (i.e. temperature/time)
      pdTdt[idx] = 0.0;

      //-----------------------------------------------------------------------
      // Stage 1:
      //-----------------------------------------------------------------------
      // Carry forward the pLocalized tag for now, alter below
      pLocalized_new[idx] = pLocalized[idx];

      // Calculate the current density and old temperature
      double rho_new = rho_0 / J;
      double temperature = pTemperature_old[idx];

      // Compute rate of change of specific volume
      double dv_dt =
        (pVolume_new[idx] - pVolume_old[idx]) / (pMass[idx] * delT);

      // Get rate of deformation and compute scalar strain rate
      Matrix3 DD = pDefRate_mid[idx];
      double trD = DD.Trace();
      Matrix3 devD = DD - Vaango::Util::Identity * (trD / 3.0);
      pStrainRate_new[idx] = Vaango::Util::sqrt_two_third * DD.Norm();

      // Deformation state for deviatoric stress (viscoelastic)
      auto defState = std::make_unique<DeformationState>();
      defState->D                    = DD;
      defState->devD                 = devD;
      defState->viscoElasticWorkRate = 0.0;

      // Set up the ModelStateBase (for t_n+1)
      auto state = std::make_unique<ModelStateBase>();
      state->strainRate          = pStrainRate_new[idx];
      state->plasticStrainRate   = pPlasticStrainRate_old[idx];
      state->plasticStrain       = pPlasticStrain_old[idx];
      state->temperature         = temperature_old;
      state->initialTemperature  = d_initialMaterialTemperature;
      state->density             = rho_new;
      state->initialDensity      = rho_0;
      state->volume              = pVolume_new[idx];
      state->initialVolume       = pMass[idx] / rho_0;
      state->bulkModulus         = bulk;
      state->initialBulkModulus  = bulk;
      state->shearModulus        = d_shear->computeShearModulus(state);
      state->initialShearModulus = shear;
      state->meltingTemp         = d_melt->computeMeltingTemp(state);
      state->initialMeltTemp     = Tm;
      state->energy              = pEnergy_old[idx];
      if (d_computeSpecificHeat) {
        double C_p          = d_Cp->computeSpecificHeat(state);
        state->specificHeat = C_p;
      } else {
        state->specificHeat = matl->getSpecificHeat();
      }

      // Get volumetric and deviatoric stress
      Matrix3 sigma_old = pStress_old[idx];
      Matrix3 devSigma_old = state->updateStressInvariants(sigma_old);
      double pressure_old = state->pressure;

      //-----------------------------------------------------------------------
      // Stage 2:
      //-----------------------------------------------------------------------
      // Assume elastic deformation to get a trial deviatoric stress
      // This is simply the previous timestep deviatoric stress plus a
      // deviatoric elastic increment based on the shear modulus supplied by
      // the strength routine in use.

      Matrix3 trialS = tensorS;

      // Calculate the equivalent stress
      // this will be removed next, it should be computed in the flow stress
      // routine
      // the flow stress routines should be passed the entire stress (not just
      // deviatoric)
      double equivStress = sqrtThreeTwo * trialS.Norm();

      // Calculate flow stress
      double flowStress =
        d_flow->computeFlowStress(state, delT, d_tol, matl, idx);
      state->yieldStress = flowStress;

      // Material has melted if flowStress <= 0.0
      bool melted  = false;
      bool plastic = false;
      if (temperature > Tm_cur || flowStress <= 0.0) {

        melted = true;
        if (d_doMelting) {
          // Set the deviatoric stress to zero
          tensorS = Vaango::Util::Zero;
        }

        d_flow->updateElastic(idx);

      } else {

        // Get the current porosity
        double porosity = pPorosity[idx];

        // Evaluate yield condition
        double traceOfTrialStress =
          3.0 * pressure + trD * (2.0 * mu_cur * delT);

        double flow_rule = d_yield->evalYieldCondition(equivStress,
                                                       flowStress,
                                                       traceOfTrialStress,
                                                       porosity,
                                                       state->yieldStress);
        // Compute the deviatoric stress
        /*
        std::cout << "flow_rule = " << flow_rule << " s_eq = " << equivStress
             << " s_flow = " << flowStress << "\n";
        */

        if (flow_rule < 0.0) {
          // Set the deviatoric stress to the trial stress
          tensorS = trialS;

          // Update the internal variables
          d_flow->updateElastic(idx);

        } else {

          plastic = true;

          double delGamma = 0.0;
          double normS    = tensorS.Norm();

          // If the material goes plastic in the first step, or
          // gammadotplus < 0 or delGamma < 0 use the Simo algorithm
          // with Newton iterations.

          // Compute Stilde using Newton iterations a la Simo
          state->plasticStrainRate = pStrainRate_new[idx];
          state->plasticStrain     = pPlasticStrain[idx];
          Matrix3 nn(0.0);
          computePlasticStateViaRadialReturn(
            trialS, delT, matl, idx, state, nn, delGamma);

          devDPlasticInc = nn * delGamma;
          tensorS =
            trialS - devDPlasticInc * (2.0 * state->shearModulus);
         

          // Update internal variables
          d_flow->updatePlastic(idx, delGamma);

        } // end of flow_rule if
      }   // end of temperature if

      // Calculate the updated hydrostatic stress
      double p =
        d_eos->computePressure(matl, state, tensorF_new, DD, delT);

      double Dkk = trD;
      /*
       // **WARNING** Produces negative Tdot
      double dTdt_isentropic = d_eos->computeIsentropicTemperatureRate(
        temperature, rho_0, rho_new, Dkk);
      pdTdt[idx] += dTdt_isentropic;
      */

      // Calculate Tdot from viscoelasticity
      double taylorQuinney = d_initialData.Chi;
      double fac           = taylorQuinney / (rho_new * state->specificHeat);
      double Tdot_VW       = defState->viscoElasticWorkRate * fac;

      pdTdt[idx] += Tdot_VW;

      double de_s = 0.;
      if (flag->d_artificialViscosity) {
        double c_bulk = sqrt(bulk / rho_new);
        p_q[idx]      = artificialBulkViscosity(Dkk, c_bulk, rho_new, dx_ave);
        de_s          = -p_q[idx] * Dkk / rho_new;
      } else {
        p_q[idx] = 0.;
        de_s     = 0.;
      }

      // Calculate Tdot due to artificial viscosity
      double Tdot_AV = de_s / state->specificHeat;
      pdTdt[idx] += Tdot_AV * include_AV_heating;

      if (pdTdt[idx] < 0.0) {
        std::cout << "dTdt = " << pdTdt[idx]
                  //<< " dTdT_isen = " << dTdt_isentropic
                  << " dTdT_plas = " << Tdot_VW << " dTdT_visc = " << Tdot_AV
                  << "\n";
      }

      Matrix3 tensorHy = Vaango::Util::Identity * p;

      // Calculate the total stress
      sigma = tensorS + tensorHy;

      // If the particle has already failed, apply various erosion algorithms
      if (flag->d_doErosion) {
        if (pLocalized[idx]) {
          if (d_allowNoTension) {
            if (p > 0.0) {
              sigma = Vaango::Util::Zero;
            } else {
              sigma = tensorHy;
            }
          }
          if (d_allowNoShear) {
            sigma = tensorHy;
          } else if (d_setStressToZero) {
            sigma = Vaango::Util::Zero;
          }
        }
      }

      //-----------------------------------------------------------------------
      // Stage 3:
      //-----------------------------------------------------------------------
      // Compute porosity/damage/temperature change
      if (!plastic) {

        // Save the updated data
        pPlasticStrain_new[idx]     = pPlasticStrain[idx];
        pPlasticStrainRate_new[idx] = 0.0;
        pDamage_new[idx]            = pDamage[idx];
        pPorosity_new[idx]          = pPorosity[idx];

      } else {

        // Update the plastic strain
        pPlasticStrain_new[idx]     = state->plasticStrain;
        pPlasticStrainRate_new[idx] = state->plasticStrainRate;

        /*
        if (pPlasticStrainRate_new[idx] > pStrainRate_new[idx]) {
          std::cout << "Patch = " << patch->getID() << " particle = " << idx
               << " edot = " << pStrainRate_new[idx]
               << " epdot = " << pPlasticStrainRate_new[idx] << "\n";
        }
        */

        // Update the porosity
        if (d_evolvePorosity) {
          pPorosity_new[idx] =
            updatePorosity(DD, delT, pPorosity[idx], state->plasticStrain);
        } else {
          pPorosity_new[idx] = pPorosity[idx];
        }

        // Calculate the updated scalar damage parameter
        if (d_evolveDamage) {
          pDamage_new[idx] =
            d_damage->computeScalarDamage(state->plasticStrainRate,
                                          sigma,
                                          temperature,
                                          delT,
                                          matl,
                                          d_tol,
                                          pDamage[idx]);
        } else {
          pDamage_new[idx] = pDamage[idx];
        }
        // Calculate rate of temperature increase due to plastic strain
        double taylorQuinney = d_initialData.Chi;
        double fac           = taylorQuinney / (rho_new * state->specificHeat);

        // Calculate Tdot (internal plastic heating rate)
        double Tdot_PW = state->yieldStress * state->plasticStrainRate * fac;

        pdTdt[idx] += Tdot_PW;
      }

      //-----------------------------------------------------------------------
      // Stage 4:
      //-----------------------------------------------------------------------
      // Find if the particle has failed/localized
      bool isLocalized = false;
      double tepla     = 0.0;

      if (flag->d_doErosion) {

        // Check 1: Look at the temperature
        if (melted)
          isLocalized = true;

        // Check 2 and 3: Look at TEPLA and stability
        else if (plastic) {

          // Check 2: Modified Tepla rule
          if (d_checkTeplaFailureCriterion) {
            tepla = pow(pPorosity_new[idx] / d_porosity.fc, 2.0) +
                    pow(pDamage_new[idx], 2.0);
            if (tepla > 1.0)
              isLocalized = true;
          }

          // Check 3: Stability criterion (only if material is plastic)
          if (d_stable->doIt() && !isLocalized) {

            // Calculate values needed for tangent modulus calculation
            state->temperature  = temperature;
            Tm_cur              = d_melt->computeMeltingTemp(state);
            state->meltingTemp  = Tm_cur;
            mu_cur              = d_shear->computeShearModulus(state);
            state->shearModulus = mu_cur;
            double sigY =
              d_flow->computeFlowStress(state, delT, d_tol, matl, idx);
            if (!(sigY > 0.0))
              isLocalized = true;
            else {
              double dsigYdep =
                d_flow->evalDerivativeWRTPlasticStrain(state, idx);
              double A = voidNucleationFactor(state->plasticStrain);

              // Calculate the elastic tangent modulus
              TangentModulusTensor Ce;
              computeElasticTangentModulus(bulk, mu_cur, Ce);

              // Calculate the elastic-plastic tangent modulus
              TangentModulusTensor Cep;
              d_yield->computeElasPlasTangentModulus(
                Ce, sigma, sigY, dsigYdep, pPorosity_new[idx], A, Cep);

              // Initialize localization direction
              Vector direction(0.0, 0.0, 0.0);
              isLocalized =
                d_stable->checkStability(sigma, DD, Cep, direction);
            }
          }
        }

        // Check 4: Look at maximum stress
        if (d_checkStressTriax) {

          // Compute eigenvalues of the stress tensor
          SymmMatrix3 stress(sigma);
          Vector eigVal(0.0, 0.0, 0.0);
          Matrix3 eigVec;
          stress.eigen(eigVal, eigVec);

          double max_stress = Max(Max(eigVal[0], eigVal[1]), eigVal[2]);
          if (max_stress > d_initialData.sigma_crit) {
            isLocalized = true;
          }
        }

        // Use erosion algorithms to treat newly localized particles
        if (isLocalized) {

          // If the localized particles fail again then set their stress to zero
          if (pLocalized[idx]) {
            pDamage_new[idx]   = 0.0;
            pPorosity_new[idx] = 0.0;
          } else {
            // set the particle localization flag to true
            pLocalized_new[idx] = 1;
            pDamage_new[idx]    = 0.0;
            pPorosity_new[idx]  = 0.0;

            // Apply various erosion algorithms
            if (d_allowNoTension) {
              if (p > 0.0) {
                sigma = Vaango::Util::Zero;
              } else {
                sigma = tensorHy;
              }
            } else if (d_allowNoShear) {
              sigma = tensorHy;
            } else if (d_setStressToZero) {
              sigma = Vaango::Util::Zero;
            }
          }
        }
      }

      //-----------------------------------------------------------------------
      // Stage 5:
      //-----------------------------------------------------------------------

      // Save the new data
      pStress_new[idx] = sigma;

      // Compute the strain energy for non-localized particles
      if (pLocalized_new[idx] == 0) {
        Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * 0.5;
        double avgVolume  = (pVolume_new[idx] + pVolume_old[idx]) * 0.5;

        double pSpecificStrainEnergy =
          (pDefRate_mid[idx](0, 0) * avgStress(0, 0) + pDefRate_mid[idx](1, 1) * avgStress(1, 1) +
           pDefRate_mid[idx](2, 2) * avgStress(2, 2) +
           2.0 * (pDefRate_mid[idx](0, 1) * avgStress(0, 1) +
                  pDefRate_mid[idx](0, 2) * avgStress(0, 2) +
                  pDefRate_mid[idx](1, 2) * avgStress(1, 2))) *
          avgVolume * delT / pMass[idx];

        pEnergy_new[idx] = pEnergy[idx] + pSpecificStrainEnergy -
                           p_q[idx] * dv_dt * delT * include_AV_heating;

        totalStrainEnergy += pSpecificStrainEnergy * pMass[idx];
      } else {
        pEnergy_new[idx] = pEnergy[idx];
      }

      // Compute wave speed at each particle, store the maximum
      double modulus = state->bulkModulus + 4.0 * state->shearModulus / 3.0;
      double c_dil = std::sqrt(modulus / rho_new);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);

    } // end particle loop

    //__________________________________
    //
    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();

    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(totalStrainEnergy), lb->StrainEnergyLabel);
    }
  }

  if (cout_EP.active())
    cout_EP << getpid() << "... End."
            << "\n";
}

////////////////////////////////////////////////////////////////////////
/*! \brief Compute Stilde, epdot, ep, and delGamma using
  Simo's approach */
////////////////////////////////////////////////////////////////////////
void
ElasticPlasticHP::computePlasticStateViaRadialReturn(const Matrix3& trialS,
                                                     const double& delT,
                                                     const MPMMaterial* matl,
                                                     const particleIndex idx,
                                                     ModelStateBase* state,
                                                     Matrix3& nn,
                                                     double& delGamma)
{
  double normTrialS = trialS.Norm();

  // Do Newton iteration to compute delGamma and updated
  // plastic strain, plastic strain rate, and yield stress
  double tolerance = std::min(delT, 1.0e-6);
  delGamma = computeDeltaGamma(delT, tolerance, normTrialS, matl, idx, state);

  nn = trialS / normTrialS;
}

////////////////////////////////////////////////////////////////////////
// Compute the quantity
//             \f$d(\gamma)/dt * \Delta T = \Delta \gamma \f$
//             using Newton iterative root finder */
////////////////////////////////////////////////////////////////////////
double
ElasticPlasticHP::computeDeltaGamma(const double& delT,
                                    const double& tolerance,
                                    const double& normTrialS,
                                    const MPMMaterial* matl,
                                    const particleIndex idx,
                                    ModelStateBase* state)
{
  // Initialize constants
  double twothird  = 2.0 / 3.0;
  double stwothird = sqrt(twothird);
  double sthreetwo = 1.0 / stwothird;
  double twomu     = 2.0 * state->shearModulus;

  // Initialize variables
  double ep            = state->plasticStrain;
  double sigma_y       = state->yieldStress;
  double deltaGamma    = state->plasticStrainRate * delT * sthreetwo;
  double deltaGammaOld = deltaGamma;
  double g             = 0.0;
  double Dg            = 1.0;

  //__________________________________
  // iterate
  int count = 0;
  do {

    ++count;

    // Compute the yield stress
    sigma_y = d_flow->computeFlowStress(state, delT, tolerance, matl, idx);

    // Compute g
    g = normTrialS - stwothird * sigma_y - twomu * deltaGamma;

    // Compute d(sigma_y)/d(epdot)
    double dsigy_depdot = d_flow->evalDerivativeWRTStrainRate(state, idx);

    // Compute d(sigma_y)/d(ep)
    double dsigy_dep = d_flow->evalDerivativeWRTPlasticStrain(state, idx);

    // Compute d(g)/d(deltaGamma)
    Dg = -twothird * (dsigy_depdot / delT + dsigy_dep) - twomu;

    // Update deltaGamma
    deltaGammaOld = deltaGamma;
    deltaGamma -= g / Dg;

    if (std::isnan(g) || std::isnan(deltaGamma)) {
      std::cout << "idx = " << idx << " iter = " << count << " g = " << g
                << " Dg = " << Dg << " deltaGamma = " << deltaGamma
                << " sigy = " << sigma_y << " dsigy/depdot = " << dsigy_depdot
                << " dsigy/dep= " << dsigy_dep
                << " epdot = " << state->plasticStrainRate
                << " ep = " << state->plasticStrain << "\n";
      throw InternalError("nans in computation", __FILE__, __LINE__);
    }

    // Update local plastic strain rate
    double stt_deltaGamma    = std::max(stwothird * deltaGamma, 0.0);
    state->plasticStrainRate = stt_deltaGamma / delT;

    // Update local plastic strain
    state->plasticStrain = ep + stt_deltaGamma;

    if (std::abs(deltaGamma - deltaGammaOld) < tolerance || count > 100)
      break;

  } while (std::abs(g) > sigma_y / 1000.);

  // Compute the yield stress
  state->yieldStress =
    d_flow->computeFlowStress(state, delT, tolerance, matl, idx);

  if (std::isnan(state->yieldStress)) {
    std::cout << "idx = " << idx << " iter = " << count
              << " sig_y = " << state->yieldStress
              << " epdot = " << state->plasticStrainRate
              << " ep = " << state->plasticStrain
              << " T = " << state->temperature << " Tm = " << state->meltingTemp
              << "\n";
  }

  return deltaGamma;
}

/*! Compute the elastic tangent modulus tensor for isotropic
    materials
    Assume: [stress] = [s11 s22 s33 s12 s23 s31]
            [strain] = [e11 e22 e33 2e12 2e23 2e31]
*/
void
ElasticPlasticHP::computeElasticTangentModulus(const double& K,
                                               const double& mu,
                                               double Ce[6][6])
{
  // Form the elastic tangent modulus tensor
  double twomu        = 2.0 * mu;
  double lambda       = K - (twomu / 3.0);
  double lambda_twomu = lambda + twomu;

  for (int ii = 0; ii < 6; ++ii) {
    for (int jj = 0; jj < 6; ++jj) {
      Ce[ii][jj] = 0.0;
    }
  }
  Ce[0][0] = lambda_twomu;
  Ce[1][1] = lambda_twomu;
  Ce[2][2] = lambda_twomu;
  Ce[3][3] = mu;
  Ce[4][4] = mu;
  Ce[5][5] = mu;
  Ce[0][1] = lambda;
  Ce[0][2] = lambda;
  Ce[1][2] = lambda;
  for (int ii = 1; ii < 3; ++ii) {
    for (int jj = 0; jj < ii; ++jj) {
      Ce[ii][jj] = Ce[jj][ii];
    }
  }
}

// Compute the elastic tangent modulus tensor for isotropic
// materials (**NOTE** can get rid of one copy operation if needed)
void
ElasticPlasticHP::computeElasticTangentModulus(double bulk,
                                               double shear,
                                               TangentModulusTensor& Ce)
{
  // Form the elastic tangent modulus tensor
  double E   = 9.0 * bulk * shear / (3.0 * bulk + shear);
  double nu  = E / (2.0 * shear) - 1.0;
  double fac = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
  double C11 = fac * (1.0 - nu);
  double C12 = fac * nu;
  FastMatrix C_6x6(6, 6);
  for (int ii = 0; ii < 6; ++ii)
    for (int jj = 0; jj < 6; ++jj)
      C_6x6(ii, jj) = 0.0;
  C_6x6(0, 0)       = C11;
  C_6x6(1, 1)       = C11;
  C_6x6(2, 2)       = C11;
  C_6x6(0, 1)       = C12;
  C_6x6(0, 2)       = C12;
  C_6x6(1, 0)       = C12;
  C_6x6(1, 2)       = C12;
  C_6x6(2, 0)       = C12;
  C_6x6(2, 1)       = C12;
  C_6x6(3, 3)       = shear;
  C_6x6(4, 4)       = shear;
  C_6x6(5, 5)       = shear;

  Ce.convertToTensorForm(C_6x6);
}

/*! Compute the elastic-plastic tangent modulus tensor for isotropic
  materials for use in the implicit stress update
  Assume: [stress] = [s11 s22 s33 s12 s23 s31]
  [strain] = [e11 e22 e33 2e12 2e23 2e31]
  Uses alogorithm for small strain plasticity (Simo 1998, p.124)
*/
void
ElasticPlasticHP::computeEPlasticTangentModulus(const double& K,
                                                const double& mu,
                                                const double& delGamma,
                                                const Matrix3& trialStress,
                                                const particleIndex idx,
                                                ModelStateBase* state,
                                                double Cep[6][6],
                                                bool consistent)
{
  double normTrialS = trialStress.Norm();
  Matrix3 n         = trialStress / normTrialS;

  // Compute theta and theta_bar
  double twomu    = 2.0 * mu;
  double dsigYdep = d_flow->evalDerivativeWRTPlasticStrain(state, idx);

  double theta = 1.0;
  if (consistent) {
    theta = 1.0 - (twomu * delGamma) / normTrialS;
  }
  double thetabar = 1.0 / (1.0 + dsigYdep / (3.0 * mu)) - (1.0 - theta);

  // Form the elastic-plastic tangent modulus tensor
  double twomu3           = twomu / 3.0;
  double twomu3theta      = twomu3 * theta;
  double kfourmu3theta    = K + 2.0 * twomu3theta;
  double twomutheta       = twomu * theta;
  double ktwomu3theta     = K - twomu3theta;
  double twomuthetabar    = twomu * thetabar;
  double twomuthetabarn11 = twomuthetabar * n(0, 0);
  double twomuthetabarn22 = twomuthetabar * n(1, 1);
  double twomuthetabarn33 = twomuthetabar * n(2, 2);
  double twomuthetabarn23 = twomuthetabar * n(1, 2);
  double twomuthetabarn31 = twomuthetabar * n(2, 0);
  double twomuthetabarn12 = twomuthetabar * n(0, 1);

  Cep[0][0] = kfourmu3theta - twomuthetabarn11 * n(0, 0);
  Cep[0][1] = ktwomu3theta - twomuthetabarn11 * n(1, 1);
  Cep[0][2] = ktwomu3theta - twomuthetabarn11 * n(2, 2);
  Cep[0][3] = -0.5 * twomuthetabarn11 * n(0, 1);
  Cep[0][4] = -0.5 * twomuthetabarn11 * n(1, 2);
  Cep[0][5] = -0.5 * twomuthetabarn11 * n(2, 0);

  Cep[1][0] = Cep[0][1];
  Cep[1][1] = kfourmu3theta - twomuthetabarn22 * n(1, 1);
  Cep[1][2] = ktwomu3theta - twomuthetabarn22 * n(2, 2);
  Cep[1][3] = -0.5 * twomuthetabarn22 * n(0, 1);
  Cep[1][4] = -0.5 * twomuthetabarn22 * n(1, 2);
  Cep[1][5] = -0.5 * twomuthetabarn22 * n(2, 0);

  Cep[2][0] = Cep[0][2];
  Cep[2][1] = Cep[1][2];
  Cep[2][2] = kfourmu3theta - twomuthetabarn33 * n(2, 2);
  Cep[2][3] = -0.5 * twomuthetabarn33 * n(0, 1);
  Cep[2][4] = -0.5 * twomuthetabarn33 * n(1, 2);
  Cep[2][5] = -0.5 * twomuthetabarn33 * n(2, 0);

  Cep[3][0] = Cep[0][3];
  Cep[3][1] = Cep[1][3];
  Cep[3][2] = Cep[2][3];
  Cep[3][3] = 0.5 * (twomutheta - twomuthetabarn12 * n(0, 1));
  Cep[3][4] = -0.5 * twomuthetabarn12 * n(1, 2);
  Cep[3][5] = -0.5 * twomuthetabarn12 * n(2, 0);

  Cep[4][0] = Cep[0][4];
  Cep[4][1] = Cep[1][4];
  Cep[4][2] = Cep[2][4];
  Cep[4][3] = Cep[3][4];
  Cep[4][4] = 0.5 * (twomutheta - twomuthetabarn23 * n(1, 2));
  Cep[4][5] = -0.5 * twomuthetabarn23 * n(2, 0);

  Cep[5][0] = Cep[0][5];
  Cep[5][1] = Cep[1][5];
  Cep[5][2] = Cep[2][5];
  Cep[5][3] = Cep[3][5];
  Cep[5][4] = Cep[4][5];
  Cep[5][5] = 0.5 * (twomutheta - twomuthetabarn31 * n(2, 0));
}

/* For copying damage parameter to SerialMPM */
void
ElasticPlasticHP::addRequiresDamageParameter(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
}

/* For copying damage parameter to SerialMPM */
void
ElasticPlasticHP::getDamageParameter(const Patch* patch,
                                     ParticleVariable<int>& damage,
                                     int matID,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, pLocalizedLabel_preReloc, pset);

  // Only update the damage variable if it hasn't been modified by a damage
  // model earlier
  for (auto particle : *pset) {
    if (damage[particle] == 0) {
      damage[particle] = pLocalized[particle];
    }
  }
}

// Update the porosity of the material
double
ElasticPlasticHP::updatePorosity(const Matrix3& D,
                                 double delT,
                                 double f,
                                 double ep)
{
  // Growth
  // Calculate trace of D
  double Dkk = D.Trace();
  Matrix3 one;
  one.Identity();
  Matrix3 eta = D - Vaango::Util::Identity * (Dkk / 3.0);

  // Calculate rate of growth
  double fdot_grow = 0.0;
  if (Dkk > 0.0)
    fdot_grow = (1.0 - f) * Dkk;

  // Nucleation
  // Calculate A
  double A = voidNucleationFactor(ep);

  // Calculate plastic strain rate
  double epdot = sqrt(eta.NormSquared() / 1.5);

  // Calculate rate of nucleation
  double fdot_nucl = A * epdot;

  // Update void volume fraction using forward euler
  double f_new = f + delT * (fdot_nucl + fdot_grow);
  // std::cout << "Porosity: D = " << D << "\n";
  // std::cout << "Porosity: eta = " << eta << "\n";
  // std::cout << "Porosity: Dkk = " << Dkk << "\n";
  // std::cout << "Porosity::fdot_gr = " << fdot_grow
  //     << " fdot_nucl = " << fdot_nucl << " f = " << f
  //     << " f_new = " << f_new << "\n";
  return f_new;
}
//______________________________________________________________________
//
// Calculate the void nucleation factor
inline double
ElasticPlasticHP::voidNucleationFactor(double ep)
{
  double temp = (ep - d_porosity.en) / d_porosity.sn;
  double A    = d_porosity.fn / (d_porosity.sn * sqrt(2.0 * M_PI)) *
             exp(-0.5 * temp * temp);
  return A;
}

/* Hardcoded specific heat computation for 4340 steel */
/*
double
ElasticPlasticHP::computeSpecificHeat(double T)
{
  // Specific heat model for 4340 steel (SI units)
  double Tc = 1040.0;
  if (T == Tc) {
    T = T - 1.0;
  }
  double Cp = 500.0;
  if (T < Tc) {
    double t = 1 - T/Tc;
    d_Cp.A = 190.14;
    d_Cp.B = 273.75;
    d_Cp.C = 418.30;
    d_Cp.n = 0.2;
    Cp = d_Cp.A - d_Cp.B*t + d_Cp.C/pow(t, d_Cp.n);
  } else {
    double t = T/Tc - 1.0;
    d_Cp.A = 465.21;
    d_Cp.B = 267.52;
    d_Cp.C = 58.16;
    d_Cp.n = 0.35;
    Cp = d_Cp.A + d_Cp.B*t + d_Cp.C/pow(t, d_Cp.n);
  }
  return Cp;

  // Specific heat model for HY-100 steel
  //return 1.0e3*(d_Cp.A + d_Cp.B*T + d_Cp.C/(T*T));
}
*/

/* For material conversion */
void
ElasticPlasticHP::allocateCMDataAddRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet* patch,
                                            MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patch);

  // Add requires local to this model
  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::NewDW, pRotationLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pStrainRateLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pPlasticStrainLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pPlasticStrainRateLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pDamageLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pPorosityLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pEnergyLabel_preReloc, matlset, gnone);
  d_flow->allocateCMDataAddRequires(task, matl, patch, lb);
}

/* For material conversion */
void
ElasticPlasticHP::allocateCMDataAdd(DataWarehouse* new_dw,
                                    ParticleSubset* addset,
                                    ParticleLabelVariableMap* newState,
                                    ParticleSubset* delset,
                                    DataWarehouse* old_dw)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  ParticleSubset::iterator n, o;

  ParticleVariable<Matrix3> pRotation;
  ParticleVariable<double> pPlasticStrain, pDamage, pPorosity, pStrainRate,
    pPlasticStrainRate, pEnergy;
  ParticleVariable<int> pLocalized;

  constParticleVariable<Matrix3> o_Rotation;
  constParticleVariable<double> o_PlasticStrain, o_Damage, o_Porosity;
  constParticleVariable<double> o_StrainRate, o_PlasticStrainRate, o_Energy;
  constParticleVariable<int> o_Localized;

  new_dw->allocateTemporary(pRotation, addset);
  new_dw->allocateTemporary(pPlasticStrain, addset);
  new_dw->allocateTemporary(pPlasticStrainRate, addset);
  new_dw->allocateTemporary(pDamage, addset);
  new_dw->allocateTemporary(pStrainRate, addset);
  new_dw->allocateTemporary(pLocalized, addset);
  new_dw->allocateTemporary(pPorosity, addset);
  new_dw->allocateTemporary(pEnergy, addset);

  new_dw->get(o_Rotation, pRotationLabel_preReloc, delset);
  new_dw->get(o_StrainRate, pStrainRateLabel_preReloc, delset);
  new_dw->get(o_PlasticStrain, pPlasticStrainLabel_preReloc, delset);
  new_dw->get(o_PlasticStrainRate, pPlasticStrainRateLabel_preReloc, delset);
  new_dw->get(o_Damage, pDamageLabel_preReloc, delset);
  new_dw->get(o_Localized, pLocalizedLabel_preReloc, delset);
  new_dw->get(o_Porosity, pPorosityLabel_preReloc, delset);
  new_dw->get(o_Energy, pEnergyLabel_preReloc, delset);

  n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pRotation[*n]          = o_Rotation[*o];
    pStrainRate[*n]        = o_StrainRate[*o];
    pPlasticStrain[*n]     = o_PlasticStrain[*o];
    pPlasticStrainRate[*n] = o_PlasticStrainRate[*o];
    pDamage[*n]            = o_Damage[*o];
    pLocalized[*n]         = o_Localized[*o];
    pPorosity[*n]          = o_Porosity[*o];
    pEnergy[*n]            = o_Energy[*o];
  }

  (*newState)[pRotationLabel]          = pRotation.clone();
  (*newState)[pStrainRateLabel]        = pStrainRate.clone();
  (*newState)[pPlasticStrainLabel]     = pPlasticStrain.clone();
  (*newState)[pPlasticStrainRateLabel] = pPlasticStrainRate.clone();
  (*newState)[pDamageLabel]            = pDamage.clone();
  (*newState)[pLocalizedLabel]         = pLocalized.clone();
  (*newState)[pPorosityLabel]          = pPorosity.clone();
  (*newState)[pEnergyLabel]            = pEnergy.clone();

  // Initialize the data for the flow model
  d_flow->allocateCMDataAdd(new_dw, addset, newState, delset, old_dw);
}

/* For material conversion */
void
ElasticPlasticHP::scheduleCheckNeedAddMPMMaterial(Task* task,
                                                  const MPMMaterial* matl,
                                                  const PatchSet*) const
{
  Ghost::GhostType gnone        = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pPlasticStrainLabel_preReloc, matlset, gnone);

  task->computes(lb->NeedAddMPMMaterialLabel);
}

/* For material conversion */
void
ElasticPlasticHP::checkNeedAddMPMMaterial(const PatchSubset* patches,
                                          const MPMMaterial* matl,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  if (cout_EP.active()) {
    cout_EP << getpid() << "checkNeedAddMPMMaterial: In : Matl = " << matl
            << " id = " << matl->getDWIndex()
            << " patch = " << (patches->get(0))->getID();
  }

  double need_add = 0.;

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int matID              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
    constParticleVariable<double> pPlasticStrain;
    new_dw->get(pPlasticStrain, pPlasticStrainLabel_preReloc, pset);

    // Loop thru particles
    ParticleSubset::iterator iter = pset->begin();
    for (; iter != pset->end(); iter++) {
      particleIndex idx = *iter;
      if (pPlasticStrain[idx] > 5.e-2) {
        need_add = -1.;
      }
    }
  }
  new_dw->put(sum_vartype(need_add), lb->NeedAddMPMMaterialLabel);
}

/* for RigidMPM */
void
ElasticPlasticHP::carryForward(const PatchSubset* patches,
                               const MPMMaterial* matl,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch   = patches->get(p);
    int matID              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    constParticleVariable<Matrix3> pRotation;
    constParticleVariable<double> pPlasticStrain, pDamage, pPorosity;
    constParticleVariable<double> pStrainRate, pPlasticStrainRate;
    constParticleVariable<int> pLocalized;

    old_dw->get(pRotation, pRotationLabel, pset);
    old_dw->get(pStrainRate, pStrainRateLabel, pset);
    old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);
    old_dw->get(pPlasticStrainRate, pPlasticStrainRateLabel, pset);
    old_dw->get(pDamage, pDamageLabel, pset);
    old_dw->get(pPorosity, pPorosityLabel, pset);
    old_dw->get(pLocalized, pLocalizedLabel, pset);

    ParticleVariable<Matrix3> pRotation_new;
    ParticleVariable<double> pPlasticStrain_new, pDamage_new;
    ParticleVariable<double> pPorosity_new, pStrainRate_new,
      pPlasticStrainRate_new;
    ParticleVariable<int> pLocalized_new;

    new_dw->allocateAndPut(pRotation_new, pRotationLabel_preReloc, pset);
    new_dw->allocateAndPut(pStrainRate_new, pStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrainRate_new, pPlasticStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(pDamage_new, pDamageLabel_preReloc, pset);
    new_dw->allocateAndPut(pPorosity_new, pPorosityLabel_preReloc, pset);
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);

    // Get the plastic strain
    d_flow->getInternalVars(pset, old_dw);
    d_flow->allocateAndPutRigid(pset, new_dw);

    d_flow->getInternalVars(pset, old_dw);
    d_flow->allocateAndPutRigid(pset, new_dw);

    for (int idx : *pset) {
      pRotation_new[idx]          = pRotation[idx];
      pStrainRate_new[idx]        = pStrainRate[idx];
      pPlasticStrain_new[idx]     = pPlasticStrain[idx];
      pPlasticStrainRate_new[idx] = pPlasticStrainRate[idx];
      pDamage_new[idx]            = pDamage[idx];
      pPorosity_new[idx]          = pPorosity[idx];
      pLocalized_new[idx]         = pLocalized[idx];
    }

    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

/* for MPMICE */
double
ElasticPlasticHP::computeRhoMicroCM(double pressure,
                                    const double p_ref,
                                    const MPMMaterial* matl,
                                    double temperature,
                                    double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double bulk     = d_initialData.Bulk;

  double p_gauge = pressure - p_ref;
  double rho_new;

  if (d_useModifiedEOS && p_gauge < 0.0) {
    double A = p_ref; // modified EOS
    double n = p_ref / bulk;
    rho_new  = rho_orig * pow(pressure / A, n);
  } else { // Standard EOS
    double p_g_over_bulk = p_gauge / bulk;
    rho_new =
      rho_orig * (p_g_over_bulk + sqrt(p_g_over_bulk * p_g_over_bulk + 1.));
  }
  return rho_new;
}

/* for MPMICE */
void
ElasticPlasticHP::computePressEOSCM(double rho_new,
                                    double& pressure,
                                    double p_ref,
                                    double& dp_drho,
                                    double& tmp,
                                    const MPMMaterial* matl,
                                    double temperature)
{
  double bulk         = d_initialData.Bulk;
  double rho_orig     = matl->getInitialDensity();
  double inv_rho_orig = 1. / rho_orig;

  if (d_useModifiedEOS && rho_new < rho_orig) {
    double A                = p_ref; // MODIFIED EOS
    double n                = bulk / p_ref;
    double rho_rat_to_the_n = pow(rho_new / rho_orig, n);
    pressure                = A * rho_rat_to_the_n;
    dp_drho                 = (bulk / rho_new) * rho_rat_to_the_n;
    tmp                     = dp_drho; // speed of sound squared
  } else {                             // STANDARD EOS
    double p_g = .5 * bulk * (rho_new * inv_rho_orig - rho_orig / rho_new);
    pressure   = p_ref + p_g;
    dp_drho    = .5 * bulk * (rho_orig / (rho_new * rho_new) + inv_rho_orig);
    tmp        = bulk / rho_new; // speed of sound squared
  }
}

/* for MPMICE */
double
ElasticPlasticHP::getCompressibility()
{
  return 1.0 / d_initialData.Bulk;
}
