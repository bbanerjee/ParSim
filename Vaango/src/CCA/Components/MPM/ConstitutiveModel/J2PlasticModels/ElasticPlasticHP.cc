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

#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/ElasticPlasticHP.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/DamageModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/DevStressModels/DevStressModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/FlowStressModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/MeltTempModels/MeltingTempModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/DeformationState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/IsoNonlinHypoelastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecHeatModels/SpecificHeatModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/StabilityModels/StabilityCheckFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldConditionFactory.h>
#include <Core/Math/Short27.h> //for Fracture

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
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
#include <Core/Exceptions/InvalidValue.h>
#include <Core/ProblemSpec/ProblemSpec.h>

//#define SUB_CYCLE_F
#undef SUB_CYCLE_F
//#define DEBUG_RETURN

using namespace Uintah;
using namespace Vaango;

static DebugStream cout_EP("EP", false);
static DebugStream cout_EP1("EP1", false);
static DebugStream cout_EP_return("EP_return", false);
static DebugStream CSTi("EPi", false);
static DebugStream CSTir("EPir", false);

ElasticPlasticHP::ElasticPlasticHP(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
  , ImplicitCM()
{
  if (cout_EP.active()) {
    cout_EP << getpid() << "Constructing ElasticPlasticHP:\n";
  }

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

  d_elastic = std::make_shared<ElasticModuli_MetalIso>(d_eos, d_shear);

  d_melt = MeltingTempModelFactory::create(ps);
  if (!d_melt) {
     std::ostringstream desc;
    desc << "ElasticPlasticHP::Error in melting temp model factory"
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_devStress = DevStressModelFactory::create(ps);
  if (!d_devStress) {
     std::ostringstream desc;
    desc << "ElasticPlasticHP::Error creating deviatoric stress model"
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  ProblemSpecP intvar_ps = ps->findBlock("internal_variable_model");
  if (!intvar_ps) {
     std::ostringstream err;
    err << "**ERROR** Please add an 'internal_variable_model' tag to the\n"
        << " 'elastic_plastic_hp' block in the input .ups file.  The\n"
        << " default type is 'metal_internal_var'.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  d_intvar = std::make_shared<Vaango::IntVar_Metal>(intvar_ps);
  if (!d_intvar) {
     std::ostringstream err;
    err << "**ERROR** An error occured while creating the internal variable \n"
         << " model. Please file a bug report.\n";
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  d_yield = Vaango::YieldConditionFactory::create(ps, d_intvar.get(), 
              const_cast<const FlowStressModel*>(d_flow));
  if (!d_yield) {
     std::ostringstream desc;
    desc << "An error occured in the YieldConditionFactory that has \n"
         << " slipped through the existing bullet proofing. Please tell \n"
         << " Biswajit.\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  d_computeSpecificHeat = false;
  ps->get("compute_specific_heat", d_computeSpecificHeat);
  d_Cp = SpecificHeatModelFactory::create(ps);

  setErosionAlgorithm();
  getInitialPorosityData(ps);
  getInitialDamageData(ps);

  if (cout_EP.active()) {
    cout_EP << getpid() << "Got all input data from ps.\n";
  }

  initializeLocalMPMLabels();

  if (cout_EP.active()) {
    cout_EP << getpid() << "Construction done\n";
  }
}

ElasticPlasticHP::ElasticPlasticHP(const ElasticPlasticHP* cm)
  : ConstitutiveModel(cm)
  , ImplicitCM(cm)
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
  d_yield  = Vaango::YieldConditionFactory::createCopy(cm->d_yield);
  d_stable = StabilityCheckFactory::createCopy(cm->d_stable);
  d_flow   = FlowStressModelFactory::createCopy(cm->d_flow);
  d_damage = DamageModelFactory::createCopy(cm->d_damage);
  d_eos    = MPMEquationOfStateFactory::createCopy(cm->d_eos);
  d_shear  = Vaango::ShearModulusModelFactory::createCopy(cm->d_shear);
  d_melt   = MeltingTempModelFactory::createCopy(cm->d_melt);
  d_eos->setBulkModulus(d_initialData.Bulk);

  d_intvar    = std::make_shared<Vaango::IntVar_Metal>(cm->d_intvar.get());
  d_elastic   = std::make_shared<ElasticModuli_MetalIso>(d_eos, d_shear);
  d_devStress = nullptr;

  initializeLocalMPMLabels();
}

ElasticPlasticHP::~ElasticPlasticHP()
{
  // Destructor
  VarLabel::destroy(pRotationLabel);
  VarLabel::destroy(pEqStrainRateLabel);
  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pEqPlasticStrainLabel);
  VarLabel::destroy(pEqPlasticStrainRateLabel);
  VarLabel::destroy(pDamageLabel);
  VarLabel::destroy(pPorosityLabel);
  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pEnergyLabel);
  VarLabel::destroy(pIntVarLabel);

  VarLabel::destroy(pRotationLabel_preReloc);
  VarLabel::destroy(pEqStrainRateLabel_preReloc);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);
  VarLabel::destroy(pEqPlasticStrainLabel_preReloc);
  VarLabel::destroy(pEqPlasticStrainRateLabel_preReloc);
  VarLabel::destroy(pDamageLabel_preReloc);
  VarLabel::destroy(pPorosityLabel_preReloc);
  VarLabel::destroy(pLocalizedLabel_preReloc);
  VarLabel::destroy(pEnergyLabel_preReloc);
  VarLabel::destroy(pIntVarLabel_preReloc);

  delete d_flow;
  delete d_yield;
  delete d_stable;
  delete d_damage;
  delete d_eos;
  delete d_shear;
  delete d_melt;
  delete d_Cp;
  delete d_devStress;
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
  d_devStress->outputProblemSpec(cm_ps);
  d_damage->outputProblemSpec(cm_ps);
  d_eos->outputProblemSpec(cm_ps);
  d_shear->outputProblemSpec(cm_ps);
  d_melt->outputProblemSpec(cm_ps);
  d_Cp->outputProblemSpec(cm_ps);
  d_intvar->outputProblemSpec(cm_ps);

  cm_ps->appendElement("evolve_porosity", d_evolvePorosity);
  cm_ps->appendElement("initial_mean_porosity", d_porosity.f0);
  cm_ps->appendElement("initial_std_porosity", d_porosity.f0_std);
  cm_ps->appendElement("critical_porosity", d_porosity.fc);
  cm_ps->appendElement("frac_nucleation", d_porosity.fn);
  cm_ps->appendElement("meanstrain_nucleation", d_porosity.en);
  cm_ps->appendElement("stddevstrain_nucleation", d_porosity.sn);
  cm_ps->appendElement("initial_porosity_distrib", d_porosity.porosityDist);

  cm_ps->appendElement("evolve_damage", d_evolveDamage);
  cm_ps->appendElement("initial_mean_scalar_damage", d_scalarDam.D0);
  cm_ps->appendElement("initial_std_scalar_damage", d_scalarDam.D0_std);
  cm_ps->appendElement("critical_scalar_damage", d_scalarDam.Dc);
  cm_ps->appendElement("initial_scalar_damage_distrib",
                       d_scalarDam.scalarDamageDist);
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
  pEqStrainRateLabel = VarLabel::create(
    "p.eqStrainRate", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel = VarLabel::create(
    "p.plasticStrain", ParticleVariable<Matrix3>::getTypeDescription());
  pEqPlasticStrainLabel = VarLabel::create(
    "p.eqPlasticStrain", ParticleVariable<double>::getTypeDescription());
  pEqPlasticStrainRateLabel = VarLabel::create(
    "p.eqPlasticStrainRate", ParticleVariable<double>::getTypeDescription());
  pDamageLabel = VarLabel::create(
    "p.damage", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel = VarLabel::create(
    "p.porosity", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pEnergyLabel = VarLabel::create(
    "p.energy", ParticleVariable<double>::getTypeDescription());
  pIntVarLabel = VarLabel::create(
    "p.intvar_metal", ParticleVariable<MetalIntVar>::getTypeDescription());

  pRotationLabel_preReloc = VarLabel::create(
    "p.rotation+", ParticleVariable<Matrix3>::getTypeDescription());
  pEqStrainRateLabel_preReloc = VarLabel::create(
    "p.eqStrainRate+", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel_preReloc = VarLabel::create(
    "p.plasticStrain+", ParticleVariable<Matrix3>::getTypeDescription());
  pEqPlasticStrainLabel_preReloc = VarLabel::create(
    "p.eqPlasticStrain+", ParticleVariable<double>::getTypeDescription());
  pEqPlasticStrainRateLabel_preReloc = VarLabel::create(
    "p.eqPlasticStrainRate+", ParticleVariable<double>::getTypeDescription());
  pDamageLabel_preReloc = VarLabel::create(
    "p.damage+", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = VarLabel::create(
    "p.porosity+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
  pEnergyLabel_preReloc = VarLabel::create(
    "p.energy+", ParticleVariable<double>::getTypeDescription());
  pIntVarLabel_preReloc = VarLabel::create(
    "p.intvar_metal+", ParticleVariable<MetalIntVar>::getTypeDescription());
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
}

/*! Compute specific heat

    double T = temperature;
    C_p = 1.0e3*(A + B*T + C/T^2)
    ** For steel **
    C_p = 1.0e3*(0.09278 + 7.454e-4*T + 12404.0/(T*T));
*/
/*
void
ElasticPlasticHP::getSpecificHeatData(ProblemSpecP& ps)
{
  d_Cp.A = 0.09278;  // Constant A (HY100)
  d_Cp.B = 7.454e-4; // Constant B (HY100)
  d_Cp.C = 12404.0;  // Constant C (HY100)
  d_Cp.n = 2.0;      // Constant n (HY100)
  ps->get("Cp_constA", d_Cp.A);
  ps->get("Cp_constB", d_Cp.B);
  ps->get("Cp_constC", d_Cp.C);
  ps->get("Cp_constn", d_Cp.n);
}
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
  // Add the local particle state data for this constitutive model.
  from.push_back(pRotationLabel);
  from.push_back(pEqStrainRateLabel);
  from.push_back(pPlasticStrainLabel);
  from.push_back(pEqPlasticStrainLabel);
  from.push_back(pEqPlasticStrainRateLabel);
  from.push_back(pDamageLabel);
  from.push_back(pPorosityLabel);
  from.push_back(pLocalizedLabel);
  from.push_back(pEnergyLabel);
  from.push_back(pIntVarLabel);

  to.push_back(pRotationLabel_preReloc);
  to.push_back(pEqStrainRateLabel_preReloc);
  to.push_back(pPlasticStrainLabel_preReloc);
  to.push_back(pEqPlasticStrainLabel_preReloc);
  to.push_back(pEqPlasticStrainRateLabel_preReloc);
  to.push_back(pDamageLabel_preReloc);
  to.push_back(pPorosityLabel_preReloc);
  to.push_back(pLocalizedLabel_preReloc);
  to.push_back(pEnergyLabel_preReloc);
  to.push_back(pIntVarLabel_preReloc);

  // Add the particle state for the flow & deviatoric stress model
  d_flow->addParticleState(from, to);
  d_devStress->addParticleState(from, to);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::addInitialComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  task->computes(pRotationLabel, matlset);
  task->computes(pEqStrainRateLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(pEqPlasticStrainLabel, matlset);
  task->computes(pEqPlasticStrainRateLabel, matlset);
  task->computes(pDamageLabel, matlset);
  task->computes(pPorosityLabel, matlset);
  task->computes(pLocalizedLabel, matlset);
  task->computes(pEnergyLabel, matlset);
  task->computes(pIntVarLabel, matlset);

  // Add internal evolution variables computed by flow & deviatoric stress model
  d_flow->addInitialComputesAndRequires(task, matl, patch);
  d_devStress->addInitialComputesAndRequires(task, matl);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::initializeCMData(const Patch* patch,
                                   const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  if (flag->d_integrator == MPMFlags::Implicit)
    initSharedDataForImplicit(patch, matl, new_dw);
  else {
    initSharedDataForExplicit(patch, matl, new_dw);
    computeStableTimestep(patch, matl, new_dw);
  }

  // Put stuff in here to initialize each particle's
  // constitutive model parameters and deformationMeasure
  // std::cout << "Initialize CM Data in ElasticPlasticHP" << "\n";
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<Matrix3> pRotation, pPlasticStrain; 
  ParticleVariable<double> pEqPlasticStrain, pDamage, pPorosity,
    pEqPlasticStrainRate, pEqStrainRate, pEnergy;
  ParticleVariable<int> pLocalized;
  ParticleVariable<MetalIntVar> pIntVar;

  new_dw->allocateAndPut(pRotation, pRotationLabel, pset);
  new_dw->allocateAndPut(pEqStrainRate, pEqStrainRateLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pEqPlasticStrainRate, pEqPlasticStrainRateLabel, pset);
  new_dw->allocateAndPut(pDamage, pDamageLabel, pset);
  new_dw->allocateAndPut(pLocalized, pLocalizedLabel, pset);
  new_dw->allocateAndPut(pPorosity, pPorosityLabel, pset);
  new_dw->allocateAndPut(pEnergy, pEnergyLabel, pset);
  new_dw->allocateAndPut(pIntVar, pIntVarLabel, pset);

  for (auto idx : *pset) {

    pRotation[idx]            = Vaango::Util::Identity;
    pEqStrainRate[idx]        = 0.0;
    pPlasticStrain[idx]       = Vaango::Util::Zero;
    pEqPlasticStrain[idx]     = 0.0;
    pEqPlasticStrainRate[idx] = 0.0;
    pDamage[idx]              = d_damage->initialize();
    pPorosity[idx]            = d_porosity.f0;
    pLocalized[idx]           = 0;
    pEnergy[idx]              = 0.;
    pIntVar[idx]              = {0.0, 0.0};
  }

  // Do some extra things if the porosity or the damage distribution
  // is not uniform.
  // ** WARNING ** Weibull distribution needs to be implemented.
  //               At present only Gaussian available.
  if (d_porosity.porosityDist != "constant") {

    // Generate a Gaussian distributed random number given the mean
    // porosity and the std.
    int patchID = patch->getID();
    int patch_div_32 = patchID / 32;
    patchID = patchID % 32;
    unsigned int unique_seed = ((12345 + patch_div_32 + 1) << patchID);
    Uintah::Gaussian gaussGen(d_porosity.f0, d_porosity.f0_std, unique_seed, 1, DBL_MAX);
    for (auto idx : *pset) {
      pPorosity[idx] = std::abs(gaussGen.rand(1.0));
    }
  }

  if (d_scalarDam.scalarDamageDist != "constant") {

    // Generate a Gaussian distributed random number given the mean
    // damage and the std.
    int patchID = patch->getID();
    int patch_div_32 = patchID / 32;
    patchID = patchID % 32;
    unsigned int unique_seed = ((612345 + patch_div_32 + 1) << patchID);
    Uintah::Gaussian gaussGen(
      d_scalarDam.D0, d_scalarDam.D0_std, unique_seed, 1, DBL_MAX);
    for (auto idx : *pset) {
      pDamage[idx] = std::abs(gaussGen.rand(1.0));
    }
  }

  // Initialize the data for the flow model
  d_flow->initializeInternalVars(pset, new_dw);

  // Deviatoric Stress Model
  d_devStress->initializeInternalVars(pset, new_dw);
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
  Vector dx     = patch->dCell();
  int matlindex = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matlindex, patch);

  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double shear = d_initialData.Shear;
  double bulk  = d_initialData.Bulk;

  for (auto idx : *pset) {

    // Compute wave speed at each particle, store the maximum
    Vector pVelocity_idx = pVelocity[idx];
    if (pMass[idx] > 0) {
      c_dil = sqrt((bulk + 4.0 * shear / 3.0) * pVolume[idx] / pMass[idx]);
    } else {
      c_dil         = 0.0;
      pVelocity_idx = Vector(0.0, 0.0, 0.0);
    }
    WaveSpeed = Vector(Max(c_dil + std::abs(pVelocity_idx.x()), WaveSpeed.x()),
                       Max(c_dil + std::abs(pVelocity_idx.y()), WaveSpeed.y()),
                       Max(c_dil + std::abs(pVelocity_idx.z()), WaveSpeed.z()));
  }

  WaveSpeed       = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
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
  if (flag->d_integrator == MPMFlags::Implicit) {
    addSharedCRForImplicitHypo(task, matlset, true);
  } else {
    addSharedCRForHypoExplicit(task, matlset, patches);
  }

  // Other constitutive model and input dependent computes and requires
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);

  task->requires(Task::OldDW, pRotationLabel, matlset, gnone);
  task->requires(Task::OldDW, pEqStrainRateLabel, matlset, gnone);
  task->requires(Task::OldDW, pPlasticStrainLabel, matlset, gnone);
  task->requires(Task::OldDW, pEqPlasticStrainLabel, matlset, gnone);
  task->requires(Task::OldDW, pEqPlasticStrainRateLabel, matlset, gnone);
  task->requires(Task::OldDW, pDamageLabel, matlset, gnone);
  task->requires(Task::OldDW, pPorosityLabel, matlset, gnone);
  task->requires(Task::OldDW, pLocalizedLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);
  task->requires(Task::OldDW, pEnergyLabel, matlset, gnone);
  task->requires(Task::OldDW, pIntVarLabel, matlset, gnone);

  task->computes(pRotationLabel_preReloc, matlset);
  task->computes(pEqStrainRateLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(pEqPlasticStrainLabel_preReloc, matlset);
  task->computes(pEqPlasticStrainRateLabel_preReloc, matlset);
  task->computes(pDamageLabel_preReloc, matlset);
  task->computes(pPorosityLabel_preReloc, matlset);
  task->computes(pLocalizedLabel_preReloc, matlset);
  task->computes(pEnergyLabel_preReloc, matlset);
  task->computes(pIntVarLabel_preReloc, matlset);

  // Add internal evolution variables computed by flow model
  d_flow->addComputesAndRequires(task, matl, patches);

  // Deviatoric stress model
  d_devStress->addComputesAndRequires(task, matl);
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

  // General stuff
  Matrix3 tensorD(0.0); // Rate of deformation
  Matrix3 tensorW(0.0); // Spin
  Matrix3 tensorF;
  tensorF.Identity(); // Deformation gradient
  Matrix3 tensorU;
  tensorU.Identity(); // Right Cauchy-Green stretch
  Matrix3 tensorR;
  tensorR.Identity();     // Rotation
  Matrix3 sigma(0.0);     // The Cauchy stress
  Matrix3 tensorEta(0.0); // Deviatoric part of tensor D
  Matrix3 tensorF_new;
  tensorF_new.Identity(); // Deformation gradient

  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double bulk         = d_initialData.Bulk;
  double shear        = d_initialData.Shear;
  double rho_0        = matl->getInitialDensity();
  double Tm           = matl->getMeltTemperature();

  double totalStrainEnergy = 0.0;

  // Loop thru patches
  for (int patchIndex = 0; patchIndex < patches->size(); patchIndex++) {
    const Patch* patch = patches->get(patchIndex);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    // std::cerr << getpid() << " patch = " << patch->getID() << "\n";
    // Get grid size
    Vector dx = patch->dCell();
    // double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
    double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;

    // Get the set of particles
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    //__________________________________
    // GET GLOBAL DATA

    // Get the deformation gradient (F)
    // Note : The deformation gradient from the old datawarehouse is no
    // longer used, but it is updated for possible use elsewhere
    constParticleVariable<Matrix3> pDefGrad;
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    // Get the particle location, particle size, particle mass, particle volume
    constParticleVariable<Point> px;
    constParticleVariable<Matrix3> pSize;
    constParticleVariable<double> pMass;
    constParticleVariable<double> pVolume;
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);

    // Get the velocity from the grid and particle velocity
    Ghost::GhostType gac = Ghost::AroundCells;
    constParticleVariable<Vector> pVelocity;
    constNCVariable<Vector> gVelocity;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(gVelocity, lb->gVelocityStarLabel, dwi, patch, gac, NGN);

    // Get the particle stress and temperature and energy
    constParticleVariable<Matrix3> pStress;
    constParticleVariable<double> pTemperature;
    old_dw->get(pStress, lb->pStressLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);

    // Get the time increment (delT)
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    constParticleVariable<Short27> pgCode;
    constNCVariable<Vector> GVelocity;
    if (flag->d_fracture) {
      new_dw->get(pgCode, lb->pgCodeLabel, pset);
      new_dw->get(GVelocity, lb->GVelocityStarLabel, dwi, patch, gac, NGN);
    }
    double include_AV_heating = 0.0;
    if (flag->d_artificialViscosityHeating) {
      include_AV_heating = 1.0;
    }

    //__________________________________
    // GET LOCAL DATA

    // Get the rotation (R) and plastic strain/rate tensors
    constParticleVariable<Matrix3> pRotation, pPlasticStrain;
    old_dw->get(pRotation, pRotationLabel, pset);
    old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);

    // Get the particle damage state
    constParticleVariable<double> pEqPlasticStrain, pDamage, pPorosity;
    constParticleVariable<double> pEqStrainRate, pEqPlasticStrainRate, pEnergy;

    old_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
    old_dw->get(pEqPlasticStrainRate, pEqPlasticStrainRateLabel, pset);
    old_dw->get(pDamage, pDamageLabel, pset);
    old_dw->get(pEqStrainRate, pEqStrainRateLabel, pset);
    old_dw->get(pPorosity, pPorosityLabel, pset);
    old_dw->get(pEnergy, pEnergyLabel, pset);

    // Get the particle localization state
    constParticleVariable<int> pLocalized;
    old_dw->get(pLocalized, pLocalizedLabel, pset);

    // Get the particle IDs, useful in case a simulation goes belly up
    constParticleVariable<long64> pParticleID;
    old_dw->get(pParticleID, lb->pParticleIDLabel, pset);

    constParticleVariable<MetalIntVar> pIntVar_old;
    old_dw->get(pIntVar_old, pIntVarLabel, pset);

    // Create and allocate arrays for storing the updated information
    // GLOBAL
    ParticleVariable<Matrix3> tensorL;
    new_dw->getModifiable(tensorL, lb->pVelGradLabel_preReloc, pset);

    ParticleVariable<Matrix3> pDefGrad_new;
    ParticleVariable<double> pVolume_deformed;
    new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->getModifiable(pVolume_deformed, lb->pVolumeLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    // LOCAL
    ParticleVariable<Matrix3> pRotation_new, pPlasticStrain_new;
    ParticleVariable<double> pEqPlasticStrain_new, pDamage_new, pPorosity_new;
    ParticleVariable<double> pEqStrainRate_new, pEqPlasticStrainRate_new;
    ParticleVariable<int> pLocalized_new;
    ParticleVariable<double> pdTdt, p_q, pEnergy_new;

    new_dw->allocateAndPut(pRotation_new, pRotationLabel_preReloc, pset);
    new_dw->allocateAndPut(pEqStrainRate_new, pEqStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrain_new, pEqPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrainRate_new, pEqPlasticStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(pDamage_new, pDamageLabel_preReloc, pset);
    new_dw->allocateAndPut(pPorosity_new, pPorosityLabel_preReloc, pset);
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pEnergy_new, pEnergyLabel_preReloc, pset);

    ParticleVariable<MetalIntVar> pIntVar_new;
    new_dw->allocateAndPut(pIntVar_new, pIntVarLabel_preReloc, pset);

    d_flow->getInternalVars(pset, old_dw);
    d_devStress->getInternalVars(pset, old_dw);

    d_flow->allocateAndPutInternalVars(pset, new_dw);
    d_devStress->allocateAndPutInternalVars(pset, new_dw);

    //______________________________________________________________________
    // Loop thru particles
    for (auto idx : *pset) {

      // Assign zero int. heating by default, modify with appropriate sources
      // This has units (in MKS) of K/s  (i.e. temperature/time)
      pdTdt[idx] = 0.0;

      //-----------------------------------------------------------------------
      // Stage 1:
      //-----------------------------------------------------------------------
      // Calculate the velocity gradient (L) from the grid velocity
      // Carry forward the pLocalized tag for now, alter below
      pLocalized_new[idx] = pLocalized[idx];

      // Update the deformation gradient tensor to its time n+1 value.
      double J = pDefGrad_new[idx].Determinant();
      if (!(J > 0.) || J > 1.e5) {
        std::cerr
          << "**ERROR** Negative (or huge) Jacobian of deformation gradient."
          << "  Deleting particle " << pParticleID[idx] << "\n";
        std::cerr << "l = " << tensorL[idx] << "\n";
        std::cerr << "F_old = " << pDefGrad[idx] << "\n";
        std::cerr << "J_old = " << pDefGrad[idx].Determinant() << "\n";
        std::cerr << "F_new = " << pDefGrad_new[idx] << "\n";
        std::cerr << "J = " << J << "\n";
        std::cerr << "Temp = " << pTemperature[idx] << "\n";
        std::cerr << "Tm = " << Tm << "\n";
        std::cerr << "DWI = " << matl->getDWIndex() << "\n";
        std::cerr << "X = " << px[idx] << "\n";
        std::cerr << "L.norm()*dt = " << tensorL[idx].Norm() * delT << "\n";
        pLocalized_new[idx] = -999;

        tensorF_new          = pDefGrad[idx];
        pDefGrad_new[idx] = pDefGrad[idx];
        tensorD              = Vaango::Util::Zero;
        tensorL[idx]         = Vaango::Util::Zero;
      }

      // Calculate the current density and deformed volume
      double rho_cur = rho_0 / J;

      // Compute rate of change of specific volume
      double Vdot =
        (pVolume_deformed[idx] - pVolume[idx]) / (pMass[idx] * delT);

      // Calculate rate of deformation tensor (D)
      tensorD = (tensorL[idx] + tensorL[idx].Transpose()) * 0.5;

      // Compute polar decomposition of F (F = RU)
      pDefGrad[idx].polarDecompositionRMB(tensorU, tensorR);

      // Rotate the total rate of deformation tensor back to the
      // material configuration
      tensorD = (tensorR.Transpose()) * (tensorD * tensorR);

      // Calculate the deviatoric part of the non-thermal part
      // of the rate of deformation tensor
      tensorEta = tensorD - Vaango::Util::Identity * (tensorD.Trace() / 3.0);

      pEqStrainRate_new[idx] = tensorD.Norm();

      // Rotate the Cauchy stress back to the
      // material configuration and calculate the deviatoric part
      sigma           = pStress[idx];
      sigma           = (tensorR.Transpose()) * (sigma * tensorR);

      // Rotate internal Cauchy stresses back to the
      // material configuration (only for viscoelasticity)
      d_devStress->rotateInternalStresses(idx, tensorR);

      // Set up the DeformationState 
      DeformationState defState;
      defState.D                    = tensorD;
      defState.devD                 = tensorEta;
      defState.viscoElasticWorkRate = 0.0;

      // Set up the ModelStateBase 
      ModelStateBase state_old;
      state_old.strainRate          = tensorD;
      state_old.eqStrainRate        = pEqStrainRate[idx];
      state_old.eqPlasticStrainRate = pEqPlasticStrainRate[idx];
      state_old.eqPlasticStrain     = pEqPlasticStrain[idx];
      state_old.temperature         = pTemperature[idx];
      state_old.initialTemperature  = d_initialMaterialTemperature;
      state_old.density             = rho_cur;
      state_old.initialDensity      = rho_0;
      state_old.volume              = pVolume[idx];
      state_old.initialVolume       = pMass[idx] / rho_0;
      state_old.bulkModulus         = bulk;
      state_old.initialBulkModulus  = bulk;
      state_old.shearModulus        = shear;
      state_old.initialShearModulus = shear;
      state_old.meltingTemp         = Tm;
      state_old.initialMeltTemp     = Tm;
      state_old.specificHeat        = matl->getSpecificHeat();
      state_old.energy              = pEnergy[idx];
      state_old.porosity            = pIntVar_old[idx].plasticPorosity;
      state_old.setStress(sigma);
      state_old.setPlasticStrain(pPlasticStrain[idx]);

      double temperature    = state_old.temperature;
      Matrix3 sigma_dev_old = state_old.devStress;

      // Get or compute the specific heat
      if (d_computeSpecificHeat) {
        double C_p          = d_Cp->computeSpecificHeat(&state_old);
        state_old.specificHeat = C_p;
      }

      // Calculate the shear modulus and the melting temperature at the
      // start of the time step and update the plasticity state
      double Tm_old = d_melt->computeMeltingTemp(&state_old);
      double mu_old = d_shear->computeShearModulus(&state_old);

      state_old.meltingTemp  = Tm_old;
      state_old.shearModulus = mu_old;

      //-----------------------------------------------------------------------
      // Stage 2:
      //-----------------------------------------------------------------------
      // Create a trial state
      ModelStateBase state_trial(state_old);
      state_trial.eqStrainRate = pEqStrainRate_new[idx];
      state_trial.volume       = pVolume_deformed[idx];

      // Create a new state to store updated results
      ModelStateBase state_new(state_trial);

      // Assume elastic deformation to get a trial stress
      double p_trial = d_eos->computePressure(matl, &state_trial, tensorF_new, tensorD, delT);
      d_devStress->computeDeviatoricStressInc(idx, &state_trial, &defState, delT);
      Matrix3 s_trial = sigma_dev_old + defState.devStressInc;
      Matrix3 sigma_trial = s_trial + Vaango::Util::Identity * p_trial;
      state_trial.setStress(sigma_trial);

      // Calculate flow stress
      double flowStress =
        d_flow->computeFlowStress(&state_trial, delT, d_tol, matl, idx);
      state_trial.yieldStress = flowStress;

      // Material has melted if flowStress <= 0.0
      bool melted  = false;
      bool plastic = false;
      Matrix3 s_new = Vaango::Util::Zero;
      if (temperature > Tm_old || flowStress <= 0.0) {

        melted = true;
        // Set the deviatoric stress to zero
        if (d_doMelting) {
          Matrix3 sigma_new = Vaango::Util::Identity * p_trial;
          state_new.setStress(sigma_new);
        }

        d_flow->updateElastic(idx);

      } else {

        // Evaluate yield condition (returns std::pair<double, Util::YieldStatus>)
        auto result = d_yield->evalYieldCondition(&state_trial);
        if (result.second == Vaango::Util::YieldStatus::IS_ELASTIC) {

          // Set the new stress to the trial stress
          state_new.setStress(sigma_trial);

          // Update the internal variables
          d_flow->updateElastic(idx);

          // Update internal Cauchy stresses (only for viscoelasticity)
          Matrix3 dp = Vaango::Util::Zero;
          d_devStress->updateInternalStresses(idx, dp, &defState, delT);

        } else {

          //  Here set to true, if all conditionals are met (immediately above)
          //  then set to false.
          plastic = true;

          // Use radial return to compute a rough estimate of stress and
          // delta Gamma = dgamma/dt * Delta t = dlambda/dt * Delta t 
          double delGamma = doApproxReturn(delT, matl, idx, 
                                           &state_old, &state_trial, &state_new);
          s_new = state_new.getStress().Deviator();

          // Update internal variables
          d_flow->updatePlastic(idx, delGamma);

          // Update internal Cauchy stresses (only for viscoelasticity)
          Matrix3 devDefRate_p = (state_new.getPlasticStrain() - 
                                  state_old.getPlasticStrain()).Deviator() / delT;
          d_devStress->updateInternalStresses(idx, devDefRate_p, &defState, delT);

        } // end of flow_rule if
      }   // end of temperature if

      double Tm_cur = d_melt->computeMeltingTemp(&state_new);
      double mu_cur = d_shear->computeShearModulus(&state_new);
      state_new.meltingTemp  = Tm_cur;
      state_new.shearModulus = mu_cur;

      // Get the elastic-plastic coupling derivatives
      IsoNonlinHypoelastic elasticityModel(d_elastic.get(), d_intvar.get());

      // Calculate the updated total stress
      double p_new =
        d_eos->computePressure(matl, &state_new, tensorF_new, tensorD, delT);
      Matrix3 sigma_new = s_new + Vaango::Util::Identity * p_new;

      /*
       // **WARNING** Produces negative Tdot
      double Dkk = tensorD.Trace();
      double dTdt_isentropic = d_eos->computeIsentropicTemperatureRate(
        temperature, rho_0, rho_cur, Dkk);
      pdTdt[idx] += dTdt_isentropic;
      */

      // Calculate Tdot from viscoelasticity
      double taylorQuinney = d_initialData.Chi;
      double fac           = taylorQuinney / (rho_cur * state_new.specificHeat);
      double Tdot_VW       = defState.viscoElasticWorkRate * fac;

      pdTdt[idx] += Tdot_VW;

      double de_s = 0.;
      if (flag->d_artificialViscosity) {
        double c_bulk = sqrt(bulk / rho_cur);
        double Dkk = tensorD.Trace();
        p_q[idx]      = artificialBulkViscosity(Dkk, c_bulk, rho_cur, dx_ave);
        de_s          = -p_q[idx] * Dkk / rho_cur;
      } else {
        p_q[idx] = 0.;
        de_s     = 0.;
      }

      // Calculate Tdot due to artificial viscosity
      double Tdot_AV = de_s / state_new.specificHeat;
      pdTdt[idx] += Tdot_AV * include_AV_heating;

      if (pdTdt[idx] < 0.0) {
        std::cout << "dTdt = " << pdTdt[idx]
                  //<< " dTdT_isen = " << dTdt_isentropic
                  << " dTdT_plas = " << Tdot_VW << " dTdT_visc = " << Tdot_AV
                  << "\n";
      }


      // If the particle has already failed, apply various erosion algorithms
      if (flag->d_doErosion) {
        if (pLocalized[idx]) {
          if (d_allowNoTension) {
            if (p_new > 0.0) {
              sigma_new = Vaango::Util::Zero;
            } else {
              sigma_new = Vaango::Util::Identity * p_new;
            }
          }
          if (d_allowNoShear) {
            sigma_new = Vaango::Util::Identity * p_new;
          } else if (d_setStressToZero) {
            sigma_new = Vaango::Util::Zero;
          }
        }
      }

      // Update state->stress
      state_new.setStress(sigma_new);

      //-----------------------------------------------------------------------
      // Stage 3:
      //-----------------------------------------------------------------------
      // Compute porosity/damage/temperature change
      if (!plastic) {

        // Save the updated data
        pIntVar_new[idx]              = pIntVar_old[idx];
        pPlasticStrain_new[idx]       = pPlasticStrain[idx];
        pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
        pEqPlasticStrainRate_new[idx] = 0.0;
        pDamage_new[idx]              = pDamage[idx];
        pPorosity_new[idx]            = pPorosity[idx];

      } else {

        // Update the plastic strain
        pIntVar_new[idx]              = {state_new.eqPlasticStrain, state_new.porosity};
        pPlasticStrain_new[idx]       = state_new.getPlasticStrain();
        pEqPlasticStrain_new[idx]     = state_new.eqPlasticStrain;
        pEqPlasticStrainRate_new[idx] = state_new.eqPlasticStrainRate;

        /*
        if (pPlasticStrainRate_new[idx] > pEqStrainRate_new[idx]) {
          std::cout << "Patch = " << patch->getID() << " particle = " << idx
               << " edot = " << pEqStrainRate_new[idx]
               << " epdot = " << pPlasticStrainRate_new[idx] << "\n";
        }
        */

        // Update the porosity
        if (d_evolvePorosity) {
          pPorosity_new[idx] = state_new.porosity;
        } else {
          pPorosity_new[idx] = pPorosity[idx];
        }

        // Calculate the updated scalar damage parameter
        if (d_evolveDamage) {
          pDamage_new[idx] =
            d_damage->computeScalarDamage(state_new.eqPlasticStrainRate,
                                          sigma_new,
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
        double fac           = taylorQuinney / (rho_cur * state_new.specificHeat);

        // Calculate Tdot (internal plastic heating rate)
        double Tdot_PW = state_new.yieldStress * state_new.eqPlasticStrainRate * fac;

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
            state_new.temperature  = temperature;
            Tm_cur                   = d_melt->computeMeltingTemp(&state_new);
            state_new.meltingTemp  = Tm_cur;
            mu_cur                   = d_shear->computeShearModulus(&state_new);
            state_new.shearModulus = mu_cur;
            state_new.yieldStress =
              d_flow->computeFlowStress(&state_new, delT, d_tol, matl, idx);
            if (!(state_new.yieldStress > 0.0))
              isLocalized = true;
            else {
  
              // Get elastic tangent modulus
              auto C_e = elasticityModel.computeElasticTangentModulus(&state_new);

              // Calculate the elastic-plastic tangent modulus
              Vaango::Tensor::Matrix6Mandel C_ep;
              Vaango::Tensor::Vector6Mandel P_vec, N_vec;
              double H;
              std::tie(C_ep, P_vec, N_vec, H) = 
                computeElasPlasTangentModulus(C_e, &state_new);

              // Initialize localization direction
              Vector direction(0.0, 0.0, 0.0);
              isLocalized =
                d_stable->checkStability(sigma_new, tensorD, 
                                         C_e, P_vec, N_vec, H,
                                         direction);
            }
          }
        }

        // Check 4: Look at maximum stress
        if (d_checkStressTriax) {

          // Compute eigenvalues of the stress tensor
          SymmMatrix3 stress(sigma_new);
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
              if (p_new > 0.0) {
                sigma_new = Vaango::Util::Zero;
              } else {
                sigma_new = Vaango::Util::Identity * p_new;
              }
            } else if (d_allowNoShear) {
              sigma_new = Vaango::Util::Identity * p_new;
            } else if (d_setStressToZero) {
              sigma_new = Vaango::Util::Zero;
            }
          }
        }
      }

      //-----------------------------------------------------------------------
      // Stage 5:
      //-----------------------------------------------------------------------

      // Rotate the stress back to the laboratory coordinates using new R
      // Compute polar decomposition of new F (F = RU)
      tensorF_new.polarDecompositionRMB(tensorU, tensorR);

      sigma_new = (tensorR * sigma_new) * (tensorR.Transpose());

      // Rotate internal Cauchy stresses back to laboratory
      // coordinates (only for viscoelasticity)
      d_devStress->rotateInternalStresses(idx, tensorR);

      // Update the kinematic variables
      pRotation_new[idx] = tensorR;

      // Save the new data
      pStress_new[idx] = sigma_new;

      // Rotate the deformation rate back to the laboratory coordinates
      tensorD = (tensorR * tensorD) * (tensorR.Transpose());

      // Compute the strain energy for non-localized particles
      if (pLocalized_new[idx] == 0) {
        Matrix3 avgStress = (pStress_new[idx] + pStress[idx]) * 0.5;
        double avgVolume  = (pVolume_deformed[idx] + pVolume[idx]) * 0.5;

        double pSpecificStrainEnergy =
          (tensorD(0, 0) * avgStress(0, 0) + tensorD(1, 1) * avgStress(1, 1) +
           tensorD(2, 2) * avgStress(2, 2) +
           2.0 * (tensorD(0, 1) * avgStress(0, 1) +
                  tensorD(0, 2) * avgStress(0, 2) +
                  tensorD(1, 2) * avgStress(1, 2))) *
          avgVolume * delT / pMass[idx];

        pEnergy_new[idx] = pEnergy[idx] + pSpecificStrainEnergy -
                           p_q[idx] * Vdot * delT * include_AV_heating;

        totalStrainEnergy += pSpecificStrainEnergy * pMass[idx];
      } else {
        pEnergy_new[idx] = pEnergy[idx];
      }

      // Compute wave speed at each particle, store the maximum
      double c_dil = sqrt((bulk + 4.0 * mu_cur / 3.0) / rho_cur);
      Vector pVel = pVelocity[idx];
      WaveSpeed   = Vector(Max(c_dil + std::abs(pVel.x()), WaveSpeed.x()),
                         Max(c_dil + std::abs(pVel.y()), WaveSpeed.y()),
                         Max(c_dil + std::abs(pVel.z()), WaveSpeed.z()));

    } // end particle loop

    //__________________________________
    //
    WaveSpeed       = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();

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
           approximate return */
////////////////////////////////////////////////////////////////////////
double
ElasticPlasticHP::doApproxReturn(const double& delT,
                                 const MPMMaterial* matl,
                                 const particleIndex idx,
                                 const ModelStateBase* state_old,
                                 const ModelStateBase* state_trial,
                                 ModelStateBase* state_new) const
{
  // Tolerance
  double tolerance = std::min(delT, 1.0e-6);

  int numSubstep = ((int) (state_trial->q / state_trial->yieldStress) + 1) * 5;
  double substep_delT = delT / (double)numSubstep;

  ModelStateBase state_trial_substep(state_trial);
  ModelStateBase state_old_substep(state_old);

  auto strain_inc = state_old->strainRate * substep_delT;
  auto strain_inc_tr = strain_inc.Trace();
  auto strain_inc_dev = strain_inc.Deviator();

  double deltaGamma = 0.0;
  for (int ii = 0; ii < numSubstep; ++ii) {
 
    double kappa = d_eos->computeBulkModulus(&state_old_substep);
    double mu = d_shear->computeShearModulus(&state_old_substep);
    state_trial_substep.bulkModulus = kappa;
    state_trial_substep.shearModulus = mu;

    auto p_trial_inc = Vaango::Util::Identity * (strain_inc_tr * kappa);
    auto s_trial_inc = strain_inc_dev * (2.0 * mu);
    auto stress_trial_substep = state_old_substep.getStress() + p_trial_inc + s_trial_inc;
    state_trial_substep.setStress(stress_trial_substep);

    #ifdef DEBUG_RETURN
      std::cout << "\n\nSubstep = " << ii << "\n";
    #endif

    double deltaGammaInc = 
     approxHardeningReturn(substep_delT, tolerance, matl, idx, 
                           &state_old_substep, &state_trial_substep, state_new);
    
    state_trial_substep = state_new;
    state_old_substep = state_new;

    deltaGamma += deltaGammaInc;
  }

  // Force the stress on to the yield surface
  double g, deltaGammaTmp, M_mag;
  Vaango::Tensor::Vector6Mandel M_vec, P_vec;
  Vaango::Tensor::Matrix6Mandel C_e;
  std::tie(g, deltaGammaTmp, M_mag, M_vec, P_vec, C_e) = 
    nonHardeningReturn(state_new, tolerance);
  auto stress_new = state_new->getStress() - 
                    Vaango::Tensor::constructMatrix3(P_vec) * deltaGamma;
  state_new->setStress(stress_new);
  
  return deltaGamma;
}

////////////////////////////////////////////////////////////////////////
// Compute the quantity
//             \f$d(\gamma)/dt * \Delta T = \Delta \gamma \f$
//             using Newton iterative root finder */
////////////////////////////////////////////////////////////////////////
double
ElasticPlasticHP::approxHardeningReturn(double delT,
                                        double tolerance,
                                        const MPMMaterial* matl,
                                        const particleIndex idx,
                                        const ModelStateBase* state_old,
                                        const ModelStateBase* state_trial,
                                        ModelStateBase* state_new) const
{
  // Compute the yield stress
  ModelStateBase state_iter(state_trial);
  state_iter.yieldStress = 
    d_flow->computeFlowStress(state_trial, delT, tolerance, matl, idx);

  // Compute the yield function to check if substep is elastic
  double f = d_yield->evalYieldCondition(state_iter.getStress(), &state_iter);

  if (f < 0.0) {

    double deltaGamma = 0.0;

    // Update stress
    state_new->setStress(state_iter.getStress());

    // Update plastic strain tensor
    state_new->setPlasticStrain(state_iter.getPlasticStrain());

    // Update plastic strain rate
    state_new->eqPlasticStrainRate = 0.0;

    // Update internal variables
    state_new->lambdaIncPlastic = deltaGamma;
    state_new->eqPlasticStrain  = state_iter.eqPlasticStrain;
    state_new->porosity         = state_iter.porosity;

    // Compute the yield stress
    state_new->yieldStress = state_iter.yieldStress;

    #ifdef DEBUG_RETURN
      std::cout << "Approx return: Elastic: f = " << f 
                << " q = " << state_iter.q << " sig_y = " << state_iter.yieldStress << "\n";
    #endif

    return deltaGamma;
  }

  // Do plastic non-hardening return
  double g, deltaGamma, M_mag;
  Vaango::Tensor::Vector6Mandel M_vec, P_vec;
  Vaango::Tensor::Matrix6Mandel C_e;
  std::tie(g, deltaGamma, M_mag, M_vec, P_vec, C_e) = 
    nonHardeningReturn(state_trial, tolerance);

  // Apply first order correction of deltaGamma (Brannon  + Leelavanichkul, 2010)
  // to allow for evolving internal variables
  // Calculate derivative wrt internal variables
  MetalIntVar f_intvar;
  d_yield->df_dintvar(state_trial, f_intvar);
  double f_eta1 = f_intvar.plasticPorosity;
  double f_eta2 = f_intvar.eqPlasticStrain;

  // Compute isotropic hardening moduli
  MetalIntVar h_eta;
  d_intvar->computeHardeningModulus(state_trial, h_eta);
  double h_eta1 = h_eta.eqPlasticStrain;
  double h_eta2 = h_eta.plasticPorosity;

  // Calculate derivative wrt plastic strain rate
  //double f_epdot = -d_flow->evalDerivativeWRTStrainRate(state_trial, idx);
  //double h_epdot = 1.0;

  // Compute H
  //double H = - (f_eta1 * h_eta1 + f_eta2 * h_eta2 + f_epdot * h_epdot) / M_mag;
  double H = - (f_eta1 * h_eta1 + f_eta2 * h_eta2) / M_mag;

  // Compute P:N and P:N + H
  double PN = P_vec.transpose() * M_vec;
  double PN_H = PN + H;

  // Apply first order correction of deltaGamma (Brannon  + Leelavanichkul, 2010)
  deltaGamma *= (PN / PN_H);

  // Update stress
  auto stress_new = state_trial->getStress() - 
                    Vaango::Tensor::constructMatrix3(P_vec) * deltaGamma;
  state_new->setStress(stress_new);

  // Elastic strain inc
  auto kappa = state_trial->bulkModulus;
  auto mu = state_trial->shearModulus;
  auto strain_inc = state_old->strainRate * delT;
  auto stress_inc = stress_new - state_old->getStress();
  auto strain_inc_elastic = 
    Vaango::Util::Identity * (1.0 / (9.0 * kappa) - 1.0 / (6.0 * mu)) * stress_inc.Trace() +
    (1.0 / (2.0 * mu)) * stress_inc;
  auto strain_inc_plastic = strain_inc - strain_inc_elastic;

  // Update plastic strain tensor
  state_new->setPlasticStrain(state_trial->getPlasticStrain() + strain_inc_plastic);

  // Update plastic strain rate
  deltaGamma                     = std::max(deltaGamma, 0.0);
  state_new->eqPlasticStrainRate = (strain_inc_plastic / delT).Norm();

  // Update internal variables
  state_new->lambdaIncPlastic = deltaGamma;
  state_new->eqPlasticStrain  = 
    d_intvar->computeInternalVariable("eqPlasticStrain", state_trial, state_new);
  state_new->porosity         = 
    d_intvar->computeInternalVariable("porosity", state_trial, state_new);

  // Compute the yield stress
  state_new->yieldStress =
    d_flow->computeFlowStress(state_new, delT, tolerance, matl, idx);

  #ifdef DEBUG_RETURN
    // Compute new g
    double g_upd = d_yield->evalYieldCondition(stress_new, state_new);

    // Compute N:delta_sigma - deltaGamma * H
    Matrix3 M = d_yield->df_dsigma(state_new);
    M /= M.Norm(); 
    double consistency = M.Contract(stress_new - state_old->getStress()) - deltaGamma * H;

    if (std::abs(g_upd) > tolerance) {
      std::cout << "Approx return: Plastic: g_old = " << f
                << " q = " << state_trial->q << " sig_y = " << state_trial->yieldStress
                << " \ng_new = " << g_upd 
                << " q = " << state_new->q << " sig_y = " << state_new->yieldStress << "\n";
      std::cout << "g_nonhard = " << g << " consistency = " << consistency 
                << " delGamma = " << deltaGamma << "\n";
      std::cout << "ep_inc = " << strain_inc_plastic << "\n deltaGamma * M = " << deltaGamma * M << "\n";
    }
  #endif

  return deltaGamma;
}

std::tuple<double, double, double, 
           Vaango::Tensor::Vector6Mandel, 
           Vaango::Tensor::Vector6Mandel, 
           Vaango::Tensor::Matrix6Mandel>
ElasticPlasticHP::nonHardeningReturn(const ModelStateBase* state_trial,
                                     double tolerance) const
{
  // Compute the projection tensor (P = C:M, C = elastic tangent, M = normal to yield surface)
  auto C_e = d_elastic->computeElasticTangentModulus(state_trial);
  auto M = d_yield->df_dsigma(state_trial);
  //auto N = M;  // Associated plasticity
  auto M_mag = M.Norm();
  M /= M_mag; 
  auto M_vec = Vaango::Tensor::constructVector6Mandel(M);
  auto P_vec = C_e * M_vec;
  auto P = Vaango::Tensor::constructMatrix3(P_vec);

  // Initialize variables
  ModelStateBase state_iter(state_trial);
  double g             = 1.0;
  double deltaGamma    = 0.0;
  double deltaGammaOld = deltaGamma;
  int count = 0;
  do {

    ++count;

    // Compute the projected stress
    auto stress_iter = state_trial->getStress() - P * deltaGamma;
    state_iter.setStress(stress_iter);

    // Compute the yield function
    g = d_yield->evalYieldCondition(stress_iter, &state_iter);

    // Compute the derivative of the yield function wrt delGamma
    auto df_dsigma = d_yield->df_dsigma(stress_iter, &state_iter);
    double Dg = - df_dsigma.Contract(P);

    // Update deltaGamma
    deltaGammaOld = deltaGamma;
    deltaGamma -= g / Dg;
    //std::cout << "iter = " << count << " g = " << g << " dg = " << Dg
    //          << "lambda_old = " << deltaGammaOld << " lambda = " << deltaGamma << "\n";

    if (std::isnan(g) || std::isnan(deltaGamma)) {
      std::cout << "iter = " << count << " g = " << g
                << " Dg = " << Dg << " deltaGamma = " << deltaGamma
                << " sigy = " << state_iter.yieldStress 
                << " epdot = " << state_iter.eqPlasticStrainRate
                << " phi = " << state_iter.porosity
                << " ep = " << state_iter.eqPlasticStrain << "\n";
      throw InternalError("Nans in computation or delta gamma < 0", __FILE__, __LINE__);
    }

    state_iter.lambdaIncPlastic = deltaGamma;

    if (std::abs(deltaGamma/deltaGammaOld - 1.0) < 1.0e-3*tolerance) break;
    if (count > 100) break;


  } while (std::abs(g) > tolerance);

  #ifdef DEBUG_RETURN
  if (deltaGamma < 0) {
    std::cout << "deltaGamma = " << deltaGamma << " iter = " << count
              << " sig_y = " << state_iter.yieldStress
              << " epdot = " << state_iter.eqPlasticStrainRate
              << " ep = " << state_iter.eqPlasticStrain
              << " T = " << state_iter.temperature 
              << " Tm = " << state_iter.meltingTemp
              << "\n";
  }
  #endif

  return std::make_tuple(g, deltaGamma, M_mag, M_vec, P_vec, C_e);
}

std::tuple<Vaango::Tensor::Matrix6Mandel,
           Vaango::Tensor::Vector6Mandel,
           Vaango::Tensor::Vector6Mandel, double>
ElasticPlasticHP::computeElasPlasTangentModulus(Vaango::Tensor::Matrix6Mandel& C_e,
                                                const ModelStateBase* state) const
{
  // Calculate the derivative of the yield function wrt sigma
  Matrix3 N = d_yield->df_dsigma(state);
  double N_mag = N.Norm();
  N /= N_mag;
  auto N_vec = Vaango::Tensor::constructVector6Mandel(N);
  auto M_vec = N_vec; // associated plasticity

  // Calculate derivative wrt internal variables
  MetalIntVar f_intvar;
  d_yield->df_dintvar(state, f_intvar);
  double f_eta1 = f_intvar.plasticPorosity;
  double f_eta2 = f_intvar.eqPlasticStrain;

  // Compute isotropic hardening moduli
  MetalIntVar h_eta;
  d_intvar->computeHardeningModulus(state, h_eta);
  double h_eta1 = h_eta.eqPlasticStrain;
  double h_eta2 = h_eta.plasticPorosity;

  // Compute P = C:M + Z
  auto P_vec = C_e * M_vec;

  // Compute H
  double H = - (f_eta1 * h_eta1 + f_eta2 * h_eta2) / N_mag;

  // Compute P:N + H
  double PN_H = P_vec.transpose() * N_vec + H;

  // Compute P(C:N)
  auto CN_vec = C_e * N_vec;
  auto P_CN = Vaango::Tensor::constructMatrix6Mandel(P_vec, CN_vec);

  // Compute C_ep
  auto C_ep = C_e - P_CN / PN_H;
  return std::make_tuple(C_ep, P_vec, N_vec, H);
}

//______________________________________________________________________
//
void
ElasticPlasticHP::computeStressTensorImplicit(const PatchSubset* patches,
                                              const MPMMaterial* matl,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw)
{
  // Constants
  int dwi              = matl->getDWIndex();
  Ghost::GhostType gac = Ghost::AroundCells;
  Matrix3 One;
  One.Identity();
  Matrix3 Zero(0.0);

  double bulk  = d_initialData.Bulk;
  double shear = d_initialData.Shear;
  double alpha = d_initialData.alpha;
  double rho_0 = matl->getInitialDensity();
  double Tm    = matl->getMeltTemperature();

  // Do thermal expansion?
  if (!flag->d_doThermalExpansion) {
    alpha = 0;
  }

  // Particle and Grid data
  delt_vartype delT;
  constParticleVariable<int> pLocalized;
  constParticleVariable<double> pMass, pVolume, pTempPrev, pTemperature,
    pEqPlasticStrain, pDamage, pPorosity, pEqStrainRate, pEqPlasticStrainRate,
    pEnergy;

  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad, pStress, pRotation, pPlasticStrain;
  constNCVariable<Vector> gDisp;

  ParticleVariable<int> pLocalized_new;
  ParticleVariable<Matrix3> pDefGrad_new, pStress_new, pRotation_new, pPlasticStrain_new;
  ParticleVariable<double> pVolume_deformed, pEqPlasticStrain_new, pDamage_new,
    pPorosity_new, pEqStrainRate_new, pEqPlasticStrainRate_new, pdTdt, pEnergy_new;

  // Local variables
  Matrix3 DispGrad(0.0); // Displacement gradient
  Matrix3 DefGrad, incDefGrad, incFFt, incFFtInv, Rotation, RightStretch;
  Matrix3 incTotalStrain(0.0), incThermalStrain(0.0), incStrain(0.0);
  Matrix3 sigma(0.0), trialStress(0.0), s_trial(0.0);
  DefGrad.Identity();
  incDefGrad.Identity();
  incFFt.Identity();
  incFFtInv.Identity();
  Rotation.Identity(), RightStretch.Identity();

  // CSTi << getpid()
  //     << "ComputeStressTensorImplicit: In : Matl = " << matl << " id = "
  //     << matl->getDWIndex() <<  " patch = "
  //     << (patches->get(0))->getID();

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    // Get the set of particles
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // GET GLOBAL DATA
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    new_dw->get(gDisp, lb->dispNewLabel, dwi, patch, gac, 1);

    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pTempPrev, lb->pTempPreviousLabel, pset);
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);
    old_dw->get(pStress, lb->pStressLabel, pset);

    // GET LOCAL DATA
    old_dw->get(pRotation, pRotationLabel, pset);
    old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);
    old_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
    old_dw->get(pEqPlasticStrainRate, pEqPlasticStrainRateLabel, pset);
    old_dw->get(pDamage, pDamageLabel, pset);
    old_dw->get(pEqStrainRate, pEqStrainRateLabel, pset);
    old_dw->get(pPorosity, pPorosityLabel, pset);
    old_dw->get(pLocalized, pLocalizedLabel, pset);
    old_dw->get(pEnergy, pEnergyLabel, pset);

    new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->getModifiable(pVolume_deformed, lb->pVolumeLabel_preReloc, pset);

    // Create and allocate arrays for storing the updated information
    // GLOBAL
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);

    // LOCAL
    new_dw->allocateAndPut(pRotation_new, pRotationLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEqStrainRate_new, pEqStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrain_new, pEqPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrainRate_new, pEqPlasticStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(pDamage_new, pDamageLabel_preReloc, pset);
    new_dw->allocateAndPut(pPorosity_new, pPorosityLabel_preReloc, pset);
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);
    new_dw->allocateAndPut(pEnergy_new, pEnergyLabel_preReloc, pset);

    // Get the plastic strain
    d_flow->getInternalVars(pset, old_dw);
    d_devStress->getInternalVars(pset, old_dw);
    d_flow->allocateAndPutInternalVars(pset, new_dw);
    d_devStress->allocateAndPutInternalVars(pset, new_dw);

    //__________________________________
    // Special case for rigid materials
    double totalStrainEnergy = 0.0;
    if (matl->getIsRigid()) {
      for (auto idx : *pset) {
        pRotation_new[idx]            = pRotation[idx];
        pEqStrainRate_new[idx]        = pEqStrainRate[idx];
        pPlasticStrain_new[idx]       = pPlasticStrain[idx];
        pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
        pEqPlasticStrainRate_new[idx] = 0.0;
        pDamage_new[idx]              = pDamage[idx];
        pPorosity_new[idx]            = pPorosity[idx];
        pLocalized_new[idx]           = pLocalized[idx];
        pStress_new[idx]              = Vaango::Util::Zero;
        pDefGrad_new[idx]             = Vaango::Util::Identity;
        pVolume_deformed[idx]         = pVolume[idx];
        pdTdt[idx]                    = 0.0;
      }

      if (flag->d_reductionVars->accStrainEnergy ||
          flag->d_reductionVars->strainEnergy) {
        new_dw->put(sum_vartype(totalStrainEnergy), lb->StrainEnergyLabel);
      }
      continue;
    }

    //__________________________________
    // Standard case for deformable materials
    // Loop thru particles
    for (auto idx : *pset) {

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx]       = 0.0;
      pEnergy_new[idx] = pEnergy[idx];

      // Update the deformation gradient
      double J = DefGrad.Determinant();

      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        std::cerr << getpid()
                  << "**ERROR** Negative Jacobian of deformation gradient"
                  << "\n";
        throw InternalError("Negative Jacobian", __FILE__, __LINE__);
      }

      // Calculate the current density and deformed volume
      double rho_cur = rho_0 / J;
      double volold  = (pMass[idx] / rho_0);

      // Compute polar decomposition of F (F = VR)
      // (**NOTE** This is being done to provide reasonable starting
      //           values for R and V if the incremental algorithm
      //           for the polar decomposition is used in the explicit
      //           calculations following an implicit calculation.)
      pDefGrad[idx].polarDecompositionRMB(RightStretch, Rotation);
      pRotation_new[idx] = Rotation;

      // Compute the current strain and strain rate
      incFFt         = incDefGrad * incDefGrad.Transpose();
      incFFtInv      = incFFt.Inverse();
      incTotalStrain = (One - incFFtInv) * 0.5;

      // Try the small strain option (BB - 06/22/05)
      // Calculate the strain
      // incTotalStrain = (DispGrad + DispGrad.Transpose())*0.5;

      pEqStrainRate_new[idx] = incTotalStrain.Norm() / delT;

      // Compute thermal strain
      double temperature = pTemperature[idx];
      double incT        = temperature - pTempPrev[idx];
      incThermalStrain   = One * (alpha * incT);
      incStrain          = incTotalStrain - incThermalStrain;

      // Compute pressure and deviatoric stress at t_n and
      // the volumetric strain and deviatoric strain increments at t_n+1
      sigma           = pStress[idx];
      double pressure = sigma.Trace() / 3.0;
      Matrix3 tensorS = sigma - One * pressure;

      // Set up the ModelStateBase
      ModelStateBase state_old;
      state_old.strainRate          = incTotalStrain / delT;
      state_old.eqStrainRate        = pEqStrainRate_new[idx];
      state_old.eqPlasticStrainRate = pEqPlasticStrainRate[idx];
      state_old.eqPlasticStrain     = pEqPlasticStrain[idx];
      state_old.pressure            = pressure;
      state_old.temperature         = temperature;
      state_old.initialTemperature  = d_initialMaterialTemperature;
      state_old.density             = rho_cur;
      state_old.initialDensity      = rho_0;
      state_old.volume              = pVolume_deformed[idx];
      state_old.initialVolume       = volold;
      state_old.bulkModulus         = bulk;
      state_old.initialBulkModulus  = bulk;
      state_old.shearModulus        = shear;
      state_old.initialShearModulus = shear;
      state_old.meltingTemp         = Tm;
      state_old.initialMeltTemp     = Tm;
      state_old.specificHeat        = matl->getSpecificHeat();
      state_old.energy              = pEnergy[idx];
      state_old.setStress(sigma);
      state_old.setPlasticStrain(pPlasticStrain[idx]);

      // Get or compute the specific heat
      if (d_computeSpecificHeat) {
        double C_p          = d_Cp->computeSpecificHeat(&state_old);
        state_old.specificHeat = C_p;
      }

      // Calculate the shear modulus and the melting temperature at the
      // start of the time step and update the plasticity state
      double Tm_cur       = d_melt->computeMeltingTemp(&state_old);
      state_old.meltingTemp  = Tm_cur;
      double mu_cur       = d_shear->computeShearModulus(&state_old);
      state_old.shearModulus = mu_cur;

      // The calculation of tensorD and tensorEta are so that the deviatoric
      // stress calculation will work.

      Matrix3 tensorD = incStrain / delT;

      // Calculate the deviatoric part of the non-thermal part
      // of the rate of deformation tensor
      Matrix3 tensorEta = tensorD - One * (tensorD.Trace() / 3.0);

      DeformationState defState;
      defState.D    = tensorD;
      defState.devD = tensorEta;

      // Assume elastic deformation to get a trial stress
      double p_trial = d_eos->computePressure(matl, &state_old, pDefGrad_new[idx], tensorD, delT);

      // This is simply the previous timestep deviatoric stress plus a
      // deviatoric elastic increment based on the shear modulus supplied by
      // the strength routine in use.
      d_devStress->computeDeviatoricStressInc(idx, &state_old, &defState, delT);
      s_trial      = tensorS + defState.devStressInc;
      trialStress = s_trial + Vaango::Util::Identity * p_trial;

      // Create a trial state
      ModelStateBase state_trial(state_old);
      state_trial.setStress(trialStress);

      // Calculate flow stress
      double flowStress =
        d_flow->computeFlowStress(&state_trial, delT, d_tol, matl, idx);
      state_trial.yieldStress = flowStress;

      // Get the current porosity
      double porosity = pPorosity[idx];
      state_trial.porosity = porosity;

      // Evaluate yield condition (returns std::pair<double, Util::YieldStatus>)
      auto result = d_yield->evalYieldCondition(&state_trial);
      if (result.second == Vaango::Util::YieldStatus::IS_ELASTIC) {

        // Save the updated data
        pStress_new[idx]              = trialStress;
        pPlasticStrain_new[idx]       = pPlasticStrain[idx];
        pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
        pEqPlasticStrainRate_new[idx] = 0.0;
        pDamage_new[idx]              = pDamage[idx];
        pPorosity_new[idx]            = pPorosity[idx];

        // Update the internal variables
        d_flow->updateElastic(idx);

        // Update internal Cauchy stresses (only for viscoelasticity)
        Matrix3 dp = Zero;
        d_devStress->updateInternalStresses(idx, dp, &defState, delT);

        // Compute stability criterion
        pLocalized_new[idx] = pLocalized[idx];

      } else {
        //__________________________________
        // Radial Return
        // Do Newton iteration to compute delGamma and updated
        // plastic strain, plastic strain rate, and yield stress
        ModelStateBase state_new(state_trial);
        double delGamma = doApproxReturn(delT, matl, idx, 
                                         &state_old, &state_trial, &state_new);
        pStress_new[idx]              = state_new.getStress();
        pPlasticStrain_new[idx]       = state_new.getPlasticStrain();
        pEqPlasticStrain_new[idx]     = state_new.eqPlasticStrain;
        pEqPlasticStrainRate_new[idx] = state_new.eqPlasticStrainRate;

        // Update internal Cauchy stresses (only for viscoelasticity)
        Matrix3 devDefRate_p = (pPlasticStrain_new[idx] - pPlasticStrain[idx]).Deviator() / delT;
        d_devStress->updateInternalStresses(idx, devDefRate_p, &defState, delT);

        // Update the porosity
        pPorosity_new[idx] = pPorosity[idx];

        if (d_evolvePorosity) {
          Matrix3 tensorD    = incStrain / delT;
          double ep          = state_trial.eqPlasticStrain;
          pPorosity_new[idx] = updatePorosity(tensorD, delT, porosity, ep);
        }

        // Calculate the updated scalar damage parameter
        if (d_evolveDamage)
          pDamage_new[idx] =
            d_damage->computeScalarDamage(state_trial.eqPlasticStrainRate,
                                          pStress_new[idx],
                                          temperature,
                                          delT,
                                          matl,
                                          d_tol,
                                          pDamage[idx]);
        else
          pDamage_new[idx] = pDamage[idx];

        // Calculate rate of temperature increase due to plastic strain
        double taylorQuinney = d_initialData.Chi;
        double fac           = taylorQuinney / (rho_cur * state_trial.specificHeat);

        // Calculate Tdot (internal plastic heating rate)
        double Tdot = state_trial.yieldStress * state_trial.eqPlasticStrainRate * fac;
        pdTdt[idx]  = Tdot;

        // No failure implemented for implcit time integration
        pLocalized_new[idx] = pLocalized[idx];

        // Update internal variables in the flow model
        d_flow->updatePlastic(idx, delGamma);
      }

      // Compute the strain energy for non-localized particles
      if (pLocalized_new[idx] == 0) {
        Matrix3 avgStress    = (pStress_new[idx] + pStress[idx]) * 0.5;
        double pStrainEnergy = (incStrain(0, 0) * avgStress(0, 0) +
                                incStrain(1, 1) * avgStress(1, 1) +
                                incStrain(2, 2) * avgStress(2, 2) +
                                2.0 * (incStrain(0, 1) * avgStress(0, 1) +
                                       incStrain(0, 2) * avgStress(0, 2) +
                                       incStrain(1, 2) * avgStress(1, 2))) *
                               pVolume_deformed[idx] * delT;
        totalStrainEnergy += pStrainEnergy;
      }
    }

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(totalStrainEnergy), lb->StrainEnergyLabel);
    }
    // delete interpolator;
  }
}
//______________________________________________________________________
//
void
ElasticPlasticHP::addComputesAndRequires(Task* task,
                                         const MPMMaterial* matl,
                                         const PatchSet* patches,
                                         const bool recurse,
                                         const bool SchedParent) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForImplicitHypo(task, matlset, true, recurse, SchedParent);

  Ghost::GhostType gnone = Ghost::None;
  if (SchedParent) {
    // For subscheduler
    task->requires(Task::ParentOldDW, lb->pTempPreviousLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pTemperatureLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, pPlasticStrainLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, pEqPlasticStrainLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, pEqPlasticStrainRateLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, pPorosityLabel, matlset, gnone);

    task->computes(pPlasticStrainLabel_preReloc, matlset);
    task->computes(pEqPlasticStrainLabel_preReloc, matlset);
    task->computes(pEqPlasticStrainRateLabel_preReloc, matlset);
    task->computes(pPorosityLabel_preReloc, matlset);
  } else {
    // For scheduleIterate
    task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pTemperatureLabel, matlset, gnone);
    task->requires(Task::OldDW, pPlasticStrainLabel, matlset, gnone);
    task->requires(Task::OldDW, pEqPlasticStrainLabel, matlset, gnone);
    task->requires(Task::OldDW, pEqPlasticStrainRateLabel, matlset, gnone);
    task->requires(Task::OldDW, pPorosityLabel, matlset, gnone);
  }

  // Add internal evolution variables computed by flow model
  d_flow->addComputesAndRequires(task, matl, patches, recurse, SchedParent);

  // Deviatoric Stress Model
  d_devStress->addComputesAndRequires(task, matl, SchedParent);
}
//______________________________________________________________________
//
void
ElasticPlasticHP::computeStressTensorImplicit(const PatchSubset* patches,
                                              const MPMMaterial* matl,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw,
                                              Solver* solver,
                                              const bool)
{
  // Constants
  Ghost::GhostType gac = Ghost::AroundCells;
  Matrix3 One;
  One.Identity();
  Matrix3 Zero(0.0);

  double bulk  = d_initialData.Bulk;
  double shear = d_initialData.Shear;
  double alpha = d_initialData.alpha;
  double rho_0 = matl->getInitialDensity();
  double Tm    = matl->getMeltTemperature();

  // Do thermal expansion?
  if (!flag->d_doThermalExpansion) {
    alpha = 0;
  }

  // Data location
  int dwi = matl->getDWIndex();
  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);

  // Particle and Grid data
  delt_vartype delT;
  constParticleVariable<double> pMass, pTempPrev, pTemperature, pEqPlasticStrain,
    pEqPlasticStrainRate, pPorosity;

  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad, pStress, pPlasticStrain;
  constNCVariable<Vector> gDisp;

  ParticleVariable<Matrix3> pDefGrad_new, pStress_new, pPlasticStrain_new;
  ParticleVariable<double> pVolume_deformed, pEqPlasticStrain_new,
    pEqPlasticStrainRate_new;

  // Local variables
  Matrix3 DispGrad(0.0); // Displacement gradient
  Matrix3 DefGrad, incDefGrad, incFFt, incFFtInv;
  Matrix3 incTotalStrain(0.0), incThermalStrain(0.0), incStrain(0.0);
  Matrix3 sigma(0.0), trialStress(0.0), s_trial(0.0);
  DefGrad.Identity();
  incDefGrad.Identity();
  incFFt.Identity();
  incFFtInv.Identity();

  // For B matrices
  double D[6][6];
  double B[6][24];
  double Bnl[3][24];
  double Kmatrix[24][24];
  int dof[24];
  double v[576];

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    // Get interpolation functions
    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    // Get patch indices for parallel solver
    IntVector lowIndex  = patch->getNodeLowIndex();
    IntVector highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    Array3<int> l2g(lowIndex, highIndex);
    solver->copyL2G(l2g, patch);

    // Get grid size
    Vector dx      = patch->dCell();
    double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

    // Get the set of particles
    ParticleSubset* pset = parent_old_dw->getParticleSubset(dwi, patch);

    // GET GLOBAL DATA
    old_dw->get(gDisp, lb->dispNewLabel, dwi, patch, gac, 1);

    parent_old_dw->get(delT, lb->delTLabel, getLevel(patches));
    parent_old_dw->get(pTempPrev, lb->pTempPreviousLabel, pset);
    parent_old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    parent_old_dw->get(px, lb->pXLabel, pset);
    parent_old_dw->get(pSize, lb->pSizeLabel, pset);
    parent_old_dw->get(pMass, lb->pMassLabel, pset);
    parent_old_dw->get(pDefGrad, lb->pDefGradLabel, pset);
    parent_old_dw->get(pStress, lb->pStressLabel, pset);

    // GET LOCAL DATA
    parent_old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);
    parent_old_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
    parent_old_dw->get(pEqPlasticStrainRate, pEqPlasticStrainRateLabel, pset);
    parent_old_dw->get(pPorosity, pPorosityLabel, pset);

    new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->getModifiable(pVolume_deformed, lb->pVolumeLabel_preReloc, pset);

    // Create and allocate arrays for storing the updated information
    // GLOBAL
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    // LOCAL
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrain_new, pEqPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrainRate_new, pEqPlasticStrainRateLabel_preReloc, pset);

    // Special case for rigid materials
    if (matl->getIsRigid()) {
      for (auto idx : *pset) {
        pPlasticStrain_new[idx]       = pPlasticStrain[idx];
        pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
        pEqPlasticStrainRate_new[idx] = 0.0;
        pStress_new[idx]              = Vaango::Util::Zero;
        pDefGrad_new[idx]             = Vaango::Util::Identity;
        pVolume_deformed[idx]         = pMass[idx] / rho_0;
      }
      continue;
    }

    // Standard case for deformable materials
    // Loop thru particles
    for (auto idx : *pset) {

      // Calculate the displacement gradient
      interpolator->findCellAndShapeDerivatives(
        px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);
      computeGradAndBmats(DispGrad, ni, d_S, oodx, gDisp, l2g, B, Bnl, dof);

      // Compute the deformation gradient increment
      incDefGrad = DispGrad + One;

      // Update the deformation gradient
      DefGrad              = incDefGrad * pDefGrad[idx];
      pDefGrad_new[idx]    = DefGrad;
      double J             = DefGrad.Determinant();

      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        std::cerr << getpid()
                  << "**ERROR** Negative Jacobian of deformation gradient"
                  << "\n";
        throw InternalError("Negative Jacobian", __FILE__, __LINE__);
      }

      // Calculate the current density and deformed volume
      double rho_cur        = rho_0 / J;
      double volold         = (pMass[idx] / rho_0);
      pVolume_deformed[idx] = volold * J;

      // Compute the current strain and strain rate
      incFFt         = incDefGrad * incDefGrad.Transpose();
      incFFtInv      = incFFt.Inverse();
      incTotalStrain = (One - incFFtInv) * 0.5;

      double pEqStrainRate_new = incTotalStrain.Norm() / delT;

      // Compute thermal strain increment
      double temperature = pTemperature[idx];
      double incT        = temperature - pTempPrev[idx];
      incThermalStrain   = One * (alpha * incT);
      incStrain          = incTotalStrain - incThermalStrain;

      // Compute pressure and deviatoric stress at t_n and
      // the volumetric strain and deviatoric strain increments at t_n+1
      sigma           = pStress[idx];
      double pressure = sigma.Trace() / 3.0;
      Matrix3 tensorS = sigma - One * pressure;

      // Set up the ModelStateBase
      ModelStateBase state_old;
      state_old.strainRate          = incTotalStrain / delT;
      state_old.eqStrainRate        = pEqStrainRate_new;
      state_old.eqPlasticStrainRate = pEqPlasticStrainRate[idx];
      state_old.eqPlasticStrain     = pEqPlasticStrain[idx];
      state_old.pressure            = pressure;
      state_old.temperature         = temperature;
      state_old.initialTemperature  = d_initialMaterialTemperature;
      state_old.density             = rho_cur;
      state_old.initialDensity      = rho_0;
      state_old.volume              = pVolume_deformed[idx];
      state_old.initialVolume       = volold;
      state_old.bulkModulus         = bulk;
      state_old.initialBulkModulus  = bulk;
      state_old.shearModulus        = shear;
      state_old.initialShearModulus = shear;
      state_old.meltingTemp         = Tm;
      state_old.initialMeltTemp     = Tm;
      state_old.specificHeat        = matl->getSpecificHeat();
      state_old.setStress(sigma);
      state_old.setPlasticStrain(pPlasticStrain[idx]);

      // Get or compute the specific heat
      if (d_computeSpecificHeat) {
        double C_p          = d_Cp->computeSpecificHeat(&state_old);
        state_old.specificHeat = C_p;
      }

      // Calculate the shear modulus and the melting temperature at the
      // start of the time step and update the plasticity state
      double Tm_cur       = d_melt->computeMeltingTemp(&state_old);
      state_old.meltingTemp  = Tm_cur;
      double mu_cur       = d_shear->computeShearModulus(&state_old);
      state_old.shearModulus = mu_cur;

      // The calculation of tensorD and tensorEta are so that the deviatoric
      // stress calculation will work.
      Matrix3 tensorD = incStrain / delT;

      // Calculate the deviatoric part of the non-thermal part
      // of the rate of deformation tensor
      Matrix3 tensorEta = tensorD - One * (tensorD.Trace() / 3.0);

      DeformationState defState;
      defState.D    = tensorD;
      defState.devD = tensorEta;

      // Assume elastic deformation to get a trial deviatoric stress
      double p_trial = d_eos->computePressure(matl, &state_old, pDefGrad_new[idx], tensorD, delT);

      // This is simply the previous timestep deviatoric stress plus a
      // deviatoric elastic increment based on the shear modulus supplied by
      // the strength routine in use.
      d_devStress->computeDeviatoricStressInc(idx, &state_old, &defState, delT);
      s_trial      = tensorS + defState.devStressInc;
      trialStress = s_trial + Vaango::Util::Identity * p_trial;

      // Create a trial state
      ModelStateBase state_trial(state_old);
      state_trial.setStress(trialStress);

      // Calculate flow stress
      double flowStress =
        d_flow->computeFlowStress(&state_trial, delT, d_tol, matl, idx);
      state_trial.yieldStress = flowStress;

      // Get the current porosity
      double porosity = pPorosity[idx];
      state_trial.porosity = porosity;

      // Evaluate yield condition (returns std::pair<double, Util::YieldStatus>)
      auto result = d_yield->evalYieldCondition(&state_trial);
      if (result.second == Vaango::Util::YieldStatus::IS_ELASTIC) {

        // Save the updated data
        pStress_new[idx]              = trialStress;
        pPlasticStrain_new[idx]       = pPlasticStrain[idx];
        pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
        pEqPlasticStrainRate_new[idx] = 0.0;

        // Update internal Cauchy stresses (only for viscoelasticity)
        Matrix3 dp = Vaango::Util::Zero;
        d_devStress->updateInternalStresses(idx, dp, &defState, delT);

        // Get elastic tangent modulus
        auto C_e = d_elastic->computeElasticTangentModulus(&state_trial);
        for (int ii = 0; ii < 6; ++ii) {
          for (int jj = 0; jj < 6; ++jj) {
            D[ii][jj] = C_e(ii, jj);
          }
        }

        /*
        computeElasticTangentModulus(bulk, shear, D);
        */

      } else {

        //__________________________________
        // Radial Return
        // Do Newton iteration to compute delGamma and updated
        // plastic strain, plastic strain rate, and yield stress
        ModelStateBase state_new(state_trial);
        [[maybe_unused]] double delGamma = doApproxReturn(delT, matl, idx, 
                                         &state_old, &state_trial, &state_new);
        pStress_new[idx]              = state_new.getStress();
        pPlasticStrain_new[idx]       = state_new.getPlasticStrain();
        pEqPlasticStrain_new[idx]     = state_new.eqPlasticStrain;
        pEqPlasticStrainRate_new[idx] = state_new.eqPlasticStrainRate;

        // Update internal Cauchy stresses (only for viscoelasticity)
        Matrix3 dp = (pPlasticStrain_new[idx] - pPlasticStrain[idx]).Deviator() / delT;
        d_devStress->updateInternalStresses(idx, dp, &defState, delT);

        // Get elastic tangent modulus
        auto C_e = d_elastic->computeElasticTangentModulus(&state_new);

        // Calculate the elastic-plastic tangent modulus
        Vaango::Tensor::Matrix6Mandel C_ep;
        Vaango::Tensor::Vector6Mandel P_vec, N_vec;
        double H;
        std::tie(C_ep, P_vec, N_vec, H) = computeElasPlasTangentModulus(C_e, &state_new);
        for (int ii = 0; ii < 6; ++ii) {
          for (int jj = 0; jj < 6; ++jj) {
            D[ii][jj] = C_ep(ii, jj);
          }
        }

        /*
        computeEPlasticTangentModulus(
          bulk, shear, delGamma, s_trial, idx, &state_trial, D, true);
        */
      }

      // Compute K matrix = Kmat + Kgeo
      computeStiffnessMatrix(
        B, Bnl, D, pStress[idx], volold, pVolume_deformed[idx], Kmatrix);

      // Assemble into global K matrix
      for (int ii = 0; ii < 24; ii++) {
        for (int jj = 0; jj < 24; jj++) {
          v[24 * ii + jj] = Kmatrix[ii][jj];
        }
      }
      solver->fillMatrix(24, dof, 24, dof, v);
    }
  }
}

//______________________________________________________________________
/*! Compute the elastic tangent modulus tensor for isotropic
    materials
    Assume: [stress] = [s11 s22 s33 s12 s23 s31]
            [strain] = [e11 e22 e33 2e12 2e23 2e31]
*/
//______________________________________________________________________
//
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
//______________________________________________________________________
//
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
//______________________________________________________________________
//
/*! Compute K matrix */
void
ElasticPlasticHP::computeStiffnessMatrix(const double B[6][24],
                                         const double Bnl[3][24],
                                         const double D[6][6],
                                         const Matrix3& sig,
                                         const double& vol_old,
                                         const double& vol_new,
                                         double Kmatrix[24][24])
{
  // Kmat = B.transpose()*D*B*volold
  double Kmat[24][24];
  BtDB(B, D, Kmat);

  // Kgeo = Bnl.transpose*sig*Bnl*volnew;
  double Kgeo[24][24];
  BnlTSigBnl(sig, Bnl, Kgeo);

  for (int ii = 0; ii < 24; ii++) {
    for (int jj = ii; jj < 24; jj++) {
      Kmatrix[ii][jj] = Kmat[ii][jj] * vol_old + Kgeo[ii][jj] * vol_new;
    }
  }
  for (int ii = 0; ii < 24; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      Kmatrix[ii][jj] = Kmatrix[jj][ii];
    }
  }
}

void
ElasticPlasticHP::BnlTSigBnl(const Matrix3& sig,
                             const double Bnl[3][24],
                             double Kgeo[24][24]) const
{
  double t1, t10, t11, t12, t13, t14, t15, t16, t17;
  double t18, t19, t2, t20, t21, t22, t23, t24, t25;
  double t26, t27, t28, t29, t3, t30, t31, t32, t33;
  double t34, t35, t36, t37, t38, t39, t4, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t5;
  double t50, t51, t52, t53, t54, t55, t56, t57, t58;
  double t59, t6, t60, t61, t62, t63, t64, t65, t66;
  double t67, t68, t69, t7, t70, t71, t72, t73, t74;
  double t75, t77, t78, t8, t81, t85, t88, t9, t90;
  double t79, t82, t83, t86, t87, t89;

  t1  = Bnl[0][0] * sig(0, 0);
  t4  = Bnl[0][0] * sig(0, 0);
  t2  = Bnl[0][0] * sig(0, 1);
  t3  = Bnl[0][0] * sig(0, 2);
  t5  = Bnl[1][1] * sig(1, 1);
  t8  = Bnl[1][1] * sig(1, 1);
  t6  = Bnl[1][1] * sig(1, 2);
  t7  = Bnl[1][1] * sig(0, 1);
  t9  = Bnl[2][2] * sig(2, 2);
  t12 = Bnl[2][2] * sig(2, 2);
  t10 = Bnl[2][2] * sig(0, 2);
  t11 = Bnl[2][2] * sig(1, 2);
  t13 = Bnl[0][3] * sig(0, 0);
  t16 = Bnl[0][3] * sig(0, 0);
  t14 = Bnl[0][3] * sig(0, 1);
  t15 = Bnl[0][3] * sig(0, 2);
  t17 = Bnl[1][4] * sig(1, 1);
  t20 = Bnl[1][4] * sig(1, 1);
  t18 = Bnl[1][4] * sig(1, 2);
  t19 = Bnl[1][4] * sig(0, 1);
  t21 = Bnl[2][5] * sig(2, 2);
  t22 = Bnl[2][5] * sig(0, 2);
  t23 = Bnl[2][5] * sig(1, 2);
  t24 = Bnl[2][5] * sig(2, 2);
  t25 = Bnl[0][6] * sig(0, 0);
  t26 = Bnl[0][6] * sig(0, 1);
  t27 = Bnl[0][6] * sig(0, 2);
  t28 = Bnl[0][6] * sig(0, 0);
  t29 = Bnl[1][7] * sig(1, 1);
  t30 = Bnl[1][7] * sig(1, 2);
  t31 = Bnl[1][7] * sig(0, 1);
  t32 = Bnl[1][7] * sig(1, 1);
  t33 = Bnl[2][8] * sig(2, 2);
  t34 = Bnl[2][8] * sig(0, 2);
  t35 = Bnl[2][8] * sig(1, 2);
  t36 = Bnl[2][8] * sig(2, 2);
  t37 = Bnl[0][9] * sig(0, 0);
  t38 = Bnl[0][9] * sig(0, 1);
  t39 = Bnl[0][9] * sig(0, 2);
  t40 = Bnl[0][9] * sig(0, 0);
  t41 = Bnl[1][10] * sig(1, 1);
  t42 = Bnl[1][10] * sig(1, 2);
  t43 = Bnl[1][10] * sig(0, 1);
  t44 = Bnl[1][10] * sig(1, 1);
  t45 = Bnl[2][11] * sig(2, 2);
  t46 = Bnl[2][11] * sig(0, 2);
  t47 = Bnl[2][11] * sig(1, 2);
  t48 = Bnl[2][11] * sig(2, 2);
  t49 = Bnl[0][12] * sig(0, 0);
  t50 = Bnl[0][12] * sig(0, 1);
  t51 = Bnl[0][12] * sig(0, 2);
  t52 = Bnl[0][12] * sig(0, 0);
  t53 = Bnl[1][13] * sig(1, 1);
  t54 = Bnl[1][13] * sig(1, 2);
  t55 = Bnl[1][13] * sig(0, 1);
  t56 = Bnl[1][13] * sig(1, 1);
  t57 = Bnl[2][14] * sig(2, 2);
  t58 = Bnl[2][14] * sig(0, 2);
  t59 = Bnl[2][14] * sig(1, 2);
  t60 = Bnl[2][14] * sig(2, 2);
  t61 = Bnl[0][15] * sig(0, 0);
  t62 = Bnl[0][15] * sig(0, 1);
  t63 = Bnl[0][15] * sig(0, 2);
  t64 = Bnl[0][15] * sig(0, 0);
  t65 = Bnl[1][16] * sig(1, 1);
  t66 = Bnl[1][16] * sig(1, 2);
  t67 = Bnl[1][16] * sig(0, 1);
  t68 = Bnl[1][16] * sig(1, 1);
  t69 = Bnl[2][17] * sig(2, 2);
  t70 = Bnl[2][17] * sig(0, 2);
  t71 = Bnl[2][17] * sig(1, 2);
  t72 = Bnl[2][17] * sig(2, 2);
  t73 = Bnl[0][18] * sig(0, 0);
  t74 = Bnl[0][18] * sig(0, 1);
  t75 = Bnl[0][18] * sig(0, 2);
  t77 = Bnl[1][19] * sig(1, 1);
  t78 = Bnl[1][19] * sig(1, 2);
  t79 = Bnl[1][19] * sig(0, 1);
  t81 = Bnl[2][20] * sig(2, 2);
  t82 = Bnl[2][20] * sig(0, 2);
  t83 = Bnl[2][20] * sig(1, 2);
  t85 = Bnl[0][21] * sig(0, 0);
  t86 = Bnl[0][21] * sig(0, 1);
  t87 = Bnl[0][21] * sig(0, 2);
  t88 = Bnl[1][22] * sig(1, 1);
  t89 = Bnl[1][22] * sig(1, 2);
  t90 = Bnl[2][23] * sig(2, 2);

  Kgeo[0][0]   = t1 * Bnl[0][0];
  Kgeo[0][1]   = t2 * Bnl[1][1];
  Kgeo[0][2]   = t3 * Bnl[2][2];
  Kgeo[0][3]   = t4 * Bnl[0][3];
  Kgeo[0][4]   = t2 * Bnl[1][4];
  Kgeo[0][5]   = t3 * Bnl[2][5];
  Kgeo[0][6]   = t4 * Bnl[0][6];
  Kgeo[0][7]   = t2 * Bnl[1][7];
  Kgeo[0][8]   = t3 * Bnl[2][8];
  Kgeo[0][9]   = t4 * Bnl[0][9];
  Kgeo[0][10]  = t2 * Bnl[1][10];
  Kgeo[0][11]  = t3 * Bnl[2][11];
  Kgeo[0][12]  = t4 * Bnl[0][12];
  Kgeo[0][13]  = t2 * Bnl[1][13];
  Kgeo[0][14]  = t3 * Bnl[2][14];
  Kgeo[0][15]  = t4 * Bnl[0][15];
  Kgeo[0][16]  = t2 * Bnl[1][16];
  Kgeo[0][17]  = t3 * Bnl[2][17];
  Kgeo[0][18]  = t4 * Bnl[0][18];
  Kgeo[0][19]  = t2 * Bnl[1][19];
  Kgeo[0][20]  = t3 * Bnl[2][20];
  Kgeo[0][21]  = t4 * Bnl[0][21];
  Kgeo[0][22]  = t2 * Bnl[1][22];
  Kgeo[0][23]  = t3 * Bnl[2][23];
  Kgeo[1][0]   = Kgeo[0][1];
  Kgeo[1][1]   = t5 * Bnl[1][1];
  Kgeo[1][2]   = t6 * Bnl[2][2];
  Kgeo[1][3]   = t7 * Bnl[0][3];
  Kgeo[1][4]   = Bnl[1][4] * t8;
  Kgeo[1][5]   = t6 * Bnl[2][5];
  Kgeo[1][6]   = t7 * Bnl[0][6];
  Kgeo[1][7]   = Bnl[1][7] * t8;
  Kgeo[1][8]   = t6 * Bnl[2][8];
  Kgeo[1][9]   = t7 * Bnl[0][9];
  Kgeo[1][10]  = Bnl[1][10] * t8;
  Kgeo[1][11]  = t6 * Bnl[2][11];
  Kgeo[1][12]  = t7 * Bnl[0][12];
  Kgeo[1][13]  = Bnl[1][13] * t8;
  Kgeo[1][14]  = t6 * Bnl[2][14];
  Kgeo[1][15]  = t7 * Bnl[0][15];
  Kgeo[1][16]  = Bnl[1][16] * t8;
  Kgeo[1][17]  = t6 * Bnl[2][17];
  Kgeo[1][18]  = t7 * Bnl[0][18];
  Kgeo[1][19]  = Bnl[1][19] * t8;
  Kgeo[1][20]  = t6 * Bnl[2][20];
  Kgeo[1][21]  = t7 * Bnl[0][21];
  Kgeo[1][22]  = Bnl[1][22] * t8;
  Kgeo[1][23]  = t6 * Bnl[2][23];
  Kgeo[2][0]   = Kgeo[0][2];
  Kgeo[2][1]   = Kgeo[1][2];
  Kgeo[2][2]   = t9 * Bnl[2][2];
  Kgeo[2][3]   = t10 * Bnl[0][3];
  Kgeo[2][4]   = Bnl[1][4] * t11;
  Kgeo[2][5]   = t12 * Bnl[2][5];
  Kgeo[2][6]   = t10 * Bnl[0][6];
  Kgeo[2][7]   = Bnl[1][7] * t11;
  Kgeo[2][8]   = t12 * Bnl[2][8];
  Kgeo[2][9]   = t10 * Bnl[0][9];
  Kgeo[2][10]  = Bnl[1][10] * t11;
  Kgeo[2][11]  = t12 * Bnl[2][11];
  Kgeo[2][12]  = t10 * Bnl[0][12];
  Kgeo[2][13]  = Bnl[1][13] * t11;
  Kgeo[2][14]  = t12 * Bnl[2][14];
  Kgeo[2][15]  = t10 * Bnl[0][15];
  Kgeo[2][16]  = Bnl[1][16] * t11;
  Kgeo[2][17]  = t12 * Bnl[2][17];
  Kgeo[2][18]  = t10 * Bnl[0][18];
  Kgeo[2][19]  = t11 * Bnl[1][19];
  Kgeo[2][20]  = t12 * Bnl[2][20];
  Kgeo[2][21]  = t10 * Bnl[0][21];
  Kgeo[2][22]  = t11 * Bnl[1][22];
  Kgeo[2][23]  = t12 * Bnl[2][23];
  Kgeo[3][0]   = Kgeo[0][3];
  Kgeo[3][1]   = Kgeo[1][3];
  Kgeo[3][2]   = Kgeo[2][3];
  Kgeo[3][3]   = t13 * Bnl[0][3];
  Kgeo[3][4]   = t14 * Bnl[1][4];
  Kgeo[3][5]   = Bnl[2][5] * t15;
  Kgeo[3][6]   = t16 * Bnl[0][6];
  Kgeo[3][7]   = t14 * Bnl[1][7];
  Kgeo[3][8]   = Bnl[2][8] * t15;
  Kgeo[3][9]   = t16 * Bnl[0][9];
  Kgeo[3][10]  = t14 * Bnl[1][10];
  Kgeo[3][11]  = Bnl[2][11] * t15;
  Kgeo[3][12]  = t16 * Bnl[0][12];
  Kgeo[3][13]  = t14 * Bnl[1][13];
  Kgeo[3][14]  = Bnl[2][14] * t15;
  Kgeo[3][15]  = t16 * Bnl[0][15];
  Kgeo[3][16]  = t14 * Bnl[1][16];
  Kgeo[3][17]  = Bnl[2][17] * t15;
  Kgeo[3][18]  = t16 * Bnl[0][18];
  Kgeo[3][19]  = t14 * Bnl[1][19];
  Kgeo[3][20]  = Bnl[2][20] * t15;
  Kgeo[3][21]  = t16 * Bnl[0][21];
  Kgeo[3][22]  = t14 * Bnl[1][22];
  Kgeo[3][23]  = Bnl[2][23] * t15;
  Kgeo[4][0]   = Kgeo[0][4];
  Kgeo[4][1]   = Kgeo[1][4];
  Kgeo[4][2]   = Kgeo[2][4];
  Kgeo[4][3]   = Kgeo[3][4];
  Kgeo[4][4]   = t17 * Bnl[1][4];
  Kgeo[4][5]   = t18 * Bnl[2][5];
  Kgeo[4][6]   = t19 * Bnl[0][6];
  Kgeo[4][7]   = t20 * Bnl[1][7];
  Kgeo[4][8]   = t18 * Bnl[2][8];
  Kgeo[4][9]   = t19 * Bnl[0][9];
  Kgeo[4][10]  = t20 * Bnl[1][10];
  Kgeo[4][11]  = t18 * Bnl[2][11];
  Kgeo[4][12]  = t19 * Bnl[0][12];
  Kgeo[4][13]  = t20 * Bnl[1][13];
  Kgeo[4][14]  = t18 * Bnl[2][14];
  Kgeo[4][15]  = t19 * Bnl[0][15];
  Kgeo[4][16]  = t20 * Bnl[1][16];
  Kgeo[4][17]  = t18 * Bnl[2][17];
  Kgeo[4][18]  = t19 * Bnl[0][18];
  Kgeo[4][19]  = t20 * Bnl[1][19];
  Kgeo[4][20]  = t18 * Bnl[2][20];
  Kgeo[4][21]  = t19 * Bnl[0][21];
  Kgeo[4][22]  = t20 * Bnl[1][22];
  Kgeo[4][23]  = t18 * Bnl[2][23];
  Kgeo[5][0]   = Kgeo[0][5];
  Kgeo[5][1]   = Kgeo[1][5];
  Kgeo[5][2]   = Kgeo[2][5];
  Kgeo[5][3]   = Kgeo[3][5];
  Kgeo[5][4]   = Kgeo[4][5];
  Kgeo[5][5]   = t21 * Bnl[2][5];
  Kgeo[5][6]   = t22 * Bnl[0][6];
  Kgeo[5][7]   = t23 * Bnl[1][7];
  Kgeo[5][8]   = t24 * Bnl[2][8];
  Kgeo[5][9]   = t22 * Bnl[0][9];
  Kgeo[5][10]  = t23 * Bnl[1][10];
  Kgeo[5][11]  = t24 * Bnl[2][11];
  Kgeo[5][12]  = t22 * Bnl[0][12];
  Kgeo[5][13]  = t23 * Bnl[1][13];
  Kgeo[5][14]  = t24 * Bnl[2][14];
  Kgeo[5][15]  = t22 * Bnl[0][15];
  Kgeo[5][16]  = t23 * Bnl[1][16];
  Kgeo[5][17]  = t24 * Bnl[2][17];
  Kgeo[5][18]  = t22 * Bnl[0][18];
  Kgeo[5][19]  = t23 * Bnl[1][19];
  Kgeo[5][20]  = t24 * Bnl[2][20];
  Kgeo[5][21]  = t22 * Bnl[0][21];
  Kgeo[5][22]  = t23 * Bnl[1][22];
  Kgeo[5][23]  = t24 * Bnl[2][23];
  Kgeo[6][0]   = Kgeo[0][6];
  Kgeo[6][1]   = Kgeo[1][6];
  Kgeo[6][2]   = Kgeo[2][6];
  Kgeo[6][3]   = Kgeo[3][6];
  Kgeo[6][4]   = Kgeo[4][6];
  Kgeo[6][5]   = Kgeo[5][6];
  Kgeo[6][6]   = t25 * Bnl[0][6];
  Kgeo[6][7]   = t26 * Bnl[1][7];
  Kgeo[6][8]   = t27 * Bnl[2][8];
  Kgeo[6][9]   = t28 * Bnl[0][9];
  Kgeo[6][10]  = t26 * Bnl[1][10];
  Kgeo[6][11]  = t27 * Bnl[2][11];
  Kgeo[6][12]  = t28 * Bnl[0][12];
  Kgeo[6][13]  = t26 * Bnl[1][13];
  Kgeo[6][14]  = t27 * Bnl[2][14];
  Kgeo[6][15]  = t28 * Bnl[0][15];
  Kgeo[6][16]  = t26 * Bnl[1][16];
  Kgeo[6][17]  = t27 * Bnl[2][17];
  Kgeo[6][18]  = t28 * Bnl[0][18];
  Kgeo[6][19]  = t26 * Bnl[1][19];
  Kgeo[6][20]  = t27 * Bnl[2][20];
  Kgeo[6][21]  = t28 * Bnl[0][21];
  Kgeo[6][22]  = t26 * Bnl[1][22];
  Kgeo[6][23]  = t27 * Bnl[2][23];
  Kgeo[7][0]   = Kgeo[0][7];
  Kgeo[7][1]   = Kgeo[1][7];
  Kgeo[7][2]   = Kgeo[2][7];
  Kgeo[7][3]   = Kgeo[3][7];
  Kgeo[7][4]   = Kgeo[4][7];
  Kgeo[7][5]   = Kgeo[5][7];
  Kgeo[7][6]   = Kgeo[6][7];
  Kgeo[7][7]   = t29 * Bnl[1][7];
  Kgeo[7][8]   = t30 * Bnl[2][8];
  Kgeo[7][9]   = t31 * Bnl[0][9];
  Kgeo[7][10]  = t32 * Bnl[1][10];
  Kgeo[7][11]  = t30 * Bnl[2][11];
  Kgeo[7][12]  = t31 * Bnl[0][12];
  Kgeo[7][13]  = t32 * Bnl[1][13];
  Kgeo[7][14]  = t30 * Bnl[2][14];
  Kgeo[7][15]  = t31 * Bnl[0][15];
  Kgeo[7][16]  = t32 * Bnl[1][16];
  Kgeo[7][17]  = t30 * Bnl[2][17];
  Kgeo[7][18]  = t31 * Bnl[0][18];
  Kgeo[7][19]  = t32 * Bnl[1][19];
  Kgeo[7][20]  = t30 * Bnl[2][20];
  Kgeo[7][21]  = t31 * Bnl[0][21];
  Kgeo[7][22]  = t32 * Bnl[1][22];
  Kgeo[7][23]  = t30 * Bnl[2][23];
  Kgeo[8][0]   = Kgeo[0][8];
  Kgeo[8][1]   = Kgeo[1][8];
  Kgeo[8][2]   = Kgeo[2][8];
  Kgeo[8][3]   = Kgeo[3][8];
  Kgeo[8][4]   = Kgeo[4][8];
  Kgeo[8][5]   = Kgeo[5][8];
  Kgeo[8][6]   = Kgeo[6][8];
  Kgeo[8][7]   = Kgeo[7][8];
  Kgeo[8][8]   = t33 * Bnl[2][8];
  Kgeo[8][9]   = t34 * Bnl[0][9];
  Kgeo[8][10]  = t35 * Bnl[1][10];
  Kgeo[8][11]  = t36 * Bnl[2][11];
  Kgeo[8][12]  = t34 * Bnl[0][12];
  Kgeo[8][13]  = t35 * Bnl[1][13];
  Kgeo[8][14]  = t36 * Bnl[2][14];
  Kgeo[8][15]  = t34 * Bnl[0][15];
  Kgeo[8][16]  = t35 * Bnl[1][16];
  Kgeo[8][17]  = t36 * Bnl[2][17];
  Kgeo[8][18]  = t34 * Bnl[0][18];
  Kgeo[8][19]  = t35 * Bnl[1][19];
  Kgeo[8][20]  = t36 * Bnl[2][20];
  Kgeo[8][21]  = t34 * Bnl[0][21];
  Kgeo[8][22]  = t35 * Bnl[1][22];
  Kgeo[8][23]  = t36 * Bnl[2][23];
  Kgeo[9][0]   = Kgeo[0][9];
  Kgeo[9][1]   = Kgeo[1][9];
  Kgeo[9][2]   = Kgeo[2][9];
  Kgeo[9][3]   = Kgeo[3][9];
  Kgeo[9][4]   = Kgeo[4][9];
  Kgeo[9][5]   = Kgeo[5][9];
  Kgeo[9][6]   = Kgeo[6][9];
  Kgeo[9][7]   = Kgeo[7][9];
  Kgeo[9][8]   = Kgeo[8][9];
  Kgeo[9][9]   = t37 * Bnl[0][9];
  Kgeo[9][10]  = t38 * Bnl[1][10];
  Kgeo[9][11]  = t39 * Bnl[2][11];
  Kgeo[9][12]  = t40 * Bnl[0][12];
  Kgeo[9][13]  = t38 * Bnl[1][13];
  Kgeo[9][14]  = t39 * Bnl[2][14];
  Kgeo[9][15]  = t40 * Bnl[0][15];
  Kgeo[9][16]  = t38 * Bnl[1][16];
  Kgeo[9][17]  = t39 * Bnl[2][17];
  Kgeo[9][18]  = t40 * Bnl[0][18];
  Kgeo[9][19]  = t38 * Bnl[1][19];
  Kgeo[9][20]  = t39 * Bnl[2][20];
  Kgeo[9][21]  = t40 * Bnl[0][21];
  Kgeo[9][22]  = t38 * Bnl[1][22];
  Kgeo[9][23]  = t39 * Bnl[2][23];
  Kgeo[10][0]  = Kgeo[0][10];
  Kgeo[10][1]  = Kgeo[1][10];
  Kgeo[10][2]  = Kgeo[2][10];
  Kgeo[10][3]  = Kgeo[3][10];
  Kgeo[10][4]  = Kgeo[4][10];
  Kgeo[10][5]  = Kgeo[5][10];
  Kgeo[10][6]  = Kgeo[6][10];
  Kgeo[10][7]  = Kgeo[7][10];
  Kgeo[10][8]  = Kgeo[8][10];
  Kgeo[10][9]  = Kgeo[9][10];
  Kgeo[10][10] = t41 * Bnl[1][10];
  Kgeo[10][11] = t42 * Bnl[2][11];
  Kgeo[10][12] = t43 * Bnl[0][12];
  Kgeo[10][13] = t44 * Bnl[1][13];
  Kgeo[10][14] = t42 * Bnl[2][14];
  Kgeo[10][15] = t43 * Bnl[0][15];
  Kgeo[10][16] = t44 * Bnl[1][16];
  Kgeo[10][17] = t42 * Bnl[2][17];
  Kgeo[10][18] = t43 * Bnl[0][18];
  Kgeo[10][19] = t44 * Bnl[1][19];
  Kgeo[10][20] = t42 * Bnl[2][20];
  Kgeo[10][21] = t43 * Bnl[0][21];
  Kgeo[10][22] = t44 * Bnl[1][22];
  Kgeo[10][23] = t42 * Bnl[2][23];
  Kgeo[11][0]  = Kgeo[0][11];
  Kgeo[11][1]  = Kgeo[1][11];
  Kgeo[11][2]  = Kgeo[2][11];
  Kgeo[11][3]  = Kgeo[3][11];
  Kgeo[11][4]  = Kgeo[4][11];
  Kgeo[11][5]  = Kgeo[5][11];
  Kgeo[11][6]  = Kgeo[6][11];
  Kgeo[11][7]  = Kgeo[7][11];
  Kgeo[11][8]  = Kgeo[8][11];
  Kgeo[11][9]  = Kgeo[9][11];
  Kgeo[11][10] = Kgeo[10][11];
  Kgeo[11][11] = t45 * Bnl[2][11];
  Kgeo[11][12] = t46 * Bnl[0][12];
  Kgeo[11][13] = t47 * Bnl[1][13];
  Kgeo[11][14] = t48 * Bnl[2][14];
  Kgeo[11][15] = t46 * Bnl[0][15];
  Kgeo[11][16] = t47 * Bnl[1][16];
  Kgeo[11][17] = t48 * Bnl[2][17];
  Kgeo[11][18] = t46 * Bnl[0][18];
  Kgeo[11][19] = t47 * Bnl[1][19];
  Kgeo[11][20] = t48 * Bnl[2][20];
  Kgeo[11][21] = t46 * Bnl[0][21];
  Kgeo[11][22] = t47 * Bnl[1][22];
  Kgeo[11][23] = t48 * Bnl[2][23];
  Kgeo[12][0]  = Kgeo[0][12];
  Kgeo[12][1]  = Kgeo[1][12];
  Kgeo[12][2]  = Kgeo[2][12];
  Kgeo[12][3]  = Kgeo[3][12];
  Kgeo[12][4]  = Kgeo[4][12];
  Kgeo[12][5]  = Kgeo[5][12];
  Kgeo[12][6]  = Kgeo[6][12];
  Kgeo[12][7]  = Kgeo[7][12];
  Kgeo[12][8]  = Kgeo[8][12];
  Kgeo[12][9]  = Kgeo[9][12];
  Kgeo[12][10] = Kgeo[10][12];
  Kgeo[12][11] = Kgeo[11][12];
  Kgeo[12][12] = t49 * Bnl[0][12];
  Kgeo[12][13] = t50 * Bnl[1][13];
  Kgeo[12][14] = t51 * Bnl[2][14];
  Kgeo[12][15] = t52 * Bnl[0][15];
  Kgeo[12][16] = t50 * Bnl[1][16];
  Kgeo[12][17] = t51 * Bnl[2][17];
  Kgeo[12][18] = t52 * Bnl[0][18];
  Kgeo[12][19] = t50 * Bnl[1][19];
  Kgeo[12][20] = t51 * Bnl[2][20];
  Kgeo[12][21] = t52 * Bnl[0][21];
  Kgeo[12][22] = t50 * Bnl[1][22];
  Kgeo[12][23] = t51 * Bnl[2][23];
  Kgeo[13][0]  = Kgeo[0][13];
  Kgeo[13][1]  = Kgeo[1][13];
  Kgeo[13][2]  = Kgeo[2][13];
  Kgeo[13][3]  = Kgeo[3][13];
  Kgeo[13][4]  = Kgeo[4][13];
  Kgeo[13][5]  = Kgeo[5][13];
  Kgeo[13][6]  = Kgeo[6][13];
  Kgeo[13][7]  = Kgeo[7][13];
  Kgeo[13][8]  = Kgeo[8][13];
  Kgeo[13][9]  = Kgeo[9][13];
  Kgeo[13][10] = Kgeo[10][13];
  Kgeo[13][11] = Kgeo[11][13];
  Kgeo[13][12] = Kgeo[12][13];
  Kgeo[13][13] = t53 * Bnl[1][13];
  Kgeo[13][14] = t54 * Bnl[2][14];
  Kgeo[13][15] = t55 * Bnl[0][15];
  Kgeo[13][16] = t56 * Bnl[1][16];
  Kgeo[13][17] = t54 * Bnl[2][17];
  Kgeo[13][18] = t55 * Bnl[0][18];
  Kgeo[13][19] = t56 * Bnl[1][19];
  Kgeo[13][20] = t54 * Bnl[2][20];
  Kgeo[13][21] = t55 * Bnl[0][21];
  Kgeo[13][22] = t56 * Bnl[1][22];
  Kgeo[13][23] = t54 * Bnl[2][23];
  Kgeo[14][0]  = Kgeo[0][14];
  Kgeo[14][1]  = Kgeo[1][14];
  Kgeo[14][2]  = Kgeo[2][14];
  Kgeo[14][3]  = Kgeo[3][14];
  Kgeo[14][4]  = Kgeo[4][14];
  Kgeo[14][5]  = Kgeo[5][14];
  Kgeo[14][6]  = Kgeo[6][14];
  Kgeo[14][7]  = Kgeo[7][14];
  Kgeo[14][8]  = Kgeo[8][14];
  Kgeo[14][9]  = Kgeo[9][14];
  Kgeo[14][10] = Kgeo[10][14];
  Kgeo[14][11] = Kgeo[11][14];
  Kgeo[14][12] = Kgeo[12][14];
  Kgeo[14][13] = Kgeo[13][14];
  Kgeo[14][14] = t57 * Bnl[2][14];
  Kgeo[14][15] = t58 * Bnl[0][15];
  Kgeo[14][16] = t59 * Bnl[1][16];
  Kgeo[14][17] = t60 * Bnl[2][17];
  Kgeo[14][18] = t58 * Bnl[0][18];
  Kgeo[14][19] = t59 * Bnl[1][19];
  Kgeo[14][20] = t60 * Bnl[2][20];
  Kgeo[14][21] = t58 * Bnl[0][21];
  Kgeo[14][22] = t59 * Bnl[1][22];
  Kgeo[14][23] = t60 * Bnl[2][23];
  Kgeo[15][0]  = Kgeo[0][15];
  Kgeo[15][1]  = Kgeo[1][15];
  Kgeo[15][2]  = Kgeo[2][15];
  Kgeo[15][3]  = Kgeo[3][15];
  Kgeo[15][4]  = Kgeo[4][15];
  Kgeo[15][5]  = Kgeo[5][15];
  Kgeo[15][6]  = Kgeo[6][15];
  Kgeo[15][7]  = Kgeo[7][15];
  Kgeo[15][8]  = Kgeo[8][15];
  Kgeo[15][9]  = Kgeo[9][15];
  Kgeo[15][10] = Kgeo[10][15];
  Kgeo[15][11] = Kgeo[11][15];
  Kgeo[15][12] = Kgeo[12][15];
  Kgeo[15][13] = Kgeo[13][15];
  Kgeo[15][14] = Kgeo[14][15];
  Kgeo[15][15] = t61 * Bnl[0][15];
  Kgeo[15][16] = t62 * Bnl[1][16];
  Kgeo[15][17] = t63 * Bnl[2][17];
  Kgeo[15][18] = t64 * Bnl[0][18];
  Kgeo[15][19] = t62 * Bnl[1][19];
  Kgeo[15][20] = t63 * Bnl[2][20];
  Kgeo[15][21] = t64 * Bnl[0][21];
  Kgeo[15][22] = t62 * Bnl[1][22];
  Kgeo[15][23] = t63 * Bnl[2][23];
  Kgeo[16][0]  = Kgeo[0][16];
  Kgeo[16][1]  = Kgeo[1][16];
  Kgeo[16][2]  = Kgeo[2][16];
  Kgeo[16][3]  = Kgeo[3][16];
  Kgeo[16][4]  = Kgeo[4][16];
  Kgeo[16][5]  = Kgeo[5][16];
  Kgeo[16][6]  = Kgeo[6][16];
  Kgeo[16][7]  = Kgeo[7][16];
  Kgeo[16][8]  = Kgeo[8][16];
  Kgeo[16][9]  = Kgeo[9][16];
  Kgeo[16][10] = Kgeo[10][16];
  Kgeo[16][11] = Kgeo[11][16];
  Kgeo[16][12] = Kgeo[12][16];
  Kgeo[16][13] = Kgeo[13][16];
  Kgeo[16][14] = Kgeo[14][16];
  Kgeo[16][15] = Kgeo[15][16];
  Kgeo[16][16] = t65 * Bnl[1][16];
  Kgeo[16][17] = t66 * Bnl[2][17];
  Kgeo[16][18] = t67 * Bnl[0][18];
  Kgeo[16][19] = t68 * Bnl[1][19];
  Kgeo[16][20] = t66 * Bnl[2][20];
  Kgeo[16][21] = t67 * Bnl[0][21];
  Kgeo[16][22] = t68 * Bnl[1][22];
  Kgeo[16][23] = t66 * Bnl[2][23];
  Kgeo[17][0]  = Kgeo[0][17];
  Kgeo[17][1]  = Kgeo[1][17];
  Kgeo[17][2]  = Kgeo[2][17];
  Kgeo[17][3]  = Kgeo[3][17];
  Kgeo[17][4]  = Kgeo[4][17];
  Kgeo[17][5]  = Kgeo[5][17];
  Kgeo[17][6]  = Kgeo[6][17];
  Kgeo[17][7]  = Kgeo[7][17];
  Kgeo[17][8]  = Kgeo[8][17];
  Kgeo[17][9]  = Kgeo[9][17];
  Kgeo[17][10] = Kgeo[10][17];
  Kgeo[17][11] = Kgeo[11][17];
  Kgeo[17][12] = Kgeo[12][17];
  Kgeo[17][13] = Kgeo[13][17];
  Kgeo[17][14] = Kgeo[14][17];
  Kgeo[17][15] = Kgeo[15][17];
  Kgeo[17][16] = Kgeo[16][17];
  Kgeo[17][17] = t69 * Bnl[2][17];
  Kgeo[17][18] = t70 * Bnl[0][18];
  Kgeo[17][19] = t71 * Bnl[1][19];
  Kgeo[17][20] = t72 * Bnl[2][20];
  Kgeo[17][21] = t70 * Bnl[0][21];
  Kgeo[17][22] = t71 * Bnl[1][22];
  Kgeo[17][23] = t72 * Bnl[2][23];
  Kgeo[18][0]  = Kgeo[0][18];
  Kgeo[18][1]  = Kgeo[1][18];
  Kgeo[18][2]  = Kgeo[2][18];
  Kgeo[18][3]  = Kgeo[3][18];
  Kgeo[18][4]  = Kgeo[4][18];
  Kgeo[18][5]  = Kgeo[5][18];
  Kgeo[18][6]  = Kgeo[6][18];
  Kgeo[18][7]  = Kgeo[7][18];
  Kgeo[18][8]  = Kgeo[8][18];
  Kgeo[18][9]  = Kgeo[9][18];
  Kgeo[18][10] = Kgeo[10][18];
  Kgeo[18][11] = Kgeo[11][18];
  Kgeo[18][12] = Kgeo[12][18];
  Kgeo[18][13] = Kgeo[13][18];
  Kgeo[18][14] = Kgeo[14][18];
  Kgeo[18][15] = Kgeo[15][18];
  Kgeo[18][16] = Kgeo[16][18];
  Kgeo[18][17] = Kgeo[17][18];
  Kgeo[18][18] = t73 * Bnl[0][18];
  Kgeo[18][19] = t74 * Bnl[1][19];
  Kgeo[18][20] = t75 * Bnl[2][20];
  Kgeo[18][21] = t73 * Bnl[0][21];
  Kgeo[18][22] = t74 * Bnl[1][22];
  Kgeo[18][23] = t75 * Bnl[2][23];
  Kgeo[19][0]  = Kgeo[0][19];
  Kgeo[19][1]  = Kgeo[1][19];
  Kgeo[19][2]  = Kgeo[2][19];
  Kgeo[19][3]  = Kgeo[3][19];
  Kgeo[19][4]  = Kgeo[4][19];
  Kgeo[19][5]  = Kgeo[5][19];
  Kgeo[19][6]  = Kgeo[6][19];
  Kgeo[19][7]  = Kgeo[7][19];
  Kgeo[19][8]  = Kgeo[8][19];
  Kgeo[19][9]  = Kgeo[9][19];
  Kgeo[19][10] = Kgeo[10][19];
  Kgeo[19][11] = Kgeo[11][19];
  Kgeo[19][12] = Kgeo[12][19];
  Kgeo[19][13] = Kgeo[13][19];
  Kgeo[19][14] = Kgeo[14][19];
  Kgeo[19][15] = Kgeo[15][19];
  Kgeo[19][16] = Kgeo[16][19];
  Kgeo[19][17] = Kgeo[17][19];
  Kgeo[19][18] = Kgeo[18][19];
  Kgeo[19][19] = t77 * Bnl[1][19];
  Kgeo[19][20] = t78 * Bnl[2][20];
  Kgeo[19][21] = t79 * Bnl[0][21];
  Kgeo[19][22] = t77 * Bnl[1][22];
  Kgeo[19][23] = t78 * Bnl[2][23];
  Kgeo[20][0]  = Kgeo[0][20];
  Kgeo[20][1]  = Kgeo[1][20];
  Kgeo[20][2]  = Kgeo[2][20];
  Kgeo[20][3]  = Kgeo[3][20];
  Kgeo[20][4]  = Kgeo[4][20];
  Kgeo[20][5]  = Kgeo[5][20];
  Kgeo[20][6]  = Kgeo[6][20];
  Kgeo[20][7]  = Kgeo[7][20];
  Kgeo[20][8]  = Kgeo[8][20];
  Kgeo[20][9]  = Kgeo[9][20];
  Kgeo[20][10] = Kgeo[10][20];
  Kgeo[20][11] = Kgeo[11][20];
  Kgeo[20][12] = Kgeo[12][20];
  Kgeo[20][13] = Kgeo[13][20];
  Kgeo[20][14] = Kgeo[14][20];
  Kgeo[20][15] = Kgeo[15][20];
  Kgeo[20][16] = Kgeo[16][20];
  Kgeo[20][17] = Kgeo[17][20];
  Kgeo[20][18] = Kgeo[18][20];
  Kgeo[20][19] = Kgeo[19][20];
  Kgeo[20][20] = t81 * Bnl[2][20];
  Kgeo[20][21] = t82 * Bnl[0][21];
  Kgeo[20][22] = t83 * Bnl[1][22];
  Kgeo[20][23] = t81 * Bnl[2][23];
  Kgeo[21][0]  = Kgeo[0][21];
  Kgeo[21][1]  = Kgeo[1][21];
  Kgeo[21][2]  = Kgeo[2][21];
  Kgeo[21][3]  = Kgeo[3][21];
  Kgeo[21][4]  = Kgeo[4][21];
  Kgeo[21][5]  = Kgeo[5][21];
  Kgeo[21][6]  = Kgeo[6][21];
  Kgeo[21][7]  = Kgeo[7][21];
  Kgeo[21][8]  = Kgeo[8][21];
  Kgeo[21][9]  = Kgeo[9][21];
  Kgeo[21][10] = Kgeo[10][21];
  Kgeo[21][11] = Kgeo[11][21];
  Kgeo[21][12] = Kgeo[12][21];
  Kgeo[21][13] = Kgeo[13][21];
  Kgeo[21][14] = Kgeo[14][21];
  Kgeo[21][15] = Kgeo[15][21];
  Kgeo[21][16] = Kgeo[16][21];
  Kgeo[21][17] = Kgeo[17][21];
  Kgeo[21][18] = Kgeo[18][21];
  Kgeo[21][19] = Kgeo[19][21];
  Kgeo[21][20] = Kgeo[20][21];
  Kgeo[21][21] = t85 * Bnl[0][21];
  Kgeo[21][22] = t86 * Bnl[1][22];
  Kgeo[21][23] = t87 * Bnl[2][23];
  Kgeo[22][0]  = Kgeo[0][22];
  Kgeo[22][1]  = Kgeo[1][22];
  Kgeo[22][2]  = Kgeo[2][22];
  Kgeo[22][3]  = Kgeo[3][22];
  Kgeo[22][4]  = Kgeo[4][22];
  Kgeo[22][5]  = Kgeo[5][22];
  Kgeo[22][6]  = Kgeo[6][22];
  Kgeo[22][7]  = Kgeo[7][22];
  Kgeo[22][8]  = Kgeo[8][22];
  Kgeo[22][9]  = Kgeo[9][22];
  Kgeo[22][10] = Kgeo[10][22];
  Kgeo[22][11] = Kgeo[11][22];
  Kgeo[22][12] = Kgeo[12][22];
  Kgeo[22][13] = Kgeo[13][22];
  Kgeo[22][14] = Kgeo[14][22];
  Kgeo[22][15] = Kgeo[15][22];
  Kgeo[22][16] = Kgeo[16][22];
  Kgeo[22][17] = Kgeo[17][22];
  Kgeo[22][18] = Kgeo[18][22];
  Kgeo[22][19] = Kgeo[19][22];
  Kgeo[22][20] = Kgeo[20][22];
  Kgeo[22][21] = Kgeo[21][22];
  Kgeo[22][22] = t88 * Bnl[1][22];
  Kgeo[22][23] = t89 * Bnl[2][23];
  Kgeo[23][0]  = Kgeo[0][23];
  Kgeo[23][1]  = Kgeo[1][23];
  Kgeo[23][2]  = Kgeo[2][23];
  Kgeo[23][3]  = Kgeo[3][23];
  Kgeo[23][4]  = Kgeo[4][23];
  Kgeo[23][5]  = Kgeo[5][23];
  Kgeo[23][6]  = Kgeo[6][23];
  Kgeo[23][7]  = Kgeo[7][23];
  Kgeo[23][8]  = Kgeo[8][23];
  Kgeo[23][9]  = Kgeo[9][23];
  Kgeo[23][10] = Kgeo[10][23];
  Kgeo[23][11] = Kgeo[11][23];
  Kgeo[23][12] = Kgeo[12][23];
  Kgeo[23][13] = Kgeo[13][23];
  Kgeo[23][14] = Kgeo[14][23];
  Kgeo[23][15] = Kgeo[15][23];
  Kgeo[23][16] = Kgeo[16][23];
  Kgeo[23][17] = Kgeo[17][23];
  Kgeo[23][18] = Kgeo[18][23];
  Kgeo[23][19] = Kgeo[19][23];
  Kgeo[23][20] = Kgeo[20][23];
  Kgeo[23][21] = Kgeo[21][23];
  Kgeo[23][22] = Kgeo[22][23];
  Kgeo[23][23] = t90 * Bnl[2][23];
}
//______________________________________________________________________
//
void
ElasticPlasticHP::carryForward(const PatchSubset* patches,
                               const MPMMaterial* matl,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch   = patches->get(p);
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    constParticleVariable<Matrix3> pRotation, pPlasticStrain;
    constParticleVariable<double> pEqPlasticStrain, pDamage, pPorosity;
    constParticleVariable<double> pEqStrainRate, pEqPlasticStrainRate;
    constParticleVariable<int> pLocalized;

    old_dw->get(pRotation, pRotationLabel, pset);
    old_dw->get(pEqStrainRate, pEqStrainRateLabel, pset);
    old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);
    old_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel, pset);
    old_dw->get(pEqPlasticStrainRate, pEqPlasticStrainRateLabel, pset);
    old_dw->get(pDamage, pDamageLabel, pset);
    old_dw->get(pPorosity, pPorosityLabel, pset);
    old_dw->get(pLocalized, pLocalizedLabel, pset);

    ParticleVariable<Matrix3> pRotation_new, pPlasticStrain_new;
    ParticleVariable<double> pEqPlasticStrain_new, pDamage_new;
    ParticleVariable<double> pPorosity_new, pEqStrainRate_new,
      pEqPlasticStrainRate_new;
    ParticleVariable<int> pLocalized_new;

    new_dw->allocateAndPut(pRotation_new, pRotationLabel_preReloc, pset);
    new_dw->allocateAndPut(pEqStrainRate_new, pEqStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrain_new, pEqPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pEqPlasticStrainRate_new, pEqPlasticStrainRateLabel_preReloc, pset);
    new_dw->allocateAndPut(pDamage_new, pDamageLabel_preReloc, pset);
    new_dw->allocateAndPut(pPorosity_new, pPorosityLabel_preReloc, pset);
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);

    // Get the plastic strain
    d_flow->getInternalVars(pset, old_dw);
    d_flow->allocateAndPutRigid(pset, new_dw);

    d_flow->getInternalVars(pset, old_dw);
    d_flow->allocateAndPutRigid(pset, new_dw);

    for (int idx : *pset) {
      pRotation_new[idx]            = pRotation[idx];
      pEqStrainRate_new[idx]        = pEqStrainRate[idx];
      pPlasticStrain_new[idx]       = pPlasticStrain[idx];
      pEqPlasticStrain_new[idx]     = pEqPlasticStrain[idx];
      pEqPlasticStrainRate_new[idx] = pEqPlasticStrainRate[idx];
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
//______________________________________________________________________
//
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
  task->requires(Task::NewDW, pEqStrainRateLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pPlasticStrainLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pEqPlasticStrainLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pEqPlasticStrainRateLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pDamageLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pPorosityLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pEnergyLabel_preReloc, matlset, gnone);
  d_flow->allocateCMDataAddRequires(task, matl, patch, lb);
}
//______________________________________________________________________
//
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

  ParticleVariable<Matrix3> pRotation, pPlasticStrain;
  ParticleVariable<double> pEqPlasticStrain, pDamage, pPorosity, pEqStrainRate,
    pEqPlasticStrainRate, pEnergy;
  ParticleVariable<int> pLocalized;

  constParticleVariable<Matrix3> o_Rotation, o_PlasticStrain;
  constParticleVariable<double> o_EqPlasticStrain, o_Damage, o_Porosity;
  constParticleVariable<double> o_EqStrainRate, o_EqPlasticStrainRate, o_Energy;
  constParticleVariable<int> o_Localized;

  new_dw->allocateTemporary(pRotation, addset);
  new_dw->allocateTemporary(pPlasticStrain, addset);
  new_dw->allocateTemporary(pEqPlasticStrain, addset);
  new_dw->allocateTemporary(pEqPlasticStrainRate, addset);
  new_dw->allocateTemporary(pDamage, addset);
  new_dw->allocateTemporary(pEqStrainRate, addset);
  new_dw->allocateTemporary(pLocalized, addset);
  new_dw->allocateTemporary(pPorosity, addset);
  new_dw->allocateTemporary(pEnergy, addset);

  new_dw->get(o_Rotation, pRotationLabel_preReloc, delset);
  new_dw->get(o_EqStrainRate, pEqStrainRateLabel_preReloc, delset);
  new_dw->get(o_PlasticStrain, pPlasticStrainLabel_preReloc, delset);
  new_dw->get(o_EqPlasticStrain, pEqPlasticStrainLabel_preReloc, delset);
  new_dw->get(o_EqPlasticStrainRate, pEqPlasticStrainRateLabel_preReloc, delset);
  new_dw->get(o_Damage, pDamageLabel_preReloc, delset);
  new_dw->get(o_Localized, pLocalizedLabel_preReloc, delset);
  new_dw->get(o_Porosity, pPorosityLabel_preReloc, delset);
  new_dw->get(o_Energy, pEnergyLabel_preReloc, delset);

  n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pRotation[*n]            = o_Rotation[*o];
    pEqStrainRate[*n]        = o_EqStrainRate[*o];
    pPlasticStrain[*n]       = o_PlasticStrain[*o];
    pEqPlasticStrain[*n]     = o_EqPlasticStrain[*o];
    pEqPlasticStrainRate[*n] = o_EqPlasticStrainRate[*o];
    pDamage[*n]              = o_Damage[*o];
    pLocalized[*n]           = o_Localized[*o];
    pPorosity[*n]            = o_Porosity[*o];
    pEnergy[*n]              = o_Energy[*o];
  }

  (*newState)[pRotationLabel]            = pRotation.clone();
  (*newState)[pEqStrainRateLabel]        = pEqStrainRate.clone();
  (*newState)[pPlasticStrainLabel]       = pPlasticStrain.clone();
  (*newState)[pEqPlasticStrainLabel]     = pEqPlasticStrain.clone();
  (*newState)[pEqPlasticStrainRateLabel] = pEqPlasticStrainRate.clone();
  (*newState)[pDamageLabel]              = pDamage.clone();
  (*newState)[pLocalizedLabel]           = pLocalized.clone();
  (*newState)[pPorosityLabel]            = pPorosity.clone();
  (*newState)[pEnergyLabel]              = pEnergy.clone();

  // Initialize the data for the flow model
  d_flow->allocateCMDataAdd(new_dw, addset, newState, delset, old_dw);
}

//______________________________________________________________________
//
void
ElasticPlasticHP::addRequiresDamageParameter(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
}
//__________________________________
//
void
ElasticPlasticHP::getDamageParameter(const Patch* patch,
                                     ParticleVariable<int>& damage,
                                     int dwi,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
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
  Matrix3 eta = D - one * (Dkk / 3.0);

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
//______________________________________________________________________
//
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
  double rho_cur;

  if (d_useModifiedEOS && p_gauge < 0.0) {
    double A = p_ref; // modified EOS
    double n = p_ref / bulk;
    rho_cur  = rho_orig * pow(pressure / A, n);
  } else { // Standard EOS
    double p_g_over_bulk = p_gauge / bulk;
    rho_cur =
      rho_orig * (p_g_over_bulk + sqrt(p_g_over_bulk * p_g_over_bulk + 1.));
  }
  return rho_cur;
}
//______________________________________________________________________
//
void
ElasticPlasticHP::computePressEOSCM(double rho_cur,
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

  if (d_useModifiedEOS && rho_cur < rho_orig) {
    double A                = p_ref; // MODIFIED EOS
    double n                = bulk / p_ref;
    double rho_rat_to_the_n = pow(rho_cur / rho_orig, n);
    pressure                = A * rho_rat_to_the_n;
    dp_drho                 = (bulk / rho_cur) * rho_rat_to_the_n;
    tmp                     = dp_drho; // speed of sound squared
  } else {                             // STANDARD EOS
    double p_g = .5 * bulk * (rho_cur * inv_rho_orig - rho_orig / rho_cur);
    pressure   = p_ref + p_g;
    dp_drho    = .5 * bulk * (rho_orig / (rho_cur * rho_cur) + inv_rho_orig);
    tmp        = bulk / rho_cur; // speed of sound squared
  }
}
//__________________________________
//
double
ElasticPlasticHP::getCompressibility()
{
  return 1.0 / d_initialData.Bulk;
}
//__________________________________
//
void
ElasticPlasticHP::scheduleCheckNeedAddMPMMaterial(Task* task,
                                                  const MPMMaterial* matl,
                                                  const PatchSet*) const
{
  Ghost::GhostType gnone        = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pEqPlasticStrainLabel_preReloc, matlset, gnone);

  task->computes(lb->NeedAddMPMMaterialLabel);
}
//__________________________________
//
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

    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<double> pEqPlasticStrain;
    new_dw->get(pEqPlasticStrain, pEqPlasticStrainLabel_preReloc, pset);

    // Loop thru particles
    for (auto idx : *pset) {
      if (pEqPlasticStrain[idx] > 5.e-2) {
        need_add = -1.;
      }
    }
  }
  new_dw->put(sum_vartype(need_add), lb->NeedAddMPMMaterialLabel);
}
