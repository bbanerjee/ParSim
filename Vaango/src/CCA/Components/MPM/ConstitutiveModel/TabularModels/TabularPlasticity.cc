/*
 * The MIT License
 *
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

// Namespace Vaango::
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularPlasticity.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldConditionFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>

// Namespace Uintah::
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>

#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>

#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>

#include <sci_values.h>

// Namespace std::
#include <cerrno>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#define CHECK_FOR_NAN
//#define CHECK_ISFINITE
//#define CHECK_PLASTIC_RATE
//#define CHECK_FOR_NAN_EXTRA
//#define WRITE_YIELD_SURF
//#define CHECK_INTERNAL_VAR_EVOLUTION
//#define DEBUG_INTERNAL_VAR_EVOLUTION
//#define DEBUG_INTERNAL_VAR_EVOLUTION_COMPUTATION
//#define CHECK_HYDROSTATIC_TENSION
//#define CHECK_TENSION_STATES
//#define CHECK_TENSION_STATES_1
//#define CHECK_DAMAGE_ALGORITHM
//#define CHECK_SUBSTEP
//#define CHECK_TRIAL_STRESS
//#define CHECK_YIELD_SURFACE_NORMAL
//#define CHECK_FLOATING_POINT_OVERFLOW
//#define DEBUG_YIELD_BISECTION_R
//#define CHECK_ELASTIC_STRAIN
//#define CHECK_RETURN_ALIGNMENT

using namespace Vaango;
using Uintah::VarLabel;
using Uintah::Matrix3;
using std::ostringstream;
using std::endl;


// MPM needs three functions to interact with ICE in MPMICE
// 1) p = f(rho) 2) rho = g(p) 3) C = 1/K(rho)
// Because the TabularPlasticity bulk modulus model does not have any closed
// form expressions for these functions, we will need spcial treatment.

// Requires the necessary input parameters CONSTRUCTORS
TabularPlasticity::TabularPlasticity(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* mpmFlags)
  : Uintah::ConstitutiveModel(mpmFlags),
    d_hydrostat(ps)
{
  //std::cout << "Creating tabular plasticity model\n";

  // Bulk and shear modulus models
  d_elastic = Vaango::ElasticModuliModelFactory::create(ps);
  if (!d_elastic) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating ElasticModuliModel."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // The capX variable is only used in tabular plasticity with cap
  // but needed here due to a design defect
  InternalVariableModel* capX = Vaango::InternalVariableModelFactory::create(ps);
  d_capX = dynamic_cast<IntVar_TabularCap*>(capX);
  if (!d_capX) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating InternalVariableModel."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Yield condition model
  d_yield = Vaango::YieldConditionFactory::create(ps, d_capX);
  if (!d_yield) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating YieldConditionModel."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Read the hydrostat for MPMICE calculations
  // (hydrostat is pressure vs total volumetric strain)
  d_hydrostat.setup();

  // Algorithmic parameters
  ps->getWithDefault("yield_surface_radius_scaling_factor",
                     d_cm.yield_scale_fac, 1.0);
  ps->getWithDefault("subcycling_characteristic_number",
                     d_cm.subcycling_characteristic_number,
                     256); // allowable subcycles

  checkInputParameters();
  initializeLocalMPMLabels();
}

void
TabularPlasticity::checkInputParameters()
{
  if (d_cm.subcycling_characteristic_number < 1) {
    ostringstream warn;
    warn << "Subcycling characteristic number should be > 1. Default = 256"
         << endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (d_cm.yield_scale_fac < 1.0 || d_cm.yield_scale_fac > 1.0e6) {
    ostringstream warn;
    warn << "Yield surface scaling factor should be between 1 and 1.0e6. "
            "Default = 1."
         << endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

TabularPlasticity::TabularPlasticity(const TabularPlasticity& cm)
  : ConstitutiveModel(cm)
{
  d_elastic = Vaango::ElasticModuliModelFactory::createCopy(cm.d_elastic);
  d_yield = Vaango::YieldConditionFactory::createCopy(cm.d_yield);
  d_hydrostat = cm.d_hydrostat;

  auto capX = Vaango::InternalVariableModelFactory::createCopy(d_capX);
  d_capX = dynamic_cast<IntVar_TabularCap*>(capX);

  // Yield surface scaling
  d_cm.yield_scale_fac = cm.d_cm.yield_scale_fac;

  // Subcycling
  d_cm.subcycling_characteristic_number =
    cm.d_cm.subcycling_characteristic_number;

  initializeLocalMPMLabels();
}

TabularPlasticity::TabularPlasticity(const TabularPlasticity* cm)
  : TabularPlasticity(*cm)
{
}

// Initialize all labels of the particle variables associated with
// TabularPlasticity.
void
TabularPlasticity::initializeLocalMPMLabels()
{
  pElasticStrainLabel = Uintah::VarLabel::create(
    "p.elasticStrain",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pElasticStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.elasticStrain+",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());

  pElasticVolStrainLabel = VarLabel::create(
    "p.elasticVolStrain", ParticleVariable<double>::getTypeDescription());
  pElasticVolStrainLabel_preReloc = VarLabel::create(
    "p.elasticVolStrain+", ParticleVariable<double>::getTypeDescription());

  pPlasticStrainLabel = Uintah::VarLabel::create(
    "p.plasticStrain",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pPlasticStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.plasticStrain+",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());

  pPlasticCumEqStrainLabel = Uintah::VarLabel::create(
    "p.plasticCumEqStrain",
    Uintah::ParticleVariable<double>::getTypeDescription());
  pPlasticCumEqStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.plasticCumEqStrain+",
    Uintah::ParticleVariable<double>::getTypeDescription());

  pPlasticVolStrainLabel = Uintah::VarLabel::create(
    "p.plasticVolStrain",
    Uintah::ParticleVariable<double>::getTypeDescription());
  pPlasticVolStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.plasticVolStrain+",
    Uintah::ParticleVariable<double>::getTypeDescription());

  pBulkModulusLabel = Uintah::VarLabel::create(
    "p.bulkModulus",
    Uintah::ParticleVariable<double>::getTypeDescription());
  pBulkModulusLabel_preReloc = Uintah::VarLabel::create(
    "p.bulkModulus+",
    Uintah::ParticleVariable<double>::getTypeDescription());

}

// DESTRUCTOR
TabularPlasticity::~TabularPlasticity()
{
  VarLabel::destroy(pElasticStrainLabel);
  VarLabel::destroy(pElasticStrainLabel_preReloc);
  VarLabel::destroy(pElasticVolStrainLabel); // Elastic Volumetric Strain
  VarLabel::destroy(pElasticVolStrainLabel_preReloc);
  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);
  VarLabel::destroy(pPlasticCumEqStrainLabel);
  VarLabel::destroy(pPlasticCumEqStrainLabel_preReloc);
  VarLabel::destroy(pPlasticVolStrainLabel);
  VarLabel::destroy(pPlasticVolStrainLabel_preReloc);
  VarLabel::destroy(pBulkModulusLabel);
  VarLabel::destroy(pBulkModulusLabel_preReloc);

  delete d_yield;
  delete d_elastic;
  delete d_capX;
}

// adds problem specification values to checkpoint data for restart
void
TabularPlasticity::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "tabular_plasticity");
  }

  d_elastic->outputProblemSpec(cm_ps);
  d_yield->outputProblemSpec(cm_ps);

  cm_ps->appendElement("yield_surface_radius_scaling_factor",
                       d_cm.yield_scale_fac);
  cm_ps->appendElement("subcycling_characteristic_number",
                       d_cm.subcycling_characteristic_number);

  d_hydrostat.outputProblemSpec(cm_ps);
}

TabularPlasticity*
TabularPlasticity::clone()
{
  return scinew TabularPlasticity(*this);
}

// When a particle is pushed from patch to patch, carry information needed for
// the particle
void
TabularPlasticity::addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to)
{
  from.push_back(pElasticStrainLabel);
  to.push_back(pElasticStrainLabel_preReloc);

  from.push_back(pElasticVolStrainLabel);
  to.push_back(pElasticVolStrainLabel_preReloc);

  from.push_back(pPlasticStrainLabel);
  to.push_back(pPlasticStrainLabel_preReloc);

  from.push_back(pPlasticCumEqStrainLabel);
  to.push_back(pPlasticCumEqStrainLabel_preReloc);

  from.push_back(pPlasticVolStrainLabel);
  to.push_back(pPlasticVolStrainLabel_preReloc);

  from.push_back(pBulkModulusLabel);
  to.push_back(pBulkModulusLabel_preReloc);

  // Add the particle state for the yield condition model
  d_yield->addParticleState(from, to);

  from.push_back(lb->pRemoveLabel);
  to.push_back(lb->pRemoveLabel_preReloc);

  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);
}

/*!------------------------------------------------------------------------*/
void
TabularPlasticity::addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patch) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  // Other constitutive model and input dependent computes and requires
  task->computes(pElasticStrainLabel, matlset);
  task->computes(pElasticVolStrainLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(pPlasticCumEqStrainLabel, matlset);
  task->computes(pPlasticVolStrainLabel, matlset);
  task->computes(pBulkModulusLabel, matlset);

  // Add yield function variablity computes
  d_yield->addInitialComputesAndRequires(task, matl, patch);
}

/*!------------------------------------------------------------------------*/
void
TabularPlasticity::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw)
{
  // Get the particles in the current patch
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  // Get the particle volume and mass
  constParticleVariable<double> pVolume, pMass;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pMass, lb->pMassLabel, pset);

  // Initialize variables for yield function parameter variability
  d_yield->initializeLocalVariables(patch, pset, new_dw, pVolume);

  // Now initialize the other variables
  ParticleVariable<double> pdTdt;
  ParticleVariable<double> pElasticVolStrain;
  ParticleVariable<double> pPlasticCumEqStrain, pPlasticVolStrain;
  ParticleVariable<double> pBulkModulus;
  ParticleVariable<Matrix3> pStress;
  ParticleVariable<Matrix3> pElasticStrain;
  ParticleVariable<Matrix3> pPlasticStrain;

  new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

  new_dw->allocateAndPut(pElasticStrain, pElasticStrainLabel, pset);
  new_dw->allocateAndPut(pElasticVolStrain, pElasticVolStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticCumEqStrain, pPlasticCumEqStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticVolStrain, pPlasticVolStrainLabel, pset);
  new_dw->allocateAndPut(pBulkModulus, pBulkModulusLabel, pset);

  ElasticModuli moduli = d_elastic->getInitialElasticModuli();
  for (const particleIndex& pidx : *pset) {
    pdTdt[pidx] = 0.0;
    pStress[pidx].set(0.0);
    pElasticStrain[pidx].set(0.0);
    pElasticVolStrain[pidx] = 0.0;
    pPlasticStrain[pidx].set(0.0);
    pPlasticCumEqStrain[pidx] = 0.0;
    pPlasticVolStrain[pidx] = 0.0;
    pBulkModulus[pidx] = moduli.bulkModulus;
  }

  // Compute timestep
  computeStableTimestep(patch, matl, new_dw);
}

// Compute stable timestep based on both the particle velocities
// and wave speed
void
TabularPlasticity::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
                             DataWarehouse* new_dw)
{
  int matID = matl->getDWIndex();

  // Compute initial elastic moduli
  ElasticModuli moduli = d_elastic->getInitialElasticModuli();
  double bulk = moduli.bulkModulus;
  double shear = moduli.shearModulus;

  // Initialize wave speed
  double c_dil = std::numeric_limits<double>::min();
  Vector dx = patch->dCell();
  Vector WaveSpeed(c_dil, c_dil, c_dil);

  // Get the particles in the current patch
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

  // Get particles mass, volume, and velocity
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<long64> pParticleID;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pParticleID, lb->pParticleIDLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  // loop over the particles in the patch
  for (const particleIndex& idx : *pset) {

    // Compute wave speed + particle velocity at each particle,
    // store the maximum
    c_dil =
      std::sqrt((bulk + Util::four_third * shear) * (pVolume[idx] / pMass[idx]));

    // std::cout << "K = " << bulk << " G = " << shear << " c_dil = " << c_dil
    // << std::endl;
    WaveSpeed =
      Vector(Max(c_dil + std::abs(pVelocity[idx].x()), WaveSpeed.x()),
             Max(c_dil + std::abs(pVelocity[idx].y()), WaveSpeed.y()),
             Max(c_dil + std::abs(pVelocity[idx].z()), WaveSpeed.z()));
  }

  // Compute the stable timestep based on maximum value of
  // "wave speed + particle velocity"
  WaveSpeed = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

/**
 * Added computes/requires for computeStressTensor
 */
void
TabularPlasticity::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pXLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pMassLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pVolumeLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pTemperatureLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, lb->pVelocityLabel, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pVolumeLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pDeformRateMidLabel, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pStressUnrotatedLabel, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pDefGradLabel_preReloc, matlset, Ghost::None);

  task->computes(lb->pStressLabel_preReloc, matlset);
  task->computes(lb->pdTdtLabel_preReloc, matlset);
  task->modifies(lb->pRemoveLabel_preReloc, matlset);

  // Add yield Function computes and requires
  d_yield->addComputesAndRequires(task, matl, patches);

  // Add internal variable computes and requires
  task->requires(Task::OldDW, pElasticStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pElasticVolStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticCumEqStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticVolStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pBulkModulusLabel, matlset, Ghost::None);

  task->computes(pElasticStrainLabel_preReloc, matlset);
  task->computes(pElasticVolStrainLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(pPlasticCumEqStrainLabel_preReloc, matlset);
  task->computes(pPlasticVolStrainLabel_preReloc, matlset);
  task->computes(pBulkModulusLabel_preReloc, matlset);
}

/**
 *  TabularPlasticity::computeStressTensor
 *  is the core of the TabularPlasticity model which computes
 *  the updated stress at the end of the current timestep along with all other
 *  required data such plastic strain and elastic strain
 */
void
TabularPlasticity::computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  // Global loop over each patch
  for (int p = 0; p < patches->size(); p++) {

    // Declare and initial value assignment for some variables
    const Patch* patch = patches->get(p);

    // Initialize wave speed
    double c_dil = std::numeric_limits<double>::min();
    Vector WaveSpeed(c_dil, c_dil, c_dil);
    Vector dx = patch->dCell();

    // Initialize strain energy
    double se = 0.0;

    // Get particle subset for the current patch
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Set up local particle variables to be read and written
    constParticleVariable<double> pEev, pEpv, pEpeq_old, pBulkModulus_old;
    constParticleVariable<Matrix3> pEe, pEp;
    old_dw->get(pEe,               pElasticStrainLabel,      pset);
    old_dw->get(pEev,              pElasticVolStrainLabel,   pset);
    old_dw->get(pEp,               pPlasticStrainLabel,      pset);
    old_dw->get(pEpv,              pPlasticVolStrainLabel,   pset);
    old_dw->get(pEpeq_old,         pPlasticCumEqStrainLabel, pset);
    old_dw->get(pBulkModulus_old,  pBulkModulusLabel,        pset);

    ParticleVariable<int> pRemove_new;
    new_dw->getModifiable(pRemove_new, lb->pRemoveLabel_preReloc, pset);

    ParticleVariable<double> pEev_new, pEpv_new, pEpeq_new, pBulkModulus_new;
    ParticleVariable<Matrix3> pEe_new, pEp_new;
    new_dw->allocateAndPut(pEe_new,     pElasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEev_new,    pElasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEp_new,     pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpv_new,    pPlasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpeq_new,   pPlasticCumEqStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pBulkModulus_new, pBulkModulusLabel_preReloc, pset);

    // Set up global particle variables to be read and written
    delt_vartype delT;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass, pVolume;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pStress_old;

    old_dw->get(delT,           lb->delTLabel, getLevel(patches));
    old_dw->get(pParticleID,    lb->pParticleIDLabel, pset);
    old_dw->get(pMass,          lb->pMassLabel, pset);
    new_dw->get(pVolume,        lb->pVolumeLabel_preReloc,  pset);
    old_dw->get(pVelocity,      lb->pVelocityLabel, pset);
    new_dw->get(pStress_old,    lb->pStressUnrotatedLabel, pset);

    // Get the particle variables computed in interpolateToParticlesAndUpdate()
    constParticleVariable<Matrix3> pDefRate_mid, pDefGrad_new;
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel,    pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> p_q, pdTdt;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(p_q,            lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt,          lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new,    lb->pStressLabel_preReloc, pset);

    // Loop over particles
    for (particleIndex& idx : *pset) {

      // No thermal effects
      pdTdt[idx] = 0.0;

      Matrix3 DD = pDefRate_mid[idx];
      Matrix3 sigma_old = pStress_old[idx];

      // Set up model state
      ModelState_Tabular state_old;
      state_old.particleID = pParticleID[idx];
      state_old.stressTensor = sigma_old;
      state_old.elasticStrainTensor = pEe[idx];
      state_old.plasticStrainTensor = pEp[idx];
      state_old.ep_cum_eq = pEpeq_old[idx];
      // std::cout << "state_old.Stress = " << state_old.stressTensor <<
      // std::endl;

      // Compute the elastic moduli at t = t_n
      computeElasticProperties(state_old);
      //std::cout << "State old: " << state_old << std::endl;

      //---------------------------------------------------------
      // Rate-independent plastic step
      // Divides the strain increment into substeps, and calls substep function
      ModelState_Tabular state_new;
      Status status = rateIndependentPlasticUpdate(
        DD, delT, idx, pParticleID[idx], state_old, state_new);

      if (status == Status::SUCCESS) {
        pStress_new[idx] =
          state_new.stressTensor; // unrotated stress at end of step
        pEe_new[idx] =
          state_new.elasticStrainTensor; // elastic strain at end of step
        pEp_new[idx] =
          state_new.plasticStrainTensor; // plastic strain at end of step
        pEpv_new[idx] =
          pEp_new[idx].Trace(); // Plastic volumetric strain at end of step
        pEpeq_new[idx] =
          state_new.ep_cum_eq; // Equivalent plastic strain at end of step

        // Elastic volumetric strain at end of step, compute from updated
        // deformation gradient.
        // H = ln(U) => tr(H) = tr(ln(U)) = ln(det(U)) = ln(sqrt(det(FT)
        // det(F))) = ln J
        pEev_new[idx] =
          log(pDefGrad_new[idx].Determinant()) - pEpv_new[idx];

        pBulkModulus_new[idx] = state_new.bulkModulus;
#ifdef CHECK_ELASTIC_STRAIN
        double pEev_integrated = pEe_new[idx].Trace();
        std::cout << "Elastic volume strain error = "
                  << (pEev_new[idx] - pEev_integrated)
                  << " Integrated = " << pEev_integrated
                  << " Defgrad-based = " << pEev_new[idx]
                  << std::endl;
#endif

      } else {

        // If the updateStressAndInternalVars function can't converge it will
        // return false.
        // This indicates substepping has failed, and the particle will be
        // deleted.
        pRemove_new[idx] = -999;
        std::cout << "** WARNING ** Bad step, deleting particle"
                  << " idx = " << idx << " particleID = " << pParticleID[idx]
                  << ":" << __FILE__ << ":" << __LINE__ << std::endl;

        pStress_new[idx] = pStress_old[idx];
        pEe_new[idx] =
          state_old.elasticStrainTensor; // elastic strain at start of step
        pEp_new[idx] =
          state_old.plasticStrainTensor; // plastic strain at start of step
        pEpv_new[idx] = pEp_new[idx].Trace();
        pEpeq_new[idx] = pEpeq_old[idx];
        pEev_new[idx] = pEe_new[idx].Trace();
        pBulkModulus_new[idx] = pBulkModulus_old[idx];
      }

      // Compute wave speed + particle velocity at each particle, store the
      // maximum
      // std::cout << "State QS new rotated";
      computeElasticProperties(state_new);
      double bulk = state_new.bulkModulus;
      double shear = state_new.shearModulus;
      double rho_cur = pMass[idx] / pVolume[idx];
      c_dil = sqrt((bulk + Util::four_third * shear) / rho_cur);
      // std::cout << "K = " << bulk << " G = " << shear << " c_dil = " << c_dil
      // << std::endl;
      WaveSpeed =
        Vector(Max(c_dil + std::abs(pVelocity[idx].x()), WaveSpeed.x()),
               Max(c_dil + std::abs(pVelocity[idx].y()), WaveSpeed.y()),
               Max(c_dil + std::abs(pVelocity[idx].z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) * Util::one_third;
        double c_bulk = sqrt(bulk / rho_cur);
        p_q[idx] = artificialBulkViscosity(DD.Trace(), c_bulk, rho_cur, dx_ave);

        #ifdef CHECK_ISFINITE
          if (!std::isfinite(p_q[idx])) {
            std::cout << "bulk = " << bulk << " rho_cur = " << rho_cur 
                      << " c_bulk = " << c_bulk 
                      << " D = " << DD << "\n";
          }
        #endif

      } else {
        p_q[idx] = 0.;
      }

      // Compute the averaged stress
      Matrix3 AvgStress = (pStress_new[idx] + pStress_old[idx]) * 0.5;

      //std::cout << "After compute stress : particle = " << idx << "\n"
      //          << " stress = " << pStress_new[idx] << "\n";

      // Compute the strain energy increment associated with the particle
      double e =
        (DD(0, 0) * AvgStress(0, 0) + DD(1, 1) * AvgStress(1, 1) +
         DD(2, 2) * AvgStress(2, 2) +
         2.0 * (DD(0, 1) * AvgStress(0, 1) + DD(0, 2) * AvgStress(0, 2) +
                DD(1, 2) * AvgStress(1, 2))) *
        pVolume[idx] * delT;

      // Accumulate the total strain energy
      // MH! Note the initialization of se needs to be fixed as it is currently
      // reset to 0
      se += e;

    } // End particle set loop

    // Compute the stable timestep based on maximum value of "wave speed +
    // particle velocity"
    WaveSpeed =
      dx / WaveSpeed; // Variable now holds critical timestep (not speed)

    double delT_new = WaveSpeed.minComponent();

    // Put the stable timestep and total strain enrgy
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }
  }
} // -----------------------------------END OF COMPUTE STRESS TENSOR FUNCTION

/**
* Function:
*   rateIndependentPlasticUpdate
*
* Purpose:
*   Divides the strain increment into substeps, and calls substep function
*   All stress values within computeStep are quasistatic.
*/
TabularPlasticity::Status
TabularPlasticity::rateIndependentPlasticUpdate(const Matrix3& D, const double& delT,
                                    particleIndex idx, long64 pParticleID,
                                    const ModelState_Tabular& state_old,
                                    ModelState_Tabular& state_new)
{
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "Rate independent update:" << std::endl;
  std::cout << " D = " << D << " delT = " << delT << std::endl;
  std::cout << "\t State old:" << state_old << std::endl;
#endif

  // Compute the strain increment
  Matrix3 strain_inc = D * delT;
  if (strain_inc.Norm() < 1.0e-30) {
    state_new = state_old;
    return Status::SUCCESS;
  }

  // Compute the trial stress
  Matrix3 stress_trial = computeTrialStress(state_old, strain_inc);

  // Set up a trial state, update the stress invariants, and compute elastic
  // properties
  ModelState_Tabular state_trial(state_old);
  state_trial.stressTensor = stress_trial;
  computeElasticProperties(state_trial);

#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "\t strain_inc = " << strain_inc << std::endl;
  std::cout << "\t State trial:" << state_trial << std::endl;
#endif

  // Determine the number of substeps (nsub) based on the magnitude of
  // the trial stress increment relative to the characteristic dimensions
  // of the yield surface.  Also compare the value of the pressure dependent
  // elastic properties at sigma_old and sigma_trial and adjust nsub if
  // there is a large change to ensure an accurate solution for nonlinear
  // elasticity even with fully elastic loading.
  int nsub = computeStepDivisions(idx, pParticleID, state_old, state_trial);

  // * Upon FAILURE *
  // Delete the particle if the number of substeps is unreasonable
  // Send ParticleDelete Flag to Host Code, Store Inputs to particle data:
  // input values for sigma_new, X_new, Zeta_new, ep_new, along with error flag
  if (nsub < 0) {
    state_new = state_old;
    proc0cout << "Step Failed: Particle idx = " << idx
              << " ID = " << pParticleID << " because nsub = " << nsub
              << std::endl;
    // bool success  = false;
    return Status::UNREASONABLE_SUBSTEPS;
  }

  // Compute a subdivided time step:
  // Loop at least once or until substepping is successful
  const int CHI_MAX = 5; // max allowed subcycle multiplier

  double dt = delT / nsub; // substep time increment

  int chi = 1; // subcycle multiplier
  double tlocal = 0.0;
  Status isSuccess = Status::SUCCESS;

  // Set up the initial states for the substeps
  ModelState_Tabular state_k_old(state_old);
  ModelState_Tabular state_k_new(state_old);
  do {

    //  Call substep function {sigma_new, ep_new, X_new, Zeta_new}
    //    = computeSubstep(D, dt, sigma_substep, ep_substep, X_substep,
    //    Zeta_substep)
    //  Repeat while substeps continue to be successful
    isSuccess = computeSubstep(D, dt, state_k_old, state_k_new);
    if (isSuccess == Status::SUCCESS) {

      tlocal += dt;

#ifdef WRITE_YIELD_SURF
      std::cout << "K = " << state_k_old.bulkModulus << std::endl;
      std::cout << "G = " << state_k_old.shearModulus << std::endl;
#endif

      state_k_old = state_k_new;

#ifdef WRITE_YIELD_SURF
      Matrix3 sig = state_k_new.stressTensor;
      std::cout << "sigma_new = np.array([[" << sig(0, 0) << "," << sig(0, 1)
                << "," << sig(0, 2) << "],[" << sig(1, 0) << "," << sig(1, 1)
                << "," << sig(1, 2) << "],[" << sig(2, 0) << "," << sig(2, 1)
                << "," << sig(2, 2) << "]])" << std::endl;
      std::cout << "plot_stress_state(K, G, sigma_trial, sigma_new, 'b')"
                << std::endl;
#endif

    } else {

      // Substepping has failed. Take tenth the timestep.
      dt /= 10.0;

      // Increase chi to keep track of the number of times the timstep has
      // been reduced
      chi += 1;
      if (chi > CHI_MAX) {
        state_new = state_k_old;
        proc0cout << "Substep failed because chi = " << chi << " > " << CHI_MAX
                  << std::endl;
        return Status::TOO_SMALL_TIMESTEP; // isSuccess = false;
      }

      proc0cout << "**WARNING** Decreasing substep time increment to " << dt
                << " because computeSubstep failed." << std::endl;
    }
#ifdef CHECK_SUBSTEP
    if (tlocal < delT) {
      std::cout << "tlocal = " << tlocal << " delT = " << delT
                << " nsub = " << nsub << std::endl;
    }
#endif
  } while (tlocal < delT);

  state_new = state_k_new;

#ifdef CHECK_INTERNAL_VAR_EVOLUTION
  // if (state_old.particleID == 3377699720593411) {
  std::cout << "rateIndependentPlasticUpdate: "
            << " ep_v_old = " << state_old.ep_v
            << " ep_v_new = " << state_new.ep_v << std::endl;
//}
#endif

  return isSuccess;
}

/**
 * Method: computeElasticProperties
 *
 * Purpose:
 *   Compute the bulk and shear modulus at a given state
 *
 * Side effects:
 *   **WARNING** Also computes stress invariants and plastic strain invariants
 */
void
TabularPlasticity::computeElasticProperties(ModelState_Tabular& state) const
{
  state.updateStressInvariants();
  state.updatePlasticStrainInvariants();
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(&state);
  state.bulkModulus = moduli.bulkModulus;
  state.shearModulus = moduli.shearModulus;
}

/**
 * Method: computeTrialStress
 * Purpose:
 *   Compute the trial stress for some increment in strain assuming linear
 * elasticity
 *   over the step.
 */
Matrix3
TabularPlasticity::computeTrialStress(const ModelState_Tabular& state_old,
                          const Matrix3& strain_inc)
{
  // Compute the trial stress
  Matrix3 stress_old = state_old.stressTensor;
  Matrix3 deps_iso = Util::Identity * (Util::one_third * strain_inc.Trace());
  Matrix3 deps_dev = strain_inc - deps_iso;
  Matrix3 stress_trial = stress_old + deps_iso * (3.0 * state_old.bulkModulus) +
                         deps_dev * (2.0 * state_old.shearModulus);
//#ifdef CHECK_TRIAL_STRESS
#ifdef CHECK_FOR_NAN
  if (std::isnan(stress_trial(0, 0))) {
    std::cout << " stress_old = " << stress_old
              << " stress_trial = " << stress_trial
              << " p_trial = " << stress_trial.Trace() / 3.0
              << " strain_inc = " << strain_inc << " deps_iso = " << deps_iso
              << " deps_dev = " << deps_dev << " K = " << state_old.bulkModulus
              << " G = " << state_old.shearModulus << std::endl;
    throw InternalError("**ERROR** Nan in compute trial stress.", __FILE__,
                        __LINE__);
  }
#endif
  //#endif

  return stress_trial;
}

/**
 * Method: computeStepDivisions
 * Purpose:
 *   Compute the number of step divisions (substeps) based on a comparison
 *   of the trial stress relative to the size of the yield surface, as well
 *   as change in elastic properties between sigma_n and sigma_trial.
 *
 * Caveat:  Uses the mean values of the yield condition parameters.
 */
int
TabularPlasticity::computeStepDivisions(particleIndex idx, long64 particleID,
                            const ModelState_Tabular& state_old,
                            const ModelState_Tabular& state_trial)
{
  // Compute change in bulk modulus:
  double bulk_old = state_old.bulkModulus;
  double bulk_trial = state_trial.bulkModulus;

  int n_bulk =
    std::max(std::ceil(std::abs(bulk_old - bulk_trial) / bulk_old), 1.0);
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "bulk_old = " << bulk_old << " bulk_trial = " << bulk_trial
            << " n_bulk = " << n_bulk << std::endl;
#endif

  // Compute trial stress increment relative to yield surface size:
  Matrix3 d_sigma = state_trial.stressTensor - state_old.stressTensor;
  auto params = d_yield->getParameters();
  double I1_size = 0.5 * (params.at("I1_max") - params.at("I1_min"));
  double J2_size = params.at("sqrtJ2_max");
  double size = std::max(I1_size, J2_size);
  size *= d_cm.yield_scale_fac;
  int n_yield = ceil(d_sigma.Norm() / size);

#ifdef CHECK_FOR_NAN_EXTRA
  proc0cout << "bulk_old = " << bulk_old << " bulk_trial = " << bulk_trial
            << " n_bulk = " << n_bulk << std::endl;
  proc0cout << "I1_max = " << params.at("I1_max") 
            << " I1_min = " << params.at("I1_min") 
            << " sqrtJ2_max = " << params.at("sqrtJ2_max") 
            << " size = " << size
            << " |dsigma| = " << d_sigma.Norm() << " n_yield = " << n_yield
            << std::endl;
#endif

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = std::max(std::max(n_bulk, n_yield), 1);
  int nmax = d_cm.subcycling_characteristic_number;

  if (nsub > nmax) {
    proc0cout << "\n **WARNING** Too many substeps needed for particle "
              << " idx = " << idx << " particle ID = " << particleID
              << std::endl;
    proc0cout << "\t" << __FILE__ << ":" << __LINE__ << std::endl;
    proc0cout << "\t State at t_n: " << state_old;

    proc0cout << "\t Trial state at t_n+1: " << state_trial;

    proc0cout << "\t Ratio of trial bulk modulus to t_n bulk modulus " << n_bulk
              << std::endl;

    proc0cout << "\t ||sig_trial - sigma_n|| " << d_sigma.Norm() << std::endl;
    proc0cout << "\t Yield surface radius in I1-space: " << size << std::endl;
    proc0cout << "\t Ratio of ||sig_trial - sigma_n|| and "
              << d_cm.yield_scale_fac << "*y.s. radius: " << n_yield
              << std::endl;
    proc0cout << "\t I1_max = " << params.at("I1_max") 
              << " I1_min = " << params.at("I1_min") 
              << std::endl;
    proc0cout << "\t I1_size = " << I1_size << " J2_size = " << J2_size
              << std::endl;

    proc0cout << "** BECAUSE** nsub = " << nsub << " > "
              << d_cm.subcycling_characteristic_number
              << " : Probably too much tension in the particle." << std::endl;
    nsub = -1;
  }
  return nsub;
}

/**
 * Method: computeSubstep
 *
 * Purpose:
 *   Computes the updated stress state for a substep that may be either
 *   elastic, plastic, or partially elastic.
 */
TabularPlasticity::Status
TabularPlasticity::computeSubstep(const Matrix3& D, const double& dt,
                      const ModelState_Tabular& state_k_old,
                      ModelState_Tabular& state_k_new)
{
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "\t D = " << D << std::endl;
  std::cout << "\t dt:" << dt << std::endl;
#endif

  // Compute the trial stress
  Matrix3 deltaEps = D * dt;
  Matrix3 stress_k_trial = computeTrialStress(state_k_old, deltaEps);

#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "\t deltaEps = " << deltaEps << std::endl;
  std::cout << "\t Stress k trial:" << stress_k_trial << std::endl;
#endif

#ifdef WRITE_YIELD_SURF
  // std::cout << "Inside computeSubstep:" << std::endl;
  std::cout << "K = " << state_k_old.bulkModulus << std::endl;
  std::cout << "G = " << state_k_old.shearModulus << std::endl;
  Matrix3 sig = stress_k_trial;
  std::cout << "sigma_trial = np.array([[" << sig(0, 0) << "," << sig(0, 1)
            << "," << sig(0, 2) << "],[" << sig(1, 0) << "," << sig(1, 1) << ","
            << sig(1, 2) << "],[" << sig(2, 0) << "," << sig(2, 1) << ","
            << sig(2, 2) << "]])" << std::endl;
  std::cout << "plot_stress_state(K, G, sigma_new, sigma_trial, 'r')"
            << std::endl;
// std::cout << "\t computeSubstep: sigma_old = " << state_k_old.stressTensor
//         << " sigma_trial = " << stress_trial
//         << " D = " << D << " dt = " << dt
//         << " deltaEps = " << deltaEps << std::endl;
#endif

  // Set up a trial state, update the stress invariants
  ModelState_Tabular state_k_trial(state_k_old);
  state_k_trial.stressTensor = stress_k_trial;

  // Compute elastic moduli at trial stress state
  // and update stress invariants
  computeElasticProperties(state_k_trial);

  // Evaluate the yield function at the trial stress:
  auto yield = d_yield->evalYieldCondition(&state_k_trial);

  // std::cout << "Has yielded ? 1 = Yes, -1 = No." << yield << std::endl;
  // std::cout << "computeSubstep:Elastic:sigma_new = " <<
  // state_k_new.stressTensor
  //          << " ep_v_trial = " << state_k_trial.ep_v
  //          << std::endl;

  // Elastic substep
  if (yield.second == Util::YieldStatus::IS_ELASTIC) {
    state_k_new = state_k_trial;
    state_k_new.elasticStrainTensor += deltaEps;

#ifdef CHECK_INTERNAL_VAR_EVOLUTION
    std::cout << "computeSubstep:Elastic:sigma_new = "
              << state_k_new.stressTensor
              << " ep_v_trial = " << state_k_trial.ep_v << std::endl;
#endif

    return Status::SUCCESS; // bool isSuccess = true;
  }

#ifdef DEBUG_YIELD_BISECTION_R
  std::cout << "before_non_hardening_return  = 1" << std::endl;
  std::cout << "I1 = " << state_k_old.I1 << std::endl;
  std::cout << "sqrt_J2 = " << state_k_old.sqrt_J2 << std::endl;
#endif
  // Elastic-plastic or fully-plastic substep
  // Compute non-hardening return to initial yield surface:
  // std::cout << "\t Doing nonHardeningReturn\n";
  Matrix3 sig_fixed(0.0); // final stress state for non-hardening return
  Matrix3 deltaEps_e_fixed(
    0.0); // increment in elastic strain for non-hardening return
  Matrix3 deltaEps_p_fixed(
    0.0); // increment in plastic strain for non-hardening return
  Status isSuccess =
    nonHardeningReturn(deltaEps, state_k_old, state_k_trial, sig_fixed,
                       deltaEps_e_fixed, deltaEps_p_fixed);
  if (isSuccess != Status::SUCCESS) {
    proc0cout << "**WARNING** nonHardeningReturn has failed." << std::endl;
    return isSuccess;
  }

  state_k_new = state_k_trial;
  state_k_new.stressTensor = sig_fixed;
  state_k_new.elasticStrainTensor += deltaEps_e_fixed;
  state_k_new.plasticStrainTensor += deltaEps_p_fixed;
  state_k_new.updateStressInvariants();
  state_k_new.updatePlasticStrainInvariants();

#ifdef DEBUG_INTERNAL_VAR_EVOLUTION
  std::cout << "computeSubstep: "
            << " ep_v_old = " << state_k_old.ep_v
            << " ep_v_new = " << state_k_new.ep_v << std::endl;
#endif

#ifdef DEBUG_YIELD_BISECTION_R
  std::cout << "after_consistency_bisection  = 1" << std::endl;
  std::cout << "I1 = " << state_k_new.I1 << std::endl;
  std::cout << "sqrt_J2 = " << state_k_new.sqrt_J2 << std::endl;
#endif

  return isSuccess;

} //===================================================================

/**
 * Method: nonHardeningReturn
 * Purpose:
 *   Computes a non-hardening return to the yield surface in the meridional
 * profile
 *   (constant Lode angle) based on the current values of the internal state
 * variables
 *   and elastic properties.  Returns the updated stress and  the increment in
 * plastic
 *   strain corresponding to this return.
 *
 *   NOTE: all values of r and z in this function are transformed!
 */
TabularPlasticity::Status
TabularPlasticity::nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                          const ModelState_Tabular& state_k_old,
                          const ModelState_Tabular& state_k_trial,
                          Uintah::Matrix3& sig_fixed,
                          Uintah::Matrix3& elasticStrain_inc_fixed,
                          Uintah::Matrix3& plasticStrain_inc_fixed)
{
  // Compute ratio of bulk and shear moduli
  double K_old = state_k_old.bulkModulus;
  double G_old = state_k_old.shearModulus;
  const double sqrt_K_over_G_old = std::sqrt(1.5 * K_old / G_old);
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << " K_old = " << K_old << " G_old = " << G_old << std::endl;
#endif

  // Save the r and z Lode coordinates for the trial stress state
  double r_trial = state_k_trial.rr;
  double z_trial = state_k_trial.zz;
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << " state_k_trial " << state_k_trial << std::endl;
#endif

  // Compute transformed r coordinates
  double rprime_trial = r_trial * sqrt_K_over_G_old;
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << " z_trial = " << z_trial
            << " r_trial = " << rprime_trial / sqrt_K_over_G_old << std::endl;
#endif

  // Find closest point
  double z_closest = 0.0, rprime_closest = 0.0;
  d_yield->getClosestPoint(&state_k_old, z_trial, rprime_trial,
                           z_closest, rprime_closest);
#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << " z_closest = " << z_closest
            << " r_closest = " << rprime_closest / sqrt_K_over_G_old
            << std::endl;
#endif

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0) * z_closest;
  double sqrtJ2_closest =
    1.0 / (sqrt_K_over_G_old * Util::sqrt_two) * rprime_closest;

#ifdef CHECK_FOR_NAN_EXTRA
  std::cout << "I1_closest = " << I1_closest 
            << " sqrtJ2_closest = " << sqrtJ2_closest << std::endl;
  std::cout << "Trial state = " << state_k_trial << std::endl;
#endif
#ifdef CHECK_HYDROSTATIC_TENSION
  if (I1_closest < 0) {
    std::cout << "I1_closest = " << I1_closest 
              << " sqrtJ2_closest = " << sqrtJ2_closest << std::endl;
    std::cout << "Trial state = " << state_k_trial << std::endl;
  }
#endif

  // Compute new stress
  Matrix3 sig_dev = state_k_trial.deviatoricStressTensor;
  if (state_k_trial.sqrt_J2 > 0.0) {
    sig_fixed = Util::one_third * I1_closest * Util::Identity +
                (sqrtJ2_closest / state_k_trial.sqrt_J2) * sig_dev;
  } else {
    sig_fixed = Util::one_third * I1_closest * Util::Identity + sig_dev;
  }

  // Compute new plastic strain increment
  //  d_ep = d_e - [C]^-1:(sigma_new-sigma_old)
  Matrix3 sig_inc = sig_fixed - state_k_old.stressTensor;
  Matrix3 sig_inc_iso = Util::one_third * sig_inc.Trace() * Util::Identity;
  Matrix3 sig_inc_dev = sig_inc - sig_inc_iso;
  elasticStrain_inc_fixed =
    sig_inc_iso * (Util::one_third / K_old) + sig_inc_dev * (0.5 / G_old);
  plasticStrain_inc_fixed = strain_inc - elasticStrain_inc_fixed;

#ifdef CHECK_ELASTIC_STRAIN
  // std::cout << "Non-hardening:\n"
  //          << "\t Delta sig = " << sig_inc << std::endl
  //          << "\t Delta Eps_e = " << elasticStrain_inc_fixed << std::endl;
  std::cout
    << "press = " << sig_fixed.Trace() / 3.0 << " K = " << K_old << " ev_e = "
    << (state_k_old.elasticStrainTensor + elasticStrain_inc_fixed).Trace()
    << std::endl;
#endif

  // Compute volumetric plastic strain and compare with p3
  Matrix3 eps_p = state_k_old.plasticStrainTensor + plasticStrain_inc_fixed;
  double ep_v = eps_p.Trace();

  if (ep_v < 0.0) {
    if (-ep_v > 10) {
      proc0cout << "**WARNING** Nonhardening return has failed because "
                << " epsbar_p_v > 10 : " << -ep_v << " > 10 " 
                << std::endl;
      proc0cout << " K_old = " << K_old << " G_old = " << G_old << std::endl;
      proc0cout << " state_k_trial " << state_k_trial << std::endl;
      proc0cout << " z_trial = " << z_trial
                << " r_trial = " << rprime_trial / sqrt_K_over_G_old
                << std::endl;
      proc0cout << " z_closest = " << z_closest
                << " r_closest = " << rprime_closest / sqrt_K_over_G_old
                << std::endl;
      proc0cout << "Delta eps = " << strain_inc << std::endl;
      proc0cout << "sig_n = " << state_k_old.stressTensor << std::endl;
      proc0cout << "sig_n+1 = " << sig_fixed << std::endl;
      proc0cout << "Delta sig = " << sig_inc << std::endl;
      proc0cout << "Delta sig_iso = " << sig_inc_iso << std::endl;
      proc0cout << "Delta sig_dev = " << sig_inc_dev << std::endl;
      proc0cout << "Delta eps_e = " << elasticStrain_inc_fixed << std::endl;
      proc0cout << "Delta eps_p = " << plasticStrain_inc_fixed << std::endl;
      proc0cout << "I1_J2_trial = [" << state_k_trial.I1 << " "
                << state_k_trial.sqrt_J2 << "];" << std::endl;
      proc0cout << "I1_J2_closest = [" << I1_closest
                << " " << sqrtJ2_closest << "];" << std::endl;
      proc0cout << "plot([I1 I1_J2_closest(1)],[sqrtJ2 I1_J2_closest(2)],'gx')"
                << ";" << std::endl;
      proc0cout << "plot([I1_J2_trial(1) I1_J2_closest(1)],[I1_J2_trial(2) "
                   "I1_J2_closest(2)],'r-')"
                << ";" << std::endl;

      return Status::TOO_LARGE_PLASTIC_STRAIN; // The plastic volume strain is too large, try again
    }
  }

#ifdef CHECK_PLASTIC_RATE
  ModelState_Tabular state_plastic_rate(state_k_old);
  state_plastic_rate.stressTensor = sig_fixed;
  state_plastic_rate.updateStressInvariants();
  Matrix3 df_dsigma = d_yield->df_dsigma(Util::Identity, &state_plastic_rate);
  df_dsigma /= df_dsigma.Norm();
  double lhs = plasticStrain_inc_fixed.Contract(df_dsigma);
  double rhs = df_dsigma.Contract(df_dsigma);
  double plastic_rate = lhs/rhs;
  if (plastic_rate < 0) {
    std::cout << "Particle = " << state_k_old.particleID << "\n";
    std::cout << "Delta eps = " << strain_inc << std::endl;
    std::cout << "Trial state = " << state_k_trial << std::endl;
    std::cout << "Delta sig = " << sig_inc << std::endl;
    std::cout << "Delta sig_iso = " << sig_inc_iso << std::endl;
    std::cout << "Delta sig_dev = " << sig_inc_dev << std::endl;
    std::cout << "Delta eps_e = " << elasticStrain_inc_fixed << std::endl;
    std::cout << "Delta eps_p = " << plasticStrain_inc_fixed << std::endl;
    std::cout << "df_dsigma = " << df_dsigma << std::endl;
    std::cout << "plastic rate = " << plastic_rate << "\n";
  }
#endif

#ifdef CHECK_YIELD_SURFACE_NORMAL
  std::cout << "Delta eps = " << strain_inc << std::endl;
  std::cout << "Trial state = " << state_k_trial << std::endl;
  std::cout << "Delta sig = " << sig_inc << std::endl;
  std::cout << "Delta sig_iso = " << sig_inc_iso << std::endl;
  std::cout << "Delta sig_dev = " << sig_inc_dev << std::endl;
  std::cout << "Delta eps_e = " << elasticStrain_inc_fixed << std::endl;
  std::cout << "Delta eps_p = " << plasticStrain_inc_fixed << std::endl;

  // Test normal to yield surface
  ModelState_Tabular state_test(state_k_old);
  state_test.stressTensor = sig_fixed;
  state_test.updateStressInvariants();

  Matrix3 df_dsigma = d_yield->df_dsigma(Identity, &state_test);
  df_dsigma /= df_dsigma.Norm();
  std::cout << "df_dsigma = " << df_dsigma << std::endl;
  std::cout << "ratio = [" << plasticStrain_inc_fixed(0, 0) / df_dsigma(0, 0)
            << "," << plasticStrain_inc_fixed(1, 1) / df_dsigma(1, 1) << ","
            << plasticStrain_inc_fixed(2, 2) / df_dsigma(2, 2) << std::endl;

  // Compute CN = C:df_dsigma
  double lambda = state_test.bulkModulus - 2.0 / 3.0 * state_test.shearModulus;
  double mu = state_test.shearModulus;
  Matrix3 CN = Identity * (lambda * df_dsigma.Trace()) + df_dsigma * (2.0 * mu);
  Matrix3 sig_diff = state_k_trial.stressTensor - sig_fixed;
  std::cout << "sig_trial = [" << state_k_trial.stressTensor << "];"
            << std::endl;
  std::cout << "sig_n+1 = [" << sig_fixed << "];" << std::endl;
  std::cout << "sig_trial - sig_n+1 = [" << sig_diff << "];" << std::endl;
  std::cout << "C_df_dsigma = [" << CN << "];" << std::endl;
  std::cout << "sig ratio = [" << sig_diff(0, 0) / CN(0, 0) << " "
            << sig_diff(1, 1) / CN(1, 1) << " " << sig_diff(2, 2) / CN(2, 2)
            << "];" << std::endl;

  // Compute a test stress to check normal
  Matrix3 sig_test = sig_fixed + df_dsigma * sig_diff(0, 0);
  ModelState_Tabular state_sig_test(state_k_old);
  state_sig_test.stressTensor = sig_test;
  state_sig_test.updateStressInvariants();
  std::cout << "I1 = " << state_sig_test.I1 << ";" << std::endl;
  std::cout << "sqrtJ2 = " << state_sig_test.sqrt_J2 << ";" << std::endl;
  std::cout << "I1_J2_trial = [" << state_k_trial.I1 << " "
            << state_k_trial.sqrt_J2 << "];" << std::endl;
  std::cout << "I1_J2_closest = [" << I1_closest
            << " " << sqrtJ2_closest << "];" << std::endl;
  std::cout << "plot([I1 I1_J2_closest(1)],[sqrtJ2 I1_J2_closest(2)],'gx-')"
            << ";" << std::endl;
  std::cout << "plot([I1_J2_trial(1) I1_J2_closest(1)],[I1_J2_trial(2) "
               "I1_J2_closest(2)],'r-')"
            << ";" << std::endl;

  // Check actual location of projected point
  Matrix3 sig_test_actual =
    state_k_trial.stressTensor - CN * (std::abs(sig_diff(0, 0) / CN(0, 0)));
  state_sig_test.stressTensor = sig_test_actual;
  state_sig_test.updateStressInvariants();
  std::cout << "I1 = " << state_sig_test.I1 << ";" << std::endl;
  std::cout << "sqrtJ2 = " << state_sig_test.sqrt_J2 << ";" << std::endl;
  std::cout << "plot([I1 I1_J2_trial(1)],[sqrtJ2 I1_J2_trial(2)],'rx')"
            << ";" << std::endl;
#endif

#ifdef CHECK_FOR_NAN
  if (std::isnan(sig_fixed(0, 0))) {
    std::cout << " K_old = " << K_old << " G_old = " << G_old << std::endl;
    std::cout << " z_trial = " << z_trial
              << " r_trial = " << rprime_trial / sqrt_K_over_G_old << std::endl;
    std::cout << " z_closest = " << z_closest
              << " r_closest = " << rprime_closest / sqrt_K_over_G_old
              << std::endl;
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << std::endl;
    std::cout << "Trial state = " << state_k_trial << std::endl;
    std::cout << "\t\t\t sig_fixed = " << sig_fixed << std::endl;
    std::cout << "\t\t\t I1_closest = " << I1_closest << std::endl;
    std::cout << "\t\t\t sqrtJ2_closest = " << sqrtJ2_closest << std::endl;
    std::cout << "\t\t\t state_k_trial.sqrt_J2 = " << state_k_trial.sqrt_J2
              << std::endl;
    std::cout << "\t\t\t sig_dev = " << sig_dev << std::endl;
    std::cout << "\t\t\t sig_inc = " << sig_inc << std::endl;
    std::cout << "\t\t\t strain_inc = " << strain_inc << std::endl;
    std::cout << "\t\t\t sig_inc_iso = " << sig_inc_iso << std::endl;
    std::cout << "\t\t\t sig_inc_dev = " << sig_inc_dev << std::endl;
    std::cout << "\t\t\t plasticStrain_inc_fixed = " << plasticStrain_inc_fixed
              << std::endl;
  }
#endif

#ifdef CHECK_HYDROSTATIC_TENSION
  if (I1_closest < 0) {
    std::cout << "\t\t\t sig_inc = " << sig_inc << std::endl;
    std::cout << "\t\t\t strain_inc = " << strain_inc << std::endl;
    std::cout << "\t\t\t sig_inc_iso = " << sig_inc_iso << std::endl;
    std::cout << "\t\t\t sig_inc_dev = " << sig_inc_dev << std::endl;
    std::cout << "\t\t\t plasticStrain_inc_fixed = " << plasticStrain_inc_fixed
              << std::endl;
  }
#endif

#ifdef CHECK_RETURN_ALIGNMENT
  ModelState_Tabular state_n(state_k_old);
  ModelState_Tabular state_np1(state_k_old);
  state_np1.stressTensor = sig_fixed;
  state_np1.updateStressInvariants();

  Matrix3 M_n = d_yield->df_dsigma(Util::Identity, &state_n);
  Matrix3 M_np1 = d_yield->df_dsigma(Util::Identity, &state_np1);
  M_n /= M_n.Norm();
  M_np1 /= M_np1.Norm();

  double angle_M_n_np1 = std::abs(M_n.Contract(M_np1) - 1.0);
  if (angle_M_n_np1 > 1.0e-6) {
    std::cout << "M_n = " << M_n << std::endl;
    std::cout << "M_np1 = " << M_np1 << std::endl;

    auto KG_dKdG_n = d_elastic->getElasticModuliAndDerivatives(&state_n);
    auto KG_dKdG_np1 = d_elastic->getElasticModuliAndDerivatives(&state_np1);
    auto KG_n = KG_dKdG_n.first;
    auto KG_np1 = KG_dKdG_np1.first;
    auto dKdG_n = KG_dKdG_n.second;
    auto dKdG_np1 = KG_dKdG_np1.second;
    auto p_n = state_n.I1/3.0;
    auto p_np1 = state_np1.I1/3.0;
    auto S_n = state_n.deviatoricStressTensor;
    auto S_np1 = state_np1.deviatoricStressTensor;
    Matrix3 Z_n = computeZMatrix(KG_n, dKdG_n, p_n, S_n, M_n);
    Matrix3 Z_np1 = computeZMatrix(KG_np1, dKdG_np1, p_np1, S_np1, M_np1);
    std::cout << "Z_n = " << Z_n << std::endl;
    std::cout << "Z_np1 = " << Z_np1 << std::endl;

    // Compute CM = C:M
    double mu_n = KG_n.shearModulus;
    double mu_np1 = KG_np1.shearModulus;
    double lambda_n = KG_n.bulkModulus - 2.0 / 3.0 * mu_n;
    double lambda_np1 = KG_np1.bulkModulus - 2.0 / 3.0 * mu_np1;
    Matrix3 CM_n = Util::Identity * (lambda_n * M_n.Trace()) + M_n * (2.0 * mu_n);
    Matrix3 CM_np1 = Util::Identity * (lambda_np1 * M_np1.Trace()) + M_np1 * (2.0 * mu_np1);
    std::cout << "CM_n = " << CM_n << std::endl;
    std::cout << "CM_np1 = " << CM_np1 << std::endl;

    // Compute P = CM + Z;
    Matrix3 P_n = CM_n + Z_n;
    Matrix3 P_np1 = CM_np1 + Z_np1;
    std::cout << "P_n = " << P_n << std::endl;
    std::cout << "P_np1 = " << P_np1 << std::endl;

    // Compute Gamma (sig_trial - sig_n+1):M_n+1/(P_n:M_n+1)
    Matrix3 sig_trial = state_k_trial.stressTensor;
    Matrix3 sig_new = sig_fixed;
    double Gamma_n = (sig_trial - sig_new).Contract(M_np1)/P_n.Contract(M_np1);
    double Gamma_np1 = (sig_trial - sig_new).Contract(M_np1)/P_np1.Contract(M_np1);
    std::cout << "Gamma_n = " << Gamma_n << " Gamma_np1 = " << Gamma_np1 << "\n";

    // Compute sigma
    Matrix3 sig_np1 = sig_trial - P_n * Gamma_n;
    std::cout << "sig_np1 = " << sig_np1 << std::endl;
    std::cout << "sig_fixed = " << sig_fixed << std::endl;

    // Compute eps_p and ev_v^p
    Matrix3 Eps_p_n = state_n.plasticStrainTensor + M_n * Gamma_n ;
    double Eps_p_v_n = state_n.ep_v + M_n.Trace() * Gamma_n;
    Matrix3 Eps_e_n = strain_inc - Eps_p_n; 
    Matrix3 Eps_p_np1 = state_n.plasticStrainTensor + M_np1 * Gamma_np1 ;
    double Eps_p_v_np1 = state_n.ep_v + M_np1.Trace() * Gamma_np1;
    Matrix3 Eps_e_np1 = strain_inc - Eps_p_np1; 
    std::cout << "Eps_p_n = " << Eps_p_n << "\n"; 
    std::cout << "eps_p = " << eps_p << "\n"; 
    std::cout << "Eps_e_n = " << Eps_e_n << "\n"; 
    std::cout << "Eps_p_np1 = " << Eps_p_np1 << "\n"; 
    std::cout << "Eps_e_np1 = " << Eps_e_np1 << "\n"; 
    std::cout << "Eps_p_v_n = " << Eps_p_v_n 
              << " Eps_p_v_np1 = " << Eps_p_v_np1 << " ep_v = " << ep_v << "\n";
  }
#endif

  return Status::SUCCESS; // isSuccess = true

} //===================================================================


/** computeZMatrix 
  *
  *  Z = -[1/K dK/dev^p p I + 1/G dG/dev^p S] tr(M)
  *
  *  where K = bulk modulus
  *        G = shear modulus
  *        ev^p = volumetric plastic strain
  *        S = deviatoric stress
  *        M = unit normal to yield surface (associated)
  *        I = Identity tensor
  */
Matrix3 
TabularPlasticity::computeZMatrix(const ElasticModuli& moduli,
                                  const ElasticModuli& derivs,
                                  double p,
                                  const Matrix3& S,
                                  const Matrix3& M) const
{
  //auto moduli_and_derivs = d_elastic->getElasticModuliAndDerivatives(&state);
  //auto moduli = moduli_and_derivs.first;
  //auto derivs = moduli_and_derivs.second;
  auto K_term = Util::Identity * (-derivs.bulkModulus/moduli.bulkModulus) * p;
  auto G_term = S * (-derivs.shearModulus/moduli.shearModulus);
  auto Z = (K_term + G_term) * M.Trace();
  return Z;
}

void
TabularPlasticity::addRequiresDamageParameter(Task* task, const MPMMaterial* matl,
                                  const PatchSet*) const
{
  // Require the damage parameter
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, lb->pRemoveLabel_preReloc, matlset,
                 Ghost::None);
}

void
TabularPlasticity::getDamageParameter(const Patch* patch, ParticleVariable<int>& damage,
                          int matID, DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  // Get the damage parameter
  ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, lb->pRemoveLabel_preReloc, pset);

  // Only update the damage variable if it hasn't been modified by a damage
  // model earlier
  for (auto particle : *pset) {
    if (damage[particle] == 0) {
      damage[particle] = pLocalized[particle];
    }
  }
}

void
TabularPlasticity::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                    DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  // Carry forward the data.
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.0), lb->StrainEnergyLabel);
    }
  }
}

// T2D: Throw exception that this is not supported
void
TabularPlasticity::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                              const bool, const bool) const
{
  std::ostringstream out;
  out << "**ERROR** NO Implicit VERSION OF addComputesAndRequires EXISTS YET FOR TabularPlasticity"
      << std::endl;
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

/*!
 * ---------------------------------------------------------------------------------------
 *  This is needed for converting from one material type to another.  The
 * functionality
 *  has been removed from the main Uintah branch.
 *  ---------------------------------------------------------------------------------------
 */
void
TabularPlasticity::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset, DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "**ERROR** Material conversion after failure DOES NOT EXIST YET for tabular plasticity.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

/*---------------------------------------------------------------------------------------
 * MPMICE Hooks
 *---------------------------------------------------------------------------------------*/
double
TabularPlasticity::computeRhoMicroCM(double pressure, const double p_ref,
                         const MPMMaterial* matl, double temperature,
                         double rho_guess)
{
  double rho_0 = matl->getInitialDensity();
  double p_gauge = pressure - p_ref;

  bool status;
  double eps_v_bar, p_bar;
  std::tie(status, eps_v_bar, p_bar) = d_hydrostat.intersect1D(-10, 10, p_gauge, p_gauge);
  double J = std::exp(-eps_v_bar);
  double rho_cur = rho_0 / J;
  if (std::isinf(rho_cur) || rho_cur < 0) {
    std::ostringstream out;
    out << "rho_0 = " << rho_0 << " rho_cur = " << rho_cur << " J = " << J << "\n";
    out << "status = " << status << " eps_v_bar = " << eps_v_bar 
        << " p_bar = " << p_bar << " p_gauge  = " << p_gauge << "\n";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  return rho_cur;
}

void
TabularPlasticity::computePressEOSCM(double rho_cur, double& pressure, double p_ref,
                         double& dp_drho, double& soundSpeedSq,
                         const MPMMaterial* matl, double temperature)
{
  double rho_0 = matl->getInitialDensity();
  double J = rho_0 / rho_cur;
  double eps_v = std::log(J);

  bool status;
  double eps_v_bar, p_bar;
  std::tie(status, eps_v_bar, p_bar) = d_hydrostat.intersect1D(-eps_v, -eps_v, -1.0e6, 1.0e10);
  double p_gauge = p_bar;
  pressure = p_ref + p_gauge;

  double eps_v_bar_hi, p_bar_hi;
  p_bar_hi = p_bar + 1.0;
  std::tie(status, eps_v_bar_hi, p_bar_hi) = d_hydrostat.intersect1D(-10, 10, p_bar_hi, p_bar_hi);

  double K = (p_bar_hi - p_bar)/(eps_v_bar_hi - eps_v_bar);
  double dJ_deps_v = J;
  double drho_dJ = -rho_0/(J * J);
  dp_drho = K/(drho_dJ * dJ_deps_v);

  
  soundSpeedSq = K / rho_cur;
  if (std::isnan(dp_drho)) {
    std::ostringstream out;
    out << "rho_0 = " << rho_0 << " rho_cur = " << rho_cur << "\n";
    out << "status = " << status << " eps_v_bar = " << eps_v_bar 
        << " p_bar = " << p_bar << " eps_v  = " << eps_v << "\n";
    out << "status = " << status << " eps_v_bar_hi = " << eps_v_bar_hi 
              << " p_bar_hi = " << p_bar_hi << " p_bar  = " << p_bar << "\n";
    out << " K = " << K << " drho_dJ = " << drho_dJ << " dp_drho = " << dp_drho
        << " c2 = " << soundSpeedSq << "\n";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }
}

double
TabularPlasticity::getCompressibility()
{
  double comp = 0.0;
  if (comp < 1.0) {
    std::ostringstream out;
    out << "The getCompressibility MPMICE hook has not been implemented for Tabular plasticity.";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }
  return comp;
}
