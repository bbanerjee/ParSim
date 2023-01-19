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

// Namespace Vaango::
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularPlasticityCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldConditionFactory.h>

// Namespace Uintah::
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
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

//#define TEST_GAMMA_H_SCALING
#ifdef TEST_GAMMA_H_SCALING
constexpr double gamma_H_factor = 1.5;
#endif

//#define TEST_K_VARIATION
#ifdef TEST_K_VARIATION
//constexpr double K_scale_factor = 2;
constexpr double K_scale_factor = 1.7;
#endif

//#define TEST_CAP_VARIATION
#ifdef TEST_CAP_VARIATION
constexpr double X_scale_factor = 1.2;
//constexpr double X_scale_factor = 0.7;
#endif

//#define TEST_CAP_TRANSLATION
#ifdef TEST_CAP_TRANSLATION
constexpr double X_trans_factor = 0.05;
#endif

#define DO_COMPUTE_SIGMA_FIXED
#define COMPUTE_PLASTIC_STRAIN_FROM_ELASTIC
//#define DO_ELASTIC_PLASTIC_COUPLING

//#define DO_CONSISTENCY_BISECTION_ELASTIC
//#define DO_CONSISTENCY_BISECTION_SIMPLIFIED
//#define DO_FIRST_ORDER_HARDENING

//#define DEBUG_WITH_PARTICLE_ID
#ifdef DEBUG_WITH_PARTICLE_ID
constexpr long64 testParticleID = 158913855488;
#endif

#define CHECK_FOR_NAN
//#define CHECK_FOR_NAN_EXTRA

//#define DEBUG_STEP_DIVISIONS
//#define DEBUG_SUBSTEP
//#define DEBUG_SIGMA_FIXED
//#define DEBUG_SIGMA_HARDENING
//#define CHECK_GAMMA_F_POSITIVE

//#define CHECK_CONSISTENCY_BISECTION_CONVERGENCE
//#define CHECK_CONSISTENCY_BISECTION_K
//#define CHECK_CONSISTENCY_BISECTION_FIXED
//#define CHECK_MODULUS_EVOLUTION
//#define DEBUG_FIRST_ORDER_HARDENING
//#define CHECK_PLASTIC_RATE
//#define WRITE_YIELD_SURF
//#define CHECK_INTERNAL_VAR_EVOLUTION
//#define DEBUG_INTERNAL_VAR_EVOLUTION
//#define DEBUG_INTERNAL_VAR_EVOLUTION_COMPUTATION
//#define CHECK_HYDROSTATIC_TENSION
//#define CHECK_TENSION_STATES
//#define CHECK_TENSION_STATES_1
//#define CHECK_DAMAGE_ALGORITHM
//#define CHECK_TRIAL_STRESS
//#define CHECK_YIELD_SURFACE_NORMAL
//#define CHECK_FLOATING_POINT_OVERFLOW
//#define DEBUG_YIELD_BISECTION_R
//#define CHECK_ELASTIC_STRAIN
//#define CHECK_RETURN_ALIGNMENT
//#define TIME_TABLE_LOOKUP
//#define TIME_SUBSTEP
//#define TIME_YIELD_PTS

using namespace Vaango;
using Uintah::VarLabel;
using Uintah::Matrix3;
using std::ostringstream;

TabularPlasticityCap::TabularPlasticityCap(Uintah::ProblemSpecP& ps,
                                           Uintah::MPMFlags* mpmFlags)
  : TabularPlasticity(ps, mpmFlags)
{
  ps->getWithDefault("consistency_bisection_tolerance",
                     d_consistency_bisection_tolerance,
                     1.0e-4);
  d_max_bisection_iterations =
    (int)std::ceil(-10.0 * std::log(d_consistency_bisection_tolerance));

  if (d_consistency_bisection_tolerance < 1.0e-16 ||
      d_consistency_bisection_tolerance > 1.0e-2) {
     std::ostringstream warn;
    warn << "Consistency bisection tolerance should be in range [1.0e-16, "
            "1.0e-2].  Default = 1.0e-4"
         << endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  ps->getWithDefault(
    "decrease_substep_at_high_curvature", d_decrease_substep, false);
}

TabularPlasticityCap::TabularPlasticityCap(const TabularPlasticityCap& cm)
  : TabularPlasticity(cm)
{
  // Consistency bisection
  d_consistency_bisection_tolerance = cm.d_consistency_bisection_tolerance;
  d_max_bisection_iterations        = cm.d_max_bisection_iterations;
  d_decrease_substep                = cm.d_decrease_substep;
}

TabularPlasticityCap::TabularPlasticityCap(const TabularPlasticityCap* cm)
  : TabularPlasticity(*cm)
{
  // Consistency bisection
  d_consistency_bisection_tolerance = cm->d_consistency_bisection_tolerance;
  d_max_bisection_iterations        = cm->d_max_bisection_iterations;
  d_decrease_substep                = cm->d_decrease_substep;
}

TabularPlasticityCap::~TabularPlasticityCap()
{
}

// adds problem specification values to checkpoint data for restart
void
TabularPlasticityCap::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "tabular_plasticity_cap");
  }

  d_elastic->outputProblemSpec(cm_ps);
  d_yield->outputProblemSpec(cm_ps);
  d_capX->outputProblemSpec(cm_ps);
  d_hydrostat.outputProblemSpec(cm_ps);

  cm_ps->appendElement("yield_surface_radius_scaling_factor",
                       d_cm.yield_scale_fac);
  cm_ps->appendElement("subcycling_characteristic_number",
                       d_cm.subcycling_characteristic_number);
  cm_ps->appendElement("consistency_bisection_tolerance",
                       d_consistency_bisection_tolerance);
  cm_ps->appendElement("decrease_substep_at_high_curvature",
                       d_decrease_substep);
}

std::unique_ptr<Uintah::ConstitutiveModel>
TabularPlasticityCap::clone()
{
  return std::make_unique<TabularPlasticityCap>(*this);
}

// When a particle is pushed from patch to patch, carry information needed for
// the particle
void
TabularPlasticityCap::addParticleState(std::vector<const VarLabel*>& from,
                                       std::vector<const VarLabel*>& to)
{
  TabularPlasticity::addParticleState(from, to);
  d_capX->addParticleState(from, to);
}

/*!------------------------------------------------------------------------*/
void
TabularPlasticityCap::addInitialComputesAndRequires(Task* task,
                                                    const MPMMaterial* matl,
                                                    const PatchSet* patch) const
{
  TabularPlasticity::addInitialComputesAndRequires(task, matl, patch);
  d_capX->addInitialComputesAndRequires(task, matl, patch);
}

/*!------------------------------------------------------------------------*/
void
TabularPlasticityCap::initializeCMData(const Patch* patch,
                                       const MPMMaterial* matl,
                                       DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  d_capX->initializeInternalVariable(pset, new_dw);
  TabularPlasticity::initializeCMData(patch, matl, new_dw);
}

// Compute stable timestep based on both the particle velocities
// and wave speed
void
TabularPlasticityCap::computeStableTimestep(const Patch* patch,
                                            const MPMMaterial* matl,
                                            DataWarehouse* new_dw)
{
  TabularPlasticity::computeStableTimestep(patch, matl, new_dw);
}

/**
 * Added computes/requires for computeStressTensor
 */
void
TabularPlasticityCap::addComputesAndRequires(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet* patches) const
{
  TabularPlasticity::addComputesAndRequires(task, matl, patches);
  d_capX->addComputesAndRequires(task, matl, patches);
}

/**
 *  TabularPlasticityCap::computeStressTensor
 *  is the core of the TabularPlasticityCap model which computes
 *  the updated stress at the end of the current timestep along with all other
 *  required data such plastic strain and elastic strain
 */
void
TabularPlasticityCap::computeStressTensor(const PatchSubset* patches,
                                          const MPMMaterial* matl,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
#ifdef TIME_YIELD_PTS
  std::chrono::time_point<std::chrono::system_clock> start_yp, end_yp;
#endif

#ifdef TIME_SUBSTEP
  std::chrono::time_point<std::chrono::system_clock> start_ss, end_ss;
#endif

  // Initialize total strain energy
  double se = 0.0;

  // Global loop over each patch
  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    // Initialize wave speed
    double c_dil = std::numeric_limits<double>::min();
    Vector waveSpeed(c_dil, c_dil, c_dil);
    Vector dx = patch->dCell();

    // Get particle subset for the current patch
    int matID            = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Set up local particle variables to be read and written
    constParticleVariable<double> pEev_old, pEpv_old, pEpeq_old;
    constParticleVariable<double> pBulkModulus_old;
    constParticleVariable<Matrix3> pEe_old, pEp_old;
    old_dw->get(pEe_old, pElasticStrainLabel, pset);
    old_dw->get(pEp_old, pPlasticStrainLabel, pset);
    old_dw->get(pEev_old, pElasticVolStrainLabel, pset);
    old_dw->get(pEpv_old, pPlasticVolStrainLabel, pset);
    old_dw->get(pEpeq_old, pPlasticCumEqStrainLabel, pset);
    old_dw->get(pBulkModulus_old, pBulkModulusLabel, pset);

    ParticleVariable<int> pRemove_new;
    new_dw->getModifiable(pRemove_new, lb->pRemoveLabel_preReloc, pset);

    ParticleVariable<double> pEev_new, pEpv_new, pEpeq_new;
    ParticleVariable<double> pBulkModulus_new;
    ParticleVariable<Matrix3> pEe_new, pEp_new;
    new_dw->allocateAndPut(pEe_new, pElasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEev_new, pElasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEp_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpv_new, pPlasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpeq_new, pPlasticCumEqStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pBulkModulus_new, pBulkModulusLabel_preReloc, pset);

    // Get and allocate the hydrostatic strength
    constParticleVariable<double> pCapX;
    ParticleVariable<double> pCapX_new;
    d_capX->getInternalVariable(pset, old_dw, pCapX);
    d_capX->allocateAndPutInternalVariable(pset, new_dw, pCapX_new);

    // Set up global particle variables to be read and written
    delt_vartype delT;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass, pVolume;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pStress_old;

    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    new_dw->get(pVolume, lb->pVolumeLabel_preReloc, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);

    // Get the particle variables computed in interpolateToParticlesAndUpdate()
    constParticleVariable<Matrix3> pDefRate_mid, pDefGrad_new;
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> p_q, pdTdt;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    // Loop over particles
    for (particleIndex idx : *pset) {

      // No thermal effects
      pdTdt[idx] = 0.0;

      Matrix3 DD        = pDefRate_mid[idx];
      Matrix3 sigma_old = pStress_old[idx];

      // Set up model state
      ModelState_TabularCap state_old;
      state_old.particleID   = pParticleID[idx];
      state_old.stressTensor = sigma_old;
      state_old.updateStressInvariants();
      state_old.elasticStrainTensor = pEe_old[idx];
      state_old.plasticStrainTensor = pEp_old[idx];
      state_old.updatePlasticStrainInvariants();
      state_old.ep_cum_eq = pEpeq_old[idx];
      state_old.capX      = pCapX[idx];
      // std::cout << "state_old.Stress = " << state_old.stressTensor <<
      // "\n";

      // Compute the elastic moduli at t = t_n
      state_old.updateStressInvariants();
      state_old.updatePlasticStrainInvariants();
      computeElasticProperties(state_old);

      // Update yield surface 
      // (**REQUIRES** elastic properties to be updated in state_old)
      #ifdef TIME_YIELD_PTS
      start_yp = std::chrono::system_clock::now();
      #endif
      Polyline yield_f_pts =
        d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_old);
      #ifdef TIME_YIELD_PTS
      end_yp = std::chrono::system_clock::now();
      #endif
      std::array<double, 3> range =
        d_yield->getYieldConditionRange(yield_f_pts);
      state_old.yield_f_pts = yield_f_pts;
      state_old.I1_max      = -range[0] * 3.0;
      state_old.I1_min      = -range[1] * 3.0;
      state_old.sqrtJ2_max  = range[2];

      // std::cout << "State old: " << state_old << "\n";

      //---------------------------------------------------------
      // Rate-independent plastic step
      // Divides the strain increment into substeps, and calls substep function
      ModelState_TabularCap state_new;
      #ifdef TIME_SUBSTEP
      start_ss = std::chrono::system_clock::now();
      #endif
      Status status = rateIndependentPlasticUpdate(
        DD, delT, idx, pParticleID[idx], state_old, state_new);
      #ifdef TIME_SUBSTEP
      end_ss = std::chrono::system_clock::now();
      #endif

      #ifdef TIME_SUBSTEP
      auto duration = std::chrono::duration<double>(end_ss - start_ss).count();
      if (duration > 0.1) {
        std::cout << "Substep : Particle = " << pParticleID[idx]
                  << " Time taken = " << duration << "\n";
        #ifdef TIME_YIELD_PTS
        std::cout << "Compute yield pts. : Particle = " << pParticleID[idx]
                  << " Time taken = "
                  << std::chrono::duration<double>(end_yp - start_yp).count()
                  << "\n";
        #endif
      }
      #endif

      if (status == Status::SUCCESS) {
        pRemove_new[idx] = 0;
        pStress_new[idx] =
          state_new.stressTensor; // unrotated stress at end of step
        pEe_new[idx] =
          state_new.elasticStrainTensor; // elastic strain at end of step
        pEp_new[idx] =
          state_new.plasticStrainTensor; // plastic strain at end of step
        pEev_new[idx] = pEe_new[idx].Trace();
        pEpv_new[idx] = pEp_new[idx].Trace();
        pEpeq_new[idx] =
          state_new.ep_cum_eq; // Equivalent plastic strain at end of step
        pCapX_new[idx]        = state_new.capX;
        pBulkModulus_new[idx] = state_new.bulkModulus;

        /*
        // Elastic volumetric strain at end of step, compute from updated
        // deformation gradient.
        // H = ln(U) => tr(H) = tr(ln(U)) = ln(det(U)) = ln(sqrt(det(FT)
        // det(F))) = ln J
        pEev_new[idx] = log(pDefGrad_new[idx].Determinant()) - pEpv_new[idx];

        #ifdef CHECK_ELASTIC_STRAIN
        double pEev_integrated = pEe_new[idx].Trace();
        std::cout << "Elastic volume strain error = "
                  << (pEev_new[idx] - pEev_integrated)
                  << " Integrated = " << pEev_integrated
                  << " Defgrad-based = " << pEev_new[idx] << "\n";
        #endif
        */

      } else {

        // If the updateStressAndInternalVars function can't converge it will
        // return false.
        // This indicates substepping has failed, and the particle will be
        // deleted.
        pRemove_new[idx] = -999;
        std::cout << "** WARNING ** Bad step, deleting particle"
                  << " idx = " << idx << " particleID = " << pParticleID[idx]
                  << ":" << __FILE__ << ":" << __LINE__ << "\n";

        pStress_new[idx] = pStress_old[idx];
        pEe_new[idx] =
          state_old.elasticStrainTensor; // elastic strain at start of step
        pEp_new[idx] =
          state_old.plasticStrainTensor; // plastic strain at start of step
        pEpv_new[idx]         = pEp_new[idx].Trace();
        pEpeq_new[idx]        = pEpeq_old[idx];
        pEev_new[idx]         = pEe_new[idx].Trace();
        pCapX_new[idx]        = state_old.capX;
        pBulkModulus_new[idx] = pBulkModulus_old[idx];
      }

      // Compute wave speed + particle velocity at each particle, store the
      // maximum
      // std::cout << "State QS new rotated";
      computeElasticProperties(state_new);
      double bulk    = state_new.bulkModulus;
      double shear   = state_new.shearModulus;
      double rho_cur = pMass[idx] / pVolume[idx];
      c_dil          = std::sqrt((bulk + Util::four_third * shear) / rho_cur);
      Vector velMax  = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed      = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) * Util::one_third;
        double c_bulk = sqrt(bulk / rho_cur);
        p_q[idx] = artificialBulkViscosity(DD.Trace(), c_bulk, rho_cur, dx_ave);
        #ifdef DEBUG_WITH_PARTICLE_ID
        if (state_new.particleID == testParticleID) {
          std::cout << "K = " << bulk << " G = " << shear
                    << " c_bulk = " << c_bulk << "\n";
          std::cout << "mass = " << pMass[idx] << " vol = " << pVolume[idx]
                    << " dx_ave = " << dx_ave << " rho_cur = " << rho_cur
                    << " tr(D) = " << DD.Trace() << "\n";
        }
        #endif
      } else {
        p_q[idx] = 0.;
      }

      // Compute the strain energy increment associated with the particle
      // and accumulate total strain energy
      Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * 0.5;
      double rateOfWork = computeRateOfWork(avgStress, DD);
      se +=  (rateOfWork * pVolume[idx] * delT);

    } // End particle set loop

    // Compute the stable timestep based on maximum value of "wave speed +
    // particle velocity"
    waveSpeed = dx / waveSpeed; 
    double delT_new = waveSpeed.minComponent();

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
TabularPlasticityCap::rateIndependentPlasticUpdate(
  const Matrix3& D,
  const double& delT,
  particleIndex idx,
  long64 pParticleID,
  const ModelState_TabularCap& state_old,
  ModelState_TabularCap& state_new)
{
  #ifdef CHECK_FOR_NAN_EXTRA
    #ifdef DEBUG_WITH_PARTICLE_ID
    if (state_old.particleID == testParticleID) {
    #endif
      std::cout << "Rate independent update:" << "\n";
      std::cout << " D = " << D << " delT = " << delT << "\n";
      std::cout << "\t State old:" << state_old << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
    }
    #endif
  #endif

  // Compute the strain increment
  Matrix3 strain_inc = D * delT;
  if (strain_inc.Norm() < 1.0e-30) {
    state_new = state_old;
    return Status::SUCCESS;
  }

  // std::cout << "delta Eps = " << strain_inc << "\n";

  // Set up a trial state, update the stress invariants, and compute elastic
  // properties
  ModelState_TabularCap state_trial(state_old);
  computeElasticProperties(state_trial, strain_inc.Trace(), 0.0);
  state_trial.updatePlasticStrainInvariants();

  // Compute the trial stress
  Matrix3 stress_trial     = computeTrialStress(state_trial, strain_inc);
  state_trial.stressTensor = stress_trial;
  state_trial.updateStressInvariants();

  #ifdef CHECK_FOR_NAN_EXTRA
    #ifdef DEBUG_WITH_PARTICLE_ID
    if (state_old.particleID == testParticleID) {
    #endif
      std::cout << "\t strain_inc = " << strain_inc << "\n";
      std::cout << "\t State trial:" << state_trial << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
    }
    #endif
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
  /*
  if (nsub < 0) {
    state_new = state_old;
    proc0cout << "Step Failed: Particle idx = " << idx
              << " ID = " << pParticleID << " because nsub = " << nsub
              << "\n";
    // bool success  = false;
    return Status::UNREASONABLE_SUBSTEPS;
  }
  */

  // Set up the initial states for the substeps
  Status status = Status::SUCCESS;
  ModelState_TabularCap state_k_old;
  ModelState_TabularCap state_k_new;

  const int CHI_MAX = 10; // max allowed number of timestep reductions
  int chi           = 0;  // subcycle timestep reduction count
  do {
    // Compute the substep time increment list
    std::vector<double> substeps(nsub, delT / static_cast<double>(nsub));

    state_k_old   = state_old;
    state_k_new   = state_old;

    // Loop through the substeps
    double tlocal = 0.0;
    for (const auto dt : substeps) {

      // Compute the substep
      status = computeSubstep(D, dt, state_k_old, state_k_new);
      if (status != Status::SUCCESS) {
        #ifdef CHECK_RETURN_ALIGNMENT
          std::cout << "dt = " << dt << " chi = " << chi << "\n";
        #endif
        nsub *= 2; // Double the number of substeps
        ++chi;     // Increment the number of dt decreases
        proc0cout << "**WARNING** Decreasing substep time increment to "
                  << dt / 2 << " because computeSubstep failed." << "\n";
        break;
      }

      state_k_old = state_k_new;
      tlocal += dt;

      #ifdef WRITE_YIELD_SURF
        std::cout << "K = " << state_k_new.bulkModulus << "\n";
        std::cout << "G = " << state_k_new.shearModulus << "\n";
        std::cout << "capX = " << state_k_new.capX << "\n";
        Matrix3 sig = state_k_new.stressTensor;
        std::cout << "sigma_new = np.array([[" << sig(0, 0) << "," << sig(0, 1)
                  << "," << sig(0, 2) << "],[" << sig(1, 0) << "," << sig(1, 1)
                  << "," << sig(1, 2) << "],[" << sig(2, 0) << "," << sig(2, 1)
                  << "," << sig(2, 2) << "]])" << "\n";
        std::cout << "plot_stress_state(K, G, sigma_trial, sigma_new, 'b')"
                  << "\n";
      #endif

      #ifdef DEBUG_SUBSTEP
        std::cout << "tlocal = " << tlocal << " delT = " << delT
                  << " nsub = " << nsub << "\n";
      #endif
    } // End substep range for loop

  } while (status != Status::SUCCESS && chi < CHI_MAX);

  if (status == Status::SUCCESS) {
    state_new = state_k_new;

    #ifdef CHECK_INTERNAL_VAR_EVOLUTION
      std::cout << "rateIndependentPlasticUpdate: "
                << " ep_v_old = " << state_old.ep_v
                << " ep_v_new = " << state_new.ep_v
                << " Xbar_old = " << -state_old.capX
                << " Xbar_new = " << -state_new.capX << "\n";
    #endif
  } else {
    if (status == Status::TOO_LARGE_YIELD_NORMAL_CHANGE) {
      state_new = state_k_new;
      status    = Status::SUCCESS;
      proc0cout << "Substep failed because chi = " << chi << " > " << CHI_MAX
                << "."
                << " Continuing with unconverged solution." << "\n";
    } else {
      state_new = state_k_old;
      proc0cout << "Substep failed because chi = " << chi << " > " << CHI_MAX
                << "."
                << " Proceeding to delete particle." << "\n";
    }
  }

  return status;
}

/**
 * Method: computeElasticProperties
 *
 * Purpose:
 *   Compute the bulk and shear modulus at a given state
 *
 */
void
TabularPlasticityCap::computeElasticProperties(
  ModelState_TabularCap& state) const
{
  #ifdef TIME_TABLE_LOOKUP
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
  #endif
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(&state);
  state.bulkModulus    = moduli.bulkModulus;
  state.shearModulus   = moduli.shearModulus;
  #ifdef TIME_TABLE_LOOKUP
    end = std::chrono::system_clock::now();
    std::cout << "Compute elastic properties : Time taken = "
              << std::chrono::duration<double>(end - start).count() << "\n";
  #endif
  #ifdef TEST_K_VARIATION
    state.bulkModulus *= K_scale_factor;
    state.shearModulus *= K_scale_factor;
  #endif
}

void
TabularPlasticityCap::computeElasticProperties(ModelState_TabularCap& state,
                                               double elasticVolStrainInc,
                                               double plasticVolStrainInc) const
{
  ModelState_TabularCap tempState = state;
  tempState.elasticStrainTensor += (Util::Identity * elasticVolStrainInc);
  tempState.plasticStrainTensor += (Util::Identity * plasticVolStrainInc);
  tempState.updatePlasticStrainInvariants();
  computeElasticProperties(tempState);
  state.bulkModulus  = tempState.bulkModulus;
  state.shearModulus = tempState.shearModulus;
}

std::tuple<double, double>
TabularPlasticityCap::computeElasticProperties(
  const ModelState_TabularCap& state,
  const Matrix3& elasticStrain,
  const Matrix3& plasticStrain) const
{
  ModelState_TabularCap temp = state;
  temp.elasticStrainTensor   = elasticStrain;
  temp.plasticStrainTensor   = plasticStrain;
  temp.updatePlasticStrainInvariants();
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(&temp);
  #ifdef TEST_K_VARIATION
    moduli.bulkModulus *= K_scale_factor;
    moduli.shearModulus *= K_scale_factor;
  #endif
  return std::make_tuple(moduli.bulkModulus, moduli.shearModulus);
}

std::tuple<double, double>
TabularPlasticityCap::computeElasticProperties(
  const ModelState_TabularCap& state) const
{
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(&state);
  #ifdef TEST_K_VARIATION
    moduli.bulkModulus *= K_scale_factor;
    moduli.shearModulus *= K_scale_factor;
  #endif
  return std::make_tuple(moduli.bulkModulus, moduli.shearModulus);
}

/**
 * Method: computeTrialStress
 * Purpose:
 *   Compute the trial stress for some increment in strain assuming linear
 * elasticity
 *   over the step.
 */
Matrix3
TabularPlasticityCap::computeTrialStress(const ModelState_TabularCap& state_old,
                                         const Matrix3& strain_inc)
{
  return TabularPlasticity::computeTrialStress(state_old, strain_inc);
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
TabularPlasticityCap::computeStepDivisions(
  particleIndex idx,
  long64 particleID,
  const ModelState_TabularCap& state_old,
  const ModelState_TabularCap& state_trial)
{
  // Compute subdivisions based on change in bulk modulus
  double bulk_old = state_old.bulkModulus;
  double bulk_trial = state_trial.bulkModulus;
  double d_bulk = std::abs(bulk_old - bulk_trial) / std::min(bulk_old, bulk_trial);
  int n_bulk = std::max(std::ceil(d_bulk * d_cm.yield_scale_fac * 1000), 1.0);
  #ifdef DEBUG_STEP_DIVISIONS
    #ifdef DEBUG_WITH_PARTICLE_ID
    if (state_old.particleID == testParticleID) {
    #endif
      std::cout << "bulk_old = " << bulk_old << " bulk_trial = " << bulk_trial
                << " d_bulk = " << d_bulk << " n_bulk = " << n_bulk << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
    }
    #endif
  #endif

  // Compute subdivisions based on trial stress increment 
  // relative to yield surface size
  double I1_old = state_old.I1;
  double I1_trial = state_trial.I1;
  double J2_old = state_old.sqrt_J2;
  double J2_trial = state_trial.sqrt_J2;
  double I1_size  = 0.5 * std::abs(state_old.I1_max - state_old.I1_min);
  double J2_size  = state_old.sqrtJ2_max;
  double ref_I1_old = (I1_old == 0.0) ? I1_size : std::abs(I1_old);
  double ref_I1_trial = (I1_trial == 0.0) ? I1_size : std::abs(I1_trial);
  double ref_J2_old = (J2_old == 0.0) ? J2_size : J2_old;
  double ref_J2_trial = (J2_trial == 0.0) ? J2_size : J2_trial;
  double ref_I1 = std::min(ref_I1_old, ref_I1_trial);
  double ref_J2 = std::min(ref_J2_old, ref_J2_trial);
  double d_I1 = std::abs(I1_old - I1_trial) / ref_I1;
  double d_J2 = std::abs(J2_old - J2_trial) / ref_J2;
  int n_yield_I1 = std::ceil(d_I1 * d_cm.yield_scale_fac * 0.1);
  int n_yield_J2 = std::ceil(d_J2 * d_cm.yield_scale_fac * 0.1);
  int n_yield = std::max(std::max(n_yield_I1, n_yield_J2), 1);

  #ifdef DEBUG_STEP_DIVISIONS
    #ifdef DEBUG_WITH_PARTICLE_ID
    if (state_old.particleID == testParticleID) {
    #endif
      proc0cout << "bulk_old = " << bulk_old << " bulk_trial = " << bulk_trial
                << " n_bulk = " << n_bulk << "\n";
      proc0cout << "I1_max = " << state_old.I1_max
                << " I1_min = " << state_old.I1_min
                << " sqrtJ2_max = " << state_old.sqrtJ2_max 
                << " I1_size = " << I1_size << " J2_size = " << J2_size
                << " d_I1 = " << d_I1 << " d_J2 = " << d_J2 << " n_yield = " << n_yield
                << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
    }
    #endif
  #endif

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = std::max(std::max(n_bulk, n_yield), 1);
  int nmax = d_cm.subcycling_characteristic_number;

  #ifdef DEBUG_STEP_DIVISIONS
  std::cout << "I1_old = " << I1_old << " I1_trial = " << I1_trial
            << " ref_I1 = " << ref_I1 << " d_I1 = " << d_I1 << "\n";
  std::cout << "n_bulk = " << n_bulk 
            << " n_yield = " << n_yield_I1 << ", " << n_yield_J2
            << " nsub = " << nsub << " nmax = " << nmax << "\n";
  #endif

  #ifdef DEBUG_STEP_DIVISIONS
  if (nsub > nmax) {
    proc0cout << "\n **WARNING** Too many substeps needed for particle "
              << " idx = " << idx << " particle ID = " << particleID << "\n";
    proc0cout << "\t" << __FILE__ << ":" << __LINE__ << "\n";
    proc0cout << "\t State at t_n: " << state_old;

    proc0cout << "\t Trial state at t_n+1: " << state_trial;

    proc0cout << "\t Ratio of trial bulk modulus to t_n bulk modulus " << n_bulk
              << "\n";

    proc0cout << "\t d_I1 " << d_I1 << " d_sqrtJ2 = " << d_J2 << "\n";
    proc0cout << "\t Yield surface radius in I1-space: " << I1_size 
              << " J2-space: " << J2_size << "\n";
    proc0cout << "\t Ratio of ||sig_trial - sigma_n|| and "
              << d_cm.yield_scale_fac << "*y.s. radius: " << n_yield
              << "\n";
    proc0cout << "\t I1_max = " << state_old.I1_max
              << " I1_min = " << state_old.I1_min << "\n";
    proc0cout << "\t I1_size = " << I1_size << " J2_size = " << J2_size
              << "\n";

    proc0cout << "** BECAUSE** nsub = " << nsub << " > "
              << d_cm.subcycling_characteristic_number
              << " : Probably too much tension in the particle." << "\n";
  }
  #endif

  nsub = std::min(nmax, nsub);

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
TabularPlasticityCap::computeSubstep(const Matrix3& D,
                                     const double& dt,
                                     const ModelState_TabularCap& state_k_old,
                                     ModelState_TabularCap& state_k_new)
{
  // Set up a trial state, update the stress invariants, and compute trial
  // stress
  ModelState_TabularCap state_k_trial(state_k_old);
  Matrix3 deltaEps       = D * dt;
  Matrix3 stress_k_trial = computeTrialStress(state_k_trial, deltaEps);

  // Update the trial stress and compute invariants
  state_k_trial.stressTensor = stress_k_trial;
  state_k_trial.updateStressInvariants();
  state_k_trial.updatePlasticStrainInvariants();
  //computeElasticProperties(state_k_trial);

  // std::cout << "\t D = " << D << "\n";
  // std::cout << "\t dt:" << dt << "\n";
  // std::cout << "\t deltaEps = " << deltaEps << "\n";

  #ifdef CHECK_FOR_NAN_EXTRA
    #ifdef DEBUG_WITH_PARTICLE_ID
    if (state_k_old.particleID == testParticleID) {
    #endif
      std::cout << "\t D = " << D << "\n";
      std::cout << "\t dt:" << dt << "\n";
      std::cout << "\t deltaEps = " << deltaEps << "\n";
      std::cout << "\t Stress k trial:" << stress_k_trial << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
    }
    #endif
  #endif

  #ifdef WRITE_YIELD_SURF
    // std::cout << "Inside computeSubstep:" << "\n";
    std::cout << "K = " << state_k_old.bulkModulus
              << "G = " << state_k_old.shearModulus << "X = " << state_k_old.capX
              << "\n";
    Matrix3 sig = stress_k_trial;
    std::cout << "sigma_trial = np.array([[" << sig(0, 0) << "," << sig(0, 1)
              << "," << sig(0, 2) << "],[" << sig(1, 0) << "," << sig(1, 1) << ","
              << sig(1, 2) << "],[" << sig(2, 0) << "," << sig(2, 1) << ","
              << sig(2, 2) << "]])" << "\n";
    std::cout << "plot_stress_state(K, G, sigma_new, sigma_trial, 'r')"
              << "\n";
    // std::cout << "\t computeSubstep: sigma_old = " << state_k_old.stressTensor
    //         << " sigma_trial = " << stress_trial
    //         << " D = " << D << " dt = " << dt
    //         << " deltaEps = " << deltaEps << "\n";
  #endif

  #ifdef CHECK_MODULUS_EVOLUTION
    std::cout << "Trial: K = " << state_k_trial.bulkModulus
              << " G = " << state_k_trial.shearModulus
              << " X = " << state_k_trial.capX << " ep_v = " << state_k_trial.ep_v
              << " ee_v = " << state_k_trial.elasticStrainTensor.Trace()
              << "\n";
  #endif

  #ifdef DEBUG_SUBSTEP
  std::cout << "p_trial = " << state_k_trial.I1/3.0 
            << " p_yield = " << state_k_trial.capX/3.0 << "\n";
  #endif

  // Evaluate the yield function at the trial stress:
  auto yield = d_yield->evalYieldCondition(&state_k_trial);

  // Elastic substep
  if (yield.second == Util::YieldStatus::IS_ELASTIC) {
    state_k_new = state_k_trial;
    state_k_new.elasticStrainTensor += deltaEps;

    #ifdef CHECK_INTERNAL_VAR_EVOLUTION
      std::cout << "computeSubstep:Elastic:sigma_new = "
                << state_k_new.stressTensor
                << " ep_v_trial = " << state_k_trial.ep_v << "\n";
    #endif

    return Status::SUCCESS; // bool isSuccess = true;
  }

  #ifdef DEBUG_YIELD_BISECTION_R
    std::cout << "before_non_hardening_return  = 1" << "\n";
    std::cout << "I1 = " << state_k_elastic.I1 << "\n";
    std::cout << "sqrt_J2 = " << state_k_elastic.sqrt_J2 << "\n";
  #endif

  // Elastic-plastic or fully-plastic substep
#ifdef DO_COMPUTE_SIGMA_FIXED

  // Compute the return to the fixed yield surface
  double H_F, Gamma_F;
  Matrix3 stress_F, P_F, N_F_norm;
  std::tie(stress_F, P_F, N_F_norm, H_F, Gamma_F) = 
    computeSigmaFixed(state_k_old, state_k_trial);

  // Compute hardening correction
  #ifdef NO_SIGMA_HARDENING
  double Gamma_H = Gamma_F;
  Matrix3 stress_H = stress_F;
  #else
  double Gamma_H;
  Matrix3 stress_H;
  std::tie(stress_H, Gamma_H) = 
    computeSigmaHardening(state_k_trial, Gamma_F, P_F, N_F_norm, H_F);
  #endif

  #ifdef COMPUTE_PLASTIC_STRAIN_FROM_ELASTIC
  // Compute elastic strain increment and then plastic strain increment 
  Matrix3 stress_inc     = stress_H - state_k_old.stressTensor;
  Matrix3 stress_inc_iso = Util::Identity * (Util::one_third * stress_inc.Trace());
  Matrix3 stress_inc_dev = stress_inc - stress_inc_iso;
  Matrix3 eps_e_inc      = stress_inc_iso * (1.0 / (3.0 * state_k_old.bulkModulus)) +
                           stress_inc_dev * (1.0 / (2.0 * state_k_old.shearModulus));
  Matrix3 eps_p_inc      = deltaEps - eps_e_inc;
  #else
  // Compute plastic strain increment and then elastic strain increment 
  auto eps_p_inc = N_F_norm * Gamma_H;
  auto eps_e_inc = deltaEps - eps_p_inc;
  #endif
  
  // Update the state
  state_k_new = state_k_trial;
  state_k_new.stressTensor = stress_H;
  state_k_new.updateStressInvariants();
  state_k_new.elasticStrainTensor = state_k_old.elasticStrainTensor + eps_e_inc;
  state_k_new.plasticStrainTensor = state_k_old.plasticStrainTensor + eps_p_inc;
  state_k_new.updatePlasticStrainInvariants();

  // Update internal variable
  state_k_new.capX = computeInternalVariable(state_k_new);

  // Update elastic properties
  computeElasticProperties(state_k_new);

  // Update yield surface
  auto yield_pts_updated =
    d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_k_new);
  state_k_new.yield_f_pts = yield_pts_updated;

  /*
  std::cout << "K_old = " << state_k_old.bulkModulus << " K_new = " << state_k_new.bulkModulus << "\n"
            << "X_old = " << state_k_old.capX << " X_new = " << state_k_new.capX << "\n"
            << "p_old = " << state_k_old.stressTensor.Trace()/3.0 
            << " p_new = " << state_k_new.stressTensor.Trace()/3.0 << "\n"
            << "eps_v_p_old = " << state_k_old.plasticStrainTensor.Trace() 
            << " eps_v_p_new = " << state_k_new.plasticStrainTensor.Trace() << "\n\n";
  std::cout << "gamma_F = " << Gamma_F << " gamma_H = " << Gamma_H << " N_F_norm = " << N_F_norm << "\n";
  std::cout << "stress_inc = " << stress_inc << "\n";
  std::cout << "deltaEps = " << deltaEps << "\n"
            << "eps_p_inc = " << eps_p_inc << "\n"
            << "eps_e_inc = " << eps_e_inc << "\n";
  std::cout << "eps_e = " << state_k_new.elasticStrainTensor << "\n"
            << "eps_p = " << state_k_new.plasticStrainTensor << "\n";
  std::cout << "capX = " << state_k_new.capX << " K = " << state_k_new.bulkModulus << "\n"
            << "sig_F = " << stress_F << "\n"
            << "sig_H = " << stress_H << "\n";
            << "sig_new = " << state_k_new.stressTensor << "\n";
  */
  
  Status status = Status::SUCCESS;

  return status;

#endif

#ifdef DO_CONSISTENCY_BISECTION_ELASTIC

  // 1) Find the purely elastic part of the timestep
  // std::cout << "\t Doing computePurelyElasticSubstep\n";
  Status status;
  double elastic_dt;
  auto state_k_elastic = state_k_old;
  std::tie(status, elastic_dt) = computePurelyElasticSubstep(
    dt, D, state_k_old, state_k_trial, state_k_elastic);
  if (status != Status::SUCCESS) {
    proc0cout << "**WARNING** computePurelyElasticSubstep has failed."
              << "\n";
    return status;
  }
  double elastic_plastic_dt = dt - elastic_dt;

  // 2) Use the state at the end of the purely elastic substep to
  //    do a elastic-plastic non-hardening return to the yield surface
  // std::cout << "\t Doing nonHardeningReturnElasticPlastic\n";
  Matrix3 deltaEps_p_nonhardening;
  auto state_k_nonhardening = state_k_elastic;
  std::tie(status, deltaEps_p_nonhardening) = nonHardeningReturnElasticPlastic(
    elastic_plastic_dt, D, state_k_elastic, state_k_nonhardening);
  if (status != Status::SUCCESS) {
    proc0cout << "**WARNING** nonHardeningReturnElasticPlastic has failed."
              << "\n";
    return status;
  }

  // 3) Use the plastic strain increment computed using the non-hardening
  //    return to estimate the final state with an evolving yield surface
  // std::cout << "\t Doing consistencyBisectionHardeningSoftening\n";
  state_k_new = state_k_nonhardening;
  status      = consistencyBisectionHardeningSoftening(elastic_plastic_dt,
                                                  D,
                                                  state_k_elastic,
                                                  state_k_nonhardening,
                                                  deltaEps_p_nonhardening,
                                                  state_k_new);
  if (status != Status::SUCCESS) {
    proc0cout
      << "**WARNING** consistencyBisectionHardeningSoftening has failed."
      << "\n";
    return status;
  }

#endif

#ifdef DO_CONSISTENCY_BISECTION_SIMPLIFIED
  // Find stress state and strain increments for non-hardening return
  Matrix3 sig_fixed(0.0);
  Matrix3 deltaEps_e_fixed(0.0);
  Matrix3 deltaEps_p_fixed(0.0);

  // std::cout << "\t Doing nonHardeningReturn\n";
  Status status = nonHardeningReturn(deltaEps,
                                     state_k_old,
                                     state_k_trial,
                                     sig_fixed,
                                     deltaEps_e_fixed,
                                     deltaEps_p_fixed);
  if (status != Status::SUCCESS) {
    proc0cout << "**WARNING** nonHardeningReturn has failed." << "\n";
    return status;
  }

  // Do "consistency bisection"
  // std::cout << "\t Doing consistencyBisection\n";
  state_k_new = state_k_old;
  status      = consistencyBisectionSimplified(deltaEps,
                                          state_k_old,
                                          state_k_trial,
                                          deltaEps_e_fixed,
                                          deltaEps_p_fixed,
                                          sig_fixed,
                                          state_k_new);
#endif

#ifdef DO_FIRST_ORDER_HARDENING
  Matrix3 sig_fixed(0.0);
  Matrix3 deltaEps_e_fixed(0.0);
  Matrix3 deltaEps_p_fixed(0.0);

  // std::cout << "\t Doing nonHardeningReturn\n";
  Status status = nonHardeningReturn(deltaEps,
                                     state_k_old,
                                     state_k_trial,
                                     sig_fixed,
                                     deltaEps_e_fixed,
                                     deltaEps_p_fixed);
  if (status != Status::SUCCESS) {
    proc0cout << "**WARNING** nonHardeningReturn has failed." << "\n";
    return status;
  }

  // Do first-order hardening update
  state_k_new = state_k_old;
  status      = firstOrderHardeningUpdate(deltaEps,
                                     state_k_old,
                                     state_k_trial,
                                     deltaEps_e_fixed,
                                     deltaEps_p_fixed,
                                     sig_fixed,
                                     state_k_new);
#endif

#ifdef DEBUG_INTERNAL_VAR_EVOLUTION
  std::cout << "computeSubstep: "
            << " ep_v_old = " << state_k_old.ep_v
            << " ep_v_new = " << state_k_new.ep_v
            << " Xbar_old = " << -state_k_old.capX
            << " Xbar_new = " << -state_k_new.capX << "\n";
#endif

#ifdef DEBUG_YIELD_BISECTION_R
  std::cout << "after_consistency_bisection  = 1" << "\n";
  std::cout << "I1 = " << state_k_new.I1 << "\n";
  std::cout << "sqrt_J2 = " << state_k_new.sqrt_J2 << "\n";
#endif

  return status;

} //===================================================================

/**
 * Compute sigma_F (stress state on the yield surface with updated normal)
 */
std::tuple<Matrix3, Matrix3, Matrix3, double, double>
TabularPlasticityCap::computeSigmaFixed(const ModelState_TabularCap& state_old, 
                                        const ModelState_TabularCap& state_trial) const
{
  // Convert the yield function points to z-rprime coordinates
  std::vector<Uintah::Point> z_r_table(state_old.yield_f_pts.size());
  double sqrtKG = std::sqrt(1.5 * state_old.bulkModulus / state_old.shearModulus);
  size_t ii = 0;
  for (const auto& pt : state_old.yield_f_pts) {
    double z = -Util::sqrt_three * pt.x();
    double r_prime = Util::sqrt_two * pt.y() * sqrtKG;
    z_r_table[ii++] = Uintah::Point(z, r_prime, 0);
  }

  // Create a KD-tree index for the yield surface polyline
  Util::PolylinePointCloud z_r_cloud(z_r_table);
  Util::PolylineKDTree z_r_index(2 /*dim*/, z_r_cloud,
                                 nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  z_r_index.buildIndex();
  
  // Make a copy of the trial state that can be updated per iteration
  ModelState_TabularCap state_trial_iter(state_trial);
  Matrix3 sig_trial = state_trial.stressTensor;

  // Make a copy of the old state that can be updated
  // This is needed for normal computation (df_dsigma)
  // The updated state contains the closest point and tangent after closest point computation
  ModelState_TabularCap state_old_for_normal = state_old;

  double Gamma_F_old = 0.0, Gamma_F = 0.0, N_c_mag;
  Matrix3 sig_c, sig_F, P_c, M_c;
  Point closest;
  Vector tangent;
  int iter = 0;

  do {
    std::tie(sig_c, closest, tangent) = closestPointInZRSpace(state_old_for_normal, state_trial_iter,
                                                              z_r_table, z_r_index);
    state_old_for_normal.closest = closest;
    state_old_for_normal.tangent = tangent;
    Matrix3 N_c = computeYieldSurfaceNormal(state_old_for_normal, sig_c);
    N_c_mag = N_c.Norm();
    Matrix3 N_c_norm = N_c / N_c_mag;
    M_c = N_c_norm;
    P_c = computeProjectionTensor(state_old_for_normal, sig_c, M_c);
    
    Gamma_F_old = Gamma_F;
    Gamma_F = computeGammaClosest(sig_c, sig_trial, N_c_norm, P_c);
    sig_F = sig_trial - P_c * Gamma_F;

    state_trial_iter.stressTensor = sig_F;
    state_trial_iter.updateStressInvariants();

    #ifdef DEBUG_SIGMA_FIXED
    std::cout << "sig_c = " << sig_c << "\n"
              << " N_c = " << N_c << "\n"
              << " P_c = " << P_c << "\n"
              << " sig_trial = " << sig_trial << "\n"
              << " sig_F = " << sig_F << "\n"
              << " Gamma_F = " << Gamma_F << "\n";
    #endif

    ++iter;
    if (iter > 20) {
      #ifdef DEBUG_SIGMA_FIXED
      std::cout << "Iter: " << iter << " Gamma_F_old = " << Gamma_F_old << " Gamma_F = " << Gamma_F << "\n";
      #endif
      break;
    }

  } while (std::abs(Gamma_F_old - Gamma_F) > 1.0e-10);

  #ifdef DEBUG_SIGMA_FIXED
  std::cout << "computing sig_F: Gamma_F = " << Gamma_F << "\n";
  #endif

  #ifdef NO_RETURN_TO_YIELD_SURFACE

  double df_deps_v = d_yield->df_depsVol(&state_trial_iter, nullptr, nullptr);
  double H_F = - df_deps_v * M_c.Trace() / N_c_mag; 
  return std::make_tuple(sig_F, P_c, M_c, H_F, Gamma_F);

  #else

  std::tie(sig_F, closest, tangent) = closestPointInZRSpace(state_old_for_normal, state_trial_iter,
                                                            z_r_table, z_r_index);
  state_trial_iter.stressTensor = sig_F;
  state_trial_iter.updateStressInvariants();
  state_old_for_normal.closest = closest;
  state_old_for_normal.tangent = tangent;
  Matrix3 N_F = computeYieldSurfaceNormal(state_old_for_normal, sig_F);
  double N_F_mag = N_F.Norm();
  Matrix3 M_F = N_F / N_F_mag;
  Matrix3 P_F = computeProjectionTensor(state_old_for_normal, sig_F, M_F);
  double df_deps_v = d_yield->df_depsVol(&state_trial_iter, nullptr, nullptr);
  double H_F = - df_deps_v * M_F.Trace() / N_F_mag; 

    #ifdef DEBUG_SIGMA_FIXED
    std::cout << "sig_F = " << sig_F << "\n";
    std::cout << "N_F = " << M_F << "\n";
    std::cout << "P_F = " << P_F << "\n";
    std::cout << "df/deps_v = " << df_deps_v << "\n";
    std::cout << "H_F = " << H_F << "\n";
    #endif

  return std::make_tuple(sig_F, P_F, M_F, H_F, Gamma_F);

  #endif
}

/**
 *  Make correction to sigma_F by accounting for hardening
 */
std::tuple<Uintah::Matrix3 , double>
TabularPlasticityCap::computeSigmaHardening(const ModelState_TabularCap& state_trial, 
                                            double Gamma_F, 
                                            const Matrix3& P_F, 
                                            const Matrix3& N_F_norm, 
                                            double H_F) const
{
  // Compute Gamma_H
  double PN = P_F.Contract(N_F_norm);
  double Gamma_H = Gamma_F * PN / (PN + H_F);

  #ifdef DEBUG_SIGMA_HARDENING
  std::cout << "PN = " << PN << " Gamma_H = " << Gamma_H << "\n";
  #endif

  #ifdef TEST_GAMMA_H_SCALING
  Gamma_H *= gamma_H_factor;
  #endif

  // Compute sig_H
  Matrix3 sig_H = state_trial.stressTensor - P_F * Gamma_H;

  #ifdef DEBUG_SIGMA_HARDENING
  std::cout << "sig_H = " << sig_H << "\n";
  #endif

  return std::make_tuple(sig_H, Gamma_H);
}

/**
 * Find closest point from the trial stress to the fixed yield surface in z-r space
 */
std::tuple<Uintah::Matrix3, Uintah::Point, Uintah::Vector>
TabularPlasticityCap::closestPointInZRSpace(const ModelState_TabularCap& state_k_old,
                                            const ModelState_TabularCap& state_k_trial) const
{
  // Compute ratio of bulk and shear moduli
  double K_old = state_k_old.bulkModulus;
  double G_old = state_k_old.shearModulus;

  const double sqrt_K_over_G_old = std::sqrt(1.5 * K_old / G_old);

  // Save the r and z Lode coordinates for the trial stress state
  double r_trial = state_k_trial.rr;
  double z_trial = state_k_trial.zz;

  // Compute transformed r coordinates
  double rprime_trial = r_trial * sqrt_K_over_G_old;

  // Find closest point and tangent at closest point
  double z_closest = 0.0, rprime_closest = 0.0;
  double z_tangent = 0.0, rprime_tangent = 0.0;
  //d_yield->getClosestPoint(
  //  &state_k_old, z_trial, rprime_trial, z_closest, rprime_closest);
  d_yield->getClosestPointAndTangent(
    &state_k_old, z_trial, rprime_trial, z_closest, rprime_closest,
    z_tangent, rprime_tangent);

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0) * z_closest;
  double sqrtJ2_closest =
    1.0 / (sqrt_K_over_G_old * Util::sqrt_two) * rprime_closest;

  #ifdef CHECK_FOR_NAN_EXTRA
    #ifdef DEBUG_WITH_PARTICLE_ID
      if (state_k_trial.particleID == testParticleID) {
    #endif
        std::cout << " K_old = " << K_old << " G_old = " << G_old << "\n";
        std::cout << " state_k_old " << state_k_old << "\n";
        std::cout << " z_trial = " << z_trial
                  << " r_trial = " << rprime_trial / sqrt_K_over_G_old << "\n";
        std::cout << " z_closest = " << z_closest
                  << " r_closest = " << rprime_closest / sqrt_K_over_G_old
                  << "\n";
        std::cout << "I1_closest = " << I1_closest
                  << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
        std::cout << "Trial state = " << state_k_trial << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
      }
    #endif
  #endif

  #ifdef CHECK_HYDROSTATIC_TENSION
    if (I1_closest < 0) {
      std::cout << "I1_closest = " << I1_closest
                << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
      std::cout << "Trial state = " << state_k_trial << "\n";
    }
  #endif

  // Compute closest point stress
  Matrix3 sig_dev = state_k_trial.deviatoricStressTensor;
  Matrix3 sig_closest = Util::one_third * I1_closest * Util::Identity;
  if (state_k_trial.sqrt_J2 > 0.0) {
    sig_closest += (sqrtJ2_closest / state_k_trial.sqrt_J2) * sig_dev;
  } 

  // Save the closest point and tangent in pbar-sqrtJ2 space
  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  double tangent_p_bar = 0.0, tangent_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(sqrt_K_over_G_old,
    z_closest, rprime_closest, closest_p_bar, closest_sqrt_J2);
  Vaango::Util::revertFromZRprime(sqrt_K_over_G_old,
    z_tangent, rprime_tangent, tangent_p_bar, tangent_sqrt_J2);
  Uintah::Point closest(closest_p_bar, closest_sqrt_J2, 0.0);
  Uintah::Vector tangent(tangent_p_bar, tangent_sqrt_J2, 0.0);

  return std::make_tuple(sig_closest, closest, tangent);
}

/**
 * Find closest point from the trial stress to the fixed yield surface in z-r space
 * given a precomputed KD-tree index
 */
std::tuple<Uintah::Matrix3, Uintah::Point, Uintah::Vector>
TabularPlasticityCap::closestPointInZRSpace(const ModelState_TabularCap& state_k_old,
                                            const ModelState_TabularCap& state_k_trial,
                                            const Polyline& z_r_table, 
                                            const Util::PolylineKDTree& z_r_index) const
{
  // Conver trial stress to z-rprime coordinates
  double K_old = state_k_old.bulkModulus;
  double G_old = state_k_old.shearModulus;
  const double sqrt_K_over_G_old = std::sqrt(1.5 * K_old / G_old);

  double z_trial = state_k_trial.zz;
  double r_trial = state_k_trial.rr;
  double rprime_trial = r_trial * sqrt_K_over_G_old;

  // Find closest point and tangent at closest point
  double z_closest = 0.0, rprime_closest = 0.0;
  double z_tangent = 0.0, rprime_tangent = 0.0;
  d_yield->getClosestPointAndTangent(
    &state_k_old, z_r_table, z_r_index, z_trial, rprime_trial, z_closest, rprime_closest,
    z_tangent, rprime_tangent);

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0) * z_closest;
  double sqrtJ2_closest =
    1.0 / (sqrt_K_over_G_old * Util::sqrt_two) * rprime_closest;

  #ifdef CHECK_FOR_NAN_EXTRA
    #ifdef DEBUG_WITH_PARTICLE_ID
      if (state_k_trial.particleID == testParticleID) {
    #endif
        std::cout << " K_old = " << K_old << " G_old = " << G_old << "\n";
        std::cout << " state_k_old " << state_k_old << "\n";
        std::cout << " z_trial = " << z_trial
                  << " r_trial = " << rprime_trial / sqrt_K_over_G_old << "\n";
        std::cout << " z_closest = " << z_closest
                  << " r_closest = " << rprime_closest / sqrt_K_over_G_old
                  << "\n";
        std::cout << "I1_closest = " << I1_closest
                  << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
        std::cout << "Trial state = " << state_k_trial << "\n";
    #ifdef DEBUG_WITH_PARTICLE_ID
      }
    #endif
  #endif

  #ifdef CHECK_HYDROSTATIC_TENSION
    if (I1_closest < 0) {
      std::cout << "I1_closest = " << I1_closest
                << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
      std::cout << "Trial state = " << state_k_trial << "\n";
    }
  #endif

  // Compute closest point stress
  Matrix3 sig_dev = state_k_trial.deviatoricStressTensor;
  Matrix3 sig_closest = Util::one_third * I1_closest * Util::Identity;
  if (state_k_trial.sqrt_J2 > 0.0) {
    sig_closest += (sqrtJ2_closest / state_k_trial.sqrt_J2) * sig_dev;
  } 

  // Save the closest point and tangent in pbar-sqrtJ2 space
  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  double tangent_p_bar = 0.0, tangent_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(sqrt_K_over_G_old,
    z_closest, rprime_closest, closest_p_bar, closest_sqrt_J2);
  Vaango::Util::revertFromZRprime(sqrt_K_over_G_old,
    z_tangent, rprime_tangent, tangent_p_bar, tangent_sqrt_J2);
  Uintah::Point closest(closest_p_bar, closest_sqrt_J2, 0.0);
  Uintah::Vector tangent(tangent_p_bar, tangent_sqrt_J2, 0.0);

  return std::make_tuple(sig_closest, closest, tangent);
}

/**
 * Find yield surface normal at the closest point
 */
Uintah::Matrix3 
TabularPlasticityCap::computeYieldSurfaceNormal(const ModelState_TabularCap& state_k_old,
                                                const Matrix3& sig_closest) const
{
  ModelState_TabularCap state(state_k_old);
  state.stressTensor = sig_closest;
  state.updateStressInvariants();

  // This version of df_dsigma assumes that the closest point and tangent
  // have already been computed and stored in the state.  The first argument
  // is a dummy that tells the code that this special version is being used.
  Matrix3 df_dsigma = d_yield->df_dsigma(Vaango::Util::Identity, &state);

  return df_dsigma;
}

/**
 * Compute projection tensor at the closest point
 */
Uintah::Matrix3 
TabularPlasticityCap::computeProjectionTensor(const ModelState_TabularCap& state_k_old, 
                                              const Matrix3& sig_closest,
                                              const Matrix3& df_dsig_closest) const
{
  ModelState_TabularCap state(state_k_old);
  state.stressTensor = sig_closest;
  state.updateStressInvariants();

  // Compute C:M
  double lambda = state.bulkModulus - 2.0 / 3.0 * state.shearModulus;
  double mu     = state.shearModulus;
  double trM    = df_dsig_closest.Trace();
  Matrix3 CM    = Util::Identity * (trM * lambda) + 
                  df_dsig_closest * (2.0 * mu);

  // Compute P = C:M (without elasticplastic coupling)
  Matrix3 P = CM;

  #ifdef DO_ELASTIC_PLASTIC_COUPLING
  std::cout << "P = CM = " << P << "\n ";

  // Compute Z (elastic-plastic coupling term)
  auto K_dK = d_elastic->getElasticModuliAndDerivatives(&state);
  auto dK_deps_v_p = K_dK.second.bulkModulus;
  auto dG_deps_v_p = K_dK.second.shearModulus;
  Matrix3 Z = Util::Identity * (state.I1 / 3.0) * (dK_deps_v_p / state.bulkModulus * trM) + 
              state.deviatoricStressTensor * (dG_deps_v_p / state.shearModulus * trM);
  Z *= (-1.0);
  std::cout << "Z = " << Z << "\n"
            << "  p = " << state.I1/3.0 
            << "  dK_eps_v_p = " << dK_deps_v_p 
            << "  K = " << state.bulkModulus << "\n";

  // Compute P = C:M + Z (with elasticplastic coupling)
  P += Z;
  std::cout << "P = CM + Z = " << P << "\n ";
  #endif

  return P; 
}

/**
 * Compute plastic strain increment (Gamma = lambda_{n+1} - \lambda_n)
 */
double 
TabularPlasticityCap::computeGammaClosest(const Matrix3& sig_c, 
                                          const Matrix3& sig_trial, 
                                          const Matrix3& N_c, 
                                          const Matrix3& P_c) const
{
  auto sig_diff = sig_trial - sig_c;
  auto numerator = sig_diff.Contract(N_c);
  auto denominator = P_c.Contract(N_c);
  auto Gamma_F = numerator / denominator;

  #ifdef CHECK_GAMMA_F_POSITIVE
  if (Gamma_F < 0) {
    std::cout << "Gamma_F = " << Gamma_F << " numer = " << numerator << " denom = " << denominator << "\n";
    std::cout << "sig_diff = " << sig_diff << "\n"
              << " sig_trial = " << sig_trial << "\n"
              << " sig_c = " << sig_c << "\n"
              << " P_c = " << P_c << "\n"
              << " N_c = " << N_c << "\n";
  }
  #endif

  return Gamma_F;
}

/**
   * Method: Compute the stress state at the end of a purely elastic
   *         substep.
   * Purpose:
   *   Find the updated stress before an elastic-plastic non-hardening
   *   stress update step
   *
   * Returns:
   *   Status  :  whether the procedure is sucessful or has failed
   *   double  :  elastic delta T
   *
*/
std::tuple<TabularPlasticity::Status, double>
TabularPlasticityCap::computePurelyElasticSubstep(
  double dt,
  const Matrix3& D,
  const ModelState_TabularCap& state_k_old,
  const ModelState_TabularCap& state_k_trial,
  ModelState_TabularCap& state_k_elastic)
{
  // Compute the elastic fraction of the time step
  bool doesIntersect;
  double t_val;
  double elastic_dt;
  Uintah::Point intersection;
  std::tie(doesIntersect, t_val, elastic_dt, intersection) =
    computeElasticDeltaT(dt, state_k_old, state_k_trial);
  std::cout << "0: doesIntersect = " << std::boolalpha << doesIntersect
            << " t_val = " << t_val << " elastic_dt = " << elastic_dt
            << " intersection = " << intersection << "\n";

  // If it doesn't intersect then the trial state is either purely
  // elastic or the initial state is just outside the yield surface.
  // The purely elastic state has been identified earlier using
  // evalYieldFunction; so the state must be on/outside the yield surface
  if (!doesIntersect) {
    elastic_dt = 0.0;
  }

  // Compute the elastic strain at the end of the substep (total and increments)
  Matrix3 deltaEps         = D * dt;
  Matrix3 elastic_deltaEps = D * elastic_dt;
  Matrix3 elasticStrain_new =
    state_k_old.elasticStrainTensor + elastic_deltaEps;
  Matrix3 plasticStrain_new = state_k_old.plasticStrainTensor;

  // std::cout << "D = " << D << "\n";
  // std::cout << "deltaEps = " << deltaEps << "\n";
  // std::cout << "delta Eps_e = " << elastic_deltaEps << "\n";
  // std::cout << "Eps_e = " << elasticStrain_new << "\n";

  // Compute stress at the end of the elastic part of the time step
  Matrix3 stress_new = computeTrialStress(state_k_old, elastic_deltaEps);

  // Compute the elastic moduli
  double K_new, G_new;
  std::tie(K_new, G_new) =
    computeElasticProperties(state_k_old, elasticStrain_new, plasticStrain_new);

  // Update the state at the end of the elastic step (no change in plastic
  // strain)
  state_k_elastic              = state_k_old;
  state_k_elastic.stressTensor = stress_new;
  state_k_elastic.updateStressInvariants();
  state_k_elastic.elasticStrainTensor = elasticStrain_new;
  state_k_elastic.bulkModulus         = K_new;
  state_k_elastic.shearModulus        = G_new;

  return std::make_tuple(Status::SUCCESS, elastic_dt);
}

/**
 * Method: nonHardeningReturnElasticPlastic
 * Purpose:
 *   Computes a non-hardening return to the yield surface in the meridional
 *   profile (constant Lode angle) based on the current values of the
 *   internal state variables and elastic properties.  Returns the updated
 *   stress.
 *   NOTE: all values of r and z in this function are transformed!
 * Inputs:
 *   elasticplastic_dt = timestep size after purely elastic update
 *   D                 = rate of deformation
 *   state_k_elastic   = state at the end of the purely elastic substep
 * Outputs:
 *   state_k_nonhardening = state at the end of the non-hardening
 *                          elastic-plastic substep
 * Returns:
 *   Status: true  = success
 *           false = failure
 *   Matrix3: plastic strain increment
 */
std::tuple<TabularPlasticity::Status, Matrix3>
TabularPlasticityCap::nonHardeningReturnElasticPlastic(
  double elasticplastic_dt,
  const Matrix3& D,
  const ModelState_TabularCap& state_k_elastic,
  ModelState_TabularCap& state_k_nonhardening)
{
  Status status = Status::SUCCESS;

  Matrix3 deltaEps = D * elasticplastic_dt;
  // std::cout << "D = " << D << "\n";
  // std::cout << "deltaEps = " << deltaEps << "\n";

  // Compute new trial stress
  Matrix3 stress_trial = computeTrialStress(state_k_elastic, deltaEps);

  // Compute ratio of bulk and shear moduli
  double K_elastic                   = state_k_elastic.bulkModulus;
  double G_elastic                   = state_k_elastic.shearModulus;
  const double sqrt_K_over_G_elastic = std::sqrt(1.5 * K_elastic / G_elastic);

  // Compute the r and z Lode coordinates for the trial stress state
  double I1_trial          = stress_trial.Trace();
  Matrix3 dev_stress_trial = stress_trial - Util::Identity * (I1_trial / 3.0);
  double J2_trial          = 0.5 * dev_stress_trial.Contract(dev_stress_trial);
  J2_trial =
    (J2_trial < 1e-16 * (I1_trial * I1_trial + J2_trial)) ? 0.0 : J2_trial;
  double sqrtJ2_trial = std::sqrt(J2_trial);
  double r_trial      = Util::sqrt_two * sqrtJ2_trial;
  double z_trial      = I1_trial / Util::sqrt_three;

  // Compute transformed r coordinates
  double rprime_trial = r_trial * sqrt_K_over_G_elastic;

  // Find closest point
  double z_closest = 0.0, rprime_closest = 0.0;
  d_yield->getClosestPoint(
    &state_k_elastic, z_trial, rprime_trial, z_closest, rprime_closest);

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0) * z_closest;
  double sqrtJ2_closest =
    1.0 / (sqrt_K_over_G_elastic * Util::sqrt_two) * rprime_closest;

#ifdef CHECK_FOR_NAN_EXTRA
#ifdef DEBUG_WITH_PARTICLE_ID
  if (state_k_elastic.particleID == testParticleID) {
#endif
    std::cout << " K_elastic = " << K_elastic << " G_elastic = " << G_elastic
              << "\n";
    std::cout << " z_trial = " << z_trial
              << " r_trial = " << rprime_trial / sqrt_K_over_G_elastic
              << "\n";
    std::cout << " z_closest = " << z_closest
              << " r_closest = " << rprime_closest / sqrt_K_over_G_elastic
              << "\n";
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "Trial state = " << state_k_elastic << "\n";
#ifdef DEBUG_WITH_PARTICLE_ID
  }
#endif
#endif

#ifdef CHECK_HYDROSTATIC_TENSION
  if (I1_closest < 0) {
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "Trial state = " << state_k_trial << "\n";
  }
#endif

  // Compute new stress (on the non-hardening yield surface)
  Matrix3 sig_fixed(0.0);
  if (J2_trial > 0.0) {
    sig_fixed = Util::one_third * I1_closest * Util::Identity +
                (sqrtJ2_closest / sqrtJ2_trial) * dev_stress_trial;
  } else {
    sig_fixed =
      Util::one_third * I1_closest * Util::Identity + dev_stress_trial;
  }
#ifdef DEBUG_WITH_PARTICLE_ID
  if (state_k_elastic.particleID == testParticleID) {
    std::cout << "sig_n = " << state_k_elastic.stressTensor
              << " sig_n+1 = " << sig_fixed << "\n";
  }
#endif

  // std::cout << "G_elastic = " << G_elastic << "\n";
  // std::cout << "I1_closest = " << I1_closest
  //          << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
  // std::cout << "Sigma_e = " << state_k_elastic.stressTensor << "\n";
  // std::cout << "Sigma_i = " << sig_fixed << "\n";

  // Compute new elastic and plastic strain increments
  //  d_ep = d_e - [C]^-1:(sigma_new-sigma_old)
  Matrix3 deltaSig     = sig_fixed - state_k_elastic.stressTensor;
  Matrix3 deltaSig_iso = Util::one_third * deltaSig.Trace() * Util::Identity;
  Matrix3 deltaSig_dev = deltaSig - deltaSig_iso;
  Matrix3 deltaEps_e   = deltaSig_iso * (Util::one_third / K_elastic) +
                       deltaSig_dev * (0.5 / G_elastic);
  Matrix3 deltaEps_p = deltaEps - deltaEps_e;

  // Update the state
  state_k_nonhardening              = state_k_elastic;
  state_k_nonhardening.stressTensor = sig_fixed;
  state_k_nonhardening.updateStressInvariants();
  state_k_nonhardening.elasticStrainTensor =
    state_k_elastic.elasticStrainTensor + deltaEps_e;
  state_k_nonhardening.plasticStrainTensor =
    state_k_elastic.plasticStrainTensor + deltaEps_p;
  state_k_nonhardening.updatePlasticStrainInvariants();

  /*
  if (deltaEps_p.NormSquared() > 5*deltaEps.NormSquared()) {
    status = Status::TOO_LARGE_YIELD_NORMAL_CHANGE;
  }
  */

  return std::make_tuple(status, deltaEps_p);
}

/**
 * Method: consistencyBisectionHardeningSoftening
 * Purpose:
 *   Find the updated stress for hardening/softening plasticity using the
 * consistency
 *   bisection algorithm
 *
 *   Returns - whether the procedure is sucessful or has failed
 */
TabularPlasticity::Status
TabularPlasticityCap::consistencyBisectionHardeningSoftening(
  double elastic_plastic_dt,
  const Matrix3& D,
  const ModelState_TabularCap& state_k_elastic,
  const ModelState_TabularCap& state_k_nonhardening,
  const Matrix3& deltaEps_p_nonhardening,
  ModelState_TabularCap& state_new)
{
  // Compute the strain increment (needed to recompute the trial stress)
  Matrix3 deltaEps    = D * elastic_plastic_dt;
  double deltaEps_vol = deltaEps.Trace();

  // Compute the volumetric part of the strain increment and the
  // non-hardening plastic volumetric strain increment
  double deltaEps_p_vol_nonhardening = deltaEps_p_nonhardening.Trace();
  double deltaEps_p_vol_min, deltaEps_p_vol_max;
  double eta_min = 0.0, eta_max = 1.0, eta_mid = 0.5;
  if (deltaEps_p_vol_nonhardening > 0) {
    deltaEps_p_vol_min = -deltaEps_p_vol_nonhardening;
    deltaEps_p_vol_max = deltaEps_p_vol_nonhardening;
  } else {
    deltaEps_p_vol_min = deltaEps_p_vol_nonhardening;
    deltaEps_p_vol_max = -deltaEps_p_vol_nonhardening;
  }

  // Compute the elastic moduli at the end of the nonhardening return state
  // (Since the total strain is used, this is a reasonable end-of-timestep
  //  estimate)
  double K_new, G_new;
  std::tie(K_new, G_new) = computeElasticProperties(state_k_nonhardening);

  // Create a new state
  ModelState_TabularCap state_k_updated(state_k_nonhardening);
  state_k_updated.bulkModulus  = K_new;
  state_k_updated.shearModulus = G_new;

  // Do bisection
  while (std::abs(eta_min - eta_max) > 1.0e-6) {
    // Compute the internal variables at the end of the nonhardening return
    // state
    double X_fixed       = computeInternalVariable(state_k_updated);
    state_k_updated.capX = X_fixed;

    // Update yield surface
    auto yield_pts_updated =
      d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_k_updated);
    state_k_updated.yield_f_pts = yield_pts_updated;

    // Eval the yield condition with the updated yield surface and the stress at
    // the end of the non-hardening update
    auto yield = d_yield->evalYieldCondition(&state_k_updated);

    if (yield.second == Util::YieldStatus::IS_ELASTIC) {
      // If the non-hardening state is on or inside the updated yield surface
      // reduce the plastic volumetric strain increment
      eta_max = eta_mid;
      eta_mid = 0.5 * (eta_min + eta_max);
    } else {
      // If the non-hardening state is outside the updated yield surface
      // increase the plastic volumetric strain increment
      eta_min = eta_mid;
      eta_mid = 0.5 * (eta_min + eta_max);
    }

    // Update the increment of plastic volumetric strain
    double deltaEps_p_vol_mid =
      (1 - eta_mid) * deltaEps_p_vol_min + eta_mid * deltaEps_p_vol_max;
    double deltaEps_e_vol_mid = deltaEps_vol - deltaEps_p_vol_mid;

    // Create a new trial state with the updated moduli
    ModelState_TabularCap state_trial_mid(state_k_elastic);
    state_trial_mid.bulkModulus  = K_new;
    state_trial_mid.shearModulus = G_new;

    // Compute the internal variables at the end of the nonhardening return
    // state
    auto status = computeInternalVariables(
      state_trial_mid, deltaEps_e_vol_mid, deltaEps_p_vol_mid);

    // Update yield surface
    auto yield_pts_mid =
      d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_trial_mid);
    state_trial_mid.yield_f_pts = yield_pts_mid;

    // Do a non-hardening return to the updated yield surface
    Matrix3 deltaEps_p_mid;
    std::tie(status, deltaEps_p_mid) = nonHardeningReturnElasticPlastic(
      elastic_plastic_dt, D, state_trial_mid, state_k_updated);
    if (std::abs(deltaEps_p_mid.Trace() - deltaEps_p_vol_mid) < 1.0e-6) {
      break;
    }
  }

  state_new = state_k_updated;

  return Status::SUCCESS;
}

/**
 * Method: nonHardeningReturn
 * Purpose:
 *   Computes a non-hardening return to the yield surface in the meridional
 *   profile (constant Lode angle) based on the current values of the
 *   internal state variables and elastic properties.  Returns the
 *   updated stress and  the increment in * plastic
 *   strain corresponding to this return.
 *
 *   NOTE: all values of r and z in this function are transformed!
 */
TabularPlasticity::Status
TabularPlasticityCap::nonHardeningReturn(
  const Uintah::Matrix3& strain_inc,
  const ModelState_TabularCap& state_k_old,
  const ModelState_TabularCap& state_k_trial,
  Uintah::Matrix3& sig_fixed,
  Uintah::Matrix3& elasticStrain_inc_fixed,
  Uintah::Matrix3& plasticStrain_inc_fixed)
{
  Status status = Status::SUCCESS;

  // Compute ratio of bulk and shear moduli
  double K_old                   = state_k_old.bulkModulus;
  double G_old                   = state_k_old.shearModulus;
  const double sqrt_K_over_G_old = std::sqrt(1.5 * K_old / G_old);

  // Save the r and z Lode coordinates for the trial stress state
  double r_trial = state_k_trial.rr;
  double z_trial = state_k_trial.zz;

  // Compute transformed r coordinates
  double rprime_trial = r_trial * sqrt_K_over_G_old;

  // Find closest point
  double z_closest = 0.0, rprime_closest = 0.0;
  d_yield->getClosestPoint(
    &state_k_old, z_trial, rprime_trial, z_closest, rprime_closest);

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0) * z_closest;
  double sqrtJ2_closest =
    1.0 / (sqrt_K_over_G_old * Util::sqrt_two) * rprime_closest;

#ifdef CHECK_FOR_NAN_EXTRA
#ifdef DEBUG_WITH_PARTICLE_ID
  if (state_k_trial.particleID == testParticleID) {
#endif
    std::cout << " K_old = " << K_old << " G_old = " << G_old << "\n";
    std::cout << " state_k_old " << state_k_old << "\n";
    std::cout << " z_trial = " << z_trial
              << " r_trial = " << rprime_trial / sqrt_K_over_G_old << "\n";
    std::cout << " z_closest = " << z_closest
              << " r_closest = " << rprime_closest / sqrt_K_over_G_old
              << "\n";
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "Trial state = " << state_k_trial << "\n";
#ifdef DEBUG_WITH_PARTICLE_ID
  }
#endif
#endif

#ifdef CHECK_HYDROSTATIC_TENSION
  if (I1_closest < 0) {
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "Trial state = " << state_k_trial << "\n";
  }
#endif

  // Compute new stress
  Matrix3 sig_dev = state_k_trial.deviatoricStressTensor;
  if (state_k_trial.sqrt_J2 > 0.0) {
    sig_fixed = Util::one_third * I1_closest * Util::Identity +
                (sqrtJ2_closest / state_k_trial.sqrt_J2) * sig_dev;
  } else {
    sig_fixed = Util::one_third * I1_closest * Util::Identity;
  }

  // Compute new plastic strain increment
  //  d_ep = d_e - [C]^-1:(sigma_new-sigma_old)
  Matrix3 sig_inc     = sig_fixed - state_k_old.stressTensor;
  Matrix3 sig_inc_iso = Util::one_third * sig_inc.Trace() * Util::Identity;
  Matrix3 sig_inc_dev = sig_inc - sig_inc_iso;
  elasticStrain_inc_fixed =
    sig_inc_iso * (Util::one_third / K_old) + sig_inc_dev * (0.5 / G_old);
  plasticStrain_inc_fixed = strain_inc - elasticStrain_inc_fixed;

#ifdef CHECK_ELASTIC_STRAIN
  // std::cout << "Non-hardening:\n"
  //          << "\t Delta sig = " << sig_inc << "\n"
  //          << "\t Delta Eps_e = " << elasticStrain_inc_fixed << "\n";
  std::cout
    << "press = " << sig_fixed.Trace() / 3.0 << " K = " << K_old << " ev_e = "
    << (state_k_old.elasticStrainTensor + elasticStrain_inc_fixed).Trace()
    << "\n";
#endif

  // Compute volumetric plastic strain and compare with p3
  Matrix3 eps_p = state_k_old.plasticStrainTensor + plasticStrain_inc_fixed;
  double ep_v   = eps_p.Trace();

  if (ep_v < 0.0) {
    if (-ep_v > 10) {
      proc0cout << "**WARNING** Nonhardening return has failed because "
                << " epsbar_p_v > 10 : " << -ep_v << " > 10 " << "\n";
      proc0cout << " K_old = " << K_old << " G_old = " << G_old << "\n";
      proc0cout << " state_k_trial " << state_k_trial << "\n";
      proc0cout << " z_trial = " << z_trial
                << " r_trial = " << rprime_trial / sqrt_K_over_G_old
                << "\n";
      proc0cout << " z_closest = " << z_closest
                << " r_closest = " << rprime_closest / sqrt_K_over_G_old
                << "\n";
      proc0cout << "Delta eps = " << strain_inc << "\n";
      proc0cout << "sig_n = " << state_k_old.stressTensor << "\n";
      proc0cout << "sig_n+1 = " << sig_fixed << "\n";
      proc0cout << "Delta sig = " << sig_inc << "\n";
      proc0cout << "Delta sig_iso = " << sig_inc_iso << "\n";
      proc0cout << "Delta sig_dev = " << sig_inc_dev << "\n";
      proc0cout << "Delta eps_e = " << elasticStrain_inc_fixed << "\n";
      proc0cout << "Delta eps_p = " << plasticStrain_inc_fixed << "\n";
      proc0cout << "I1_J2_trial = [" << state_k_trial.I1 << " "
                << state_k_trial.sqrt_J2 << "];" << "\n";
      proc0cout << "I1_J2_closest = [" << I1_closest << " " << sqrtJ2_closest
                << "];" << "\n";
      proc0cout << "plot([I1 I1_J2_closest(1)],[sqrtJ2 I1_J2_closest(2)],'gx')"
                << ";" << "\n";
      proc0cout << "plot([I1_J2_trial(1) I1_J2_closest(1)],[I1_J2_trial(2) "
                   "I1_J2_closest(2)],'r-')"
                << ";" << "\n";

      status = Status::TOO_LARGE_PLASTIC_STRAIN; // The plastic volume strain is
                                                 // too large, try again
    }
  }

#ifdef CHECK_PLASTIC_RATE
  ModelState_TabularCap state_plastic_rate(state_k_old);
  state_plastic_rate.stressTensor = sig_fixed;
  state_plastic_rate.updateStressInvariants();
  Matrix3 df_dsigma = d_yield->df_dsigma(Util::Identity, &state_plastic_rate);
  df_dsigma /= df_dsigma.Norm();
  double lhs          = plasticStrain_inc_fixed.Contract(df_dsigma);
  double rhs          = df_dsigma.Contract(df_dsigma);
  double plastic_rate = lhs / rhs;
  if (plastic_rate < 0) {
    std::cout << "Particle = " << state_k_old.particleID << "\n";
    std::cout << "Delta eps = " << strain_inc << "\n";
    std::cout << "Trial state = " << state_k_trial << "\n";
    std::cout << "Delta sig = " << sig_inc << "\n";
    std::cout << "Delta sig_iso = " << sig_inc_iso << "\n";
    std::cout << "Delta sig_dev = " << sig_inc_dev << "\n";
    std::cout << "Delta eps_e = " << elasticStrain_inc_fixed << "\n";
    std::cout << "Delta eps_p = " << plasticStrain_inc_fixed << "\n";
    std::cout << "df_dsigma = " << df_dsigma << "\n";
    std::cout << "plastic rate = " << plastic_rate << "\n";
  }
#endif

#ifdef CHECK_YIELD_SURFACE_NORMAL
  std::cout << "Delta eps = " << strain_inc << "\n";
  std::cout << "Trial state = " << state_k_trial << "\n";
  std::cout << "Delta sig = " << sig_inc << "\n";
  std::cout << "Delta sig_iso = " << sig_inc_iso << "\n";
  std::cout << "Delta sig_dev = " << sig_inc_dev << "\n";
  std::cout << "Delta eps_e = " << elasticStrain_inc_fixed << "\n";
  std::cout << "Delta eps_p = " << plasticStrain_inc_fixed << "\n";

  // Test normal to yield surface
  ModelState_TabularCap state_test(state_k_old);
  state_test.stressTensor = sig_fixed;
  state_test.updateStressInvariants();

  Matrix3 df_dsigma = d_yield->df_dsigma(Util::Identity, &state_test);
  df_dsigma /= df_dsigma.Norm();
  std::cout << "df_dsigma = " << df_dsigma << "\n";
  std::cout << "ratio = [" << plasticStrain_inc_fixed(0, 0) / df_dsigma(0, 0)
            << "," << plasticStrain_inc_fixed(1, 1) / df_dsigma(1, 1) << ","
            << plasticStrain_inc_fixed(2, 2) / df_dsigma(2, 2) << "\n";

  // Compute CN = C:df_dsigma
  double lambda = state_test.bulkModulus - 2.0 / 3.0 * state_test.shearModulus;
  double mu     = state_test.shearModulus;
  Matrix3 CN =
    Util::Identity * (lambda * df_dsigma.Trace()) + df_dsigma * (2.0 * mu);
  Matrix3 sig_diff = state_k_trial.stressTensor - sig_fixed;
  std::cout << "sig_trial = [" << state_k_trial.stressTensor << "];"
            << "\n";
  std::cout << "sig_n+1 = [" << sig_fixed << "];" << "\n";
  std::cout << "sig_trial - sig_n+1 = [" << sig_diff << "];" << "\n";
  std::cout << "C_df_dsigma = [" << CN << "];" << "\n";
  std::cout << "sig ratio = [" << sig_diff(0, 0) / CN(0, 0) << " "
            << sig_diff(1, 1) / CN(1, 1) << " " << sig_diff(2, 2) / CN(2, 2)
            << "];" << "\n";

  // Compute a test stress to check normal
  Matrix3 sig_test = sig_fixed + df_dsigma * sig_diff(0, 0);
  ModelState_TabularCap state_sig_test(state_k_old);
  state_sig_test.stressTensor = sig_test;
  state_sig_test.updateStressInvariants();
  std::cout << "I1 = " << state_sig_test.I1 << ";" << "\n";
  std::cout << "sqrtJ2 = " << state_sig_test.sqrt_J2 << ";" << "\n";
  std::cout << "I1_J2_trial = [" << state_k_trial.I1 << " "
            << state_k_trial.sqrt_J2 << "];" << "\n";
  std::cout << "I1_J2_closest = [" << I1_closest << " " << sqrtJ2_closest
            << "];" << "\n";
  std::cout << "plot([I1 I1_J2_closest(1)],[sqrtJ2 I1_J2_closest(2)],'gx-')"
            << ";" << "\n";
  std::cout << "plot([I1_J2_trial(1) I1_J2_closest(1)],[I1_J2_trial(2) "
               "I1_J2_closest(2)],'r-')"
            << ";" << "\n";

  // Check actual location of projected point
  Matrix3 sig_test_actual =
    state_k_trial.stressTensor - CN * (std::abs(sig_diff(0, 0) / CN(0, 0)));
  state_sig_test.stressTensor = sig_test_actual;
  state_sig_test.updateStressInvariants();
  std::cout << "I1 = " << state_sig_test.I1 << ";" << "\n";
  std::cout << "sqrtJ2 = " << state_sig_test.sqrt_J2 << ";" << "\n";
  std::cout << "plot([I1 I1_J2_trial(1)],[sqrtJ2 I1_J2_trial(2)],'rx')"
            << ";" << "\n";
#endif

#ifdef CHECK_FOR_NAN
  if (std::isnan(sig_fixed(0, 0))) {
    std::cout << " K_old = " << K_old << " G_old = " << G_old << "\n";
    std::cout << " z_trial = " << z_trial
              << " r_trial = " << rprime_trial / sqrt_K_over_G_old << "\n";
    std::cout << " z_closest = " << z_closest
              << " r_closest = " << rprime_closest / sqrt_K_over_G_old
              << "\n";
    std::cout << "I1_closest = " << I1_closest
              << " sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "Trial state = " << state_k_trial << "\n";
    std::cout << "\t\t\t sig_fixed = " << sig_fixed << "\n";
    std::cout << "\t\t\t I1_closest = " << I1_closest << "\n";
    std::cout << "\t\t\t sqrtJ2_closest = " << sqrtJ2_closest << "\n";
    std::cout << "\t\t\t state_k_trial.sqrt_J2 = " << state_k_trial.sqrt_J2
              << "\n";
    std::cout << "\t\t\t sig_dev = " << sig_dev << "\n";
    std::cout << "\t\t\t sig_inc = " << sig_inc << "\n";
    std::cout << "\t\t\t strain_inc = " << strain_inc << "\n";
    std::cout << "\t\t\t sig_inc_iso = " << sig_inc_iso << "\n";
    std::cout << "\t\t\t sig_inc_dev = " << sig_inc_dev << "\n";
    std::cout << "\t\t\t plasticStrain_inc_fixed = " << plasticStrain_inc_fixed
              << "\n";
  }
#endif

#ifdef CHECK_HYDROSTATIC_TENSION
  if (I1_closest < 0) {
    std::cout << "\t\t\t sig_inc = " << sig_inc << "\n";
    std::cout << "\t\t\t strain_inc = " << strain_inc << "\n";
    std::cout << "\t\t\t sig_inc_iso = " << sig_inc_iso << "\n";
    std::cout << "\t\t\t sig_inc_dev = " << sig_inc_dev << "\n";
    std::cout << "\t\t\t plasticStrain_inc_fixed = " << plasticStrain_inc_fixed
              << "\n";
  }
#endif

  return status; // isSuccess = true

} //===================================================================

/**
 * Method: firstOrderHardeningUpdate
 * Purpose:
 *   Find the updated stress for hardening plasticity using a first-order
 *   update based on the velocity of the yield surface
 *
 *   Returns whether the procedure is sucessful or has failed
 */
TabularPlasticity::Status
TabularPlasticityCap::firstOrderHardeningUpdate(
  const Matrix3& deltaEps_new,
  const ModelState_TabularCap& state_k_old,
  const ModelState_TabularCap& state_k_trial,
  const Matrix3& deltaEps_e_fixed,
  const Matrix3& deltaEps_p_fixed,
  const Matrix3& sig_fixed,
  ModelState_TabularCap& state_k_new)
{
  Status status = Status::SUCCESS;

  ModelState_TabularCap state_n(state_k_old);
  ModelState_TabularCap state_np1(state_k_old);
  state_np1.stressTensor = sig_fixed;
  state_np1.updateStressInvariants();

  Matrix3 M_n   = d_yield->df_dsigma(Util::Identity, &state_n);
  Matrix3 M_np1 = d_yield->df_dsigma(Util::Identity, &state_np1);
  M_n /= M_n.Norm();
  double norm_M_np1 = M_np1.Norm();
  M_np1 /= norm_M_np1;
  double angle_M_n_np1 = std::abs(M_n.Contract(M_np1) - 1.0);

  if (angle_M_n_np1 > 1.0e-6) {
    status = Status::TOO_LARGE_YIELD_NORMAL_CHANGE;

#ifdef CHECK_RETURN_ALIGNMENT
    std::cout << "M_n = " << M_n << "\n";
    std::cout << "M_np1 = " << M_np1 << "\n";
    std::cout << "angle =" << angle_M_n_np1 << "\n";

    auto KG_dKdG_n   = d_elastic->getElasticModuliAndDerivatives(&state_n);
    auto KG_dKdG_np1 = d_elastic->getElasticModuliAndDerivatives(&state_np1);
    auto KG_n        = KG_dKdG_n.first;
    auto KG_np1      = KG_dKdG_np1.first;
    auto dKdG_n      = KG_dKdG_n.second;
    auto dKdG_np1    = KG_dKdG_np1.second;
    auto p_n         = state_n.I1 / 3.0;
    auto p_np1       = state_np1.I1 / 3.0;
    auto S_n         = state_n.deviatoricStressTensor;
    auto S_np1       = state_np1.deviatoricStressTensor;
    Matrix3 Z_n      = computeZMatrix(KG_n, dKdG_n, p_n, S_n, M_n);
    Matrix3 Z_np1    = computeZMatrix(KG_np1, dKdG_np1, p_np1, S_np1, M_np1);
    std::cout << "Z_n = " << Z_n << "\n";
    std::cout << "Z_np1 = " << Z_np1 << "\n";

    // Compute CM = C:M
    double mu_n       = KG_n.shearModulus;
    double mu_np1     = KG_np1.shearModulus;
    double lambda_n   = KG_n.bulkModulus - 2.0 / 3.0 * mu_n;
    double lambda_np1 = KG_np1.bulkModulus - 2.0 / 3.0 * mu_np1;
    Matrix3 CM_n =
      Util::Identity * (lambda_n * M_n.Trace()) + M_n * (2.0 * mu_n);
    Matrix3 CM_np1 =
      Util::Identity * (lambda_np1 * M_np1.Trace()) + M_np1 * (2.0 * mu_np1);
    std::cout << "CM_n = " << CM_n << "\n";
    std::cout << "CM_np1 = " << CM_np1 << "\n";

    // Compute P = CM + Z;
    Matrix3 P_n   = CM_n + Z_n;
    Matrix3 P_np1 = CM_np1 + Z_np1;
    std::cout << "P_n = " << P_n << "\n";
    std::cout << "P_np1 = " << P_np1 << "\n";

    // Compute Gamma (sig_trial - sig_n+1):M_n+1/(P_n:M_n+1)
    Matrix3 sig_trial = state_k_trial.stressTensor;
    Matrix3 sig_new   = sig_fixed;
    double Gamma_n =
      (sig_trial - sig_new).Contract(M_np1) / P_n.Contract(M_np1);
    double Gamma_np1 =
      (sig_trial - sig_new).Contract(M_np1) / P_np1.Contract(M_np1);
    std::cout << "Gamma_n = " << Gamma_n << " Gamma_np1 = " << Gamma_np1
              << "\n";

    // Compute sigma
    Matrix3 sig_np1 = sig_trial - P_n * Gamma_n;
    std::cout << "sig_np1 = " << sig_np1 << "\n";
    std::cout << "sig_fixed = " << sig_fixed << "\n";

    // Compute volumetric plastic strain
    Matrix3 eps_p = state_n.plasticStrainTensor + deltaEps_p_fixed;
    double ep_v   = eps_p.Trace();

    // Compute eps_p and ev_v^p
    Matrix3 Eps_p_n    = state_n.plasticStrainTensor + M_n * Gamma_n;
    double Eps_p_v_n   = state_n.ep_v + M_n.Trace() * Gamma_n;
    Matrix3 Eps_e_n    = deltaEps_new - Eps_p_n;
    Matrix3 Eps_p_np1  = state_n.plasticStrainTensor + M_np1 * Gamma_np1;
    double Eps_p_v_np1 = state_n.ep_v + M_np1.Trace() * Gamma_np1;
    Matrix3 Eps_e_np1  = deltaEps_new - Eps_p_np1;
    std::cout << "Eps_p_n = " << Eps_p_n << "\n";
    std::cout << "eps_p = " << eps_p << "\n";
    std::cout << "Eps_e_n = " << Eps_e_n << "\n";
    std::cout << "Eps_p_np1 = " << Eps_p_np1 << "\n";
    std::cout << "Eps_e_np1 = " << Eps_e_np1 << "\n";
    std::cout << "Eps_p_v_n = " << Eps_p_v_n << " Eps_p_v_np1 = " << Eps_p_v_np1
              << " ep_v = " << ep_v << "\n";
#endif
  } else {

    // Compute elastic moduli at sigma^F
    auto KG_dKdG_np1 = d_elastic->getElasticModuliAndDerivatives(&state_np1);
    auto KG_np1      = KG_dKdG_np1.first;
    auto dKdG_np1    = KG_dKdG_np1.second;
    auto p_np1       = state_np1.I1 / 3.0;
    auto S_np1       = state_np1.deviatoricStressTensor;
    Matrix3 Z_np1    = computeZMatrix(KG_np1, dKdG_np1, p_np1, S_np1, M_np1);
#ifdef DEBUG_FIRST_ORDER_HARDENING
    std::cout << "Z_np1 = " << Z_np1 << "\n";
#endif

    // Compute CM = C:M
    double mu_np1     = KG_np1.shearModulus;
    double lambda_np1 = KG_np1.bulkModulus - 2.0 / 3.0 * mu_np1;
    Matrix3 CM_np1 =
      Util::Identity * (lambda_np1 * M_np1.Trace()) + M_np1 * (2.0 * mu_np1);
#ifdef DEBUG_FIRST_ORDER_HARDENING
    std::cout << "CM_np1 = " << CM_np1 << "\n";
#endif

    // Compute P = CM + Z;
    Matrix3 P_np1 = CM_np1 + Z_np1;
#ifdef DEBUG_FIRST_ORDER_HARDENING
    std::cout << "P_np1 = " << P_np1 << "\n";
#endif

    // Compute PN = P:N (note that we are using associated plasticity)
    double PN_np1 = P_np1.Contract(M_np1);

    // Compute Gamma (sig_trial - sig_n+1):M_n+1/(P_n:M_n+1)
    Matrix3 sig_trial = state_k_trial.stressTensor;
    double Gamma_np1 =
      (sig_trial - sig_fixed).Contract(M_np1) / P_np1.Contract(M_np1);
#ifdef DEBUG_FIRST_ORDER_HARDENING
    std::cout << " Gamma_np1 = " << Gamma_np1 << "\n";
    std::cout << " sig_trial = " << sig_trial << "\n";
    std::cout << " sig_fixed = " << sig_fixed << "\n";
    std::cout << " M_np1 = " << M_np1 << "\n";
#endif

    // Compute unscaled hardening modulus
    double df_dep_v = d_yield->df_depsVol(&state_np1, nullptr, nullptr);
    double H        = -(df_dep_v * M_np1.Trace()) / norm_M_np1;

    // Compute Gamma^H
    double Gamma_H = (Gamma_np1 * PN_np1) / (PN_np1 + H);

    // Compute updated stress
    Matrix3 sig_new = sig_trial - P_np1 * Gamma_H;
#ifdef DEBUG_FIRST_ORDER_HARDENING
    std::cout << "sig_new = " << sig_new << "\n";
#endif

    // Compute plastic strain increments
    Matrix3 delta_eps_p  = M_np1 * Gamma_H;
    double delta_eps_p_v = delta_eps_p.Trace();

    // Compute the elastic strain increment
    Matrix3 delta_eps_e  = deltaEps_new - delta_eps_p;
    double delta_eps_e_v = delta_eps_e.Trace();

    // Recompute internal variables
    status = computeInternalVariables(state_np1, delta_eps_e_v, delta_eps_p_v);
    if (status != Status::SUCCESS) {
      state_k_new = state_k_old;
      proc0cout << "computeInternalVariables has failed." << "\n";
      return status;
    }

    // Update the new state
    state_k_new              = state_np1;
    state_k_new.stressTensor = sig_new;
    state_k_new.updateStressInvariants();
    state_k_new.elasticStrainTensor += delta_eps_e;
    state_k_new.plasticStrainTensor += delta_eps_p;
    state_k_new.updatePlasticStrainInvariants();
    state_k_new.ep_v += delta_eps_p_v;
    computeElasticProperties(state_k_new);

    // Update the cumulative equivalent plastic strain
    Uintah::Matrix3 delta_eps_p_dev =
      delta_eps_p - Util::Identity * (delta_eps_p_v / 3.0);
    state_k_new.ep_cum_eq +=
      Util::sqrt_two_third *
      std::sqrt(delta_eps_p_dev.Contract(delta_eps_p_dev));

#ifdef DEBUG_INTERNAL_VAR_EVOLUTION
    std::cout << "firstOrderHardening: " << "\n"
              << "\t state_old = " << state_k_old << "\n"
              << "\t state_new = " << state_k_new << "\n";
#endif

    status = Status::SUCCESS;
  }

  return status;
}

/**
 * Method: consistencyBisectionSimplified
 * Purpose:
 *   Find the updated stress for hardening plasticity using the consistency
 *   bisection algorithm
 *   Returns whether the procedure is sucessful or has failed
 */
TabularPlasticity::Status
TabularPlasticityCap::consistencyBisectionSimplified(
  const Matrix3& deltaEps_new,
  const ModelState_TabularCap& state_k_old,
  const ModelState_TabularCap& state_k_trial,
  const Matrix3& deltaEps_e_fixed,
  const Matrix3& deltaEps_p_fixed,
  const Matrix3& sig_fixed,
  ModelState_TabularCap& state_k_new)
{
  Status status = Status::SUCCESS;

  ModelState_TabularCap state_n(state_k_old);
  ModelState_TabularCap state_np1(state_k_old);
  state_np1.stressTensor = sig_fixed;
  state_np1.updateStressInvariants();

  if (d_decrease_substep) {
    Matrix3 M_n   = d_yield->df_dsigma(Util::Identity, &state_n);
    Matrix3 M_np1 = d_yield->df_dsigma(Util::Identity, &state_np1);
    M_n /= M_n.Norm();
    double norm_M_np1 = M_np1.Norm();
    M_np1 /= norm_M_np1;
    double angle_M_n_np1 = std::abs(M_n.Contract(M_np1) - 1.0);

    if (angle_M_n_np1 > 1.0e-6) {
      status = Status::TOO_LARGE_YIELD_NORMAL_CHANGE;

#ifdef CHECK_RETURN_ALIGNMENT
      std::cout << "M_n = " << M_n << "\n";
      std::cout << "M_np1 = " << M_np1 << "\n";
      std::cout << "angle =" << angle_M_n_np1 << "\n";

      auto KG_dKdG_n   = d_elastic->getElasticModuliAndDerivatives(&state_n);
      auto KG_dKdG_np1 = d_elastic->getElasticModuliAndDerivatives(&state_np1);
      auto KG_n        = KG_dKdG_n.first;
      auto KG_np1      = KG_dKdG_np1.first;
      auto dKdG_n      = KG_dKdG_n.second;
      auto dKdG_np1    = KG_dKdG_np1.second;
      auto p_n         = state_n.I1 / 3.0;
      auto p_np1       = state_np1.I1 / 3.0;
      auto S_n         = state_n.deviatoricStressTensor;
      auto S_np1       = state_np1.deviatoricStressTensor;
      Matrix3 Z_n      = computeZMatrix(KG_n, dKdG_n, p_n, S_n, M_n);
      Matrix3 Z_np1    = computeZMatrix(KG_np1, dKdG_np1, p_np1, S_np1, M_np1);
      std::cout << "Z_n = " << Z_n << "\n";
      std::cout << "Z_np1 = " << Z_np1 << "\n";

      // Compute CM = C:M
      double mu_n       = KG_n.shearModulus;
      double mu_np1     = KG_np1.shearModulus;
      double lambda_n   = KG_n.bulkModulus - 2.0 / 3.0 * mu_n;
      double lambda_np1 = KG_np1.bulkModulus - 2.0 / 3.0 * mu_np1;
      Matrix3 CM_n =
        Util::Identity * (lambda_n * M_n.Trace()) + M_n * (2.0 * mu_n);
      Matrix3 CM_np1 =
        Util::Identity * (lambda_np1 * M_np1.Trace()) + M_np1 * (2.0 * mu_np1);
      std::cout << "CM_n = " << CM_n << "\n";
      std::cout << "CM_np1 = " << CM_np1 << "\n";

      // Compute P = CM + Z;
      Matrix3 P_n   = CM_n + Z_n;
      Matrix3 P_np1 = CM_np1 + Z_np1;
      std::cout << "P_n = " << P_n << "\n";
      std::cout << "P_np1 = " << P_np1 << "\n";

      // Compute Gamma (sig_trial - sig_n+1):M_n+1/(P_n:M_n+1)
      Matrix3 sig_trial = state_k_trial.stressTensor;
      Matrix3 sig_new   = sig_fixed;
      double Gamma_n =
        (sig_trial - sig_new).Contract(M_np1) / P_n.Contract(M_np1);
      double Gamma_np1 =
        (sig_trial - sig_new).Contract(M_np1) / P_np1.Contract(M_np1);
      std::cout << "Gamma_n = " << Gamma_n << " Gamma_np1 = " << Gamma_np1
                << "\n";

      // Compute sigma
      Matrix3 sig_np1 = sig_trial - P_n * Gamma_n;
      std::cout << "sig_np1 = " << sig_np1 << "\n";
      std::cout << "sig_fixed = " << sig_fixed << "\n";

      // Compute volumetric plastic strain
      Matrix3 eps_p = state_n.plasticStrainTensor + deltaEps_p_fixed;
      double ep_v   = eps_p.Trace();

      // Compute eps_p and ev_v^p
      Matrix3 Eps_p_n    = state_n.plasticStrainTensor + M_n * Gamma_n;
      double Eps_p_v_n   = state_n.ep_v + M_n.Trace() * Gamma_n;
      Matrix3 Eps_e_n    = deltaEps_new - Eps_p_n;
      Matrix3 Eps_p_np1  = state_n.plasticStrainTensor + M_np1 * Gamma_np1;
      double Eps_p_v_np1 = state_n.ep_v + M_np1.Trace() * Gamma_np1;
      Matrix3 Eps_e_np1  = deltaEps_new - Eps_p_np1;
      std::cout << "Eps_p_n = " << Eps_p_n << "\n";
      std::cout << "eps_p = " << eps_p << "\n";
      std::cout << "Eps_e_n = " << Eps_e_n << "\n";
      std::cout << "Eps_p_np1 = " << Eps_p_np1 << "\n";
      std::cout << "Eps_e_np1 = " << Eps_e_np1 << "\n";
      std::cout << "Eps_p_v_n = " << Eps_p_v_n
                << " Eps_p_v_np1 = " << Eps_p_v_np1 << " ep_v = " << ep_v
                << "\n";
#endif

      return status;
    }
  }

  // bisection convergence tolerance on eta (if changed, change imax)
  const double TOLERANCE = d_consistency_bisection_tolerance;
  // imax = ceil(-10.0*log(TOL)); // Update this if TOL changes
  const int IMAX = d_max_bisection_iterations;

  // Get the old state
  Matrix3 sig_old   = state_k_old.stressTensor;
  Matrix3 eps_e_old = state_k_old.elasticStrainTensor;
  Matrix3 eps_p_old = state_k_old.plasticStrainTensor;

  // Get the fixed non-hardening return state and compute invariants
  //double deltaEps_e_v_fixed = deltaEps_e_fixed.Trace();
  double deltaEps_p_v_fixed = deltaEps_p_fixed.Trace();

  // Create a state for the fixed non-hardening yield surface state
  // and update only the stress and plastic strain
  ModelState_TabularCap state_k_fixed(state_k_old);
  state_k_fixed.stressTensor = sig_fixed;
  state_k_fixed.updateStressInvariants();
  state_k_fixed.elasticStrainTensor = eps_e_old + deltaEps_e_fixed;
  state_k_fixed.plasticStrainTensor = eps_p_old + deltaEps_p_fixed;
  state_k_fixed.updatePlasticStrainInvariants();
  computeElasticProperties(state_k_fixed);
  //computeElasticProperties(
  //  state_k_fixed, deltaEps_e_v_fixed, deltaEps_p_v_fixed);

  // Initialize the new consistently updated state
  ModelState_TabularCap state_trial_mid(state_k_old);
  Matrix3 sig_fixed_mid       = sig_fixed;
  Matrix3 deltaEps_e_mid      = deltaEps_e_fixed;
  Matrix3 deltaEps_p_mid      = deltaEps_p_fixed;
  double deltaEps_p_v_mid_new = 0.5 * deltaEps_p_v_fixed;
  double deltaEps_p_v_mid_old = deltaEps_p_v_fixed;

#ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
  double norm_deltaEps_p_fixed = deltaEps_p_fixed.Norm();
  double norm_deltaEps_p_mid   = norm_deltaEps_p_fixed;
#endif

  // Start loop
  int ii        = 1;
  double eta_lo = 0.0, eta_hi = 1.0, eta_mid = 0.5;

  while (std::abs(eta_hi - eta_lo) > TOLERANCE) {

#ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
    if (state_k_fixed.particleID == 124554051588) {
      std::cout << "Enter: consistency_iter = " << ii << std::setprecision(16)
                << " eta_hi = " << eta_hi << " eta_lo = " << eta_lo
                << " eta_mid = " << eta_mid << "\n";
      std::cout << std::setprecision(16)
                << " delta_eps_p_v_fixed = " << deltaEps_p_v_fixed
                << " delta_eps_p_v_mid_old = " << deltaEps_p_v_mid_old
                << " delta_eps_p_v_mid_new = " << deltaEps_p_v_mid_new
                << " ratio - 1 = "
                << (deltaEps_p_v_mid_new / deltaEps_p_v_mid_old - 1.0)
                << "\n";
    }
#endif

    // Reset the local trial state
    state_trial_mid = state_k_old;

    // Compute the volumetric plastic strain at eta = eta_mid
    double deltaEps_p_v_mid = eta_mid * deltaEps_p_v_fixed;
    double deltaEps_e_v_mid = deltaEps_new.Trace() - deltaEps_p_v_mid;

#ifdef CHECK_CONSISTENCY_BISECTION_K
    std::cout << "Before: K = " << state_trial_mid.bulkModulus
              << " G = " << state_trial_mid.shearModulus
              << " X = " << state_trial_mid.capX
              << " del eps_e_v = " << deltaEps_e_v_mid
              << " del eps_p_v = " << deltaEps_p_v_mid << "\n";
#endif

    // Update the elastic moduli
    ModelState_TabularCap state_temp(state_trial_mid);
    state_temp.stressTensor = sig_fixed_mid;
    state_temp.plasticStrainTensor += deltaEps_p_v_mid;
    state_temp.elasticStrainTensor += deltaEps_p_v_mid;
    state_temp.updateStressInvariants();
    state_temp.updatePlasticStrainInvariants();
    computeElasticProperties(state_temp);
    state_trial_mid.bulkModulus  = state_temp.bulkModulus;
    state_trial_mid.shearModulus = state_temp.shearModulus;
    
    // Update the internal variables at eta = eta_mid in the local trial state
    status = computeInternalVariables(
      state_trial_mid, deltaEps_e_v_mid, deltaEps_p_v_mid);
    if (status != Status::SUCCESS) {
      state_k_new = state_k_old;
      proc0cout << "computeInternalVariables has failed." << "\n";
      return status;
    }

#ifdef CHECK_CONSISTENCY_BISECTION_K
    std::cout << "After: K = " << state_trial_mid.bulkModulus
              << " G = " << state_trial_mid.shearModulus
              << " X = " << state_trial_mid.capX << "\n";
#endif

    // Update yield surface
    Polyline yield_f_pts =
      d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_trial_mid);
    state_trial_mid.yield_f_pts = yield_f_pts;

    // Update the trial stress
    auto stress_trial = computeTrialStress(state_trial_mid, deltaEps_new);
    state_trial_mid.stressTensor = stress_trial;
    state_trial_mid.updateStressInvariants();

    // Test the yield condition to check whether the yield surface moves beyond
    // the trial stress state when the internal variables are changed.
    // If the yield surface is too big, the plastic strain is reduced
    // by bisecting <eta> and the loop is repeated.
    auto yield = d_yield->evalYieldCondition(&state_trial_mid);

    // If the local trial state is inside the updated yield surface the yield
    // condition evaluates to "elastic".  We need to reduce the size of the
    // yield surface by decreasing the plastic strain increment.
    // Elastic or on yield surface
    if (yield.second == Util::YieldStatus::IS_ELASTIC) {
      eta_hi  = eta_mid;
      eta_mid = 0.5 * (eta_lo + eta_hi);

      deltaEps_p_v_mid_old = deltaEps_p_v_mid;
      deltaEps_p_v_mid_new = eta_mid * deltaEps_p_v_fixed;

#ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
      if (state_k_fixed.particleID == 124554051588) {
        std::cout << "Elastic: consistency_iter = " << ii
                  << " eta_hi = " << eta_hi << " eta_lo = " << eta_lo
                  << " eta_mid = " << eta_mid << "\n";
        std::cout << " delta_eps_p_v_mid_old = " << deltaEps_p_v_mid_old
                  << " delta_eps_p_v_mid_new = " << deltaEps_p_v_mid_new
                  << " ratio - 1 = "
                  << (deltaEps_p_v_mid_new / deltaEps_p_v_mid_old - 1)
                  << "\n";
      }
#endif
      ii++;
      continue;
    }

    // At this point, state_trial_mid contains the trial stress, the plastic
    // strain at the beginning of the timestep, the updated elastic moduli,
    // and the updated values of the internal  variables
    // The yield surface depends only on X.  We will compute the updated
    // location of the yield surface based on the updated internal variables
    // and do a non-hardening return to that yield surface.
    ModelState_TabularCap state_old_mid(state_k_old);
    state_old_mid.bulkModulus  = state_trial_mid.bulkModulus;
    state_old_mid.shearModulus = state_trial_mid.shearModulus;
    state_old_mid.capX         = state_trial_mid.capX;
    state_old_mid.yield_f_pts  = state_trial_mid.yield_f_pts;
    state_old_mid.updateStressInvariants();
    state_old_mid.updatePlasticStrainInvariants();

    status = nonHardeningReturn(deltaEps_new,
                                state_old_mid,
                                state_trial_mid,
                                sig_fixed_mid,
                                deltaEps_e_mid,
                                deltaEps_p_mid);
    if (status != Status::SUCCESS) {
      proc0cout << "nonHardeningReturn inside consistencyBisection failed."
                << "\n";
      return status;
    }

#ifdef CHECK_CONSISTENCY_BISECTION_FIXED
    std::cout << "Fixed: \n"
              << " sig_fixed_diff =  " << sig_fixed_mid - sig_fixed << "\n"
              << " deltaEps_e_diff = " << deltaEps_e_mid - deltaEps_e_fixed
              << "\n"
              << " deltaEps_p_diff = " << deltaEps_p_mid - deltaEps_p_fixed
              << "\n";
#endif

    // Update volumetric plastic strain
    deltaEps_p_v_mid_old = deltaEps_p_v_mid_new;
    deltaEps_p_v_mid_new = deltaEps_p_mid.Trace();

    // Check whether the isotropic component of the return has changed sign, as
    // this would indicate that the cap apex has moved past the trial stress,
    // indicating too much plastic strain in the return.
    Matrix3 sig_trial           = state_trial_mid.stressTensor;
    double diff_trial_fixed_new = (sig_trial - sig_fixed_mid).Trace();
    double diff_trial_fixed     = (sig_trial - sig_fixed).Trace();
    if (std::signbit(diff_trial_fixed_new) != std::signbit(diff_trial_fixed)) {
      eta_hi  = eta_mid;
      eta_mid = 0.5 * (eta_lo + eta_hi);
    } else {
      // if (norm_deltaEps_p_fixed_new > eta_mid*norm_deltaEps_p_fixed) {
      if (std::abs(deltaEps_p_v_mid_new) > std::abs(deltaEps_p_v_mid_old)) {
        eta_lo  = eta_mid;
        eta_mid = 0.5 * (eta_lo + eta_hi);
      } else {
        eta_hi  = eta_mid;
        eta_mid = 0.5 * (eta_lo + eta_hi);
      }
    }

// Compare magnitude of plastic strain with prior update
#ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
    norm_deltaEps_p_mid   = deltaEps_p_mid.Norm();
    norm_deltaEps_p_fixed = eta_mid * deltaEps_p_fixed.Norm();
    if (state_k_fixed.particleID == 124554051588) {
      std::cout << "Exit: \n"
                << " eta_mid = " << eta_mid
                << " eta_mid*||deltaEps_p_fixed|| = "
                << eta_mid * norm_deltaEps_p_fixed
                << " ||deltaEps_p_mid_new|| = " << norm_deltaEps_p_mid
                << " ratio = "
                << eta_mid * norm_deltaEps_p_fixed / norm_deltaEps_p_mid
                << "\n";
      std::cout << " delta_eps_p_v_mid_old = " << deltaEps_p_v_mid_old
                << " delta_eps_p_v_mid_new = " << deltaEps_p_v_mid_new
                << " ratio - 1 = "
                << (deltaEps_p_v_mid_new / deltaEps_p_v_mid_old - 1)
                << "\n";
    }
#endif

    // Increment i and check
    ii++;
    if (ii > IMAX) {
      state_k_new = state_k_old;
      proc0cout << "Consistency bisection has failed because ii > IMAX." << ii
                << " > " << IMAX << "\n";
      return Status::TOO_MANY_CONSISTENCY_ITERATIONS; // bool isSuccess = false;
    }

  } // end  while (std::abs(eta_hi - eta_lo) > TOLERANCE);

  // Update the stress and plastic strain of the new state +  the elastic moduli
  // and capX
  state_k_new.stressTensor = sig_fixed_mid;
  state_k_new.updateStressInvariants();
  state_k_new.elasticStrainTensor = eps_e_old + deltaEps_e_mid;
  state_k_new.plasticStrainTensor = eps_p_old + deltaEps_p_mid;
  state_k_new.updatePlasticStrainInvariants();
  state_k_new.bulkModulus  = state_trial_mid.bulkModulus;
  state_k_new.shearModulus = state_trial_mid.shearModulus;
  state_k_new.capX         = state_trial_mid.capX;

#ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
  if (state_k_fixed.particleID == 124554051588) {
    std::cout << "Final: K = " << state_k_new.bulkModulus
              << " G = " << state_k_new.shearModulus
              << " X = " << state_k_new.capX << "\n";
    std::cout << "K_before_consistency = " << state_k_old.bulkModulus
              << "\n";
    std::cout << "K_after_consistency = " << state_k_new.bulkModulus
              << "\n";
    std::cout << "I1_before_consistency = " << state_k_old.stressTensor.Trace()
              << "\n";
    std::cout << "I1_after_consistency = " << state_k_new.stressTensor.Trace()
              << "\n";
  }
#endif

  // Update the cumulative equivalent plastic strain
  double deltaEps_p_v = deltaEps_p_mid.Trace();
  Uintah::Matrix3 deltaEps_p_dev =
    deltaEps_p_mid - Util::Identity * (deltaEps_p_v / 3.0);
  state_k_new.ep_cum_eq =
    state_k_old.ep_cum_eq +
    std::sqrt(2.0 / 3.0 * deltaEps_p_dev.Contract(deltaEps_p_dev));

#ifdef DEBUG_INTERNAL_VAR_EVOLUTION
  std::cout << "consistencyBisection: " << "\n"
            << "\t state_old = " << state_k_old << "\n"
            << "\t state_new = " << state_k_new << "\n";
#endif

  // Return success = true
  return Status::SUCCESS; // bool isSuccess = true;
}

/**
 * Method: computeInternalVariables
 * Purpose:
 *   Update an old state with new values of internal variables given the old
 *   state and an increment in volumetric plastic strain
 */
TabularPlasticity::Status
TabularPlasticityCap::computeInternalVariables(ModelState_TabularCap& state,
                                               const double& delta_eps_e_v,
                                               const double& delta_eps_p_v)
{
  ModelState_TabularCap tempState = state;
  tempState.elasticStrainTensor += (Util::Identity * delta_eps_e_v);
  tempState.plasticStrainTensor += (Util::Identity * delta_eps_p_v);
  tempState.updatePlasticStrainInvariants();
#ifdef TEST_CAP_TRANSLATION
  tempState.ep_v -= X_trans_factor;
#endif

  // Update the hydrostatic compressive strength
  double X_new = d_capX->computeInternalVariable("capX", nullptr, &tempState);

  // Update the state with new values of the internal variables
  state.capX = X_new;
#ifdef TEST_CAP_VARIATION
  state.capX *= X_scale_factor;
#endif

  return Status::SUCCESS;
}

double
TabularPlasticityCap::computeInternalVariable(
  const ModelState_TabularCap& state,
  const Matrix3& elasticStrain,
  const Matrix3& plasticStrain) const
{
  ModelState_TabularCap tempState = state;
  tempState.elasticStrainTensor   = elasticStrain;
  tempState.plasticStrainTensor   = plasticStrain;
  tempState.updatePlasticStrainInvariants();

#ifdef TEST_CAP_TRANSLATION
  tempState.ep_v -= X_trans_factor;
#endif

  // Update the hydrostatic compressive strength
  double X_new = d_capX->computeInternalVariable("capX", nullptr, &tempState);

#ifdef TEST_CAP_VARIATION
  X_new *= X_scale_factor;
#endif

  return X_new;
}

double
TabularPlasticityCap::computeInternalVariable(
  const ModelState_TabularCap& state) const
{
  ModelState_TabularCap tempState = state;

#ifdef TEST_CAP_TRANSLATION
  tempState.ep_v -= X_trans_factor;
#endif

  // Update the hydrostatic compressive strength
  double X_new = d_capX->computeInternalVariable("capX", nullptr, &tempState);

#ifdef TEST_CAP_VARIATION
  X_new *= X_scale_factor;
#endif

  return X_new;
}

/**
 * Method: computeElasticDeltaT
 * Purpose:
 *   Compute the fraction of the time step that is elastic
 * Inputs:
 *   double       - Total delta t
 *   ModelState   - Old state
 *   ModelState   - Trial state
 * Returns:
 *   bool         - true if success
 *                  false if failure
 *   double       - t value of intersection of line segment joining
 *                  old state and trial state with yield surface
 *                  in pbar-sqrtJ2 space
 *   double       - elastic delta t
 *   Point        - point of intersection in pbar-sqrtJ2 space
 */
std::tuple<bool, double, double, Uintah::Point>
TabularPlasticityCap::computeElasticDeltaT(
  double dt,
  const ModelState_TabularCap& state_old,
  const ModelState_TabularCap& state_trial)
{
  // Set up the line segment
  Uintah::Point seg_start(-state_old.I1 / 3.0, state_old.sqrt_J2, 0);
  Uintah::Point seg_end(-state_trial.I1 / 3.0, state_trial.sqrt_J2, 0);

  // Set up the polyline
  Polyline polyline = state_old.yield_f_pts;

  // Compute intersection
  bool status;
  double t_val;
  Uintah::Point intersection;
  std::tie(status, t_val, intersection) =
    Vaango::Util::intersectionPointBSpline(polyline, seg_start, seg_end);

  // Compute elastic delta t
  double elastic_dt = t_val * dt;
  if (!status) {
    /*
    if (t_val > 10 || t_val < -10) {
      std::ostringstream out;
      out << "seg_start = " << seg_start << ";\n";
      out << "seg_end = " << seg_end << ";\n";
      out << "polyline = [";
      std::copy(polyline.begin(), polyline.end(),
                std::ostream_iterator<Point>(out, " "));
      out << "];\n";

      out << "status = " << status << " t = " << t_val
                << " elastic_dt = " << elastic_dt << "\n";
      out << "intersection = " << intersection << ";\n";
      //std::cout << out.str();
      throw InvalidValue(out.str(), __FILE__, __LINE__);
    }
    */
    return std::make_tuple(false, t_val, elastic_dt, intersection);
  }

  return std::make_tuple(true, t_val, elastic_dt, intersection);
}

void
TabularPlasticityCap::addRequiresDamageParameter(Task* task,
                                                 const MPMMaterial* matl,
                                                 const PatchSet*) const
{
  // Require the damage parameter
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, lb->pRemoveLabel_preReloc, matlset, Ghost::None);
}

void
TabularPlasticityCap::getDamageParameter(const Patch* patch,
                                         ParticleVariable<int>& damage,
                                         int matID,
                                         DataWarehouse* old_dw,
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
TabularPlasticityCap::carryForward(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  // Carry forward the data.
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch   = patches->get(p);
    int matID            = matl->getDWIndex();
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
TabularPlasticityCap::addComputesAndRequires(Task*,
                                             const MPMMaterial*,
                                             const PatchSet*,
                                             const bool,
                                             const bool) const
{
  std::cout << "NO Implicit VERSION OF addComputesAndRequires EXISTS YET FOR "
               "TabularPlasticityCap"
            << endl;
}

/*!
 * ---------------------------------------------------------------------------------------
 *  This is needed for converting from one material type to another.  The
 * functionality
 *  has been removed from the main Uintah branch.
 *  ---------------------------------------------------------------------------------------
 */
void
TabularPlasticityCap::allocateCMDataAdd(DataWarehouse* new_dw,
                                        ParticleSubset* addset,
                                        ParticleLabelVariableMap* newState,
                                        ParticleSubset* delset,
                                        DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "Material conversion after failure not implemented for "
         "TabularPlasticityCap.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}
