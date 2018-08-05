/*
 * The MIT License
 *
 * Copyright (c) 2015-2018 Parresia Research Limited, New Zealand
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
#include <CCA/Components/MPM/ConstitutiveModel/TabularPlasticityCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldConditionFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVariableModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>

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
#define CHECK_PLASTIC_RATE
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
//#define TIME_TABLE_LOOKUP
//#define TIME_SUBSTEP
//#define TIME_YIELD_PTS

using namespace Vaango;
using Uintah::VarLabel;
using Uintah::Matrix3;
using std::ostringstream;
using std::endl;

TabularPlasticityCap::TabularPlasticityCap(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* mpmFlags)
  : TabularPlasticity(ps, mpmFlags)
{
  d_capX = Vaango::InternalVariableModelFactory::create(ps);
  if (!d_capX) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating InternalVariableModel."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  ps->getWithDefault("consistency_bisection_tolerance",
                     d_consistency_bisection_tolerance, 1.0e-4);
  d_max_bisection_iterations =
    (int)std::ceil(-10.0 * std::log(d_consistency_bisection_tolerance));

  if (d_consistency_bisection_tolerance < 1.0e-16 ||
      d_consistency_bisection_tolerance > 1.0e-2) {
    ostringstream warn;
    warn << "Consistency bisection tolerance should be in range [1.0e-16, "
            "1.0e-2].  Default = 1.0e-4"
         << endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

TabularPlasticityCap::TabularPlasticityCap(const TabularPlasticityCap& cm)
  : TabularPlasticity(cm)
{
  d_capX = Vaango::InternalVariableModelFactory::createCopy(d_capX);

  // Consistency bisection
  d_consistency_bisection_tolerance = cm.d_consistency_bisection_tolerance;
  d_max_bisection_iterations = cm.d_max_bisection_iterations;

}

TabularPlasticityCap::TabularPlasticityCap(const TabularPlasticityCap* cm)
  : TabularPlasticity(*cm)
{
  d_capX = Vaango::InternalVariableModelFactory::createCopy(cm->d_capX);

  // Consistency bisection
  d_consistency_bisection_tolerance = cm->d_consistency_bisection_tolerance;
  d_max_bisection_iterations = cm->d_max_bisection_iterations;
}

TabularPlasticityCap::~TabularPlasticityCap()
{
  delete d_capX;
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

  cm_ps->appendElement("yield_surface_radius_scaling_factor",
                       d_cm.yield_scale_fac);
  cm_ps->appendElement("consistency_bisection_tolerance",
                       d_consistency_bisection_tolerance);
  cm_ps->appendElement("subcycling_characteristic_number",
                       d_cm.subcycling_characteristic_number);
}

TabularPlasticityCap*
TabularPlasticityCap::clone()
{
  return scinew TabularPlasticityCap(*this);
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
TabularPlasticityCap::addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patch) const
{
  TabularPlasticity::addInitialComputesAndRequires(task, matl, patch);
  d_capX->addInitialComputesAndRequires(task, matl, patch);
}

/*!------------------------------------------------------------------------*/
void
TabularPlasticityCap::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  d_capX->initializeInternalVariable(pset, new_dw);
  TabularPlasticity::initializeCMData(patch, matl, new_dw);
}

// Compute stable timestep based on both the particle velocities
// and wave speed
void
TabularPlasticityCap::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
                                            DataWarehouse* new_dw)
{
  TabularPlasticity::computeStableTimestep(patch, matl, new_dw);
}

/**
 * Added computes/requires for computeStressTensor
 */
void
TabularPlasticityCap::addComputesAndRequires(Task* task, const MPMMaterial* matl,
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
TabularPlasticityCap::computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  // Global loop over each patch
  for (int p = 0; p < patches->size(); p++) {

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
    constParticleVariable<int> pRemove;
    constParticleVariable<double> pEev, pEpv, pEpeq_old;
    constParticleVariable<Matrix3> pEe, pEp;
    old_dw->get(pEe,       pElasticStrainLabel, pset);
    old_dw->get(pEev,      pElasticVolStrainLabel, pset);
    old_dw->get(pEp,       pPlasticStrainLabel, pset);
    old_dw->get(pEpv,      pPlasticVolStrainLabel, pset);
    old_dw->get(pEpeq_old, pPlasticCumEqStrainLabel, pset);
    old_dw->get(pRemove,   pRemoveLabel, pset);

    ParticleVariable<int>    pRemove_new;
    ParticleVariable<double> pEev_new, pEpv_new, pEpeq_new;
    ParticleVariable<Matrix3> pEe_new, pEp_new;
    new_dw->allocateAndPut(pEe_new,     pElasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEev_new,    pElasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEp_new,     pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpv_new,    pPlasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pEpeq_new,   pPlasticCumEqStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pRemove_new, pRemoveLabel_preReloc, pset);

    // Get and allocate the hydrostatic strength
    constParticleVariable<double> pCapX;
    ParticleVariable<double> pCapX_new;
    d_capX->getInternalVariable(pset, old_dw, pCapX);
    d_capX->allocateAndPutInternalVariable(pset, new_dw, pCapX_new);

    // Set up global particle variables to be read and written
    delt_vartype delT;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pDefGrad, pStress_old;

    old_dw->get(delT,           lb->delTLabel, getLevel(patches));
    old_dw->get(pMass,          lb->pMassLabel, pset);
    old_dw->get(pParticleID,    lb->pParticleIDLabel, pset);
    old_dw->get(pVelocity,      lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad,       lb->pDefGradLabel, pset);
    old_dw->get(pStress_old,    lb->pStressLabel, pset);

    // Get the particle variables computed in interpolateToParticlesAndUpdate()
    constParticleVariable<double> pVolume;
    constParticleVariable<Matrix3> pVelGrad_new, pDefGrad_new;
    new_dw->get(pVolume,      lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);
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

      // Compute the symmetric part of the velocity gradient
      // std::cout << "DefGrad = " << pDefGrad_new[idx] << std::endl;
      // std::cout << "VelGrad = " << pVelGrad_new[idx] << std::endl;
      Matrix3 DD = (pVelGrad_new[idx] + pVelGrad_new[idx].Transpose()) * .5;

      // Use polar decomposition to compute the rotation and stretch tensors
      Matrix3 FF = pDefGrad[idx];
      Matrix3 RR, UU;
      FF.polarDecompositionRMB(UU, RR);

      // Compute the unrotated symmetric part of the velocity gradient
      DD = (RR.Transpose()) * (DD * RR);

#ifdef CHECK_FOR_NAN
      if (std::isnan(DD(0, 0))) {
        std::cout << " L_new = " << pVelGrad_new[idx]
                  << " F_new = " << pDefGrad_new[idx] << " F = " << FF
                  << " R = " << RR << " U = " << UU << " D = " << DD
                  << " delT = " << delT << std::endl;
        throw InternalError("**ERROR** Zero or Nan in rate of deformation",
                             __FILE__, __LINE__);
      }
#endif

      // Compute the unrotated stress at the start of the current timestep
      Matrix3 sigma_old = (RR.Transpose()) * (pStress_old[idx] * RR);
      // std::cout << "pStress_old = " << pStress_old[idx] << std::endl
      //           << "sigma_old = " << sigma_old << std::endl;

      // Set up model state
      ModelState_TabularCap state_old;
      state_old.particleID = pParticleID[idx];
      state_old.stressTensor = sigma_old;
      state_old.elasticStrainTensor = pEe[idx];
      state_old.plasticStrainTensor = pEp[idx];
      state_old.ep_cum_eq = pEpeq_old[idx];
      state_old.capX = pCapX[idx];
      // std::cout << "state_old.Stress = " << state_old.stressTensor <<
      // std::endl;

      // Update yield surface
      #ifdef TIME_YIELD_PTS
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
      #endif
      Polyline yield_f_pts = d_yield->computeYieldSurfacePolylinePbarSqrtJ2(&state_old);
      #ifdef TIME_YIELD_PTS
        end = std::chrono::system_clock::now();
        std::cout << "Compute yield pts. : Time taken = " 
                  << std::chrono::duration<double>(end-start).count() << std::endl;
      #endif
      std::array<double,3> range = d_yield->getYieldConditionRange(yield_f_pts);
      state_old.yield_f_pts = yield_f_pts;
      state_old.I1_max = -range[0]*3.0;
      state_old.I1_min = -range[1]*3.0;
      state_old.sqrtJ2_max = range[2];

      // Compute the elastic moduli at t = t_n
      computeElasticProperties(state_old);
      //std::cout << "State old: " << state_old << std::endl;

      //---------------------------------------------------------
      // Rate-independent plastic step
      // Divides the strain increment into substeps, and calls substep function
      ModelState_TabularCap state_new;
      #ifdef TIME_SUBSTEP
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
      #endif
      Status status = rateIndependentPlasticUpdate(
        DD, delT, idx, pParticleID[idx], state_old, state_new);
      #ifdef TIME_SUBSTEP
        end = std::chrono::system_clock::now();
        std::cout << "Substep : Time taken = " 
                  << std::chrono::duration<double>(end-start).count() << std::endl;
      #endif

      if (status == Status::SUCCESS) {
        pRemove_new[idx] = 0;
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
        pCapX_new[idx] = state_new.capX;

        // Elastic volumetric strain at end of step, compute from updated
        // deformation gradient.
        // H = ln(U) => tr(H) = tr(ln(U)) = ln(det(U)) = ln(sqrt(det(FT)
        // det(F))) = ln J
        pEev_new[idx] =
          log(pDefGrad_new[idx].Determinant()) - pEpv_new[idx];

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
        pCapX_new[idx] = state_old.capX;
      }

      //---------------------------------------------------------
      // Use polar decomposition to compute the rotation and stretch tensors.
      // These checks prevent
      // failure of the polar decomposition algorithm if [F_new] has some
      // extreme values.
      Matrix3 FF_new = pDefGrad_new[idx];
      double Fmax_new = FF_new.MaxAbsElem();
      double JJ_new = FF_new.Determinant();
      if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
        pRemove_new[idx] = -999;
        proc0cout << "Deformation gradient component unphysical: [F] = " << FF
                  << std::endl;
        proc0cout << "Resetting [F]=[I] for this step and deleting particle"
                  << " idx = " << idx << " particleID = " << pParticleID[idx]
                  << std::endl;
        Util::Identity.polarDecompositionRMB(UU, RR);
      } else {
        FF_new.polarDecompositionRMB(UU, RR);
      }

      // Compute the rotated dynamic and quasistatic stress at the end of the
      // current timestep
      pStress_new[idx] = (RR * pStress_new[idx]) * (RR.Transpose());
      // std::cout << "pStress_new = " << pStress_new[idx]
      //          << "pStressQS_new = " << pStressQS_new[idx] << std::endl;

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
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) * Util::one_third;
        double c_bulk = sqrt(bulk / rho_cur);
        p_q[idx] = artificialBulkViscosity(DD.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      // Compute the averaged stress
      Matrix3 AvgStress = (pStress_new[idx] + pStress_old[idx]) * 0.5;

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
TabularPlasticityCap::rateIndependentPlasticUpdate(const Matrix3& D, const double& delT,
                                    particleIndex idx, long64 pParticleID,
                                    const ModelState_TabularCap& state_old,
                                    ModelState_TabularCap& state_new)
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
  ModelState_TabularCap state_trial(state_old);
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

  // Set up the initial states for the substeps
  Status status = Status::SUCCESS;
  ModelState_TabularCap state_k_old;
  ModelState_TabularCap state_k_new;

  const int CHI_MAX = 10; // max allowed number of timestep reductions
  int chi = 0;            // subcycle timestep reduction count
  do {
    // Compute the substep time increment list 
    std::vector<double> substeps(nsub, delT/ static_cast<double>(nsub));

    // Loop through the substeps
    state_k_old = state_old;
    state_k_new = state_old;
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
         proc0cout << "**WARNING** Decreasing substep time increment to " << dt/2
                   << " because computeSubstep failed." << std::endl;
         break;
      }

      state_k_old = state_k_new;
      tlocal += dt;

      #ifdef WRITE_YIELD_SURF
        std::cout << "K = " << state_k_new.bulkModulus << std::endl;
        std::cout << "G = " << state_k_new.shearModulus << std::endl;
        std::cout << "capX = " << state_k_new.capX << std::endl;
        Matrix3 sig = state_k_new.stressTensor;
        std::cout << "sigma_new = np.array([[" << sig(0, 0) << "," << sig(0, 1)
                  << "," << sig(0, 2) << "],[" << sig(1, 0) << "," << sig(1, 1)
                  << "," << sig(1, 2) << "],[" << sig(2, 0) << "," << sig(2, 1)
                  << "," << sig(2, 2) << "]])" << std::endl;
        std::cout << "plot_stress_state(K, G, sigma_trial, sigma_new, 'b')"
                  << std::endl;
      #endif

      #ifdef CHECK_SUBSTEP
        std::cout << "tlocal = " << tlocal << " delT = " << delT
                  << " nsub = " << nsub << std::endl;
      #endif
    } // End substep range for loop

  } while (status != Status::SUCCESS && chi < CHI_MAX);

  if (status == Status::SUCCESS) {
    state_new = state_k_new;

    #ifdef CHECK_INTERNAL_VAR_EVOLUTION
      // if (state_old.particleID == 3377699720593411) {
      std::cout << "rateIndependentPlasticUpdate: "
                << " ep_v_old = " << state_old.ep_v
                << " ep_v_new = " << state_new.ep_v 
                << " Xbar_old = " << -state_old.capX
                << " Xbar_new = " << -state_new.capX
                << std::endl;
      //}
    #endif
  } else {
    if (status == Status::TOO_LARGE_YIELD_NORMAL_CHANGE) {
      state_new = state_k_new;
      status = Status::SUCCESS;
      proc0cout << "Substep failed because chi = " << chi << " > " << CHI_MAX << "."
                << " Continuing with unconverged solution."
                << std::endl;
    } else {
      state_new = state_k_old;
      proc0cout << "Substep failed because chi = " << chi << " > " << CHI_MAX << "."
                << " Proceeding to delete particle."
                << std::endl;
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
 * Side effects:
 *   **WARNING** Also computes stress invariants and plastic strain invariants
 */
void
TabularPlasticityCap::computeElasticProperties(ModelState_TabularCap& state)
{
  #ifdef TIME_TABLE_LOOKUP
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
  #endif
  TabularPlasticity::computeElasticProperties(state);
  #ifdef TIME_TABLE_LOOKUP
    end = std::chrono::system_clock::now();
    std::cout << "Compute elastic properties : Time taken = " 
              << std::chrono::duration<double>(end-start).count() << std::endl;
  #endif
}

/**
 * Method: computeElasticProperties
 *
 * Purpose:
 *   Compute the bulk and shear modulus at a given state
 */
void
TabularPlasticityCap::computeElasticProperties(ModelState_TabularCap& state,
                                               double plasticVolStrainInc)
{
  ModelState_TabularCap tempState = state;
  tempState.ep_v += plasticVolStrainInc;
  TabularPlasticity::computeElasticProperties(tempState);
  state.bulkModulus = tempState.bulkModulus;
  state.shearModulus = tempState.shearModulus;
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
TabularPlasticityCap::computeStepDivisions(particleIndex idx, long64 particleID,
                            const ModelState_TabularCap& state_old,
                            const ModelState_TabularCap& state_trial)
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
  double I1_size = 0.5 * (state_old.I1_max - state_old.I1_min);
  double J2_size = state_old.sqrtJ2_max;
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
    proc0cout << "\t I1_max = " << state_old.I1_max 
              << " I1_min = " << state_old.I1_min 
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
TabularPlasticityCap::computeSubstep(const Matrix3& D, const double& dt,
                      const ModelState_TabularCap& state_k_old,
                      ModelState_TabularCap& state_k_new)
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
  ModelState_TabularCap state_k_trial(state_k_old);
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
  Status status =
    nonHardeningReturn(deltaEps, state_k_old, state_k_trial, sig_fixed,
                       deltaEps_e_fixed, deltaEps_p_fixed);
  if (status != Status::SUCCESS) {
    proc0cout << "**WARNING** nonHardeningReturn has failed." << std::endl;
    return status;
  }

  // Do first-order hardening update
  state_k_new = state_k_old;
  status = firstOrderHardeningUpdate(
    deltaEps, state_k_old, state_k_trial, deltaEps_e_fixed, deltaEps_p_fixed,
    sig_fixed, state_k_new);

  // Do "consistency bisection"
  // std::cout << "\t Doing consistencyBisection\n";
  /*
  state_k_new = state_k_old;
  status = consistencyBisectionSimplified(
    deltaEps, state_k_old, state_k_trial, deltaEps_e_fixed, deltaEps_p_fixed,
    sig_fixed, state_k_new);
  */

  #ifdef DEBUG_INTERNAL_VAR_EVOLUTION
    std::cout << "computeSubstep: "
              << " ep_v_old = " << state_k_old.ep_v
              << " ep_v_new = " << state_k_new.ep_v 
              << " Xbar_old = " << -state_k_old.capX
              << " Xbar_new = " << -state_k_new.capX
              << std::endl;
  #endif

#ifdef DEBUG_YIELD_BISECTION_R
  std::cout << "after_consistency_bisection  = 1" << std::endl;
  std::cout << "I1 = " << state_k_new.I1 << std::endl;
  std::cout << "sqrt_J2 = " << state_k_new.sqrt_J2 << std::endl;
#endif

  return status;

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
TabularPlasticityCap::nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                          const ModelState_TabularCap& state_k_old,
                          const ModelState_TabularCap& state_k_trial,
                          Uintah::Matrix3& sig_fixed,
                          Uintah::Matrix3& elasticStrain_inc_fixed,
                          Uintah::Matrix3& plasticStrain_inc_fixed)
{

  Status status = Status::SUCCESS;

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

      status = Status::TOO_LARGE_PLASTIC_STRAIN; // The plastic volume strain is too large, try again
    }
  }

#ifdef CHECK_PLASTIC_RATE
  ModelState_TabularCap state_plastic_rate(state_k_old);
  state_plastic_rate.stressTensor = sig_fixed;
  state_plastic_rate.updateStressInvariants();
  Matrix3 df_dsigma;
  d_yield->eval_df_dsigma(Util::Identity, &state_plastic_rate, df_dsigma);
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
  ModelState_TabularCap state_test(state_k_old);
  state_test.stressTensor = sig_fixed;
  state_test.updateStressInvariants();

  Matrix3 df_dsigma;
  d_yield->eval_df_dsigma(Util::Identity, &state_test, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  std::cout << "df_dsigma = " << df_dsigma << std::endl;
  std::cout << "ratio = [" << plasticStrain_inc_fixed(0, 0) / df_dsigma(0, 0)
            << "," << plasticStrain_inc_fixed(1, 1) / df_dsigma(1, 1) << ","
            << plasticStrain_inc_fixed(2, 2) / df_dsigma(2, 2) << std::endl;

  // Compute CN = C:df_dsigma
  double lambda = state_test.bulkModulus - 2.0 / 3.0 * state_test.shearModulus;
  double mu = state_test.shearModulus;
  Matrix3 CN = Util::Identity * (lambda * df_dsigma.Trace()) + df_dsigma * (2.0 * mu);
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
  ModelState_TabularCap state_sig_test(state_k_old);
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
TabularPlasticityCap::firstOrderHardeningUpdate(const Matrix3& deltaEps_new,
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

  Matrix3 M_n, M_np1;
  d_yield->eval_df_dsigma(Util::Identity, &state_n, M_n);
  d_yield->eval_df_dsigma(Util::Identity, &state_np1, M_np1);
  M_n /= M_n.Norm();
  double norm_M_np1 = M_np1.Norm();
  M_np1 /= norm_M_np1;
  double angle_M_n_np1 = std::abs(M_n.Contract(M_np1) - 1.0);

  if (angle_M_n_np1 > 1.0e-6) {
    status = Status::TOO_LARGE_YIELD_NORMAL_CHANGE;

    std::cout << "M_n = " << M_n << std::endl;
    std::cout << "M_np1 = " << M_np1 << std::endl;
    std::cout << "angle =" << angle_M_n_np1 << "\n";

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

    // Compute volumetric plastic strain
    Matrix3 eps_p = state_n.plasticStrainTensor + deltaEps_p_fixed;
    double ep_v = eps_p.Trace();

    // Compute eps_p and ev_v^p
    Matrix3 Eps_p_n = state_n.plasticStrainTensor + M_n * Gamma_n ;
    double Eps_p_v_n = state_n.ep_v + M_n.Trace() * Gamma_n;
    Matrix3 Eps_e_n = deltaEps_new - Eps_p_n; 
    Matrix3 Eps_p_np1 = state_n.plasticStrainTensor + M_np1 * Gamma_np1 ;
    double Eps_p_v_np1 = state_n.ep_v + M_np1.Trace() * Gamma_np1;
    Matrix3 Eps_e_np1 = deltaEps_new - Eps_p_np1; 
    std::cout << "Eps_p_n = " << Eps_p_n << "\n"; 
    std::cout << "eps_p = " << eps_p << "\n"; 
    std::cout << "Eps_e_n = " << Eps_e_n << "\n"; 
    std::cout << "Eps_p_np1 = " << Eps_p_np1 << "\n"; 
    std::cout << "Eps_e_np1 = " << Eps_e_np1 << "\n"; 
    std::cout << "Eps_p_v_n = " << Eps_p_v_n 
              << " Eps_p_v_np1 = " << Eps_p_v_np1 << " ep_v = " << ep_v << "\n";
  } else {

    // Compute elastic moduli at sigma^F
    auto KG_dKdG_np1 = d_elastic->getElasticModuliAndDerivatives(&state_np1);
    auto KG_np1 = KG_dKdG_np1.first;
    auto dKdG_np1 = KG_dKdG_np1.second;
    auto p_np1 = state_np1.I1/3.0;
    auto S_np1 = state_np1.deviatoricStressTensor;
    Matrix3 Z_np1 = computeZMatrix(KG_np1, dKdG_np1, p_np1, S_np1, M_np1);
    //std::cout << "Z_np1 = " << Z_np1 << std::endl;

    // Compute CM = C:M
    double mu_np1 = KG_np1.shearModulus;
    double lambda_np1 = KG_np1.bulkModulus - 2.0 / 3.0 * mu_np1;
    Matrix3 CM_np1 = Util::Identity * (lambda_np1 * M_np1.Trace()) + M_np1 * (2.0 * mu_np1);
    //std::cout << "CM_np1 = " << CM_np1 << std::endl;

    // Compute P = CM + Z;
    Matrix3 P_np1 = CM_np1 + Z_np1;
    //std::cout << "P_np1 = " << P_np1 << std::endl;

    // Compute Gamma (sig_trial - sig_n+1):M_n+1/(P_n:M_n+1)
    Matrix3 sig_trial = state_k_trial.stressTensor;
    double Gamma_np1 = (sig_trial - sig_fixed).Contract(M_np1)/P_np1.Contract(M_np1);
    //std::cout << " Gamma_np1 = " << Gamma_np1 << "\n";

    // Compute PN = P:N (note that we are using associated plasticity)
    double PN_np1 = P_np1.Contract(M_np1);

    // Compute unscaled hardening modulus
    double df_dep_v = d_yield->computeVolStrainDerivOfYieldFunction(&state_np1,
      nullptr, nullptr, d_capX);
    double H = -(df_dep_v * M_np1.Trace())/norm_M_np1;

    // Compute Gamma^H
    double Gamma_H = (Gamma_np1*PN_np1)/(PN_np1 + H);

    // Compute updated stress
    Matrix3 sig_new = sig_trial - P_np1 * Gamma_H;
    //std::cout << "sig_new = " << sig_new << std::endl;

    // Compute plastic strain increments
    Matrix3 delta_eps_p = M_np1 * Gamma_H;
    double delta_eps_p_v = delta_eps_p.Trace();

    // Compute the elastic strain increment
    Matrix3 delta_eps_e = deltaEps_new - delta_eps_p; 

    // Recompute internal variables
    status = computeInternalVariables(state_np1, delta_eps_p_v);
    if (status != Status::SUCCESS) {
      state_k_new = state_k_old;
      proc0cout << "computeInternalVariables has failed." << std::endl;
      return status;
    }

    // Update the new state
    state_k_new = state_np1;
    state_k_new.stressTensor = sig_new;
    state_k_new.elasticStrainTensor += delta_eps_e;
    state_k_new.plasticStrainTensor += delta_eps_p;
    state_k_new.ep_v += delta_eps_p_v;
    computeElasticProperties(state_k_new);

    // Update the cumulative equivalent plastic strain
    Uintah::Matrix3 delta_eps_p_dev =
      delta_eps_p - Util::Identity * (delta_eps_p_v / 3.0);
    state_k_new.ep_cum_eq += 
      Util::sqrt_two_third * std::sqrt(delta_eps_p_dev.Contract(delta_eps_p_dev));

    #ifdef DEBUG_INTERNAL_VAR_EVOLUTION
      std::cout << "consistencyBisection: " << std::endl
                << "\t state_old = " << state_k_old << std::endl
                << "\t state_new = " << state_k_new << std::endl;
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
TabularPlasticityCap::consistencyBisectionSimplified(const Matrix3& deltaEps_new,
                                      const ModelState_TabularCap& state_k_old,
                                      const ModelState_TabularCap& state_k_trial,
                                      const Matrix3& deltaEps_e_fixed,
                                      const Matrix3& deltaEps_p_fixed,
                                      const Matrix3& sig_fixed,
                                      ModelState_TabularCap& state_k_new)
{
  // bisection convergence tolerance on eta (if changed, change imax)
  const double TOLERANCE = d_consistency_bisection_tolerance;
  // imax = ceil(-10.0*log(TOL)); // Update this if TOL changes
  const int IMAX = d_max_bisection_iterations;

  // Get the old state
  Matrix3 sig_old = state_k_old.stressTensor;
  Matrix3 eps_e_old = state_k_old.elasticStrainTensor;
  Matrix3 eps_p_old = state_k_old.plasticStrainTensor;

  // Get the fixed non-hardening return state and compute invariants
  double deltaEps_p_v_fixed = deltaEps_p_fixed.Trace();

  // Create a state for the fixed non-hardening yield surface state
  // and update only the stress and plastic strain
  ModelState_TabularCap state_k_fixed(state_k_old);
  state_k_fixed.stressTensor = sig_fixed;
  state_k_fixed.elasticStrainTensor = eps_e_old + deltaEps_e_fixed;
  state_k_fixed.plasticStrainTensor = eps_p_old + deltaEps_p_fixed;

  // Initialize the new consistently updated state
  Matrix3 sig_fixed_new = sig_fixed;
  Matrix3 deltaEps_e_fixed_new = deltaEps_e_fixed;
  Matrix3 deltaEps_p_fixed_new = deltaEps_p_fixed;
  double deltaEps_p_v_fixed_new = deltaEps_p_v_fixed;

  #ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
    double norm_deltaEps_p_fixed = deltaEps_p_fixed.Norm();
    double norm_deltaEps_p_fixed_new = norm_deltaEps_p_fixed;
  #endif

  // Set up a local trial state
  ModelState_TabularCap state_trial_local(state_k_trial);

  // Start loop
  int ii = 1;
  double eta_lo = 0.0, eta_hi = 1.0, eta_mid = 0.5;
  Status status = Status::SUCCESS;

  while (std::abs(eta_hi - eta_lo) > TOLERANCE) {

    #ifdef DEBUG_YIELD_BISECTION_R
      std::cout << "consistency_iter = " << ii << std::endl;
      std::cout << "eta_hi = " << eta_hi << std::endl;
      std::cout << "eta_lo = " << eta_lo << std::endl;
    #endif

    // Reset the local trial state
    state_trial_local = state_k_trial;

    // Compute the volumetric plastic strain at eta = eta_mid
    eta_mid = 0.5 * (eta_lo + eta_hi);
    double deltaEps_p_v_mid = eta_mid * deltaEps_p_v_fixed;

    // Update the elastic moduli at eta = eta_mid in the local trial state
    // and update stress invariants
    computeElasticProperties(state_trial_local, deltaEps_p_v_mid);

    // Update the internal variables at eta = eta_mid in the local trial state
    status = computeInternalVariables(state_trial_local, deltaEps_p_v_mid);
    if (status != Status::SUCCESS) {
      state_k_new = state_k_old;
      proc0cout << "computeInternalVariables has failed." << std::endl;
      return status;
    }

    // Test the yield condition to check whether the yield surface moves beyond
    // the trial stress state when the internal variables are changed.
    // If the yield surface is too big, the plastic strain is reduced
    // by bisecting <eta> and the loop is repeated.
    auto yield = d_yield->evalYieldCondition(&state_trial_local);

    // If the local trial state is inside the updated yield surface the yield
    // condition evaluates to "elastic".  We need to reduce the size of the
    // yield surface by decreasing the plastic strain increment.
    // Elastic or on yield surface
    if (yield.second == Util::YieldStatus::IS_ELASTIC) {
      eta_hi = eta_mid;
      ii++;
      continue;
    }

    #ifdef CHECK_RETURN_ALIGNMENT
     // std::cout << "consistency_iter = " << ii << std::endl;
     // std::cout << "eta_hi = " << eta_hi << std::endl;
     // std::cout << "eta_lo = " << eta_lo << std::endl;
    #endif

    // At this point, state_trial_local contains the trial stress, the plastic
    // strain at the beginning of the timestep, the updated elastic moduli,
    // and the updated values of the internal  variables
    // The yield surface depends only on X.  We will compute the updated
    // location of the yield surface based on the updated internal variables
    // and do a non-hardening return to that yield surface.
    ModelState_TabularCap state_k_updated(state_k_old);
    state_k_updated.bulkModulus = state_trial_local.bulkModulus;
    state_k_updated.shearModulus = state_trial_local.shearModulus;
    state_k_updated.capX = state_trial_local.capX;
    state_k_updated.updateStressInvariants();
    status = nonHardeningReturn(deltaEps_new, state_k_updated,
                                   state_trial_local, sig_fixed_new,
                                   deltaEps_e_fixed_new, deltaEps_p_fixed_new);
    if (status != Status::SUCCESS) {
      proc0cout << "nonHardeningReturn inside consistencyBisection failed."
                << std::endl;
      return status;
    }

    // Check whether the isotropic component of the return has changed sign, as
    // this would indicate that the cap apex has moved past the trial stress,
    // indicating too much plastic strain in the return.
    Matrix3 sig_trial = state_trial_local.stressTensor;
    double diff_trial_fixed_new = (sig_trial - sig_fixed_new).Trace();
    double diff_trial_fixed = (sig_trial - sig_fixed).Trace();
    if (std::signbit(diff_trial_fixed_new) != std::signbit(diff_trial_fixed)) {
      eta_hi = eta_mid;
      ii++;
      continue;
    }

    // Compare magnitude of plastic strain with prior update
    #ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
      norm_deltaEps_p_fixed_new = deltaEps_p_fixed_new.Norm();
      norm_deltaEps_p_fixed = eta_mid * deltaEps_p_fixed.Norm();
    #endif

    deltaEps_p_v_fixed_new = deltaEps_p_fixed_new.Trace();
    deltaEps_p_v_fixed = eta_mid * deltaEps_p_fixed.Trace();

    #ifdef CHECK_CONSISTENCY_BISECTION_CONVERGENCE
      std::cout << "eta_mid = " << eta_mid << " eta_mid*||deltaEps_p_fixed|| = "
                << eta_mid * norm_deltaEps_p_fixed
                << " ||deltaEps_p_fixed_new|| = " << norm_deltaEps_p_fixed_new
                << " ratio = "
                << eta_mid * norm_deltaEps_p_fixed / norm_deltaEps_p_fixed_new
                << std::endl;
      std::cout << " delta_eps_p_v_mid = " << deltaEps_p_v_mid
                << " delta_eps_p_v_fixed_new = " << deltaEps_p_v_fixed_new
                << " ratio = " << deltaEps_p_v_mid / deltaEps_p_fixed_new.Trace()
                << std::endl;
    #endif

    // if (norm_deltaEps_p_fixed_new > eta_mid*norm_deltaEps_p_fixed) {
    if (std::abs(deltaEps_p_v_fixed_new) >
        eta_mid * std::abs(deltaEps_p_v_fixed)) {
      eta_lo = eta_mid;
    } else {
      eta_hi = eta_mid;
    }

    // Increment i and check
    ii++;
    if (ii > IMAX) {
      state_k_new = state_k_old;
      proc0cout << "Consistency bisection has failed because ii > IMAX." << ii
                << " > " << IMAX << std::endl;
      return Status::TOO_MANY_CONSISTENCY_ITERATIONS; // bool isSuccess = false;
    }

  } // end  while (std::abs(eta_hi - eta_lo) > TOLERANCE);

  // Set the new state to the old trial state
  // The volumetric strain may not have converged so recompute internal
  // variables
  state_k_new = state_k_trial;
  computeElasticProperties(state_k_new, deltaEps_p_v_fixed_new);
  status = computeInternalVariables(state_k_new, deltaEps_p_v_fixed_new);
  if (status != Status::SUCCESS) {
    state_k_new = state_k_old;
    proc0cout << "computeInternalVariables has failed." << std::endl;
    return status;
  }

  // Update the stress and plastic strain of the new state +  the elastic moduli
  state_k_new.stressTensor = sig_fixed_new;
  state_k_new.elasticStrainTensor = eps_e_old + deltaEps_e_fixed_new;
  state_k_new.plasticStrainTensor = eps_p_old + deltaEps_p_fixed_new;
  computeElasticProperties(state_k_new);

  #ifdef DEBUG_YIELD_BISECTION_R
    std::cout << "pbar_w_before_consistency_3 = " << 3.0 * state_k_old.pbar_w
              << std::endl;
    std::cout << "pbar_w_after_consistency_3 = " << 3.0 * state_k_new.pbar_w
              << std::endl;
    std::cout << "K_before_consistency = " << state_k_old.bulkModulus
              << std::endl;
    std::cout << "K_after_consistency = " << state_k_new.bulkModulus << std::endl;
    std::cout << "I1_before_consistency = " << state_k_old.stressTensor.Trace()
              << std::endl;
    std::cout << "I1_after_consistency = " << state_k_new.stressTensor.Trace()
              << std::endl;
    std::cout << "I1_eff_before_consistency = " << state_k_old.I1_eff
              << std::endl;
    std::cout << "I1_eff_after_consistency = " << state_k_new.I1_eff << std::endl;
  #endif

  // Update the cumulative equivalent plastic strain
  double deltaEps_p_v = deltaEps_p_fixed_new.Trace();
  Uintah::Matrix3 deltaEps_p_dev =
    deltaEps_p_fixed_new - Util::Identity * (deltaEps_p_v / 3.0);
  state_k_new.ep_cum_eq =
    state_k_old.ep_cum_eq +
    std::sqrt(2.0 / 3.0 * deltaEps_p_dev.Contract(deltaEps_p_dev));

  #ifdef DEBUG_INTERNAL_VAR_EVOLUTION
    std::cout << "consistencyBisection: " << std::endl
              << "\t state_old = " << state_k_old << std::endl
              << "\t state_new = " << state_k_new << std::endl;
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
                                               const double& delta_eps_p_v)
{
  ModelState_TabularCap tempState = state;
  tempState.ep_v += delta_eps_p_v;

  // Update the hydrostatic compressive strength
  double X_new = d_capX->computeInternalVariable(&tempState);

  // Update the state with new values of the internal variables
  state.capX = X_new;

  return Status::SUCCESS;
}

void
TabularPlasticityCap::addRequiresDamageParameter(Task* task, const MPMMaterial* matl,
                                  const PatchSet*) const
{
  // Require the damage parameter
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pRemoveLabel_preReloc, matlset,
                 Ghost::None);
}

void
TabularPlasticityCap::getDamageParameter(const Patch* patch, ParticleVariable<int>& damage,
                          int matID, DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  // Get the damage parameter
  ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, pRemoveLabel_preReloc, pset);

  // Loop over the particle in the current patch.
  for (int& iter : *pset) {
    damage[iter] = pLocalized[iter];
  }
}

void
TabularPlasticityCap::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
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
TabularPlasticityCap::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                              const bool, const bool) const
{
  std::cout
    << "NO Implicit VERSION OF addComputesAndRequires EXISTS YET FOR TabularPlasticityCap"
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
TabularPlasticityCap::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset, DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "Material conversion after failure not implemented for TabularPlasticityCap.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

/*---------------------------------------------------------------------------------------
 * MPMICE Hooks
 *---------------------------------------------------------------------------------------*/
double
TabularPlasticityCap::computeRhoMicroCM(double pressure, const double p_ref,
                         const MPMMaterial* matl, double temperature,
                         double rho_guess)
{
  double rho_cur = 0.0;
  if (rho_cur < 1.0) {
    std::ostringstream out;
    out << "MPMICE hooks not implemented for TabularPlasticityCap.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  return rho_cur;
}

void
TabularPlasticityCap::computePressEOSCM(double rho_cur, double& pressure, double p_ref,
                         double& dp_drho, double& soundSpeedSq,
                         const MPMMaterial* matl, double temperature)
{
  std::ostringstream out;
  out << "MPMICE hooks not implemented for TabularPlasticityCap.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

double
TabularPlasticityCap::getCompressibility()
{
  double comp = 0.0;
  if (comp < 1.0) {
    std::ostringstream out;
    out << "MPMICE hooks not implemented for TabularPlasticityCap.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return comp;
}
