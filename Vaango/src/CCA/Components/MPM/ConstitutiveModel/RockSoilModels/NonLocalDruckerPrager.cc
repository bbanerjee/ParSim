/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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
// #include </usr/include/valgrind/callgrind.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/NonLocalDruckerPrager.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Weibull.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sci_values.h>
#include <string>

using namespace Uintah;

NonLocalDruckerPrager::NonLocalDruckerPrager(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{

  ps->require("alpha", d_initialData.alpha);
  ps->require("alpha_p", d_initialData.alpha_p);
  ps->require("k_o", d_initialData.k_o);
  ps->require("bulk_modulus", d_initialData.bulk_modulus);
  ps->require("shear_modulus", d_initialData.shear_modulus);
  ps->getWithDefault("l_nonlocal", d_initialData.l_nonlocal, 0.0);
  ps->getWithDefault("h_local", d_initialData.h_local, 0.0);
  ps->getWithDefault("h_nonlocal", d_initialData.h_nonlocal, 0.0);
  ps->getWithDefault("minimum_yield_stress",
                     d_initialData.minimum_yield_stress,
                     0.0);
  ps->getWithDefault("initial_xstress", d_initialData.initial_xstress, 0.0);
  ps->getWithDefault("initial_ystress", d_initialData.initial_ystress, 0.0);
  ps->getWithDefault("initial_zstress", d_initialData.initial_zstress, 0.0);
  ps->getWithDefault("hardening_type", d_initialData.hardening_type, 0);
  ps->get("k_o_dist", wdist.WeibDist);
  WeibullParser(wdist);

  initializeLocalMPMLabels();
}

NonLocalDruckerPrager::NonLocalDruckerPrager(const NonLocalDruckerPrager* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.alpha                = cm->d_initialData.alpha;
  d_initialData.alpha_p              = cm->d_initialData.alpha_p;
  d_initialData.k_o                  = cm->d_initialData.k_o;
  d_initialData.bulk_modulus         = cm->d_initialData.bulk_modulus;
  d_initialData.shear_modulus        = cm->d_initialData.shear_modulus;
  d_initialData.l_nonlocal           = cm->d_initialData.l_nonlocal;
  d_initialData.h_local              = cm->d_initialData.h_local;
  d_initialData.h_nonlocal           = cm->d_initialData.h_nonlocal;
  d_initialData.minimum_yield_stress = cm->d_initialData.minimum_yield_stress;
  d_initialData.initial_xstress      = cm->d_initialData.initial_xstress;
  d_initialData.initial_ystress      = cm->d_initialData.initial_ystress;
  d_initialData.initial_zstress      = cm->d_initialData.initial_zstress;
  d_initialData.hardening_type       = cm->d_initialData.hardening_type;

  wdist.WeibMed    = cm->wdist.WeibMed;
  wdist.WeibMod    = cm->wdist.WeibMod;
  wdist.WeibRefVol = cm->wdist.WeibRefVol;
  wdist.WeibSeed   = cm->wdist.WeibSeed;
  wdist.Perturb    = cm->wdist.Perturb;
  wdist.WeibDist   = cm->wdist.WeibDist;

  initializeLocalMPMLabels();
}

NonLocalDruckerPrager::~NonLocalDruckerPrager()
{
  VarLabel::destroy(etaLabel);
  VarLabel::destroy(etaLabel_preReloc);
  VarLabel::destroy(eta_nlLabel);
  VarLabel::destroy(eta_nlLabel_preReloc);
  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);
  VarLabel::destroy(k_o_distLabel);
  VarLabel::destroy(k_o_distLabel_preReloc);
}

void
NonLocalDruckerPrager::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "nonlocal_drucker_prager");
  }

  cm_ps->appendElement("alpha", d_initialData.alpha);
  cm_ps->appendElement("alpha_p", d_initialData.alpha_p);
  cm_ps->appendElement("k_o", d_initialData.k_o);
  cm_ps->appendElement("bulk_modulus", d_initialData.bulk_modulus);
  cm_ps->appendElement("shear_modulus", d_initialData.shear_modulus);
  cm_ps->appendElement("l_nonlocal", d_initialData.l_nonlocal);
  cm_ps->appendElement("h_local", d_initialData.h_local);
  cm_ps->appendElement("h_nonlocal", d_initialData.h_nonlocal);
  cm_ps->appendElement("minimum_yield_stress",
                       d_initialData.minimum_yield_stress);
  cm_ps->appendElement("initial_xstress", d_initialData.initial_xstress);
  cm_ps->appendElement("initial_ystress", d_initialData.initial_ystress);
  cm_ps->appendElement("initial_zstress", d_initialData.initial_zstress);
  cm_ps->appendElement("hardening_type", d_initialData.hardening_type);
  cm_ps->appendElement("k_o_Perturb", wdist.Perturb);
  cm_ps->appendElement("k_o_Med", wdist.WeibMed);
  cm_ps->appendElement("k_o_Mod", wdist.WeibMod);
  cm_ps->appendElement("k_o_RefVol", wdist.WeibRefVol);
  cm_ps->appendElement("k_o_Seed", wdist.WeibSeed);
  cm_ps->appendElement("k_o_dist", wdist.WeibDist);
}

std::unique_ptr<ConstitutiveModel>
NonLocalDruckerPrager::clone()
{
  return std::make_unique<NonLocalDruckerPrager>(*this);
}

void
NonLocalDruckerPrager::initializeCMData(const Patch* patch,
                                        const MPMMaterial* matl,
                                        DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.

  initSharedDataForExplicit(patch, matl, new_dw);
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  Uintah::Weibull weibGen(wdist.WeibMed,
                          wdist.WeibMod,
                          wdist.WeibRefVol,
                          wdist.WeibSeed,
                          wdist.WeibMod);
  std::cout << "Weibull Variables for PEAKI1I: (initialize CMData)\n"
            << "Median:            " << wdist.WeibMed
            << "\nModulus:         " << wdist.WeibMod
            << "\nReference Vol:   " << wdist.WeibRefVol
            << "\nSeed:            " << wdist.WeibSeed
            << "\nPerturb?:        " << wdist.Perturb << std::endl;

  constParticleVariable<double> pVolume;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);

  ParticleVariable<double> eta, eta_nl, pPlasticStrain;
  ParticleVariable<double> k_o_dist;
  new_dw->allocateAndPut(eta, etaLabel, pset);
  new_dw->allocateAndPut(eta_nl, eta_nlLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);
  new_dw->allocateAndPut(k_o_dist, k_o_distLabel, pset);
  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    eta[*iter]            = 0.0;
    eta_nl[*iter]         = 0.0;
    pPlasticStrain[*iter] = 0.0;
    if (wdist.Perturb) {
      k_o_dist[*iter] = weibGen.rand(pVolume[*iter]);
    }
  }

  computeStableTimestep(patch, matl, new_dw);
}

///////////////////////////////////////////////////////////////////////////
/*! Allocate data required during the conversion of failed particles
    from one material to another */
///////////////////////////////////////////////////////////////////////////
void
NonLocalDruckerPrager::allocateCMDataAddRequires(Task* task,
                                                 const MPMMaterial* matl,
                                                 const PatchSet* patches,
                                                 [[maybe_unused]] MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
  task->needs(Task::NewDW, etaLabel_preReloc, matlset, Ghost::None);
  task->needs(Task::NewDW, eta_nlLabel_preReloc, matlset, Ghost::None);
  task->needs(Task::NewDW, k_o_distLabel_preReloc, matlset, Ghost::None);
}

void
NonLocalDruckerPrager::allocateCMDataAdd([[maybe_unused]] DataWarehouse* new_dw,
                                         [[maybe_unused]] ParticleSubset* addset,
                                         [[maybe_unused]] ParticleLabelVariableMap* newState,
                                         [[maybe_unused]] ParticleSubset* delset,
                                         DataWarehouse*)
{
}

void
NonLocalDruckerPrager::computeStableTimestep(const Patch* patch,
                                             const MPMMaterial* matl,
                                             DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi   = matl->getDWIndex();
  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;
  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);
  double bulk  = d_initialData.bulk_modulus;
  double shear = d_initialData.shear_modulus;

  for (int idx : *pset) {
    // Compute wave speed + particle velocity at each particle,
    // store the maximum

    c_dil     = sqrt((bulk + 4.0 * shear / 3.0) * pVolume[idx] / pMass[idx]);
    WaveSpeed = Vector(Max(c_dil + fabs(pVelocity[idx].x()), WaveSpeed.x()),
                       Max(c_dil + fabs(pVelocity[idx].y()), WaveSpeed.y()),
                       Max(c_dil + fabs(pVelocity[idx].z()), WaveSpeed.z()));
  }
  WaveSpeed       = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  if (delT_new < 1.e-12) {
    new_dw->put(delt_vartype(DBL_MAX), lb->delTLabel, patch->getLevel());
  } else {
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
  }
}

void
NonLocalDruckerPrager::computeStressTensor(const PatchSubset* patches,
                                           const MPMMaterial* matl,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw)
{

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Matrix3 Identity, D;
    double J;
    Identity.Identity();
    double c_dil = 0.0, se = 0.0;
    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    Vector dx = patch->dCell();
    // double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    int dwi = matl->getDWIndex();

    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<Matrix3> pDefGrad_new, velGrad;
    constParticleVariable<Matrix3> pDefGrad;
    constParticleVariable<Matrix3> stress_old;
    ParticleVariable<Matrix3> stress_new;
    constParticleVariable<Point> px;
    constParticleVariable<double> pMass;
    constParticleVariable<double> pVolume;
    ParticleVariable<double> p_q;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pSize;
    ParticleVariable<double> pdTdt;
    constParticleVariable<double> eta_old;
    constParticleVariable<double> eta_nl_old;
    ParticleVariable<double> eta_new;
    ParticleVariable<double> eta_nl_new;
    constParticleVariable<double> pPlasticStrain;
    ParticleVariable<double> pPlasticStrain_new;
    constParticleVariable<double> k_o_dist;
    ParticleVariable<double> k_o_dist_new;

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);
    new_dw->allocateAndPut(pPlasticStrain_new,
                           pPlasticStrainLabel_preReloc,
                           pset);
    Ghost::GhostType gac = Ghost::AroundCells;
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);
    old_dw->get(stress_old, lb->pStressLabel, pset);
    old_dw->get(eta_old, etaLabel, pset);
    old_dw->get(eta_nl_old, eta_nlLabel, pset);
    old_dw->get(k_o_dist, k_o_distLabel, pset);

    new_dw->get(pVolume, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(velGrad, lb->pVelGradLabel_preReloc, pset);

    new_dw->allocateAndPut(stress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);

    new_dw->allocateAndPut(eta_new, etaLabel_preReloc, pset);
    new_dw->allocateAndPut(eta_nl_new, eta_nlLabel_preReloc, pset);
    new_dw->allocateAndPut(k_o_dist_new, k_o_distLabel_preReloc, pset);

    k_o_dist_new.copyData(k_o_dist);

    ParticleVariable<Matrix3> rotation, trial_stress;
    ParticleVariable<double> f_trial, pdlambda, rho_cur;
    ParticleVariable<int> softened;
    ParticleVariable<double> k_o;
    new_dw->allocateTemporary(k_o, pset);
    new_dw->allocateTemporary(rotation, pset);
    new_dw->allocateTemporary(pdlambda, pset);
    new_dw->allocateTemporary(trial_stress, pset);
    new_dw->allocateTemporary(f_trial, pset);
    new_dw->allocateTemporary(rho_cur, pset);
    new_dw->allocateTemporary(softened, pset);

    const double alpha          = d_initialData.alpha;
    const double alpha_p        = d_initialData.alpha_p;
    const double k_o_const      = d_initialData.k_o;
    const double bulk           = d_initialData.bulk_modulus;
    const double shear          = d_initialData.shear_modulus;
    const double l_nonlocal     = d_initialData.l_nonlocal;
    double h_local              = d_initialData.h_local;
    double h_nonlocal           = d_initialData.h_nonlocal;
    double minimum_yield_stress = d_initialData.minimum_yield_stress;
    double initial_xstress      = d_initialData.initial_xstress;
    double initial_ystress      = d_initialData.initial_ystress;
    double initial_zstress      = d_initialData.initial_zstress;
    int hardening_type          = d_initialData.hardening_type;

    // assemble the initial stress tensor
    Matrix3 initial_stress(0.0);
    initial_stress.set(0, 0, initial_xstress);
    initial_stress.set(1, 1, initial_ystress);
    initial_stress.set(2, 2, initial_zstress);

    // create node data for the plastic multiplier field
    NCVariable<double> gdlambda, gmat;
    int NGhost = (int)std::ceil(l_nonlocal / dx.maxComponent());
    new_dw->allocateTemporary(gdlambda, patch, gac, NGhost);
    new_dw->allocateTemporary(gmat, patch, gac, NGhost);
    gdlambda.initialize(0.0);
    gmat.initialize(0.0);
    double rho_orig = matl->getInitialDensity();
    Matrix3 tensorL(0.0);

    for (int idx : *pset) {
      // cout<<"BEGIN eta_nl_old["<<idx<<"]= "<<eta_nl_old[idx]<<endl;
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      J = pDefGrad_new[idx].Determinant();
      // Update particle volumes
      rho_cur[idx] = rho_orig / J;

      // Compute the rate of deformation tensor
      Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose()) * .5;
      // NEED TO FIND R
      Matrix3 tensorR, tensorU;
      pDefGrad_new[idx].polarDecompositionRMB(tensorU, tensorR);
      // Compute the Jacobian of the deformation gradient

      rotation[idx] = tensorR;

      // This is the previous timestep Cauchy stress
      // unrotated tensorSig=R^T*pstress*R
      Matrix3 unrotated_stress =
        (tensorR.Transpose()) * ((stress_old[idx] + initial_stress) * tensorR);
      D = (tensorR.Transpose()) * (D * tensorR);
      if (wdist.Perturb) {
        k_o[idx] = k_o_dist[idx];
      } else {
        k_o[idx]          = k_o_const;
        k_o_dist_new[idx] = k_o_const;
      }
      double lame      = bulk - 2.0 / 3.0 * shear;
      double eta_in    = eta_old[idx];
      double eta_nl_in = eta_nl_old[idx];
      trial_stress[idx] =
        unrotated_stress +
        (Identity * lame * (D.Trace() * delT) + D * delT * 2.0 * shear);

      // evaluate yield function:
      f_trial[idx]            = YieldFunction(trial_stress[idx],
                                   alpha,
                                   k_o[idx],
                                   eta_in,
                                   eta_nl_in,
                                   hardening_type);
      eta_new[idx]            = eta_old[idx];
      eta_nl_new[idx]         = eta_nl_old[idx];
      pPlasticStrain_new[idx] = pPlasticStrain[idx];
      pdlambda[idx]           = 0.0;
      softened[idx]           = 0;
      if (f_trial[idx] < 0.0) {
        stress_new[idx] = trial_stress[idx];
      } else {

        // define frequently used constants:

        double I1_trial, J2_trial;
        Matrix3 S_trial;
        Matrix3 S_old;
        double J2_old, I1_old;
        computeInvariants(trial_stress[idx], S_trial, I1_trial, J2_trial);
        computeInvariants(stress_old[idx], S_old, I1_old, J2_old);

        Matrix3 A, M;
        // comput M and N using the new stress estimate:
        // M = Identity*alpha_p +
        // S_trial*(1.0/(sqrt_two*sqrt_J2_trial))*(1.0-alpha_p);
        Matrix3 Strial_hat = S_trial * (1.0 / sqrt(2.0 * J2_trial));
        M = (Identity * alpha_p + Strial_hat * (sqrt(2.0) / 2.0)) *
            (1.0 / sqrt(0.5 + 3.0 * alpha_p * alpha_p));
        A                  = (Identity * lame * (M.Trace()) + M * 2.0 * shear);
        double dlambda_new = 0.0;

        double current_yield_strength;
        // do a local iteration in order to provide a good initial guess for the
        // nonlocal iteration
        double factorp = sqrt(0.5 + 3.0 * alpha_p * alpha_p);
        dlambda_new =
          (f_trial[idx]) /
          (shear / factorp + alpha * 9.0 * bulk * alpha_p / factorp + h_local);

        current_yield_strength =
          k_o[idx] + (eta_old[idx] + h_local * dlambda_new);

        pdlambda[idx] = dlambda_new;

        // update the stress using the the last estimate for dlambda
        current_yield_strength =
          k_o[idx] + (eta_old[idx] + h_local * dlambda_new);
        if (current_yield_strength < minimum_yield_stress) {
          eta_new[idx] = eta_old[idx];
        } else {
          eta_new[idx] = eta_old[idx] + h_local * pdlambda[idx];
        }

        // compute a new stress estimate
        stress_new[idx] = trial_stress[idx] - A * pdlambda[idx];

        // see if we are beyond the vertex
        double J2_new, I1_new;
        Matrix3 S_new;
        computeInvariants(stress_new[idx], S_new, I1_new, J2_new);
        if (current_yield_strength < 0) {
          std::cout << "zero yield strength detected" << std::endl;
          // just set deviator to zero (von-Mise ONLY!!)
          stress_new[idx] = Identity * (1.0 / 3.0) * trial_stress[idx].Trace();
        } else if (alpha * I1_new > current_yield_strength) {
          // just put the stress on the vertex
          std::cout << "stress on yield vertex" << std::endl;
          stress_new[idx] =
            Identity * (1.0 / (3.0 * alpha)) * current_yield_strength;
        }

        [[maybe_unused]] double f_new;
        if (current_yield_strength < minimum_yield_stress) {
          f_new = YieldFunction(stress_new[idx],
                                alpha,
                                minimum_yield_stress,
                                0.0,
                                0.0,
                                hardening_type);
        } else {
          f_new = YieldFunction(stress_new[idx],
                                alpha,
                                k_o[idx],
                                eta_new[idx],
                                eta_nl_new[idx],
                                hardening_type);
        }
        // this is just generating an initial estimate, so it doesn't need to
        // actually hit the yield surface
        /*if (abs(f_new)>10.0){
          cerr<<"ERROR!  did not return to yield surface"<<endl;
          cerr<<"eta_new[idx]= "<<eta_new[idx]<<endl;
          cerr<<"sqrt(J2_new)= "<<sqrt(J2_new)<<endl;
          cerr<<"I1_new= "<<I1_new<<endl;
          cerr<<"current_yield_strength= "<<current_yield_strength<<endl;
          }*/
        // cerr<<"final increment in plastic strain is"<<dlambda_new<<endl;
      }

    } // end loop over particles

    int n8or27 = flag->d_8or27;
    // first interpolate the initial guess for plastic multipliers to the grid
    // loop over plastic particles
    for (int idx : *pset) {
      // interpolate plastic multiplier to the grid
      if (!flag->d_axisymmetric) {
        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(px[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad[idx]);
      } else { // axi-symmetric kinematics
        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(px[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad[idx]);
      }
      IntVector node;
      for (int k = 0; k < n8or27; k++) {
        node = ni[k];
        if (patch->containsNode(node)) {
          gmat[node] += S[k];
          gdlambda[node] += pdlambda[idx] * S[k];
        }
      }

    } // end loop over particles

    bool done = false;

    // now begin iteration
    int q = 0;
    while (!done) {

      q++;
      bool soften_elastic = false;
      double tolerance    = 1.e-14;
      double error_max    = 0.0;
      double error        = 0.0;
      double dlambda_max  = 0.0;
      double dlambda_nl   = 0.0;
      double V_alpha      = 0.0;

      // loop over particles
      for (int idx : *pset) {
        if (f_trial[idx] > 0.0 || softened[idx] == 1) {
          dlambda_nl = 0.0;
          // use the local estimate
          double dlambda_old = pdlambda[idx];
          // evaluate nonlocal average
          // if l_nonlocal is less than the cell spacing, just stick with the
          // local value:
          EvaluateNonLocalAverage(dlambda_nl,
                                  V_alpha,
                                  pdlambda,
                                  px,
                                  gdlambda,
                                  gmat,
                                  patch,
                                  idx,
                                  l_nonlocal);
          eta_nl_new[idx] = eta_nl_old[idx] + h_nonlocal * dlambda_nl;
          if (softened[idx] == 1) {
            f_trial[idx] = YieldFunction(trial_stress[idx],
                                         alpha,
                                         k_o[idx],
                                         eta_old[idx],
                                         eta_nl_new[idx],
                                         hardening_type);
          }
          if (f_trial[idx] > 0) {
            double factorp = sqrt(0.5 + 3.0 * alpha_p * alpha_p);
            pdlambda[idx]  = (f_trial[idx] - h_nonlocal * dlambda_nl) /
                            (shear / factorp +
                             alpha * 9.0 * bulk * alpha_p / factorp + h_local);

          } else {
            pdlambda[idx] = 0.0;
          }
          // check if we have over-softened
          double current_yield_strength =
            k_o[idx] + (eta_old[idx] + h_local * pdlambda[idx]) +
            (eta_nl_old[idx] + h_nonlocal * dlambda_nl);
          if (current_yield_strength < minimum_yield_stress) {

            eta_new[idx]    = eta_old[idx];
            eta_nl_new[idx] = eta_nl_old[idx];

            double f_trial_min = YieldFunction(trial_stress[idx],
                                               alpha,
                                               minimum_yield_stress,
                                               0.0,
                                               0.0,
                                               hardening_type);
            double factorp     = sqrt(0.5 + 3.0 * alpha_p * alpha_p);
            pdlambda[idx] =
              (f_trial_min) /
              (shear / factorp + alpha * 9.0 * bulk * alpha_p / factorp);

          } else {
            eta_new[idx] = eta_old[idx] + h_local * pdlambda[idx];

            eta_nl_new[idx] = eta_nl_old[idx] + h_nonlocal * dlambda_nl;
          }

          // compute the return derection based upon the trial stress
          Matrix3 S_trial;
          double I1_trial, J2_trial;
          computeInvariants(trial_stress[idx], S_trial, I1_trial, J2_trial);
          Matrix3 Strial_hat = S_trial * (1.0 / sqrt(2.0 * J2_trial));
          double factorp     = sqrt(0.5 + 3.0 * alpha_p * alpha_p);
          Matrix3 M = (Identity * alpha_p + Strial_hat * (sqrt(2.0) / 2.0)) *
                      (1.0 / factorp);
          double lame = bulk - 2.0 / 3.0 * shear;
          Matrix3 A   = (Identity * lame * (M.Trace()) + M * 2.0 * shear);
          // compute a new estimate for the stress
          stress_new[idx] = trial_stress[idx] - A * pdlambda[idx];

          [[maybe_unused]] double f_new;
          current_yield_strength = k_o[idx] +
                                   (eta_old[idx] + h_local * pdlambda[idx]) +
                                   (eta_nl_old[idx] + h_nonlocal * dlambda_nl);
          if (current_yield_strength < minimum_yield_stress) {
            f_new = YieldFunction(stress_new[idx],
                                  alpha,
                                  minimum_yield_stress,
                                  0.0,
                                  0.0,
                                  hardening_type);
          } else {
            f_new = YieldFunction(stress_new[idx],
                                  alpha,
                                  k_o[idx],
                                  eta_new[idx],
                                  eta_nl_new[idx],
                                  hardening_type);
          }

          // cerr<<"yield function value after nonlocal iteration "<<q-1<<" is
          // "<<f_new<<endl;

          // compute the iteration error
          error = abs(pdlambda[idx] - dlambda_old);

          // find out if this is the maximum value for dlambda or error
          if (pdlambda[idx] > dlambda_max) {
            dlambda_max = pdlambda[idx];
          }
          if (error > error_max) {
            error_max = error;
          }

          // get the interpolation data
          if (!flag->d_axisymmetric) {
            // Get the node indices that surround the cell
            interpolator->findCellAndWeights(px[idx],
                                             ni,
                                             S,
                                             pSize[idx],
                                             pDefGrad[idx]);
          } else { // axi-symmetric kinematics
            // Get the node indices that surround the cell
            interpolator->findCellAndWeights(px[idx],
                                             ni,
                                             S,
                                             pSize[idx],
                                             pDefGrad[idx]);
          }

          // replace the old value of gdlambda with the interpolated values of
          // the new estimate
          // this results in a Gauss-Seidell iteration rather than a Jacobi
          for (int k = 0; k < n8or27; k++) {
            IntVector node;
            node = ni[k];
            if (patch->containsNode(node)) {
              gdlambda[node] -= dlambda_old * S[k];
              gdlambda[node] += pdlambda[idx] * S[k];
            }
          }
          // done with plastic part

        } else {
          // update the nonlocal isv for the elastic particles
          EvaluateNonLocalAverage(dlambda_nl,
                                  V_alpha,
                                  pdlambda,
                                  px,
                                  gdlambda,
                                  gmat,
                                  patch,
                                  idx,
                                  l_nonlocal);
          eta_nl_new[idx] = eta_nl_old[idx] + h_nonlocal * dlambda_nl;
          // now evaluate the yield function with this value for eta_nl
          double f_max = YieldFunction(stress_new[idx],
                                       alpha,
                                       k_o[idx],
                                       eta_new[idx],
                                       eta_nl_new[idx],
                                       hardening_type);
          if (f_max > 0.0) {
            softened[idx]  = 1;
            soften_elastic = true;
            std::cout << "elastic softening" << std::endl;
          }

        } // end elastic part

      } // end loop particles

      // check if the error in the plastic multiplier is small and that
      // there are no softened elastic particles
      std::cout << "maximum iteration error after iteration " << q << " is "
                << error_max << std::endl;
      if (error_max < tolerance && !soften_elastic) {
        done = true;
      }

    } // end iterative loop

    // final loop over all particles
    for (int idx : *pset) {
      stress_new[idx] =
        (rotation[idx] * stress_new[idx]) * (rotation[idx].Transpose()) -
        initial_stress;

      pPlasticStrain_new[idx] = pPlasticStrain[idx] + pdlambda[idx];

      // Compute wave speed + particle velocity at each particle,
      // store the maximum

      c_dil     = sqrt((bulk + (4.0 / 3.0) * shear) / (rho_cur[idx]));
      WaveSpeed = Vector(Max(c_dil + fabs(pVelocity[idx].x()), WaveSpeed.x()),
                         Max(c_dil + fabs(pVelocity[idx].y()), WaveSpeed.y()),
                         Max(c_dil + fabs(pVelocity[idx].z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulk / rho_cur[idx]);
        p_q[idx] =
          artificialBulkViscosity(D.Trace(), c_bulk, rho_cur[idx], dx_ave);
      } else {
        p_q[idx] = 0.;
      }
      Matrix3 AvgStress = (stress_new[idx] + stress_old[idx]) * .5;

      double e = (D(0, 0) * AvgStress(0, 0) + D(1, 1) * AvgStress(1, 1) +
                  D(2, 2) * AvgStress(2, 2) +
                  2. * (D(0, 1) * AvgStress(0, 1) + D(0, 2) * AvgStress(0, 2) +
                        D(1, 2) * AvgStress(1, 2))) *
                 pVolume[idx] * delT;

      se += e;

    } // end loop over particles

    WaveSpeed       = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();

    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }

    // delete interpolator;

  } // end loop over patches
}

void
NonLocalDruckerPrager::EvaluateNonLocalAverage(
  double& dlambda_nl,
  double& V_alpha,
  ParticleVariable<double>& pdlambda,
  constParticleVariable<Point>& px,
  NCVariable<double>& gdlambda,
  NCVariable<double>& gmat,
  const Patch*& patch,
  particleIndex& idx,
  const double& l_nonlocal)
{
  Vector dx = patch->dCell();
  std::vector<double> l_nonlocal_index(3);
  l_nonlocal_index[0] = l_nonlocal / dx.x();
  l_nonlocal_index[1] = l_nonlocal / dx.y();
  l_nonlocal_index[2] = l_nonlocal / dx.z();
  if (std::min(std::min(dx[0], dx[1]), dx[2]) > 2.0 * l_nonlocal) {
    dlambda_nl = pdlambda[idx];

    V_alpha = 1;
  } else {
    Point px_index = patch->getLevel()->positionToIndex(px[idx]);
    // loop over grid nodes within the nonlocal support domain:
    double pxx = px[idx].x();
    double pxy = px[idx].y();
    double pxz = px[idx].z();
    Point nl_upper_corner(pxx + l_nonlocal, pxy + l_nonlocal, pxz + l_nonlocal);
    Point nl_lower_corner(pxx - l_nonlocal, pxy - l_nonlocal, pxz - l_nonlocal);
    Point nl_upper_index = patch->getLevel()->positionToIndex(nl_upper_corner);
    Point nl_lower_index = patch->getLevel()->positionToIndex(nl_lower_corner);
    /*
    int xmax_index = Floor(nl_upper_index.x());
    int xmin_index = std::ceil(nl_lower_index.x());
    int ymax_index = Floor(nl_upper_index.y());
    int ymin_index = std::ceil(nl_lower_index.y());
    int zmax_index = Floor(nl_upper_index.z());
    int zmin_index = std::ceil(nl_lower_index.z());
    */
    int xmax_index = (int)std::ceil(nl_upper_index.x());
    int xmin_index = (int)std::floor(nl_lower_index.x());
    int ymax_index = (int)std::ceil(nl_upper_index.y());
    int ymin_index = (int)std::floor(nl_lower_index.y());
    int zmax_index = (int)std::ceil(nl_upper_index.z());
    int zmin_index = (int)std::floor(nl_lower_index.z());
    dlambda_nl     = 0.0;
    V_alpha        = 0.0;
    for (int i = xmin_index; i < xmax_index + 1; i++) {
      for (int j = ymin_index; j < ymax_index + 1; j++) {
        for (int k = zmin_index; k < zmax_index + 1; k++) {
          IntVector index_ijk(i, j, k);
          if (patch->containsNode(index_ijk)) {
            Point x_ijk(i, j, k);
            double alpha_ijk =
              alpha_nl(px_index, x_ijk, l_nonlocal_index) * gmat[index_ijk];
            V_alpha += alpha_ijk;

            if (gmat[index_ijk] > 1e-10) {
              dlambda_nl += alpha_ijk * gdlambda[index_ijk] / gmat[index_ijk];
            }
          }
        } // end z-loop
      }   // end y-loop
    }     // end x-loop

  } // end loop over nodes in nonlocal support

  dlambda_nl = dlambda_nl / V_alpha;
}

double
NonLocalDruckerPrager::alpha_nl(const Point& x,
                                Point& s,
                                const vector<double>& l_nl)
{
  if (l_nl[0] < 1e-14) {
    return 1.0;
  }
  double pi    = 3.14159265358979323;
  double k     = pow(6.0 * sqrt(pi), (1.0 / 3.0));
  double diffx = x.x() - s.x();
  double diffy = x.y() - s.y();
  double diffz = x.z() - s.z();
  double arg   = k * (diffx / l_nl[0] + diffy / l_nl[1] + diffz / l_nl[2]);
  return exp(-pow(arg, 2));
}

void
NonLocalDruckerPrager::computeInvariants(Matrix3& stress,
                                         Matrix3& S,
                                         double& I1,
                                         double& J2)
{
  Matrix3 Identity;
  Identity.Identity();
  I1 = stress.Trace();
  S  = stress - Identity * (1.0 / 3.0) * I1;
  J2 = 0.5 * S.Contract(S);
}

void
NonLocalDruckerPrager::computeInvariants(const Matrix3& stress,
                                         Matrix3& S,
                                         double& I1,
                                         double& J2)
{
  Matrix3 Identity;
  Identity.Identity();
  I1 = stress.Trace();
  S  = stress - Identity * (1.0 / 3.0) * I1;
  J2 = 0.5 * S.Contract(S);
}

double
NonLocalDruckerPrager::YieldFunction(Matrix3& stress,
                                     const double& alpha,
                                     double& k_o,
                                     double& eta,
                                     double& eta_nl,
                                     [[maybe_unused]] const int& hardening_type)
{

  Matrix3 S;
  double I1, J2;
  computeInvariants(stress, S, I1, J2);
  return sqrt(J2) + alpha * I1 - k_o - eta - eta_nl;
}

double
NonLocalDruckerPrager::YieldFunction(const Matrix3& stress,
                                     const double& alpha,
                                     double& k_o,
                                     const double& eta,
                                     const double& eta_nl,
                                     [[maybe_unused]] const int& hardening_type)
{

  Matrix3 S;
  double I1, J2;
  computeInvariants(stress, S, I1, J2);
  return sqrt(J2) + alpha * I1 - k_o - eta - eta_nl;
}

double
NonLocalDruckerPrager::YieldFunction(Matrix3& stress,
                                     const double& alpha,
                                     double& k_o,
                                     const double& eta,
                                     const double& eta_nl,
                                     [[maybe_unused]] const int& hardening_type)
{

  Matrix3 S;
  double I1, J2;
  computeInvariants(stress, S, I1, J2);
  return sqrt(J2) + alpha * I1 - k_o - eta - eta_nl;
}

void
NonLocalDruckerPrager::carryForward(const PatchSubset* patches,
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
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

void
NonLocalDruckerPrager::addParticleState(std::vector<const VarLabel*>& from,
                                        std::vector<const VarLabel*>& to)
{

  from.push_back(etaLabel);
  from.push_back(eta_nlLabel);
  from.push_back(pPlasticStrainLabel);
  from.push_back(k_o_distLabel);
  to.push_back(etaLabel_preReloc);
  to.push_back(eta_nlLabel_preReloc);
  to.push_back(pPlasticStrainLabel_preReloc);
  to.push_back(k_o_distLabel_preReloc);
}

void
NonLocalDruckerPrager::addInitialComputesAndRequires(Task* task,
                                                     const MPMMaterial* matl,
                                                     const PatchSet*) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  // Other constitutive model and input dependent computes and requires
  task->computes(etaLabel, matlset);
  task->computes(eta_nlLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(k_o_distLabel, matlset);
}

void
NonLocalDruckerPrager::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl,
                                              const PatchSet* patches) const
{

  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);
  task->needs(Task::OldDW, etaLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, eta_nlLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pPlasticStrainLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, k_o_distLabel, matlset, Ghost::None);
  task->computes(etaLabel_preReloc, matlset);
  task->computes(eta_nlLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(k_o_distLabel_preReloc, matlset);
}

void
NonLocalDruckerPrager::addComputesAndRequires(Task*,
                                              const MPMMaterial*,
                                              const PatchSet*,
                                              const bool,
                                              const bool) const
{
}

double
NonLocalDruckerPrager::computeRhoMicroCM(double pressure,
                                         const double p_ref,
                                         const MPMMaterial* matl,
                                         [[maybe_unused]] double temperature,
                                         [[maybe_unused]] double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge  = pressure - p_ref;
  double rho_cur;
  double bulk = d_initialData.bulk_modulus;

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;

#if 1
  std::cout
    << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR NonLocalDruckerPrager"
    << std::endl;
#endif
}

void
NonLocalDruckerPrager::computePressEOSCM(double rho_cur,
                                         double& pressure,
                                         double p_ref,
                                         double& dp_drho,
                                         double& tmp,
                                         const MPMMaterial* matl,
                                         [[maybe_unused]] double temperature)
{

  double bulk     = d_initialData.bulk_modulus;
  double shear    = d_initialData.shear_modulus;
  double rho_orig = matl->getInitialDensity();

  double p_g = .5 * bulk * (rho_cur / rho_orig - rho_orig / rho_cur);
  pressure   = p_ref + p_g;
  dp_drho    = .5 * bulk * (rho_orig / (rho_cur * rho_cur) + 1. / rho_orig);
  tmp        = (bulk + 4. * shear / 3.) / rho_cur; // speed of sound squared

  std::cout
    << "NO VERSION OF computePressEOSCM EXISTS YET FOR NonLocalDruckerPrager"
    << std::endl;
}

double
NonLocalDruckerPrager::getCompressibility()
{
  std::cout
    << "NO VERSION OF computePressEOSCM EXISTS YET FOR NonLocalDruckerPrager"
    << std::endl;
  return 1.0;
}

void
NonLocalDruckerPrager::initializeLocalMPMLabels()
{

  etaLabel =
    VarLabel::create("eta", ParticleVariable<double>::getTypeDescription());
  eta_nlLabel =
    VarLabel::create("eta_nl", ParticleVariable<double>::getTypeDescription());
  etaLabel_preReloc =
    VarLabel::create("eta+", ParticleVariable<double>::getTypeDescription());
  eta_nlLabel_preReloc =
    VarLabel::create("eta_nl+", ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel =
    VarLabel::create("p.plasticStrain",
                     ParticleVariable<double>::getTypeDescription());
  pPlasticStrainLabel_preReloc =
    VarLabel::create("p.plasticStrain+",
                     ParticleVariable<double>::getTypeDescription());
  k_o_distLabel =
    VarLabel::create("k_o_dist",
                     ParticleVariable<double>::getTypeDescription());
  k_o_distLabel_preReloc =
    VarLabel::create("k_o_dist+",
                     ParticleVariable<double>::getTypeDescription());
}

// Weibull input parser that accepts a structure of input
// parameters defined as:
//
// bool Perturb        'True' for perturbed parameter
// double WeibMed       Medain distrib. value OR const value
//                         depending on bool Perturb
// double WeibMod       Weibull modulus
// double WeibRefVol    Reference Volume
// int    WeibSeed      Seed for random number generator
// std::string WeibDist  String for Distribution
//
// the string 'WeibDist' accepts strings of the following form
// when a perturbed value is desired:
//
// --Distribution--|-Median-|-Modulus-|-Reference Vol -|- Seed -|
// "    weibull,      45e6,      4,        0.0001,          0"
//
// or simply a number if no perturbed value is desired.

void
NonLocalDruckerPrager::WeibullParser(WeibParameters& iP)
{

  // Remove all unneeded characters
  // only remaining are alphanumeric '.' and ','
  for (int i = iP.WeibDist.length() - 1; i >= 0; i--) {
    iP.WeibDist[i] = tolower(iP.WeibDist[i]);
    if (!isalnum(iP.WeibDist[i]) && iP.WeibDist[i] != '.' &&
        iP.WeibDist[i] != ',' && iP.WeibDist[i] != '-' &&
        iP.WeibDist[i] != EOF) {
      iP.WeibDist.erase(i, 1);
    }
  } // End for

  if (iP.WeibDist.substr(0, 4) == "weib") {
    iP.Perturb = true;
  } else {
    iP.Perturb = false;
  }

  // ######
  // If perturbation is NOT desired
  // ######
  if (!iP.Perturb) {
    bool escape        = false;
    int num_of_e       = 0;
    int num_of_periods = 0;
    for (char i : iP.WeibDist) {
      if (i != '.' && i != 'e' && i != '-' && !isdigit(i)) {
        escape = true;
      }

      if (i == 'e') {
        num_of_e += 1;
      }

      if (i == '.') {
        num_of_periods += 1;
      }

      if (num_of_e > 1 || num_of_periods > 1 || escape) {
        std::cerr << "\n\nERROR:\nInput value cannot be parsed. Please\n"
                     "check your input values.\n"
                  << std::endl;
        exit(1);
      }
    } // end for(int i = 0;....)

    if (escape) {
      exit(1);
    }

    iP.WeibMed = atof(iP.WeibDist.c_str());
  }

  // ######
  // If perturbation IS desired
  // ######
  if (iP.Perturb) {
    int weibValues[4];
    int weibValuesCounter = 0;

    for (unsigned int r = 0; r < iP.WeibDist.length(); r++) {
      if (iP.WeibDist[r] == ',') {
        weibValues[weibValuesCounter] = r;
        weibValuesCounter += 1;
      } // end if(iP.WeibDist[r] == ',')
    }   // end for(int r = 0; ...... )

    if (weibValuesCounter != 4) {
      std::cerr << "\n\nERROR:\nWeibull perturbed input string must contain\n"
                   "exactly 4 commas. Verify that your input string is\n"
                   "of the form 'weibull, 45e6, 4, 0.001, 1'.\n"
                << std::endl;
      exit(1);
    } // end if(weibValuesCounter != 4)

    std::string weibMedian;
    std::string weibModulus;
    std::string weibRefVol;
    std::string weibSeed;

    weibMedian =
      iP.WeibDist.substr(weibValues[0] + 1, weibValues[1] - weibValues[0] - 1);
    weibModulus =
      iP.WeibDist.substr(weibValues[1] + 1, weibValues[2] - weibValues[1] - 1);
    weibRefVol =
      iP.WeibDist.substr(weibValues[2] + 1, weibValues[3] - weibValues[2] - 1);
    weibSeed = iP.WeibDist.substr(weibValues[3] + 1);

    iP.WeibMed    = atof(weibMedian.c_str());
    iP.WeibMod    = atof(weibModulus.c_str());
    iP.WeibRefVol = atof(weibRefVol.c_str());
    iP.WeibSeed   = atoi(weibSeed.c_str());

  } // End if (iP.Perturb)
}
