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

#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/MWViscoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h> // just added
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Short27.h> //for Fracture
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using std::cerr;
using namespace Uintah;

MWViscoElastic::MWViscoElastic(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("e_shear_modulus", d_initialData.E_Shear);
  ps->require("e_bulk_modulus", d_initialData.E_Bulk);
  ps->require("ve_shear_modulus", d_initialData.VE_Shear);
  ps->require("ve_bulk_modulus", d_initialData.VE_Bulk);
  ps->require("ve_volumetric_viscosity", d_initialData.V_Viscosity);
  ps->require("ve_deviatoric_viscosity", d_initialData.D_Viscosity);

  pStress_eLabel = VarLabel::create(
    "p.stress_e", ParticleVariable<Matrix3>::getTypeDescription());
  pStress_ve_vLabel = VarLabel::create(
    "p.stress_ve_v", ParticleVariable<double>::getTypeDescription());
  pStress_ve_dLabel = VarLabel::create(
    "p.stress_ve_d", ParticleVariable<Matrix3>::getTypeDescription());
  pStress_e_vLabel = VarLabel::create(
    "p.stress_e_v", ParticleVariable<double>::getTypeDescription());
  pStress_e_dLabel = VarLabel::create(
    "p.stress_e_d", ParticleVariable<Matrix3>::getTypeDescription());

  pStress_eLabel_preReloc = VarLabel::create(
    "p.stress_e+", ParticleVariable<Matrix3>::getTypeDescription());
  pStress_ve_vLabel_preReloc = VarLabel::create(
    "p.stress_ve_v+", ParticleVariable<double>::getTypeDescription());
  pStress_ve_dLabel_preReloc = VarLabel::create(
    "p.stress_ve_d+", ParticleVariable<Matrix3>::getTypeDescription());
  pStress_e_vLabel_preReloc = VarLabel::create(
    "p.stress_e_v+", ParticleVariable<double>::getTypeDescription());
  pStress_e_dLabel_preReloc = VarLabel::create(
    "p.stress_e_d+", ParticleVariable<Matrix3>::getTypeDescription());
}

MWViscoElastic::MWViscoElastic(const MWViscoElastic* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.E_Shear = cm->d_initialData.E_Shear;
  d_initialData.E_Bulk = cm->d_initialData.E_Bulk;
  d_initialData.VE_Shear = cm->d_initialData.VE_Shear;
  d_initialData.VE_Bulk = cm->d_initialData.VE_Bulk;
  d_initialData.V_Viscosity = cm->d_initialData.V_Viscosity;
  d_initialData.D_Viscosity = cm->d_initialData.D_Viscosity;
}

MWViscoElastic::~MWViscoElastic()
{
  VarLabel::destroy(pStress_eLabel);
  VarLabel::destroy(pStress_ve_vLabel);
  VarLabel::destroy(pStress_ve_dLabel);
  VarLabel::destroy(pStress_e_vLabel);
  VarLabel::destroy(pStress_e_dLabel);
  VarLabel::destroy(pStress_eLabel_preReloc);
  VarLabel::destroy(pStress_ve_vLabel_preReloc);
  VarLabel::destroy(pStress_ve_dLabel_preReloc);
  VarLabel::destroy(pStress_e_vLabel_preReloc);
  VarLabel::destroy(pStress_e_dLabel_preReloc);
}

void
MWViscoElastic::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "mw_visco_elastic");
  }

  cm_ps->appendElement("e_shear_modulus", d_initialData.E_Shear);
  cm_ps->appendElement("e_bulk_modulus", d_initialData.E_Bulk);
  cm_ps->appendElement("ve_shear_modulus", d_initialData.VE_Shear);
  cm_ps->appendElement("ve_bulk_modulus", d_initialData.VE_Bulk);
  cm_ps->appendElement("ve_volumetric_viscosity", d_initialData.V_Viscosity);
  cm_ps->appendElement("ve_deviatoric_viscosity", d_initialData.D_Viscosity);
}

std::unique_ptr<ConstitutiveModel>
MWViscoElastic::clone()
{
  return std::make_unique<MWViscoElastic>(*this);
}

void
MWViscoElastic::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                                 DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  // Put stuff in here to initialize each particle's
  // constitutive model parameters and deformationMeasure
  Matrix3 Identity, zero(0.);
  Identity.Identity();

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<Matrix3> pstress_e, pstress_ve_d, pstress_e_d;
  ParticleVariable<double> pstress_ve_v, pstress_e_v;
  new_dw->allocateAndPut(pstress_e, pStress_eLabel, pset);
  new_dw->allocateAndPut(pstress_ve_v, pStress_ve_vLabel, pset);
  new_dw->allocateAndPut(pstress_ve_d, pStress_ve_dLabel, pset);
  new_dw->allocateAndPut(pstress_e_v, pStress_e_vLabel, pset);
  new_dw->allocateAndPut(pstress_e_d, pStress_e_dLabel, pset);

  for (int& iter : *pset) {
    pstress_e[iter] = zero;
    pstress_ve_v[iter] = 0.0;
    pstress_ve_d[iter] = zero;
    pstress_e_v[iter] = 0.0;
    pstress_e_d[iter] = zero;
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
MWViscoElastic::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                          const PatchSet* patches,
                                          [[maybe_unused]] MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);

  // Add requires local to this model
  Ghost::GhostType gnone = Ghost::None;
  task->needs(Task::NewDW, pStress_ve_vLabel, matlset, gnone);
  task->needs(Task::NewDW, pStress_ve_dLabel, matlset, gnone);
  task->needs(Task::NewDW, pStress_e_vLabel, matlset, gnone);
  task->needs(Task::NewDW, pStress_e_dLabel, matlset, gnone);
  task->needs(Task::NewDW, pStress_eLabel, matlset, gnone);
}

void
MWViscoElastic::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                                  ParticleLabelVariableMap* newState,
                                  ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  ParticleVariable<Matrix3> pstress_e, pstress_ve_d, pstress_e_d;
  ParticleVariable<double> pstress_ve_v, pstress_e_v;

  constParticleVariable<Matrix3> o_stress_e, o_stress_ve_d, o_stress_e_d;
  constParticleVariable<double> o_stress_ve_v, o_stress_e_v;

  new_dw->allocateTemporary(pstress_e, addset);
  new_dw->allocateTemporary(pstress_ve_v, addset);
  new_dw->allocateTemporary(pstress_ve_d, addset);
  new_dw->allocateTemporary(pstress_e_v, addset);
  new_dw->allocateTemporary(pstress_e_d, addset);

  new_dw->get(o_stress_e, pStress_eLabel, delset);
  new_dw->get(o_stress_ve_v, pStress_ve_vLabel, delset);
  new_dw->get(o_stress_ve_d, pStress_ve_dLabel, delset);
  new_dw->get(o_stress_e_v, pStress_e_vLabel, delset);
  new_dw->get(o_stress_e_d, pStress_e_dLabel, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pstress_e[*n] = o_stress_e[*o];
    pstress_ve_v[*n] = o_stress_ve_v[*o];
    pstress_ve_d[*n] = o_stress_ve_d[*o];
    pstress_e_v[*n] = o_stress_e_v[*o];
    pstress_e_d[*n] = o_stress_e_d[*o];
  }

  (*newState)[pStress_eLabel] = pstress_e.clone();
  (*newState)[pStress_ve_vLabel] = pstress_ve_v.clone();
  (*newState)[pStress_ve_dLabel] = pstress_ve_d.clone();
  (*newState)[pStress_e_vLabel] = pstress_e_v.clone();
  (*newState)[pStress_e_dLabel] = pstress_e_d.clone();
}

void
MWViscoElastic::addParticleState(std::vector<const VarLabel*>& from,
                                 std::vector<const VarLabel*>& to)
{

  // Add the local particle state data for this constitutive model.
  from.push_back(pStress_eLabel);
  from.push_back(pStress_ve_vLabel);
  from.push_back(pStress_ve_dLabel);
  from.push_back(pStress_e_vLabel);
  from.push_back(pStress_e_dLabel);

  to.push_back(pStress_eLabel_preReloc);
  to.push_back(pStress_ve_vLabel_preReloc);
  to.push_back(pStress_ve_dLabel_preReloc);
  to.push_back(pStress_e_vLabel_preReloc);
  to.push_back(pStress_e_dLabel_preReloc);
}

void
MWViscoElastic::computeStableTimestep(const Patch* patch,
                                      const MPMMaterial* matl,
                                      DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double e_shear = d_initialData.E_Shear;
  double e_bulk = d_initialData.E_Bulk;
  double ve_shear = d_initialData.VE_Shear;
  double ve_bulk = d_initialData.VE_Bulk;
  double bulk = e_bulk + ve_bulk;
  double shear = e_shear + ve_shear;

  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    c_dil = sqrt((bulk + 4. * shear / 3.) * pVolume[idx] / pMass[idx]);
    WaveSpeed = Vector(Max(c_dil + fabs(pVelocity[idx].x()), WaveSpeed.x()),
                       Max(c_dil + fabs(pVelocity[idx].y()), WaveSpeed.y()),
                       Max(c_dil + fabs(pVelocity[idx].z()), WaveSpeed.z()));
  }
  WaveSpeed = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
MWViscoElastic::computeStressTensor(const PatchSubset* patches,
                                    const MPMMaterial* matl,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    double se = 0.0, ve = 0.0;
    const Patch* patch = patches->get(p);
    Matrix3 pDefGradInc, Identity, zero(0.), One(1.);
    double c_dil = 0.0;
    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);
    double onethird = (1.0 / 3.0);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    Identity.Identity();

    Vector dx = patch->dCell();
    // double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    int dwi = matl->getDWIndex();
    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<Point> px;
    constParticleVariable<Matrix3> pDefGrad, pstress_e;
    constParticleVariable<Matrix3> pstress_ve_d, pstress_e_d;
    ParticleVariable<Matrix3> pstress_e_new, pstress_ve_new;
    ParticleVariable<Matrix3> pstress_ve_d_new, pstress_e_d_new;
    ParticleVariable<Matrix3> pstress_new;
    constParticleVariable<double> pMass, ptemperature;
    constParticleVariable<double> pstress_ve_v, pstress_e_v;
    ParticleVariable<double> pstress_ve_v_new, pstress_e_v_new;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pSize;
    ParticleVariable<double> pdTdt, p_q;

    constParticleVariable<Matrix3> pDefGrad_new, velGrad;
    constParticleVariable<double> pVolume_new;
    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(velGrad, lb->pVelGradLabel_preReloc, pset);

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    // Ghost::GhostType  gac   = Ghost::AroundCells;

    new_dw->allocateTemporary(pstress_ve_new, pset);
    new_dw->allocateAndPut(pstress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pstress_e_new, pStress_eLabel_preReloc, pset);
    new_dw->allocateAndPut(pstress_ve_v_new, pStress_ve_vLabel_preReloc, pset);
    new_dw->allocateAndPut(pstress_ve_d_new, pStress_ve_dLabel_preReloc, pset);
    new_dw->allocateAndPut(pstress_e_v_new, pStress_e_vLabel_preReloc, pset);
    new_dw->allocateAndPut(pstress_e_d_new, pStress_e_dLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);

    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pstress_e, pStress_eLabel, pset);
    old_dw->get(pstress_ve_v, pStress_ve_vLabel, pset);
    old_dw->get(pstress_ve_d, pStress_ve_dLabel, pset);
    old_dw->get(pstress_e_v, pStress_e_vLabel, pset);
    old_dw->get(pstress_e_d, pStress_e_dLabel, pset);
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(ptemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    double e_shear = d_initialData.E_Shear;
    double e_bulk = d_initialData.E_Bulk;
    double ve_shear = d_initialData.VE_Shear;
    double ve_bulk = d_initialData.VE_Bulk;
    double v_viscosity = d_initialData.V_Viscosity;
    double d_viscosity = d_initialData.D_Viscosity;
    double bulk = e_bulk + ve_bulk;
    double shear = e_shear + ve_shear;
    double rho_orig = matl->getInitialDensity();

    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Calculate rate of deformation D, and deviatoric rate DPrime

      Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose()) * .5;
      Matrix3 DPrime = D - Identity * onethird * D.Trace();

      // standard solid element:
      pstress_e_new[idx] =
        pstress_e[idx] +
        (DPrime * 2. * e_shear + Identity * e_bulk * D.Trace()) * delT;

      pstress_ve_v_new[idx] =
        (-onethird * D.Trace() * v_viscosity +
         pstress_ve_v[idx] * (v_viscosity / 3 / ve_bulk / delT)) /
        (1 + v_viscosity / 3 / ve_bulk / delT);

      pstress_ve_d_new[idx] =
        (DPrime * d_viscosity +
         pstress_ve_d[idx] * (d_viscosity / 2 / ve_shear / delT)) /
        (1 + d_viscosity / 2 / ve_shear / delT);

      pstress_ve_new[idx] =
        pstress_ve_d_new[idx] - Identity * pstress_ve_v_new[idx];

      pstress_new[idx] = pstress_e_new[idx] + pstress_ve_new[idx];

      // get the volumetric part of the deformation
      double J = pDefGrad[idx].Determinant();

      // Compute the strain energy for all the particles
      pstress_e_v_new[idx] = pstress_e_v[idx] - D.Trace() * e_bulk * delT;

      pstress_e_d_new[idx] = pstress_e_d[idx] + DPrime * 2 * e_shear * delT;

      double p = pstress_e_v_new[idx] + pstress_ve_v_new[idx];

      double ee = (pstress_ve_d_new[idx].NormSquared() / 4 / ve_shear +
                   pstress_e_d_new[idx].NormSquared() / 4 / e_shear +
                   p * p / 2 / (ve_bulk + e_bulk)) *
                  pVolume_new[idx];

      ve = ve +
           (pstress_ve_d_new[idx].NormSquared() / d_viscosity) *
             pVolume_new[idx] * delT;

      double e = ee + ve;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pVelocity_idx = pVelocity[idx];
      c_dil = sqrt((bulk + 4. * shear / 3.) * pVolume_new[idx] / pMass[idx]);
      WaveSpeed = Vector(Max(c_dil + fabs(pVelocity_idx.x()), WaveSpeed.x()),
                         Max(c_dil + fabs(pVelocity_idx.y()), WaveSpeed.y()),
                         Max(c_dil + fabs(pVelocity_idx.z()), WaveSpeed.z()));
      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double rho_cur = rho_orig / J;
        double c_bulk = sqrt(bulk / rho_cur);
        Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose()) * 0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    WaveSpeed = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }

    //delete interpolator;
  }
}

void
MWViscoElastic::addInitialComputesAndRequires(Task* task,
                                              const MPMMaterial* matl,
                                              const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pStress_ve_vLabel, matlset);
  task->computes(pStress_ve_dLabel, matlset);
  task->computes(pStress_e_vLabel, matlset);
  task->computes(pStress_e_dLabel, matlset);
  task->computes(pStress_eLabel, matlset);
}

void
MWViscoElastic::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                       const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);

  // Add requires local to this model
  Ghost::GhostType gnone = Ghost::None;
  task->computes(pStress_ve_vLabel_preReloc, matlset);
  task->computes(pStress_ve_dLabel_preReloc, matlset);
  task->computes(pStress_e_vLabel_preReloc, matlset);
  task->computes(pStress_e_dLabel_preReloc, matlset);
  task->computes(pStress_eLabel_preReloc, matlset);

  task->needs(Task::OldDW, pStress_ve_vLabel, matlset, gnone);
  task->needs(Task::OldDW, pStress_ve_dLabel, matlset, gnone);
  task->needs(Task::OldDW, pStress_e_vLabel, matlset, gnone);
  task->needs(Task::OldDW, pStress_e_dLabel, matlset, gnone);
  task->needs(Task::OldDW, pStress_eLabel, matlset, gnone);
}

void
MWViscoElastic::addComputesAndRequires(Task*, const MPMMaterial*,
                                       const PatchSet*, const bool,
                                       const bool) const
{
}

double
MWViscoElastic::computeRhoMicroCM(double pressure, const double p_ref,
                                  const MPMMaterial* matl, [[maybe_unused]] double temperature,
                                  [[maybe_unused]] double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  // double p_ref=101325.0;
  double p_gauge = pressure - p_ref;
  double rho_cur;
  double e_bulk = d_initialData.E_Bulk;
  double ve_bulk = d_initialData.VE_Bulk;
  double bulk = e_bulk + ve_bulk;

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;

#if 0
  std::cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR MWViscoElastic"
       << std::endl;
#endif
}

void
MWViscoElastic::computePressEOSCM(const double rho_cur, double& pressure,
                                  const double p_ref, double& dp_drho,
                                  double& tmp, const MPMMaterial* matl,
                                  [[maybe_unused]] double temperature)
{

  double e_shear = d_initialData.E_Shear;
  double e_bulk = d_initialData.E_Bulk;
  double ve_shear = d_initialData.VE_Shear;
  double ve_bulk = d_initialData.VE_Bulk;
  double bulk = e_bulk + ve_bulk;
  double shear = e_shear + ve_shear;
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = bulk * rho_orig / (rho_cur * rho_cur);
  tmp = sqrt((bulk + 4. * shear / 3.) / rho_cur); // speed of sound squared

#if 0
  std::cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR MWViscoElastic"
       << std::endl;
#endif
}

double
MWViscoElastic::getCompressibility()
{
  return 1.0 / (d_initialData.E_Bulk + d_initialData.VE_Bulk);
}

namespace Uintah {

#if 0
static MPI_Datatype makeMPI_CMData()
{
   ASSERTEQ(sizeof(MWViscoElastic::StateData), sizeof(double)*2);
   MPI_Datatype mpitype;
   MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &mpitype);
   MPI_Type_commit(&mpitype);
   return mpitype;
}

const TypeDescription* fun_getTypeDescription(MWViscoElastic::StateData*)
{
   static TypeDescription* td = 0;
   if(!td){
      td = scinew
        TypeDescription(TypeDescription::Type::Other,
                        "MWViscoElastic::StateData", true, &makeMPI_CMData);
   }
   return td;
}
#endif

} // End namespace Uintah
