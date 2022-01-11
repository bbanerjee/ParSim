/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/IdealGasMP.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
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

IdealGasMP::IdealGasMP(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("gamma", d_initialData.gamma);
  ps->require("specific_heat", d_initialData.cv);
  ps->getWithDefault("reference_pressure", d_initialData.Pref, 101325.0);
  ps->getWithDefault("reference_temperature", d_initialData.Tref, 298.0);
}

IdealGasMP::IdealGasMP(const IdealGasMP* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.gamma = cm->d_initialData.gamma;
  d_initialData.cv = cm->d_initialData.cv;
  d_initialData.Pref = cm->d_initialData.Pref;
  d_initialData.Tref = cm->d_initialData.Tref;
}

void
IdealGasMP::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "ideal_gas");
  }

  cm_ps->appendElement("gamma", d_initialData.gamma);
  cm_ps->appendElement("specific_heat", d_initialData.cv);
  cm_ps->appendElement("reference_pressure", d_initialData.Pref);
  cm_ps->appendElement("reference_temperature", d_initialData.Tref);
}

IdealGasMP*
IdealGasMP::clone()
{
  return scinew IdealGasMP(*this);
}

void
IdealGasMP::addParticleState(std::vector<const VarLabel*>&,
                             std::vector<const VarLabel*>&)
{
}

void
IdealGasMP::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                             DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);
}

void
IdealGasMP::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
                                  DataWarehouse* new_dw)
{
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pMass, pVolume, pTemperature;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pTemperature, lb->pTemperatureLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double gamma = d_initialData.gamma;
  double pRef = d_initialData.Pref;

  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (int idx : *pset) {
    double rho = pMass[idx] / pVolume[idx];
    double dp_dJ = gamma * pRef;
    double c_dil = std::sqrt(dp_dJ / rho) * 10.0;
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
IdealGasMP::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                   const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);
}

void
IdealGasMP::computeStressTensor(const PatchSubset* patches,
                                const MPMMaterial* matl, DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  Matrix3 Identity; Identity.Identity();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  double gamma = d_initialData.gamma;
  double cv = d_initialData.cv;
  double rho_orig = matl->getInitialDensity();

  for (int pp = 0; pp < patches->size(); pp++) {

    const Patch* patch = patches->get(pp);
    Vector dx = patch->dCell();

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pTemperature, pVolume_new;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad_mid, pDefGrad_old, pDefGrad_new;
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pVelocity,    lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    new_dw->get(pVolume_new,  lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad_mid, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    double strainEnergy = 0.;
    for (int idx : *pset) {

      double J_new = pDefGrad_new[idx].Determinant();
      double eps_v = (J_new > 1.0) ? 0.0 : -std::log(J_new);
      double pressure_new = d_initialData.Pref * (std::exp(gamma * eps_v) - 1.0);
      pressure_new = (pressure_new > 0.0) ? pressure_new : 0.0;
      pStress_new[idx] = Identity * (-pressure_new);

      double J_old = pDefGrad_old[idx].Determinant();
      double rho_new = rho_orig / J_new;

      strainEnergy += (-pressure_new * pVolume_new[idx]);
      double dT = (1.0 - J_new/J_old) * (pressure_new / (rho_new * cv));
      pdTdt[idx] = dT / delT;
      //std::cout << "J = " << J_new <<  " p = " << pressure_new 
      //          << " T = " << pTemperature[idx] << " dT/dt = " << pdTdt[idx] << "\n";

      p_q[idx] = 0.;
      if (flag->d_artificialViscosity) {
        Matrix3 D = (pVelGrad_mid[idx] + pVelGrad_mid[idx].Transpose()) * 0.5;
        double DTrace = D.Trace();
        if (DTrace < 0.) {
          double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
          p_q[idx] = 2.5 * 2.5 * dx_ave * dx_ave * rho_new * DTrace * DTrace;
        } else {
          p_q[idx] = 0.;
        }
      }
 
      double dp_dJ = gamma * pressure_new / J_new;
      double c_dil = std::sqrt(dp_dJ / rho_new) * 10.0;
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);
    }

    waveSpeed = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
IdealGasMP::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                                   const bool, const bool) const
{
}

void
IdealGasMP::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                      const PatchSet* patches, MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
}

void
IdealGasMP::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                              ParticleLabelVariableMap* newState,
                              ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
}

// The "CM" versions use the pressure-volume relationship of the CNH model
double
IdealGasMP::computeRhoMicroCM(double press, const double, const MPMMaterial*,
                              double Temp, double rho_guess)
{
  double gamma = d_initialData.gamma;
  double cv = d_initialData.cv;
  return press / ((gamma - 1.0) * cv * Temp);
}

void
IdealGasMP::computePressEOSCM(double rhoM, double& pressure, double,
                              double& dp_drho, double& tmp, const MPMMaterial*,
                              double Temp)
{
  double gamma = d_initialData.gamma;
  double cv = d_initialData.cv;

  pressure = (gamma - 1.0) * rhoM * cv * Temp;
  dp_drho = (gamma - 1.0) * cv * Temp;
  double dp_de = (gamma - 1.0) * rhoM;
  tmp = dp_drho + dp_de * pressure / (rhoM * rhoM); // C^2
}

double
IdealGasMP::getCompressibility()
{
  return 1.0 / d_initialData.Pref;
}

namespace Uintah {
} // End namespace Uintah
