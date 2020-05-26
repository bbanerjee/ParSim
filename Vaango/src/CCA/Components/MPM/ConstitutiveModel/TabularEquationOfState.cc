/*
 * The MIT License
 *
 * Copyright (c) 2017-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularEquationOfState.h>
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
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <sci_values.h>

using namespace Vaango;

using ProblemSpecP = Uintah::ProblemSpecP;
using Task = Uintah::Task;
using Patch = Uintah::Patch;
using PatchSet = Uintah::PatchSet;
using PatchSubset = Uintah::PatchSubset;
using MaterialSubset = Uintah::MaterialSubset;
using ParticleSubset = Uintah::ParticleSubset;
using constParticleVariableBase = Uintah::constParticleVariableBase;
using ParticleVariableBase = Uintah::ParticleVariableBase;
using ParticleLabelVariableMap = Uintah::ParticleLabelVariableMap;
using DataWarehouse = Uintah::DataWarehouse;
using VarLabel = Uintah::VarLabel;
using MPMFlags = Uintah::MPMFlags;
using MPMLabel = Uintah::MPMLabel;
using MPMMaterial = Uintah::MPMMaterial;
using delt_vartype = Uintah::delt_vartype;
using sum_vartype = Uintah::sum_vartype;
using Vector = Uintah::Vector;
using Matrix3 = Uintah::Matrix3;

TabularEquationOfState::TabularEquationOfState(ProblemSpecP& ps,
                                               MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
  , d_table(ps)
{
  d_table.setup();
}

TabularEquationOfState::TabularEquationOfState(const TabularEquationOfState* cm)
  : ConstitutiveModel(cm)
{
  d_table = cm->d_table;
}

TabularEquationOfState::~TabularEquationOfState()
= default;

void
TabularEquationOfState::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "tabular_eos");
  }
  d_table.outputProblemSpec(cm_ps);
}

TabularEquationOfState*
TabularEquationOfState::clone()
{
  return scinew TabularEquationOfState(*this);
}

void
TabularEquationOfState::initializeCMData(const Patch* patch,
                                         const MPMMaterial* matl,
                                         DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  computeStableTimestep(patch, matl, new_dw);
}

///////////////////////////////////////////////////////////////////////////
/*! Allocate data required during the conversion of failed particles
    from one material to another */
///////////////////////////////////////////////////////////////////////////
void
TabularEquationOfState::allocateCMDataAddRequires(Task* task,
                                                  const MPMMaterial* matl,
                                                  const PatchSet* patches,
                                                  MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
}

void
TabularEquationOfState::allocateCMDataAdd(DataWarehouse* new_dw,
                                          ParticleSubset* addset,
                                          ParticleLabelVariableMap* newState,
                                          ParticleSubset* delset,
                                          DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
}

// This is only called for the initial timestep - all other timesteps
// are computed as a side-effect of computeStressTensor
void
TabularEquationOfState::computeStableTimestep(const Patch* patch,
                                              const MPMMaterial* matl,
                                              DataWarehouse* new_dw)
{
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();
  double rho_0 = matl->getInitialDensity();

  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  Uintah::constParticleVariable<double> pMass, pVolume;
  Uintah::constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  // Compute wave speed + particle velocity at each particle,
  // store the maximum
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (const auto& pidx : *pset) {
    double rho = pMass[pidx]/pVolume[pidx];
    std::cout << "Density = " << rho << "\n";
    double bulkModulus = computeBulkModulus(rho_0, rho);
    double c_bulk = std::sqrt(bulkModulus/rho);
    std::cout << "K = " << bulkModulus 
              << " c_p = " << c_bulk << "\n";
    WaveSpeed = Vector(Uintah::Max(c_bulk + std::abs(pVelocity[pidx].x()), WaveSpeed.x()),
                       Uintah::Max(c_bulk + std::abs(pVelocity[pidx].y()), WaveSpeed.y()),
                       Uintah::Max(c_bulk + std::abs(pVelocity[pidx].z()), WaveSpeed.z()));
  }
  WaveSpeed = dx / WaveSpeed;
  double delT = WaveSpeed.minComponent();
  std::cout << "delT = " << delT << "\n";
  new_dw->put(delt_vartype(delT), lb->delTLabel, patch->getLevel());
}

double
TabularEquationOfState::computeBulkModulus(const double& rho_0,
                                           const double& rho) const
{
  double epsilon = 1.0e-6;
  double rho_over_rho0 = rho/rho_0;
  DoubleVec1D pressure_lo = d_table.interpolate<1>({{rho_over_rho0 - epsilon}});
  DoubleVec1D pressure_hi = d_table.interpolate<1>({{rho_over_rho0 + epsilon}});
  return rho*(pressure_hi[0] - pressure_lo[0])/(2*epsilon);
}

void
TabularEquationOfState::computeStressTensor(const PatchSubset* patches,
                                            const MPMMaterial* matl,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw)
{
  Matrix3 Identity;
  Identity.Identity();
  double rho_0 = matl->getInitialDensity();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    Vector dx = patch->dCell();
    int matID = matl->getDWIndex();

    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    Uintah::constParticleVariable<double> pMass;
    Uintah::constParticleVariable<Vector> pVelocity;
    old_dw->get(pMass,     lb->pMassLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);

    Uintah::constParticleVariable<double> pVolume;
    Uintah::constParticleVariable<Matrix3> pVelGrad;
    new_dw->get(pVolume,  lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);

    Uintah::ParticleVariable<Matrix3> pStress;
    Uintah::ParticleVariable<double> pdTdt, p_q;
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt,   lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q,     lb->p_qLabel_preReloc, pset);

    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);
    double se = 0;
    for (auto& pidx : *pset) {

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[pidx] = 0.0;

      double rho = pMass[pidx]/pVolume[pidx];
      DoubleVec1D pressure = d_table.interpolate<1>({{rho/rho_0}});
      pStress[pidx] = Identity * (-pressure[0]);

      // Compute wave speed + particle velocity at each particle,
      // store the maximum
      double bulkModulus = computeBulkModulus(rho_0, rho);
      double c_bulk = std::sqrt(bulkModulus/rho);
      std::cout << "Density = " << rho << " K = " << bulkModulus
                << " c_p = " << c_bulk << "\n";
      WaveSpeed = Vector(Uintah::Max(c_bulk + fabs(pVelocity[pidx].x()), WaveSpeed.x()),
                         Uintah::Max(c_bulk + fabs(pVelocity[pidx].y()), WaveSpeed.y()),
                         Uintah::Max(c_bulk + fabs(pVelocity[pidx].z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dxAvg = (dx.x() + dx.y() + dx.z()) / 3.0;
        Matrix3 D = (pVelGrad[pidx] + pVelGrad[pidx].Transpose()) * 0.5;
        p_q[pidx] = artificialBulkViscosity(D.Trace(), c_bulk, rho, dxAvg);
      } else {
        p_q[pidx] = 0.;
      }

      // Compute the strain energy for all the particles (0.5*p*eps_v)
      double e = 0.5*pressure[0]*std::log(rho_0/rho);
      se += e;
    } // end loop over particles

    WaveSpeed = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }
  }
}

void
TabularEquationOfState::carryForward(const PatchSubset* patches,
                                     const MPMMaterial* matl,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
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
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

void
TabularEquationOfState::addParticleState(std::vector<const VarLabel*>& from,
                                         std::vector<const VarLabel*>& to)
{
}

void
TabularEquationOfState::addComputesAndRequires(Task* task,
                                               const MPMMaterial* matl,
                                               const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);
}

void
TabularEquationOfState::addComputesAndRequires(Task*, const MPMMaterial*,
                                               const PatchSet*, const bool,
                                               const bool) const
{
}

double
TabularEquationOfState::computeRhoMicroCM(double /*pressure*/,
                                          const double /*p_ref*/,
                                          const MPMMaterial* /*matl*/,
                                          double temperature, double rho_guess)
{
#if 0
  double rho_0 = matl->getInitialDensity();
  double bulk = d_initialData.Bulk;

  double p_gauge = pressure - p_ref;
  double rho_cur;

  rho_cur = rho_0*(p_gauge/bulk + sqrt((p_gauge/bulk)*(p_gauge/bulk) +1));
#endif

  std::cout
    << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR TabularEquationOfState"
    << std::endl;

  double rho_cur = 0.;

  return rho_cur;
}

void
TabularEquationOfState::computePressEOSCM(double /*rho_cur*/,
                                          double& /*pressure*/,
                                          double /*p_ref*/, double& /*dp_drho*/,
                                          double& /*tmp*/,
                                          const MPMMaterial* /*matl*/,
                                          double temperature)
{
#if 0
  double bulk = d_initialData.Bulk;
  double shear = d_initialData.Shear;
  double rho_0 = matl->getInitialDensity();

  double p_g = .5*bulk*(rho_cur/rho_0 - rho_0/rho_cur);
  pressure = p_ref + p_g;
  dp_drho  = .5*bulk*(rho_0/(rho_cur*rho_cur) + 1./rho_0);
  tmp = (bulk + 4.*shear/3.)/rho_cur;  // speed of sound squared
#endif

  std::cout
    << "NO VERSION OF computePressEOSCM EXISTS YET FOR TabularEquationOfState"
    << std::endl;
}

double
TabularEquationOfState::getCompressibility()
{
  std::cout
    << "NO VERSION OF getCompressibility EXISTS YET FOR TabularEquationOfState"
    << std::endl;
  return 1.0;
}
