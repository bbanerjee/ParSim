/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/MaterialModels/IsotropicElasticNeoHookeanStateModel.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype

#include <Core/Util/DebugStream.h>

#include <limits>                           // for std::numeric_limits

using namespace Vaango;

using Uintah::ProblemSpecP;
using Uintah::Task;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::DataWarehouse;
using Uintah::ParticleSubset;
using Uintah::Ghost;
using Uintah::VarLabel;
using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::particleIndex;
using Uintah::delt_vartype;
using Uintah::Matrix3;
using Uintah::Vector;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDNeoHookeanDoing:+,PDNeoHookeanDebug:+".....
//  bash     : export SCI_DEBUG="PDNeoHookeanDoing:+,PDNeoHookeanDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDNeoHookeanDoing", false);
static DebugStream dbg("PDNeoHookeanDebug", false);
static DebugStream dbg_extra("PDNeoHookeanDebugExtra", false);

IsotropicElasticNeoHookeanStateModel::IsotropicElasticNeoHookeanStateModel(ProblemSpecP& ps,
                                                                           PeridynamicsFlags* flags)
  : PeridynamicsMaterialModel(flags)
{
  // Get either the Young's modulus and Poisson's ratio, or bulk and shear moduli
  double youngModulus = -1.0, poissonRatio = -1.0;
  if (ps->get("young_modulus", youngModulus)) {
    ps->require("poisson_ratio", poissonRatio);
    d_cm.bulkModulus = youngModulus/(3.0*(1.0-2.0*poissonRatio));
    d_cm.shearModulus = youngModulus/(2.0*(1.0+poissonRatio));
  } else {
    ps->require("bulk_modulus", d_cm.bulkModulus);
    ps->require("shear_modulus", d_cm.shearModulus);
  }
 
}

IsotropicElasticNeoHookeanStateModel::IsotropicElasticNeoHookeanStateModel(const IsotropicElasticNeoHookeanStateModel* cm)
  : PeridynamicsMaterialModel(cm)
{
  d_cm.bulkModulus = cm->d_cm.bulkModulus;
  d_cm.shearModulus = cm->d_cm.shearModulus;
}

// Make a clone of the constitutive model
IsotropicElasticNeoHookeanStateModel* 
IsotropicElasticNeoHookeanStateModel::clone()
{
  return scinew IsotropicElasticNeoHookeanStateModel(*this);
}

IsotropicElasticNeoHookeanStateModel::~IsotropicElasticNeoHookeanStateModel()
{
}

void 
IsotropicElasticNeoHookeanStateModel::outputProblemSpec(ProblemSpecP& ps,
                                                        bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("material_model");
    cm_ps->setAttribute("type", "elastic_neo_hookean_state");
  }
  cm_ps->appendElement("bulk_modulus", d_cm.bulkModulus);
  cm_ps->appendElement("shear_modulus", d_cm.shearModulus);
}

/*! Identify the variabless to be used in the initialization task */
void 
IsotropicElasticNeoHookeanStateModel::addInitialComputesAndRequires(Task* task,
                                                                    const PeridynamicsMaterial* material,
                                                                    const PatchSet* patches) const
{
  // Identify this material
  const Uintah::MaterialSubset* matlset = material->thisMaterial();

  // Add compute flags for the initialization of the stress
  task->computes(d_label->pPK1StressLabel, matlset);
}

/*! Initialize the variables used in the CM */
void 
IsotropicElasticNeoHookeanStateModel::initialize(const Patch* patch,
                                                 const PeridynamicsMaterial* matl,
                                                 DataWarehouse* new_dw)
{
  // Get the set of particles of this material type in the current patch  
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  // Allocate for saving
  ParticleVariable<Matrix3> pStress, pPK1Stress;
  new_dw->allocateAndPut(pStress, d_label->pStressLabel, pset);
  new_dw->allocateAndPut(pPK1Stress, d_label->pPK1StressLabel, pset);

  // Initialize the stress to zero (for now)
  Matrix3 zero(0.0);
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
    pStress[*iter] = zero; 
    pPK1Stress[*iter] = zero; 
  }

  // Compute a stable time step
  computeStableTimestep(patch, matl, new_dw);
}

/* Compute a stable initial timestep */
void
IsotropicElasticNeoHookeanStateModel::computeStableTimestep(const Patch* patch,
                                                            const PeridynamicsMaterial* matl,
                                                            DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matlIndex = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass,     d_label->pMassLabel,     pset);
  new_dw->get(pVolume,   d_label->pVolumeLabel,   pset);
  new_dw->get(pVelocity, d_label->pVelocityLabel, pset);

  double speed_of_sound = 0.0;
  Vector waveSpeed(std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min());

  double kappa = d_cm.bulkModulus;
  double mu = d_cm.shearModulus;
  double pWaveModulus = kappa + mu*(4.0/3.0);

  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

     particleIndex idx = *iter;

     // Compute wave speed at each particle, store the maximum
     Vector vel(0.0, 0.0, 0.0);
     if (pMass[idx] > 0.0) {
       speed_of_sound = std::sqrt(pWaveModulus*pVolume[idx]/pMass[idx]);
       vel[0] = speed_of_sound + std::abs(pVelocity[idx].x());
       vel[1] = speed_of_sound + std::abs(pVelocity[idx].y());
       vel[2] = speed_of_sound + std::abs(pVelocity[idx].z());
     } else {
       speed_of_sound = 0.0;
     }
     waveSpeed = Vector(std::max(vel.x(), waveSpeed.x()),
                        std::max(vel.y(), waveSpeed.y()),
                        std::max(vel.z(), waveSpeed.z()));
  }

  waveSpeed = dx/waveSpeed;
  double delT_new = waveSpeed.minComponent();
  if(delT_new < 1.e-12) {
    new_dw->put(Uintah::delt_vartype(std::numeric_limits<double>::max()), d_label->delTLabel, patch->getLevel());
  } else {
    new_dw->put(Uintah::delt_vartype(delT_new), d_label->delTLabel, patch->getLevel());
  }
}

void 
IsotropicElasticNeoHookeanStateModel::addComputesAndRequires(Task* task, 
                                                             const PeridynamicsMaterial* matl,
                                                             const PatchSet* patches) const
{
  // Constants
  Ghost::GhostType gnone = Ghost::None;
  //Ghost::GhostType gac   = Ghost::AroundCells;

  // Get the current material
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();

  // List the variables needed for this task to execute
  task->requires(Task::OldDW, d_label->delTLabel,              matlset, gnone);
  task->requires(Task::OldDW, d_label->pMassLabel,             matlset, gnone);
  task->requires(Task::OldDW, d_label->pVelocityLabel,         matlset, gnone);
  task->requires(Task::NewDW, d_label->pDefGradLabel_preReloc, matlset, gnone);

  // List the variables computed by this task
  task->computes(d_label->pVolumeLabel_preReloc,    matlset);
  task->computes(d_label->pStressLabel_preReloc,    matlset);
  task->computes(d_label->pPK1StressLabel_preReloc, matlset);
}

void 
IsotropicElasticNeoHookeanStateModel::computeStressTensor(const PatchSubset* patches,
                                                          const PeridynamicsMaterial* matl,
                                                          DataWarehouse* old_dw,
                                                          DataWarehouse* new_dw)
{
  // Set up constants
  Matrix3 One; One.Identity();

  // Get the timestep size
  Uintah::delt_vartype delT;
  old_dw->get(delT, d_label->delTLabel, getLevel(patches));
  
  // Loop through patches
  for (int p = 0; p < patches->size(); p++) {

    // Get the current patch
    const Patch* patch = patches->get(p);

    // Get the material index
    int matlIndex = matl->getDWIndex();
 
    // Get the particle subset for this material
    ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle variables needed
    constParticleVariable<double> pMass;
    old_dw->get(pMass, d_label->pMassLabel, pset);

    constParticleVariable<Vector> pVelocity_old;
    old_dw->get(pVelocity_old, d_label->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pDefGrad_new;
    new_dw->get(pDefGrad_new, d_label->pDefGradLabel_preReloc, pset);

    // Initialize the variables to be updated
    ParticleVariable<double> pVolume_new;
    new_dw->allocateAndPut(pVolume_new, d_label->pVolumeLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStress_new, pPK1Stress_new;
    new_dw->allocateAndPut(pStress_new, d_label->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pPK1Stress_new, d_label->pPK1StressLabel_preReloc, pset);

    // Loop through particles
    double rho_0 = matl->getInitialDensity();
    double kappa = d_cm.bulkModulus;
    double mu = d_cm.shearModulus;
    double pWaveModulus = kappa + mu*(4.0/3.0);
    Vector waveSpeed(std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min());
    dbg_extra << "rho0 = " << rho_0 << " kappa = " << kappa << " mu = " << mu 
              << " pM = " << pWaveModulus << " no. particles = " << pset->numParticles() <<  std::endl;

    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
      
      // Get particle index
      particleIndex idx = *iter;

      // Compute J = det(F)
      Matrix3 FF =  pDefGrad_new[idx];
      Matrix3 FFT =  FF.Transpose();
      double J = FF.Determinant();
      //Matrix3 BB = FF*FFT;
      //std::cout << "F = " << FF << std::endl
      //          << "FT = " << FFT << std::endl
      //          << "J = " << J << std::endl
      //          << "BB = " << BB << std::endl;
      //double a = 1.7e-16;
      //double b = 2.1e-17;
      //double c = a*b;
      //std::cout << "a = " << a << "b = " << b << "c = " << c << std::endl;
 
      dbg_extra << "\t Def Grad = " << pDefGrad_new[idx] << " J = " << J << std::endl;

      // Compute Bbar = J^{-2/3} (F F^T)  and dev(Bbar) = Bbar - 1/3 Tr(Bbar) I
      Matrix3 Bbar = (FF*FFT)*std::pow(J, -2.0/3.0);
      Matrix3 BbarDev = Bbar - One*(Bbar.Trace()/3.0); 

      // Computes stress
      double pressure = -d_cm.bulkModulus*(J-1.0);
      pStress_new[idx] = One*pressure + BbarDev*(d_cm.shearModulus/J);

      // Compute PK1 stress
      pPK1Stress_new[idx] = pStress_new[idx]*(FFT*J);

      if (dbg.active()) {
        dbg << "IsoNeoHookean:" << std::endl
            << "\t DefGrad = " << pDefGrad_new[idx] << " J = " << J << std::endl
            << "\t Bbar = " << Bbar 
            << "\t BbarDev = " << BbarDev << std::endl
            << "\t p = " << pressure << " sig = " << pStress_new[idx]
            << "\t PK1 = " << pPK1Stress_new[idx] << std::endl;
      }

      // Update the particle volume
      double rho_new = rho_0/J;
      pVolume_new[idx] = pMass[idx]/rho_new;
      
      // Compute the wavespeed at each particle and store the maximum
      double speed_of_sound = std::sqrt(pWaveModulus/rho_new);
      Vector vel(0.0, 0.0, 0.0);
      vel[0] = speed_of_sound + std::abs(pVelocity_old[idx].x());
      vel[1] = speed_of_sound + std::abs(pVelocity_old[idx].y());
      vel[2] = speed_of_sound + std::abs(pVelocity_old[idx].z());
      waveSpeed = Vector(std::max(vel.x(), waveSpeed.x()),
                         std::max(vel.y(), waveSpeed.y()),
                         std::max(vel.z(), waveSpeed.z()));

      dbg_extra << " rho_new = " << rho_new << " speed_of sound = " << speed_of_sound
                << " wave speed = " << waveSpeed << std::endl;

    } // end particles loop

    // Find the grid spacing and update deltT
    Vector dx = patch->dCell();
    waveSpeed = dx/waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), d_label->delTLabel, patch->getLevel());

    dbg_extra << " dx = " << dx << " delt = " << waveSpeed << " delT_new = " << delT_new << std::endl;

  } // end patches loop
}

// Register the permanent particle variable states specific to this material
void 
IsotropicElasticNeoHookeanStateModel::addParticleState(std::vector<const VarLabel*>& ,
                                                       std::vector<const VarLabel*>& )
{
}


