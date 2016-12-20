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

#include <CCA/Components/Peridynamics/MaterialModels/PolarOrthotropicLinearElasticStateModel.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Math/Matrix3Rotation.h>

#include <iostream>
#include <limits>                           // for std::numeric_limits

using namespace Vaango;

using Uintah::Point;
using Uintah::Vector;
using Uintah::Matrix3;
using Uintah::SymmMatrix6;
using Uintah::delt_vartype;
using Uintah::ProblemSetupException;
using Uintah::ProblemSpecP;
using Uintah::ParticleSubset;
using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::Ghost;
using Uintah::MaterialSubset;
using Uintah::Task;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::DataWarehouse;

PolarOrthotropicLinearElasticStateModel::PolarOrthotropicLinearElasticStateModel(ProblemSpecP& ps,
                                                                                 PeridynamicsFlags* flags)
  : PeridynamicsMaterialModel(flags)
{
  // Get the axis of symmetry of the cylinder 
  ps->require("symmetry_axis_top", d_cm.top);
  ps->require("symmetry_axis_bottom", d_cm.bottom);
  if ((d_cm.top-d_cm.bottom).length2() == 0.0) {
    std::ostringstream out;
    out << "The axis of symmetry is a point. ";
    out << "Please check the values for the top and bottom points of the axis in the input file";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Get the elastic moduli with respect to the axis of symmetry 
  // 1 <-> r, 2 <-> theta, 3 <-> z
  ps->require("E_r", d_cm.Er);
  ps->require("E_theta", d_cm.Etheta);
  ps->require("E_z", d_cm.Ez);
  ps->require("nu_theta_r", d_cm.nuthetar);
  ps->require("nu_z_r", d_cm.nuzr);
  ps->require("nu_z_theta", d_cm.nuztheta);
  ps->require("G_theta_z", d_cm.Gthetaz);
  ps->require("G_z_r", d_cm.Gzr);
  ps->require("G_r_theta", d_cm.Grtheta);

  // Compute the compliance matrix
  SymmMatrix6 complianceMatrix;
  complianceMatrix(0,0) = 1.0/d_cm.Er;
  complianceMatrix(0,1) = -d_cm.nuthetar/d_cm.Etheta;
  complianceMatrix(0,2) = -d_cm.nuzr/d_cm.Ez;
  complianceMatrix(1,1) = 1.0/d_cm.Etheta;
  complianceMatrix(1,2) = -d_cm.nuztheta/d_cm.Ez;
  complianceMatrix(2,2) = 1.0/d_cm.Ez;
  complianceMatrix(3,3) = 1.0/d_cm.Gthetaz;
  complianceMatrix(4,4) = 1.0/d_cm.Gzr;
  complianceMatrix(5,5) = 1.0/d_cm.Grtheta;

  complianceMatrix(1,0) = complianceMatrix(0,1);
  complianceMatrix(2,0) = complianceMatrix(0,2);
  complianceMatrix(2,1) = complianceMatrix(1,2);
  
  // Check that everything is consistent
  if ((complianceMatrix(0,0) < 0.0) || (complianceMatrix(1,1) < 0.0) || (complianceMatrix(2,2) < 0.0) ||
      (complianceMatrix(3,3) < 0.0) || (complianceMatrix(4,4) < 0.0) || (complianceMatrix(5,5) < 0.0)) {
    std::ostringstream out;
    out << "The compliance matrix has negative diagonal components";
    out << "Please check the values in the input file to make sure the Young's and shear moduli are positive.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  double delta2 = complianceMatrix(0,0)*complianceMatrix(1,1) - complianceMatrix(0,1)*complianceMatrix(0,1);
  if (delta2 < 0.0) {
    std::ostringstream out;
    out << "Compliance matrix submatrix has negative determinant: S11 S22 - S12^2 < 0.";
    out << "Please check the values in the input file to make sure the input data are correct.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  double delta3 = delta2*complianceMatrix(2,2) - complianceMatrix(0,0)*complianceMatrix(1,2)*complianceMatrix(1,2) +
    2.0*complianceMatrix(0,1)*complianceMatrix(1,2)*complianceMatrix(0,2) - 
    complianceMatrix(1,1)*complianceMatrix(0,2)*complianceMatrix(0,2);
  if (delta3 < 0.0) {
    std::ostringstream out;
    out << "Compliance matrix submatrix has negative determinant: ";
    out << "  (S11 S22 - S12^2)S33 - S11 S23^2 + 2 S12 S23 S13 - S22 S13^2 < 0.";
    out << "Please check the values in the input file to make sure the input data are correct.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  
  // Compute stiffness matrix
  complianceMatrix.inverse(d_cm.stiffnessMatrix);
}

PolarOrthotropicLinearElasticStateModel::PolarOrthotropicLinearElasticStateModel(const PolarOrthotropicLinearElasticStateModel* cm)
  : PeridynamicsMaterialModel(cm)
{
  d_cm.top = cm->d_cm.top;
  d_cm.bottom = cm->d_cm.bottom;
  d_cm.Er = cm->d_cm.Er;
  d_cm.Etheta = cm->d_cm.Etheta;
  d_cm.Ez = cm->d_cm.Ez;
  d_cm.nuthetar = cm->d_cm.nuthetar;
  d_cm.nuzr = cm->d_cm.nuzr;
  d_cm.nuztheta = cm->d_cm.nuztheta;
  d_cm.Gthetaz = cm->d_cm.Gthetaz;
  d_cm.Gzr = cm->d_cm.Gzr;
  d_cm.Grtheta = cm->d_cm.Grtheta;
  d_cm.stiffnessMatrix = cm->d_cm.stiffnessMatrix;
}

// Make a clone of the constitutive model
PolarOrthotropicLinearElasticStateModel* 
PolarOrthotropicLinearElasticStateModel::clone()
{
  return scinew PolarOrthotropicLinearElasticStateModel(*this);
}

PolarOrthotropicLinearElasticStateModel::~PolarOrthotropicLinearElasticStateModel()
{
}

void 
PolarOrthotropicLinearElasticStateModel::outputProblemSpec(ProblemSpecP& ps,
                                                           bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("material_model");
    cm_ps->setAttribute("type", "polar_orthotropic_linear_elastic_state");
  }
  cm_ps->appendElement("symmetry_axis_top", d_cm.top);
  cm_ps->appendElement("symmetry_axis_bottom", d_cm.bottom);

  // Get the elastic moduli with respect to the axis of symmetry 
  // 1 <-> r, 2 <-> theta, 3 <-> z
  cm_ps->appendElement("E_r", d_cm.Er);
  cm_ps->appendElement("E_theta", d_cm.Etheta);
  cm_ps->appendElement("E_z", d_cm.Ez);
  cm_ps->appendElement("nu_theta_r", d_cm.nuthetar);
  cm_ps->appendElement("nu_z_r", d_cm.nuzr);
  cm_ps->appendElement("nu_z_theta", d_cm.nuztheta);
  cm_ps->appendElement("G_theta_z", d_cm.Gthetaz);
  cm_ps->appendElement("G_z_r", d_cm.Gzr);
  cm_ps->appendElement("G_r_heta", d_cm.Grtheta);
}

/*! Identify the variabless to be used in the initialization task */
void 
PolarOrthotropicLinearElasticStateModel::addInitialComputesAndRequires(Task* task,
                                                                    const PeridynamicsMaterial* material,
                                                                    const PatchSet* patches) const
{
  // Identify this material
  const MaterialSubset* matlset = material->thisMaterial();

  // Add compute flags for the initialization of the stress
  task->computes(d_label->pPK1StressLabel, matlset);
}

/*! Initialize the variables used in the CM */
void 
PolarOrthotropicLinearElasticStateModel::initialize(const Patch* patch,
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

  // Compute an initial stable timestep
  computeStableTimestep(patch, matl, new_dw);
}

/* Compute a stable initial timestep */
void
PolarOrthotropicLinearElasticStateModel::computeStableTimestep(const Uintah::Patch* patch,
                                                               const PeridynamicsMaterial* matl,
                                                               Uintah::DataWarehouse* new_dw)
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

  double longitudinal_modulus = std::max(std::max(d_cm.stiffnessMatrix(1,1),d_cm.stiffnessMatrix(2,2)),
                                         d_cm.stiffnessMatrix(3,3));

  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

     Uintah::particleIndex idx = *iter;

     // Compute wave speed at each particle, store the maximum
     Vector vel(0.0, 0.0, 0.0);
     if (pMass[idx] > 0.0) {
       speed_of_sound = std::sqrt(longitudinal_modulus*pVolume[idx]/pMass[idx]);
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
    new_dw->put(delt_vartype(std::numeric_limits<double>::max()), d_label->delTLabel, patch->getLevel());
  } else {
    new_dw->put(delt_vartype(delT_new), d_label->delTLabel, patch->getLevel());
  }
}

void 
PolarOrthotropicLinearElasticStateModel::addComputesAndRequires(Task* task, 
                                                             const PeridynamicsMaterial* matl,
                                                             const PatchSet* patches) const
{
  // Constants
  Ghost::GhostType gnone = Ghost::None;
  //Ghost::GhostType gac   = Ghost::AroundCells;

  // Get the current material
  const MaterialSubset* matlset = matl->thisMaterial();

  // List the variables needed for this task to execute
  task->requires(Task::OldDW, d_label->delTLabel,              matlset, gnone);
  task->requires(Task::OldDW, d_label->pPositionLabel,         matlset, gnone);
  task->requires(Task::OldDW, d_label->pMassLabel,             matlset, gnone);
  task->requires(Task::OldDW, d_label->pVelocityLabel,         matlset, gnone);
  task->requires(Task::OldDW, d_label->pDefGradLabel,          matlset, gnone);
  task->requires(Task::NewDW, d_label->pDefGradLabel_preReloc, matlset, gnone);
  task->requires(Task::OldDW, d_label->pStressLabel,           matlset, gnone);

  // List the variables computed by this task
  task->computes(d_label->pVolumeLabel_preReloc, matlset);
  task->computes(d_label->pStressLabel_preReloc, matlset);
  task->computes(d_label->pPK1StressLabel_preReloc, matlset);
}

void 
PolarOrthotropicLinearElasticStateModel::computeStressTensor(const PatchSubset* patches,
                                                             const PeridynamicsMaterial* matl,
                                                             DataWarehouse* old_dw,
                                                             DataWarehouse* new_dw)
{
  // Set up constants
  Matrix3 One; One.Identity();

  // Get the timestep size
  delt_vartype delT;
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

    constParticleVariable<Point> pPosition_old;
    old_dw->get(pPosition_old, d_label->pPositionLabel, pset);

    constParticleVariable<Vector> pVelocity_old;
    old_dw->get(pVelocity_old, d_label->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new;
    old_dw->get(pDefGrad_old, d_label->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, d_label->pDefGradLabel_preReloc, pset);

    constParticleVariable<Matrix3> pStress_old;
    old_dw->get(pStress_old, d_label->pStressLabel, pset);

    // Initialize the variables to be updated
    ParticleVariable<double> pVolume_new;
    new_dw->allocateAndPut(pVolume_new, d_label->pVolumeLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStress_new, pPK1Stress_new;
    new_dw->allocateAndPut(pStress_new, d_label->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pPK1Stress_new, d_label->pPK1StressLabel_preReloc, pset);

    // Loop through particles
    double rho_0 = matl->getInitialDensity();
    double longitudinalModulus = std::max(std::max(d_cm.stiffnessMatrix(1,1),d_cm.stiffnessMatrix(2,2)),
                                                   d_cm.stiffnessMatrix(3,3));
    Vector waveSpeed(std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min());

    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
      
      // Get particle index
      Uintah::particleIndex idx = *iter;

      // Compute the polar decomposition of the deformation gradient (F = RU)
      Matrix3 FF = pDefGrad_new[idx];
      Matrix3 RR, UU;
      FF.polarDecompositionRMB(UU, RR);

      // Compute the rate of deformation (d)
      // 1) Estimate the material time derivative of the deformation gradient (Forward Euler)   
      // 2) Compute F^{-1}
      // 3) Compute the velocity gradient l = Fdot.Finv
      // 4) Compute the rate of deformation d = 1/2(l + l^T)
      Matrix3 Fdot = (FF - pDefGrad_old[idx])*(1.0/delT);
      Matrix3 Finv = FF.Inverse();
      Matrix3 ll = Fdot*Finv;
      Matrix3 dd = (ll + ll.Transpose())*0.5;

      // Unrotate the stress and the rate of deformation (sig_rot = R^T sig R, d_rot = R^T d R)
      Matrix3 stress_old_unrotated = (RR.Transpose())*(pStress_old[idx]*RR);
      Matrix3 d_unrotated = (RR.Transpose())*(dd*RR);
      
      // Compute stress
      // This is the operation dsigma_rot/dt
      // 1) Express the stress and rate of deformation components in a rectangular coordinate system aligned with the
      //    axis of cylindrical anisotropy (assuming that the default state is that the global 3-axis is 
      //    aligned with the z-axis of the cylinder)
      // 2) Convert the stress and rate of deformation components from rectangular to cylindrical
      // 3) Update the stress using the constitutive relation
      // 4) Convert the stress components from cylindrical to rectangular 
      // 5) Express the stress components in a coordinate system aligned with the global coordinate system

      // Step 1:
      Vector axis_e3(0.0, 0.0, 1.0);
      Vector axis_ez = d_cm.top - d_cm.bottom;
      axis_ez.normalize();
      double angle = std::acos(Uintah::Dot(axis_e3,axis_ez));
      Vector rot_axis = Uintah::Cross(axis_e3, axis_ez); 
      rot_axis.normalize();

      Matrix3 stress_zaligned, d_zaligned;
      Vaango::Matrix9d QQ;
      Vaango::formRotationMatrix(angle, rot_axis, QQ);
      Vaango::rotateMatrix(QQ, stress_old_unrotated, stress_zaligned);
      Vaango::rotateMatrix(QQ, d_unrotated, d_zaligned);

      // Step 2:
      //   a) Project the particle on to the r-theta plane assuming that the bottom of the
      //      axis vector is the origin.  
      //   b) Compute (r,theta) for the particle
      //   c) Transform stress and rate of deformation to cylindrical coordinates
      Matrix3 nn(axis_ez, axis_ez);
      Matrix3 projMatrix = One - nn;
      Vector particleLoc = pPosition_old[idx] - d_cm.bottom;
      Vector particleProj = projMatrix*particleLoc;
      // The radial position will be need if there is r-variability
      //double rr = (particleProj - d_cm.bottom).length();
      Vector axis_e1(1.0, 0.0, 0.0);
      Vector axisLoc = axis_e1 - d_cm.bottom;
      Vector axisProj = projMatrix*axisLoc;
      double theta = std::acos(Uintah::Dot((particleProj - d_cm.bottom).normal(),
                                          (axisProj - d_cm.bottom).normal()));
      double cc = std::cos(theta);
      double ss = std::sin(theta);
      Matrix3 Transform(cc, ss, 0.0, -ss, cc, 0.0, 0.0, 0.0, 1.0);
      Matrix3 stress_cyl = (Transform*stress_zaligned)*Transform.Transpose();
      Matrix3 d_cyl = (Transform*d_zaligned)*Transform.Transpose();

      // Step 3:
      Matrix3 stress_cyl_new = stress_cyl + d_cm.stiffnessMatrix*(d_cyl*delT);

      // Step 4:     
      Matrix3 stress_rect_new = Transform.Transpose()*(stress_cyl_new*Transform);
      
      // Step 5:
      Matrix3 stress_global_new;
      Vaango::rotateMatrix(QQ.transpose(), stress_rect_new, stress_global_new);

      // Rotate the stress back (sig = R sigma_rot R^T)
      Matrix3 stress_new_rotated = (RR*stress_global_new)*(RR.Transpose());
      pStress_new[idx] = stress_new_rotated;

      // Compute PK1 stress
      double J = pDefGrad_new[idx].Determinant();
      pPK1Stress_new[idx] = pStress_new[idx]*(pDefGrad_new[idx].Transpose()*J);
 
      // Update the particle volume
      double rho_new = rho_0/J;
      pVolume_new[idx] = pMass[idx]/rho_new;
      
      // Compute the wavespeed at each particle and store the maximum
      double speed_of_sound = std::sqrt(longitudinalModulus/rho_new);
      Vector vel(0.0, 0.0, 0.0);
      vel[0] = speed_of_sound + std::abs(pVelocity_old[idx].x());
      vel[1] = speed_of_sound + std::abs(pVelocity_old[idx].y());
      vel[2] = speed_of_sound + std::abs(pVelocity_old[idx].z());
      waveSpeed = Vector(std::max(vel.x(), waveSpeed.x()),
                         std::max(vel.y(), waveSpeed.y()),
                         std::max(vel.z(), waveSpeed.z()));

    } // end particles loop

    // Find the grid spacing and update deltT
    Vector dx = patch->dCell();
    waveSpeed = dx/waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), d_label->delTLabel, patch->getLevel());

  } // end patches loop
}

// Register the permanent particle state associated with this material
void 
PolarOrthotropicLinearElasticStateModel::addParticleState(std::vector<const Uintah::VarLabel*>& ,
                                                          std::vector<const Uintah::VarLabel*>& )
{
}


