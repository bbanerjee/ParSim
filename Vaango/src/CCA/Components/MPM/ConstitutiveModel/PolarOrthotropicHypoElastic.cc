/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/PolarOrthotropicHypoElastic.h>

#include <CCA/Components/MPM/ConstitutiveModel/Constants.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <Core/Labels/MPMLabel.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Variables/VarTypes.h> // for delt_vartype

#include <Core/Math/Matrix3Rotation.h>
#include <Core/Util/DebugStream.h>

#include <iostream>
#include <limits> // for std::numeric_limits

using Uintah::Point;
using Uintah::Vector;
using Uintah::Matrix3;
using Uintah::SymmMatrix6;
using Uintah::delt_vartype;
using Uintah::ProblemSetupException;
using Uintah::InternalError;
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
using Uintah::MPMFlags;
using Uintah::PolarOrthotropicHypoElastic;

static Uintah::DebugStream dbg("Polar", false);
static Uintah::DebugStream dbg_extra("PolarExtra", false);

PolarOrthotropicHypoElastic::PolarOrthotropicHypoElastic(ProblemSpecP& ps,
                                                         MPMFlags* flags)
  : ConstitutiveModel(flags)
{
  // Get the axis of symmetry of the cylinder
  ps->require("symmetry_axis_top", d_cm.top);
  ps->require("symmetry_axis_bottom", d_cm.bottom);
  if ((d_cm.top - d_cm.bottom).length2() == 0.0) {
    std::ostringstream out;
    out << "The axis of symmetry is a point. ";
    out << "Please check the values for the top and bottom points of the axis "
           "in the input file";
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
  complianceMatrix(0, 0) = 1.0 / d_cm.Er;
  complianceMatrix(0, 1) = -d_cm.nuthetar / d_cm.Etheta;
  complianceMatrix(0, 2) = -d_cm.nuzr / d_cm.Ez;
  complianceMatrix(1, 1) = 1.0 / d_cm.Etheta;
  complianceMatrix(1, 2) = -d_cm.nuztheta / d_cm.Ez;
  complianceMatrix(2, 2) = 1.0 / d_cm.Ez;
  complianceMatrix(3, 3) = 1.0 / d_cm.Gthetaz;
  complianceMatrix(4, 4) = 1.0 / d_cm.Gzr;
  complianceMatrix(5, 5) = 1.0 / d_cm.Grtheta;

  complianceMatrix(1, 0) = complianceMatrix(0, 1);
  complianceMatrix(2, 0) = complianceMatrix(0, 2);
  complianceMatrix(2, 1) = complianceMatrix(1, 2);

  // Check that everything is consistent
  if ((complianceMatrix(0, 0) < 0.0) || (complianceMatrix(1, 1) < 0.0) ||
      (complianceMatrix(2, 2) < 0.0) || (complianceMatrix(3, 3) < 0.0) ||
      (complianceMatrix(4, 4) < 0.0) || (complianceMatrix(5, 5) < 0.0)) {
    std::ostringstream out;
    out << "The compliance matrix has negative diagonal components";
    out << "Please check the values in the input file to make sure the Young's "
           "and shear moduli are positive.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  double delta2 = complianceMatrix(0, 0) * complianceMatrix(1, 1) -
                  complianceMatrix(0, 1) * complianceMatrix(0, 1);
  if (delta2 < 0.0) {
    std::ostringstream out;
    out << "Compliance matrix submatrix has negative determinant: S11 S22 - "
           "S12^2 < 0.";
    out << "Please check the values in the input file to make sure the input "
           "data are correct.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  double delta3 =
    delta2 * complianceMatrix(2, 2) -
    complianceMatrix(0, 0) * complianceMatrix(1, 2) * complianceMatrix(1, 2) +
    2.0 * complianceMatrix(0, 1) * complianceMatrix(1, 2) *
      complianceMatrix(0, 2) -
    complianceMatrix(1, 1) * complianceMatrix(0, 2) * complianceMatrix(0, 2);
  if (delta3 < 0.0) {
    std::ostringstream out;
    out << "Compliance matrix submatrix has negative determinant: ";
    out
      << "  (S11 S22 - S12^2)S33 - S11 S23^2 + 2 S12 S23 S13 - S22 S13^2 < 0.";
    out << "Please check the values in the input file to make sure the input "
           "data are correct.";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Compute stiffness matrix
  complianceMatrix.inverse(d_cm.stiffnessMatrix);

  // Set up the variables for the (r, theta, z) coordinates of the particles
  pRCoordLabel = VarLabel::create(
    "p.RCoord", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel = VarLabel::create(
    "p.ThetaCoord", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel = VarLabel::create(
    "p.ZCoord", ParticleVariable<double>::getTypeDescription());
  pRCoordLabel_preReloc = VarLabel::create(
    "p.RCoord+", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel_preReloc = VarLabel::create(
    "p.ThetaCoord+", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel_preReloc = VarLabel::create(
    "p.ZCoord+", ParticleVariable<double>::getTypeDescription());
}

PolarOrthotropicHypoElastic::PolarOrthotropicHypoElastic(
  const PolarOrthotropicHypoElastic* cm)
  : ConstitutiveModel(cm)
{
  d_cm.top             = cm->d_cm.top;
  d_cm.bottom          = cm->d_cm.bottom;
  d_cm.Er              = cm->d_cm.Er;
  d_cm.Etheta          = cm->d_cm.Etheta;
  d_cm.Ez              = cm->d_cm.Ez;
  d_cm.nuthetar        = cm->d_cm.nuthetar;
  d_cm.nuzr            = cm->d_cm.nuzr;
  d_cm.nuztheta        = cm->d_cm.nuztheta;
  d_cm.Gthetaz         = cm->d_cm.Gthetaz;
  d_cm.Gzr             = cm->d_cm.Gzr;
  d_cm.Grtheta         = cm->d_cm.Grtheta;
  d_cm.stiffnessMatrix = cm->d_cm.stiffnessMatrix;

  // Set up the variables for the (r, theta, z) coordinates of the particles
  pRCoordLabel = VarLabel::create(
    "p.RCoord", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel = VarLabel::create(
    "p.ThetaCoord", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel = VarLabel::create(
    "p.ZCoord", ParticleVariable<double>::getTypeDescription());
  pRCoordLabel_preReloc = VarLabel::create(
    "p.RCoord+", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel_preReloc = VarLabel::create(
    "p.ThetaCoord+", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel_preReloc = VarLabel::create(
    "p.ZCoord+", ParticleVariable<double>::getTypeDescription());
}

// Make a clone of the constitutive model
PolarOrthotropicHypoElastic*
PolarOrthotropicHypoElastic::clone()
{
  return scinew PolarOrthotropicHypoElastic(*this);
}

PolarOrthotropicHypoElastic::~PolarOrthotropicHypoElastic()
{
  VarLabel::destroy(pRCoordLabel);
  VarLabel::destroy(pThetaCoordLabel);
  VarLabel::destroy(pZCoordLabel);
  VarLabel::destroy(pRCoordLabel_preReloc);
  VarLabel::destroy(pThetaCoordLabel_preReloc);
  VarLabel::destroy(pZCoordLabel_preReloc);
}

void
PolarOrthotropicHypoElastic::outputProblemSpec(ProblemSpecP& ps,
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

// Register the permanent particle state associated with this material
void
PolarOrthotropicHypoElastic::addParticleState(
  std::vector<const Uintah::VarLabel*>& from,
  std::vector<const Uintah::VarLabel*>& to)
{
  // This is an INCREMENTAL model. Needs polar decomp R to be saved.
  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);

  from.push_back(pRCoordLabel);
  from.push_back(pThetaCoordLabel);
  from.push_back(pZCoordLabel);
  to.push_back(pRCoordLabel_preReloc);
  to.push_back(pThetaCoordLabel_preReloc);
  to.push_back(pZCoordLabel_preReloc);
}

/*! Identify the variabless to be used in the initialization task */
void
PolarOrthotropicHypoElastic::addInitialComputesAndRequires(
  Task* task,
  const MPMMaterial* matl,
  const PatchSet*) const
{
  // **NOTE** The initialization is done in ConstitutiveModel.cc
  //   for all particle variables other than the ones defined here

  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pRCoordLabel, matlset);
  task->computes(pThetaCoordLabel, matlset);
  task->computes(pZCoordLabel, matlset);
}

/*! Initialize the variables used in the CM */
void
PolarOrthotropicHypoElastic::initializeCMData(const Patch* patch,
                                              const MPMMaterial* matl,
                                              DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  // Set up the (r, theta, z) coordinates of each particle.  These are
  // material coordinates - not geometrical coordinates
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  constParticleVariable<Point> pPosition;
  new_dw->get(pPosition, lb->pXLabel, pset);

  ParticleVariable<double> pRCoord;
  ParticleVariable<double> pThetaCoord;
  ParticleVariable<double> pZCoord;
  new_dw->allocateAndPut(pRCoord, pRCoordLabel, pset);
  new_dw->allocateAndPut(pThetaCoord, pThetaCoordLabel, pset);
  new_dw->allocateAndPut(pZCoord, pZCoordLabel, pset);

  //  a) Project the particle on to the r-theta plane assuming that the bottom
  //  of the
  //     axis vector is the origin.
  //  b) Compute (r,theta,z) for the particle
  Vector globalAxisE1(1.0, 0.0, 0.0);
  Vector globalAxisE3(0.0, 0.0, 1.0);
  Vector cylAxisEz = d_cm.top - d_cm.bottom;
  cylAxisEz.normalize();
  Matrix3 nn(cylAxisEz, cylAxisEz);
  Matrix3 pp = Vaango::Util::Identity - nn;

  for (int idx : *pset) {

    // Find the cylindrical z-coord of particle
    Vector cylTopTran = d_cm.top - d_cm.bottom;
    Vector pPosTran   = pPosition[idx] - d_cm.bottom;
    Vector pPosProjEz = nn * pPosTran;
    double tt =
      Uintah::Dot(pPosProjEz, pPosProjEz) / Uintah::Dot(cylTopTran, pPosProjEz);
    double cylZMax = cylTopTran.length();
    double zz      = tt * cylZMax;

    // Find the cylindrical r-coord of particle
    Vector pPosRVecProj = pp * pPosTran;
    double rr           = pPosRVecProj.length();

    dbg_extra << " pp = " << pp << " pPosTran = " << pPosTran
              << " pPosRVecProj = " << pPosRVecProj << std::endl;

    // Find the cylindrical theta-coord of particle
    //  (between 0 and pi)
    double theta = 0.0;
    if (rr > 0.0) {
      Vector pPosThetaVecProj = pp * globalAxisE1;
      dbg_extra << " pp = " << pp << " globalAxisE1 = " << globalAxisE1
                << " pPosThetaVecProj = " << pPosThetaVecProj << std::endl;
      theta = std::acos(
        Uintah::Dot(pPosRVecProj.normal(), pPosThetaVecProj.normal()));
    }

    pRCoord[idx]     = rr;
    pThetaCoord[idx] = theta;
    pZCoord[idx]     = zz;

    dbg_extra << "Particle " << idx << " (x,y,z) = " << pPosition[idx]
              << " (r,theta,z) = " << rr << ", " << theta << ", " << zz
              << std::endl;
  }

  // Compute an initial stable timestep
  computeStableTimestep(patch, matl, new_dw);
}

/* Compute a stable initial timestep */
void
PolarOrthotropicHypoElastic::computeStableTimestep(
  const Uintah::Patch* patch,
  const MPMMaterial* matl,
  Uintah::DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx     = patch->dCell();
  int matlIndex = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double longitudinal_modulus =
    std::max(std::max(d_cm.stiffnessMatrix(1, 1), d_cm.stiffnessMatrix(2, 2)),
             d_cm.stiffnessMatrix(3, 3));

  double speed_of_sound = 0.0;
  Vector waveSpeed(std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min());

  for (int idx : *pset) {

    // Compute wave speed at each particle, store the maximum
    Vector vel(0.0, 0.0, 0.0);
    if (pMass[idx] > 0.0) {
      speed_of_sound =
        std::sqrt(longitudinal_modulus * pVolume[idx] / pMass[idx]);
    } else {
      speed_of_sound = 0.0;
    }
    Vector velMax = pVelocity[idx].cwiseAbs() + speed_of_sound;
    waveSpeed     = Max(velMax, waveSpeed);
  }

  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  if (delT_new < 1.e-12) {
    new_dw->put(delt_vartype(std::numeric_limits<double>::max()),
                lb->delTLabel,
                patch->getLevel());
  } else {
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
  }
}

void
PolarOrthotropicHypoElastic::addComputesAndRequires(
  Task* task,
  const MPMMaterial* matl,
  const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);

  // Only the local computes and requires
  task->requires(Task::OldDW, pRCoordLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pThetaCoordLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pZCoordLabel, matlset, Ghost::None);
  task->computes(pRCoordLabel_preReloc, matlset);
  task->computes(pThetaCoordLabel_preReloc, matlset);
  task->computes(pZCoordLabel_preReloc, matlset);
}

void
PolarOrthotropicHypoElastic::computeStressTensor(const PatchSubset* patches,
                                                 const MPMMaterial* matl,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw)
{
  double rho_0 = matl->getInitialDensity();
  double longitudinalModulus =
    std::max(std::max(d_cm.stiffnessMatrix(1, 1), d_cm.stiffnessMatrix(2, 2)),
             d_cm.stiffnessMatrix(3, 3));

  // Get the timestep size
  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  // Loop through patches
  for (int p = 0; p < patches->size(); p++) {

    // Get the current patch
    const Patch* patch = patches->get(p);

    // Get the material index
    int matlIndex = matl->getDWIndex();

    // Get the particle subset for this material
    ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle variables needed
    constParticleVariable<double> pRCoord, pThetaCoord, pZCoord;
    old_dw->get(pRCoord, pRCoordLabel, pset);
    old_dw->get(pThetaCoord, pThetaCoordLabel, pset);
    old_dw->get(pZCoord, pZCoordLabel, pset);

    constParticleVariable<Vector> pVelocity_old;
    old_dw->get(pVelocity_old, lb->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pDefRate_mid, pDefGrad_new, pStress_old;
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);

    // Initialize the variables to be updated
    ParticleVariable<double> pdTdt, pQ;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(pQ, lb->p_qLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    //  Copy the material coordinates
    ParticleVariable<double> pRCoord_new, pThetaCoord_new, pZCoord_new;
    new_dw->allocateAndPut(pRCoord_new, pRCoordLabel_preReloc, pset);
    new_dw->allocateAndPut(pThetaCoord_new, pThetaCoordLabel_preReloc, pset);
    new_dw->allocateAndPut(pZCoord_new, pZCoordLabel_preReloc, pset);

    pRCoord_new.copyData(pRCoord);
    pThetaCoord_new.copyData(pThetaCoord);
    pZCoord_new.copyData(pZCoord);

    // Loop through particles
    Vector waveSpeed(std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min(),
                     std::numeric_limits<double>::min());

    for (int idx : *pset) {

      // Get particle index
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // No artificial viscosity
      pQ[idx] = 0.0;

      Matrix3 stress_old_unrotated = pStress_old[idx];
      Matrix3 d_unrotated          = pDefRate_mid[idx];

      dbg << " sig_rot = " << stress_old_unrotated << std::endl
          << " d_rot = " << d_unrotated << std::endl;

      // Compute stress
      // This is the operation dsigma_rot/dt
      // 1) Express the stress and rate of deformation components in a
      // rectangular coordinate system aligned with the
      //    axis of cylindrical anisotropy (assuming that the default state is
      //    that the global 3-axis is
      //    aligned with the z-axis of the cylinder)
      // 2) Convert the stress and rate of deformation components from
      // rectangular to cylindrical
      // 3) Update the stress using the constitutive relation
      // 4) Convert the stress components from cylindrical to rectangular
      // 5) Express the stress components in a coordinate system aligned with
      // the global coordinate system

      // Step 1:
      Vector axis_e3(0.0, 0.0, 1.0);
      Vector axis_ez = d_cm.top - d_cm.bottom;
      axis_ez.normalize();
      double angle    = std::acos(Uintah::Dot(axis_e3, axis_ez));
      Vector rot_axis = Uintah::Cross(axis_e3, axis_ez);
      rot_axis.normalize();

      dbg << " Rot angle = " << angle << " rot axis = " << rot_axis
          << std::endl;

      Matrix3 stress_zaligned, d_zaligned;
      Vaango::Matrix9d QQ;
      Vaango::formRotationMatrix(angle, rot_axis, QQ);
      dbg << " Rotation matrix = " << QQ << std::endl;

      Vaango::rotateMatrix(QQ, stress_old_unrotated, stress_zaligned);
      dbg << " Stress z aligned = " << stress_zaligned << std::endl;

      Vaango::rotateMatrix(QQ, d_unrotated, d_zaligned);
      dbg << " d z aligned = " << d_zaligned << std::endl;

      // Step 2: Transform stress and rate of deformation to cylindrical
      // coordinates
      double cc = std::cos(pThetaCoord[idx]);
      double ss = std::sin(pThetaCoord[idx]);
      Matrix3 Transform(cc, ss, 0.0, -ss, cc, 0.0, 0.0, 0.0, 1.0);
      Matrix3 stress_cyl =
        (Transform * stress_zaligned) * Transform.Transpose();
      dbg << " Stress cyl = " << stress_cyl << std::endl;

      Matrix3 d_cyl = (Transform * d_zaligned) * Transform.Transpose();
      dbg << " d cyl = " << d_cyl << std::endl;

      // Step 3:
      Matrix3 stress_cyl_new =
        stress_cyl + d_cm.stiffnessMatrix * (d_cyl * delT);

      // Step 4:
      Matrix3 stress_rect_new =
        Transform.Transpose() * (stress_cyl_new * Transform);

      // Step 5:
      Matrix3 stress_global_new;
      Vaango::rotateMatrix(QQ.transpose(), stress_rect_new, stress_global_new);
      pStress_new[idx]           = stress_global_new;

      // Compute the wavespeed at each particle and store the maximum
      double J              = pDefGrad_new[idx].Determinant();
      double rho_new        = rho_0 / J;
      double speed_of_sound = std::sqrt(longitudinalModulus / rho_new);
      Vector velMax = pVelocity_old[idx].cwiseAbs() + speed_of_sound;
      waveSpeed     = Max(velMax, waveSpeed);

    } // end particles loop

    // Find the grid spacing and update deltT
    Vector dx       = patch->dCell();
    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

  } // end patches loop
}

// Set up computes and requires for implicit time integration.
//        @todo:  This task has not been implemented yet.
void
PolarOrthotropicHypoElastic::addComputesAndRequires(Task* task,
                                                    const MPMMaterial* matl,
                                                    const PatchSet* patches,
                                                    const bool recursion,
                                                    const bool schedPar) const
{
  std::ostringstream out;
  out << "**ERROR** Implicit time integration not implemented yet for";
  out << " the polar orthotropic hypoelastic material model.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Set up task variables for situations where particles are moved to
//        another material type
void
PolarOrthotropicHypoElastic::allocateCMDataAddRequires(Task* task,
                                                       const MPMMaterial* matl,
                                                       const PatchSet* patch,
                                                       MPMLabel* lb) const
{
  std::ostringstream out;
  out << "**ERROR** Conversion to another material cannot be use in "
         "conjunction with ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to model material "
         "conversion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Allocate variables for situations where particles are to be
//        transformed into a different type of material
void
PolarOrthotropicHypoElastic::allocateCMDataAdd(
  DataWarehouse* new_dw,
  ParticleSubset* subset,
  ParticleLabelVariableMap* newState,
  ParticleSubset* delset,
  DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "**ERROR** Conversion to another material cannot be use in "
         "conjunction with ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to model material "
         "conversion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Carry forward CM data for RigidMPM
void
PolarOrthotropicHypoElastic::carryForward(const PatchSubset* patches,
                                          const MPMMaterial* matl,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  std::ostringstream out;
  out << "**ERROR** RigigMPM cannot be use in conjunction with ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use RigidMPM.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Compute grid cell microscopic density for MPMICE calculations.
//        @todo:  This task has not been implemented yet.
double
PolarOrthotropicHypoElastic::computeRhoMicroCM(double pressure,
                                               const double p_ref,
                                               const MPMMaterial* matl,
                                               double temperature,
                                               double rho_guess)
{
  std::ostringstream out;
  out << "**ERROR** No computation of rho_micro is available for ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use MPMICE.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Compute grid cell pressure using an equation of state
//        for MPMICE calculations.
//        @todo:  This task has not been implemented yet.
void
PolarOrthotropicHypoElastic::computePressEOSCM(double rho_m,
                                               double& press_eos,
                                               double p_ref,
                                               double& dp_drho,
                                               double& ss_new,
                                               const MPMMaterial* matl,
                                               double temperature)
{
  std::ostringstream out;
  out << "**ERROR** No computation of pressure EOS is available for ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use MPMICE.";
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Get the compressiblity (inverse of bulk modulus)
//        for MPMICE calculations.
//        @todo:  This task has not been implemented yet.
double
PolarOrthotropicHypoElastic::getCompressibility()
{
  std::ostringstream out;
  out << "**ERROR** No computation of compressibility is available for ";
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use MPMICE.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

