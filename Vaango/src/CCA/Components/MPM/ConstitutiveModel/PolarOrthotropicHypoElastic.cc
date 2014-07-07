#include <CCA/Components/MPM/ConstitutiveModel/PolarOrthotropicHypoElastic.h>

#include <Core/Labels/MPMLabel.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InternalError.h>

#include <Core/Math/Matrix3Rotation.h>

#include <iostream>
#include <limits>                           // for std::numeric_limits

using SCIRun::Point;
using SCIRun::Vector;
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

PolarOrthotropicHypoElastic::PolarOrthotropicHypoElastic(ProblemSpecP& ps,
                                                         MPMFlags* flags)
  : ConstitutiveModel(flags)
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

  // Set up the variables for the (r, theta, z) coordinates of the particles
  pRCoordLabel = VarLabel::create("p.RCoord", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel = VarLabel::create("p.ThetaCoord", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel = VarLabel::create("p.ZCoord", ParticleVariable<double>::getTypeDescription());
  pRCoordLabel_preReloc = VarLabel::create("p.RCoord+", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel_preReloc = VarLabel::create("p.ThetaCoord+", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel_preReloc = VarLabel::create("p.ZCoord+", ParticleVariable<double>::getTypeDescription());
}

PolarOrthotropicHypoElastic::PolarOrthotropicHypoElastic(const PolarOrthotropicHypoElastic* cm)
  : ConstitutiveModel(cm)
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

  // Set up the variables for the (r, theta, z) coordinates of the particles
  pRCoordLabel = VarLabel::create("p.RCoord", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel = VarLabel::create("p.ThetaCoord", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel = VarLabel::create("p.ZCoord", ParticleVariable<double>::getTypeDescription());
  pRCoordLabel_preReloc = VarLabel::create("p.RCoord+", ParticleVariable<double>::getTypeDescription());
  pThetaCoordLabel_preReloc = VarLabel::create("p.ThetaCoord+", ParticleVariable<double>::getTypeDescription());
  pZCoordLabel_preReloc = VarLabel::create("p.ZCoord+", ParticleVariable<double>::getTypeDescription());
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

/*! Identify the variabless to be used in the initialization task */
void 
PolarOrthotropicHypoElastic::addInitialComputesAndRequires(Task* task,
                                                           const MPMMaterial* matl,
                                                           const PatchSet* ) const
{
  // **NOTE** The initialization is done in ConstitutiveModel.cc
  //   for all particle variables other than the ones definde here

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

  // Set up constants
  Matrix3 One; One.Identity();

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

  //  a) Project the particle on to the r-theta plane assuming that the bottom of the
  //     axis vector is the origin.  
  //  b) Compute (r,theta,z) for the particle
  Vector globalAxisE1(1.0, 0.0, 0.0);
  Vector globalAxisE3(0.0, 0.0, 1.0);
  Vector cylAxisEz = d_cm.top - d_cm.bottom;
  cylAxisEz.normalize();
  Matrix3 nn(cylAxisEz, cylAxisEz);
  Matrix3 pp = One - nn;

  for (auto iter = pset->begin(); iter != pset->end(); iter++) {

    particleIndex idx = *iter;

    // Find the cylindrical z-coord of particle
    Vector cylTopTran = d_cm.top - d_cm.bottom;
    Vector pPosTran = pPosition[idx] - d_cm.bottom;
    Vector pPosProjEz = nn*pPosTran;
    double tt = SCIRun::Dot(pPosProjEz, pPosProjEz)/SCIRun::Dot(cylTopTran, pPosProjEz);
    double cylZMax = cylTopTran.length();
    double zz = tt*cylZMax;

    // Find the cylindrical r-coord of particle
    Vector pPosRVecProj = pp*pPosTran;
    double rr = pPosRVecProj.length();

    // Find the cylindrical theta-coord of particle
    Vector pPosThetaVecProj = pp*globalAxisE1;
    double theta = std::acos(SCIRun::Dot(pPosRVecProj.normal(),
                                         pPosThetaVecProj.normal()));

    pRCoord[idx] = rr;
    pThetaCoord[idx] = theta;
    pZCoord[idx] = zz;
   
    std::cout << "Particle " << idx << " (x,y,z) = " << pPosition[idx]
             << " (r,theta,z) = " << rr << ", " << theta << ", " << zz << std::endl;
  }

  // Compute an initial stable timestep
  computeStableTimestep(patch, matl, new_dw);
}

/* Compute a stable initial timestep */
void
PolarOrthotropicHypoElastic::computeStableTimestep(const Uintah::Patch* patch,
                                                   const MPMMaterial* matl,
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

  new_dw->get(pMass,     lb->pMassLabel,     pset);
  new_dw->get(pVolume,   lb->pVolumeLabel,   pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

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
    new_dw->put(delt_vartype(std::numeric_limits<double>::max()), lb->delTLabel, patch->getLevel());
  } else {
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
  }
}

void 
PolarOrthotropicHypoElastic::addComputesAndRequires(Task* task, 
                                                    const MPMMaterial* matl,
                                                    const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit 
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

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
  // Set up constants
  Matrix3 One; One.Identity();

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
    constParticleVariable<double> pMass;
    old_dw->get(pMass, lb->pMassLabel, pset);

    constParticleVariable<Point> pPosition_old;
    old_dw->get(pPosition_old, lb->pXLabel, pset);

    constParticleVariable<double> pRCoord, pThetaCoord, pZCoord;
    old_dw->get(pRCoord, pRCoordLabel, pset);
    old_dw->get(pThetaCoord, pThetaCoordLabel, pset);
    old_dw->get(pZCoord, pZCoordLabel, pset);

    constParticleVariable<Vector> pVelocity_old;
    old_dw->get(pVelocity_old, lb->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new;
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    constParticleVariable<Matrix3> pStress_old;
    old_dw->get(pStress_old, lb->pStressLabel, pset);

    // Initialize the variables to be updated
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    ParticleVariable<double> pdTdt;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);

    ParticleVariable<double> pQ;
    new_dw->allocateAndPut(pQ, lb->p_qLabel_preReloc, pset);

    //  Copy the material coordinates
    ParticleVariable<double> pRCoord_new, pThetaCoord_new, pZCoord_new;
    new_dw->allocateAndPut(pRCoord_new, pRCoordLabel_preReloc, pset);
    new_dw->allocateAndPut(pThetaCoord_new, pThetaCoordLabel_preReloc, pset);
    new_dw->allocateAndPut(pZCoord_new, pZCoordLabel_preReloc, pset);

    pRCoord_new.copyData(pRCoord);
    pThetaCoord_new.copyData(pThetaCoord);
    pZCoord_new.copyData(pZCoord);

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

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // No artificial viscosity
      pQ[idx] = 0.0;

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
      double angle = std::acos(SCIRun::Dot(axis_e3,axis_ez));
      Vector rot_axis = SCIRun::Cross(axis_e3, axis_ez); 
      rot_axis.normalize();

      Matrix3 stress_zaligned, d_zaligned;
      Vaango::Matrix9d QQ;
      Vaango::formRotationMatrix(angle, rot_axis, QQ);
      Vaango::rotateMatrix(QQ, stress_old_unrotated, stress_zaligned);
      Vaango::rotateMatrix(QQ, d_unrotated, d_zaligned);

      // Step 2: Transform stress and rate of deformation to cylindrical coordinates
      double cc = std::cos(pThetaCoord[idx]);
      double ss = std::sin(pThetaCoord[idx]);
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

      // Compute the wavespeed at each particle and store the maximum
      double J = pDefGrad_new[idx].Determinant();
      double rho_new = rho_0/J;
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
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

  } // end patches loop
}

// Register the permanent particle state associated with this material
void 
PolarOrthotropicHypoElastic::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                              std::vector<const Uintah::VarLabel*>& to)
{
  from.push_back(pRCoordLabel);
  from.push_back(pThetaCoordLabel);
  from.push_back(pZCoordLabel);
  to.push_back(pRCoordLabel_preReloc);
  to.push_back(pThetaCoordLabel_preReloc);
  to.push_back(pZCoordLabel_preReloc);
}

// Set up computes and requires for implicit time integration.
//        @todo:  This task has not been implemented yet. 
void 
PolarOrthotropicHypoElastic::addComputesAndRequires(Task* task,
                                                    const MPMMaterial* matl,
                                                    const PatchSet* patches,
                                                    const bool recursion) const
{
  std::ostringstream out;
  out << "**ERROR** Implicit time integration not implemented yet for" ;
  out << " the polar orthotropic hypoelastic material model."; 
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
  out << "**ERROR** No computation of rho_micro is available for " ;
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
  out << "**ERROR** No computation of pressure EOS is available for " ;
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
  out << "**ERROR** No computation of compressibility is available for " ;
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use MPMICE."; 
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Carry forward CM data for RigidMPM 
void 
PolarOrthotropicHypoElastic::carryForward(const PatchSubset* patches,
                                          const MPMMaterial* matl,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  std::ostringstream out;
  out << "**ERROR** RigigMPM cannot be use in conjunction with " ;
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to use RigidMPM."; 
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
  out << "**ERROR** Conversion to another material cannot be use in conjunction with " ;
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to model material conversion."; 
  throw InternalError(out.str(), __FILE__, __LINE__);
}

// Allocate variables for situations where particles are to be 
//        transformed into a different type of material 
void 
PolarOrthotropicHypoElastic::allocateCMDataAdd(DataWarehouse* new_dw,
                                               ParticleSubset* subset,
                                               map<const VarLabel*, ParticleVariableBase*>* newState,
                                               ParticleSubset* delset,
                                               DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "**ERROR** Conversion to another material cannot be use in conjunction with " ;
  out << " the polar orthotropic hypoelastic material model. " << std::endl;
  out << " Please choose another material model if you wish to model material conversion."; 
  throw InternalError(out.str(), __FILE__, __LINE__);
}


