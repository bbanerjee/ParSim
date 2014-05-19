#include <CCA/Components/Peridynamics/MaterialModels/PolarOrthotropicLinearElasticStateModel.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype
#include <Core/Exceptions/ProblemSetupException.h>

#include <iostream>

using namespace Vaango;

using Uintah::Matrix3;
using Uintah::SymmMatrix6;
using Uintah::ProblemSetupException;
using Uintah::ProblemSpecP;
using Uintah::ParticleSubset;
using Uintah::ParticleVariable;
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
  ps->require("G_r_heta", d_cm.Grtheta);

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
  task->computes(d_label->pStressLabel, matlset);
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
  ParticleVariable<Matrix3> pStress;
  new_dw->allocateAndPut(pStress, d_label->pStressLabel, pset);

  // Initialize the stress to zero (for now)
  Matrix3 zero(0.0);
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
    pStress[*iter] = zero; 
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
  task->requires(Task::OldDW, d_label->pVolumeLabel,           matlset, gnone);
  task->requires(Task::OldDW, d_label->pDefGradLabel,          matlset, gnone);
  task->requires(Task::NewDW, d_label->pDefGradLabel_preReloc, matlset, gnone);
  task->requires(Task::OldDW, d_label->pStressLabel,           matlset, gnone);

  // List the variables computed by this task
  task->computes(d_label->pStressLabel_preReloc, matlset);
}

void 
PolarOrthotropicLinearElasticStateModel::computeStress(const PatchSubset* patches,
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
    Uintah::constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new;
    old_dw->get(pDefGrad_old, d_label->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, d_label->pDefGradLabel_preReloc, pset);

    Uintah::constParticleVariable<Matrix3> pStress_old;
    old_dw->get(pStress_old, d_label->pStressLabel, pset);

    // Initialize the variables to be updated
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, d_label->pStressLabel_preReloc, pset);

    // Loop through particles
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
      Matrix3 dd_unrotated = (RR.Transpose())*(dd*RR);
      
      // Compute stress
      // This is the operation dsigma_rot/dt
      // 1) Express the stress and rate of deformation components in a rectangular coordinate system aligned with the
      //    axis of cylindrical anisotropy
      // 2) Convert the stress and rate of deformation components from rectangular to cylindrical
      // 3) Update the stress using the constitutive relation
      // 4) Convert the stress components from cylindrical to rectangular 
      // 5) Express the stress components in a coordinate system aligned with the global coordinate system
      double pressure = 0.0;
      pStress_new[idx] = One*pressure;

      // Rotate the stress back (sig = R sigma_rot R^T)
      Matrix3 pStress = pStress_new[idx];
      Matrix3 stress_new_rotated = (RR*pStress)*(RR.Transpose());
     
    }

  } // end patches loop
}

// Register the permanent particle state associated with this material
void 
PolarOrthotropicLinearElasticStateModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                                       std::vector<const Uintah::VarLabel*>& to)
{
  // These are common to all models and will have to be moved up
  from.push_back(d_label->pDefGradLabel);
  to.push_back(d_label->pDefGradLabel_preReloc);

  from.push_back(d_label->pStressLabel);
  to.push_back(d_label->pStressLabel_preReloc);
}


