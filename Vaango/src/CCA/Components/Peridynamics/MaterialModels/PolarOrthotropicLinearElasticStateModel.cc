#include <CCA/Components/Peridynamics/MaterialModels/PolarOrthotropicLinearElasticStateModel.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype

#include <iostream>

using namespace Vaango;

using Uintah::Matrix3;

PolarOrthotropicLinearElasticStateModel::PolarOrthotropicLinearElasticStateModel(Uintah::ProblemSpecP& ps,
                                                                                 PeridynamicsFlags* flags)
  : PeridynamicsMaterialModel(flags)
{
  // Get the axis of symmetry of the cylinder 
  ps->require("symmetry_axis_top", d_top);
  ps->require("symmetry_axis_bottom", d_bottom);
  if ((d_top-d_bottom).length2() == 0.0) {
    std::ostringstream out;
    out << "The axis of symmetry is a point. ";
    out << "Please check the values for the top and bottom points of the axis in the input file";
    SCI_THROW(ProblemSetupException(out.str(), __FILE__, __LINE__));
  }

  // Get the elastic moduli with respect to the axis of symmetry 
  // 1 <-> r, 2 <-> theta, 3 <-> z
  ps->require("E_r", d_Er);
  ps->require("E_theta", d_Etheta);
  ps->require("E_z", d_Ez);
  ps->require("nu_theta_r", d_nuThetaR);
  ps->require("nu_z_r", d_nuZR);
  ps->require("nu_z_theta", d_nuZTheta);
  ps->require("G_theta_z", d_GThetaZ);
  ps->require("G_z_r", d_GZR);
  ps->require("G_r_heta", d_GRTheta);

  // Compute the compliance matrix
  
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

PolarOrthotropicLinearElasticStateModel::PolarOrthotropicLinearElasticStateModel(const PolarOrthotropicLinearElasticStateModel* cm)
  : PeridynamicsMaterialModel(cm)
{
  d_cm.bulkModulus = cm->d_cm.bulkModulus;
  d_cm.shearModulus = cm->d_cm.shearModulus;
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
PolarOrthotropicLinearElasticStateModel::outputProblemSpec(Uintah::ProblemSpecP& ps,
                                                        bool output_cm_tag)
{
  Uintah::ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("material_model");
    cm_ps->setAttribute("type", "elastic_neo_hookean_state");
  }
  cm_ps->appendElement("bulk_modulus", d_cm.bulkModulus);
  cm_ps->appendElement("shear_modulus", d_cm.shearModulus);
}

/*! Identify the variabless to be used in the initialization task */
void 
PolarOrthotropicLinearElasticStateModel::addInitialComputesAndRequires(Uintah::Task* task,
                                                                    const PeridynamicsMaterial* material,
                                                                    const Uintah::PatchSet* patches) const
{
  // Identify this material
  const Uintah::MaterialSubset* matlset = material->thisMaterial();

  // Add compute flags for the initialization of the stress
  task->computes(d_label->pStressLabel, matlset);
}

/*! Initialize the variables used in the CM */
void 
PolarOrthotropicLinearElasticStateModel::initialize(const Uintah::Patch* patch,
                                                 const PeridynamicsMaterial* matl,
                                                 Uintah::DataWarehouse* new_dw)
{
  // Get the set of particles of this material type in the current patch  
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  // Allocate for saving
  Uintah::ParticleVariable<Uintah::Matrix3> pStress;
  new_dw->allocateAndPut(pStress, d_label->pStressLabel, pset);

  // Initialize the stress to zero (for now)
  Uintah::Matrix3 zero(0.0);
  for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
    pStress[*iter] = zero; 
  }
}

void 
PolarOrthotropicLinearElasticStateModel::addComputesAndRequires(Uintah::Task* task, 
                                                             const PeridynamicsMaterial* matl,
                                                             const Uintah::PatchSet* patches) const
{
  // Constants
  Uintah::Ghost::GhostType gnone = Uintah::Ghost::None;
  //Uintah::Ghost::GhostType gac   = Uintah::Ghost::AroundCells;

  // Get the current material
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();

  // List the variables needed for this task to execute
  task->requires(Uintah::Task::OldDW, d_label->delTLabel,            matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pPositionLabel,       matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pMassLabel,           matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pVolumeLabel,         matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pDefGradLabel,        matlset, gnone);

  // List the variables computed by this task
  task->computes(d_label->pStressLabel, matlset);
}

void 
PolarOrthotropicLinearElasticStateModel::computeStress(const Uintah::PatchSubset* patches,
                                                    const PeridynamicsMaterial* matl,
                                                    Uintah::DataWarehouse* old_dw,
                                                    Uintah::DataWarehouse* new_dw)
{
  // Set up constants
  Uintah::Matrix3 One; One.Identity();

  // Get the timestep size
  Uintah::delt_vartype delT;
  old_dw->get(delT, d_label->delTLabel, getLevel(patches));
  
  // Loop through patches
  for (int p = 0; p < patches->size(); p++) {

    // Get the current patch
    const Uintah::Patch* patch = patches->get(p);

    // Set up variables used to compute stress
    Uintah::Matrix3 defGrad(1.0);

    // Get the material index
    int matlIndex = matl->getDWIndex();
 
    // Get the particle subset for this material
    Uintah::ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle variables needed
    Uintah::constParticleVariable<Uintah::Matrix3> pDefGrad;
    old_dw->get(pDefGrad, d_label->pDefGradLabel, pset);

    // Initialize the variables to be updated
    Uintah::ParticleVariable<Uintah::Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, d_label->pStressLabel_preReloc, pset);

    // Loop through particles
    for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
      
      // Get particle index
      Uintah::particleIndex idx = *iter;

      // Compute J = det(F)
      double J = pDefGrad[idx].Determinant();

      // Compute Bbar = J^{-2/3} (F F^T)  and dev(Bbar) = Bbar - 1/3 Tr(Bbar) I
      Uintah::Matrix3 Bbar = pDefGrad[idx]*(pDefGrad[idx].Transpose())*std::pow(J, -2.0/3.0);
      Uintah::Matrix3 BbarDev = Bbar - One*(Bbar.Trace()/3.0); 

      // Computes stress
      double pressure = -d_cm.bulkModulus*(J-1.0);
      pStress_new[idx] = One*pressure + BbarDev*(d_cm.shearModulus/J);
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


