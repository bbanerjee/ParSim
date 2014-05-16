#include <CCA/Components/Peridynamics/MaterialModels/LinearElasticStateModel.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype

using namespace Vaango;

LinearElasticStateModel::LinearElasticStateModel(Uintah::ProblemSpecP& ps,
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

LinearElasticStateModel::LinearElasticStateModel(const LinearElasticStateModel* cm)
  : PeridynamicsMaterialModel(cm)
{
  d_cm.bulkModulus = cm->d_cm.bulkModulus;
  d_cm.shearModulus = cm->d_cm.shearModulus;
}

// Make a clone of the constitutive model
LinearElasticStateModel* 
LinearElasticStateModel::clone()
{
  return scinew LinearElasticStateModel(*this);
}

LinearElasticStateModel::~LinearElasticStateModel()
{
}

void 
LinearElasticStateModel::outputProblemSpec(Uintah::ProblemSpecP& ps,
                                          bool output_cm_tag)
{
  Uintah::ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("material_model");
    cm_ps->setAttribute("type", "linear_elastic_state_based");
  }
  cm_ps->appendElement("bulk_modulus", d_cm.bulkModulus);
  cm_ps->appendElement("shear_modulus", d_cm.shearModulus);
}

/*! Identify the variabless to be used in the initialization task */
void 
LinearElasticStateModel::addInitialComputesAndRequires(Uintah::Task* task,
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
LinearElasticStateModel::initialize(const Uintah::Patch* patch,
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
LinearElasticStateModel::addComputesAndRequires(Uintah::Task* task, 
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
LinearElasticStateModel::computeStress(const Uintah::PatchSubset* patches,
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
LinearElasticStateModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                          std::vector<const Uintah::VarLabel*>& to)
{
  // These are common to all models and will have to be moved up
  from.push_back(d_label->pDefGradLabel);
  to.push_back(d_label->pDefGradLabel_preReloc);

  from.push_back(d_label->pStressLabel);
  to.push_back(d_label->pStressLabel_preReloc);
}


