#include <CCA/Components/Peridynamics/MaterialModels/IsotropicElasticNeoHookeanStateModel.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarTypes.h>   // for delt_vartype

#include <limits>                           // for std::numeric_limits

using namespace Vaango;

IsotropicElasticNeoHookeanStateModel::IsotropicElasticNeoHookeanStateModel(Uintah::ProblemSpecP& ps,
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
IsotropicElasticNeoHookeanStateModel::outputProblemSpec(Uintah::ProblemSpecP& ps,
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
IsotropicElasticNeoHookeanStateModel::addInitialComputesAndRequires(Uintah::Task* task,
                                                                    const PeridynamicsMaterial* material,
                                                                    const Uintah::PatchSet* patches) const
{
  // Identify this material
  const Uintah::MaterialSubset* matlset = material->thisMaterial();

  // Add compute flags for the initialization of the stress
  task->computes(d_label->pPK1StressLabel, matlset);
}

/*! Initialize the variables used in the CM */
void 
IsotropicElasticNeoHookeanStateModel::initialize(const Uintah::Patch* patch,
                                                 const PeridynamicsMaterial* matl,
                                                 Uintah::DataWarehouse* new_dw)
{
  // Get the set of particles of this material type in the current patch  
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  // Allocate for saving
  Uintah::ParticleVariable<Uintah::Matrix3> pStress, pPK1Stress;
  new_dw->allocateAndPut(pStress, d_label->pStressLabel, pset);
  new_dw->allocateAndPut(pPK1Stress, d_label->pPK1StressLabel, pset);

  // Initialize the stress to zero (for now)
  Uintah::Matrix3 zero(0.0);
  for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
    pStress[*iter] = zero; 
    pPK1Stress[*iter] = zero; 
  }

  // Compute a stable time step
  computeStableTimestep(patch, matl, new_dw);
}

/* Compute a stable initial timestep */
void
IsotropicElasticNeoHookeanStateModel::computeStableTimestep(const Uintah::Patch* patch,
                                                            const PeridynamicsMaterial* matl,
                                                            Uintah::DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  SCIRun::Vector dx = patch->dCell();
  int matlIndex = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);
  Uintah::constParticleVariable<double> pMass, pVolume;
  Uintah::constParticleVariable<SCIRun::Vector> pVelocity;

  new_dw->get(pMass,     d_label->pMassLabel,     pset);
  new_dw->get(pVolume,   d_label->pVolumeLabel,   pset);
  new_dw->get(pVelocity, d_label->pVelocityLabel, pset);

  double speed_of_sound = 0.0;
  SCIRun::Vector waveSpeed(std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min(),
                   std::numeric_limits<double>::min());

  double kappa = d_cm.bulkModulus;
  double mu = d_cm.shearModulus;
  double pWaveModulus = kappa + mu*(4.0/3.0);

  for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

     Uintah::particleIndex idx = *iter;

     // Compute wave speed at each particle, store the maximum
     SCIRun::Vector vel(0.0, 0.0, 0.0);
     if (pMass[idx] > 0.0) {
       speed_of_sound = std::sqrt(pWaveModulus*pVolume[idx]/pMass[idx]);
       vel[0] = speed_of_sound + std::abs(pVelocity[idx].x());
       vel[1] = speed_of_sound + std::abs(pVelocity[idx].y());
       vel[2] = speed_of_sound + std::abs(pVelocity[idx].z());
     } else {
       speed_of_sound = 0.0;
     }
     waveSpeed = SCIRun::Vector(std::max(vel.x(), waveSpeed.x()),
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
IsotropicElasticNeoHookeanStateModel::addComputesAndRequires(Uintah::Task* task, 
                                                             const PeridynamicsMaterial* matl,
                                                             const Uintah::PatchSet* patches) const
{
  // Constants
  Uintah::Ghost::GhostType gnone = Uintah::Ghost::None;
  //Uintah::Ghost::GhostType gac   = Uintah::Ghost::AroundCells;

  // Get the current material
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();

  // List the variables needed for this task to execute
  task->requires(Uintah::Task::OldDW, d_label->delTLabel,              matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pPositionLabel,         matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pMassLabel,             matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_label->pVolumeLabel,           matlset, gnone);
  task->requires(Uintah::Task::NewDW, d_label->pDefGradLabel_preReloc, matlset, gnone);

  // List the variables computed by this task
  task->computes(d_label->pStressLabel_preReloc,    matlset);
  task->computes(d_label->pPK1StressLabel_preReloc, matlset);
}

void 
IsotropicElasticNeoHookeanStateModel::computeStressTensor(const Uintah::PatchSubset* patches,
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

    // Get the material index
    int matlIndex = matl->getDWIndex();
 
    // Get the particle subset for this material
    Uintah::ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle variables needed
    Uintah::constParticleVariable<Uintah::Matrix3> pDefGrad_new;
    new_dw->get(pDefGrad_new, d_label->pDefGradLabel_preReloc, pset);

    // Initialize the variables to be updated
    Uintah::ParticleVariable<Uintah::Matrix3> pStress_new, pPK1Stress_new;
    new_dw->allocateAndPut(pStress_new, d_label->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pPK1Stress_new, d_label->pPK1StressLabel_preReloc, pset);

    // Loop through particles
    for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
      
      // Get particle index
      Uintah::particleIndex idx = *iter;

      // Compute J = det(F)
      double J = pDefGrad_new[idx].Determinant();

      // Compute Bbar = J^{-2/3} (F F^T)  and dev(Bbar) = Bbar - 1/3 Tr(Bbar) I
      Uintah::Matrix3 Bbar = pDefGrad_new[idx]*(pDefGrad_new[idx].Transpose())*std::pow(J, -2.0/3.0);
      Uintah::Matrix3 BbarDev = Bbar - One*(Bbar.Trace()/3.0); 

      // Computes stress
      double pressure = -d_cm.bulkModulus*(J-1.0);
      pStress_new[idx] = One*pressure + BbarDev*(d_cm.shearModulus/J);

      // Compute PK1 stress
      pPK1Stress_new[idx] = pStress_new[idx]*(pDefGrad_new[idx].Transpose()*J);
    }

  } // end patches loop
}

// Register the permanent particle state associated with this material
void 
IsotropicElasticNeoHookeanStateModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                                       std::vector<const Uintah::VarLabel*>& to)
{
  // These are common to all models and will have to be moved up
  from.push_back(d_label->pDefGradLabel);
  to.push_back(d_label->pDefGradLabel_preReloc);

  from.push_back(d_label->pStressLabel);
  to.push_back(d_label->pStressLabel_preReloc);

  from.push_back(d_label->pPK1StressLabel);
  to.push_back(d_label->pPK1StressLabel_preReloc);
}


