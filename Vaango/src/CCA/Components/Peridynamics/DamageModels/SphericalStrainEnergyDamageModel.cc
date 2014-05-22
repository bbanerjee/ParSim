#include <CCA/Components/Peridynamics/DamageModels/SphericalStrainEnergyDamageModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <cmath>
#include <iostream>

using namespace Vaango;

using Uintah::ProblemSpecP;
using Uintah::MaterialSubset;
using Uintah::Task;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::DataWarehouse;
using Uintah::ParticleSubset;
using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::Matrix3;
using Uintah::NeighborList;
using Uintah::NeighborConnectivity;
using Uintah::NeighborBondEnergy;
using Uintah::Ghost;
using SCIRun::Point;
using Uintah::particleIndex;
using Uintah::long64;

SphericalStrainEnergyDamageModel::SphericalStrainEnergyDamageModel(ProblemSpecP& ps,
                                                                   PeridynamicsLabel* labels,
                                                                   PeridynamicsFlags* flags)
  : PeridynamicsDamageModel(labels, flags)
{

  // Get the mode I critical strain energy release rate
  ps->require("G_Ic", d_GIc);
}

SphericalStrainEnergyDamageModel::SphericalStrainEnergyDamageModel(const SphericalStrainEnergyDamageModel* cm)
  : PeridynamicsDamageModel(cm)
{
  d_GIc = cm->d_GIc;
}

SphericalStrainEnergyDamageModel* 
SphericalStrainEnergyDamageModel::clone()
{
  return scinew SphericalStrainEnergyDamageModel(*this);
}

SphericalStrainEnergyDamageModel::~SphericalStrainEnergyDamageModel()
{
}

void 
SphericalStrainEnergyDamageModel::outputProblemSpec(ProblemSpecP& ps,
                                                    bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("damage_model");
    cm_ps->setAttribute("type", "spherical_strain_energy");
  }
  cm_ps->appendElement("G_Ic", d_GIc);
}

void 
SphericalStrainEnergyDamageModel::addInitialComputesAndRequires(Task* task ,
                                                                const PeridynamicsMaterial* material,
                                                                const PatchSet*) const
{
  // Get material
  const MaterialSubset* matl = material->thisMaterial();

  // Initial task will initialize the damage tensor
  task->computes(d_label->pDamageLabel, matl);
}

void 
SphericalStrainEnergyDamageModel::initialize(const Patch* patch,
                                             const PeridynamicsMaterial* matl,
                                             DataWarehouse* new_dw)
{
 // Get the set of particles of this material type in the current patch  
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  
  // Allocate for saving
  ParticleVariable<Matrix3> pDamage;
  new_dw->allocateAndPut(pDamage, d_label->pDamageLabel, pset);
  
  // Initialize the damage to zero (for now)
  Matrix3 zero(0.0);
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {
    pDamage[*iter] = zero;
  }
}

void 
SphericalStrainEnergyDamageModel::addComputesAndRequires(Task* task, 
                                                         const PeridynamicsMaterial* matl,
                                                         const PatchSet* patches) const
{
  // Constants
  Ghost::GhostType gac = Ghost::AroundCells;
  
  // Get the current material
  const MaterialSubset* matlset = matl->thisMaterial();
  
  // List the variables needed for this task to execute
  int numGhostCells = 3;
  task->requires(Task::OldDW, d_label->pPositionLabel,           matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pParticleIDLabel,         matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pHorizonLabel,            matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborListLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborConnLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborCountLabel,      matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborBondEnergyLabel, matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pStressLabel_preReloc,    matlset, gac, numGhostCells);

  // List the variables computed by this task
  task->computes(d_label->pNeighborBondEnergyLabel_preReloc, matlset);
  task->computes(d_label->pNeighborConnLabel_preReloc,       matlset);
  task->computes(d_label->pDamageLabel_preReloc,             matlset);
}

void 
SphericalStrainEnergyDamageModel::computeDamageTensor(const PatchSubset* patches,
                                                      const PeridynamicsMaterial* matl,
                                                      DataWarehouse* old_dw,
                                                      DataWarehouse* new_dw)
{

  // Loop through patches
  for (int pp = 0; pp < patches->size(); pp++) {

    // Get the current patch
    const Patch* patch = patches->get(pp);

    // Get the index of current material in the data warehouse
    int matlIndex = matl->getDWIndex();

    // Get the particle subset for this material in this patch
    ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle data needed for damage computation
    constParticleVariable<Point> pPosition;
    old_dw->get(pPosition, d_label->pPositionLabel, pset);

    constParticleVariable<long64> pParticleID;
    old_dw->get(pParticleID, d_label->pParticleIDLabel, pset);

    constParticleVariable<double> pHorizon;
    old_dw->get(pHorizon, d_label->pHorizonLabel, pset);

    constParticleVariable<Matrix3> pStress_new;
    new_dw->get(pStress_new, d_label->pStressLabel_preReloc, pset);

    constParticleVariable<int> pNeighborCount;
    old_dw->get(pNeighborCount, d_label->pNeighborCountLabel, pset);

    constParticleVariable<NeighborList> pNeighborList;
    old_dw->get(pNeighborList, d_label->pNeighborListLabel, pset);

    constParticleVariable<NeighborConnectivity> pNeighborConn;
    old_dw->get(pNeighborConn, d_label->pNeighborConnLabel, pset);

    constParticleVariable<NeighborBondEnergy> pNeighborBondEnergy;
    old_dw->get(pNeighborBondEnergy, d_label->pNeighborBondEnergyLabel, pset);

    // Initialize variables that will be updated in this task
    ParticleVariable<NeighborConnectivity> pNeighborConn_new;
    new_dw->allocateAndPut(pNeighborConn_new, d_label->pNeighborConnLabel_preReloc, pset);

    ParticleVariable<NeighborBondEnergy> pNeighborBondEnergy_new;
    new_dw->allocateAndPut(pNeighborBondEnergy_new, d_label->pNeighborBondEnergyLabel_preReloc, pset);

    ParticleVariable<Matrix3> pDamage_new;;
    new_dw->allocateAndPut(pDamage_new, d_label->pDamageLabel_preReloc, pset);

    // Loop through particles
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

      // Get particle index
      particleIndex idx = *iter;
    }
  }  // End patch loop
}


void 
SphericalStrainEnergyDamageModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                                   std::vector<const Uintah::VarLabel*>& to)
{
}
