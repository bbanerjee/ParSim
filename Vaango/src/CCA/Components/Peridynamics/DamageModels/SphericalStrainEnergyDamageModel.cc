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
using Uintah::particleIndex;
using Uintah::long64;
using Uintah::ParticleID;
using SCIRun::Vector;
using SCIRun::Point;

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
  int numGhostCells = d_flags->d_numCellsInHorizon;
  
  // Get the current material
  const MaterialSubset* matlset = matl->thisMaterial();
  
  // List the variables needed for this task to execute
  task->requires(Task::OldDW, d_label->pDisplacementLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pParticleIDLabel,         matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pHorizonLabel,            matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborListLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborConnLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborCountLabel,      matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborBondEnergyLabel, matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pPositionLabel_preReloc,       matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pDisplacementLabel_preReloc,   matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pPK1StressLabel_preReloc,      matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pShapeTensorInvLabel_preReloc, matlset, gac, numGhostCells);

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
  // Assume influence function is always 1 (**TODO** Compute it using a factory.)
  double pInfluence = 1.0;

  // Loop through patches
  for (int pp = 0; pp < patches->size(); pp++) {

    // Get the current patch
    const Patch* patch = patches->get(pp);

    // Get the index of current material in the data warehouse
    int matlIndex = matl->getDWIndex();

    // Get the particle subset for this material in this patch
    ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle subset for the neighbor particles of the same material in this patch + ghost regions
    ParticleSubset* familySet = old_dw->getParticleSubset(matlIndex, patch, Ghost::AroundCells,
                                                          d_flags->d_numCellsInHorizon,
                                                          d_label->pPositionLabel);

    // Create a map that takes particle IDs to the array index in the larger particle subset
    Uintah::ParticleIDMap familyIdMap;
    old_dw->createParticleIDMap(familySet, d_label->pParticleIDLabel, familyIdMap);

    // Get the particle data needed for damage computation
    constParticleVariable<Point> pPosition_new, pPosition_new_family;
    new_dw->get(pPosition_new,        d_label->pPositionLabel_preReloc, pset);
    new_dw->get(pPosition_new_family, d_label->pPositionLabel_preReloc, familySet);

    constParticleVariable<Vector> pDisplacement_old, pDisplacement_old_family;
    old_dw->get(pDisplacement_old,        d_label->pDisplacementLabel, pset);
    old_dw->get(pDisplacement_old_family, d_label->pDisplacementLabel, familySet);

    constParticleVariable<Vector> pDisplacement_new, pDisplacement_new_family;
    new_dw->get(pDisplacement_new,        d_label->pDisplacementLabel_preReloc, pset);
    new_dw->get(pDisplacement_new_family, d_label->pDisplacementLabel_preReloc, familySet);

    constParticleVariable<Matrix3> pPK1Stress_new, pPK1Stress_new_family;
    new_dw->get(pPK1Stress_new,        d_label->pPK1StressLabel_preReloc, pset);
    new_dw->get(pPK1Stress_new_family, d_label->pPK1StressLabel_preReloc, familySet);

    constParticleVariable<Matrix3> pShapeInv_new, pShapeInv_new_family;
    new_dw->get(pShapeInv_new,        d_label->pShapeTensorInvLabel_preReloc, pset);
    new_dw->get(pShapeInv_new_family, d_label->pShapeTensorInvLabel_preReloc, familySet);

    constParticleVariable<double> pHorizon;
    old_dw->get(pHorizon, d_label->pHorizonLabel, pset);

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

      // Find the position and displacement of the current particle
      Point cur_pos = pPosition_new[idx];
      Vector cur_disp_old = pDisplacement_old[idx];
      Vector cur_disp = pDisplacement_new[idx];
      Point cur_ref_pos = cur_pos - cur_disp;

      // Find the first PK stress and shape tensor inverse
      Matrix3 cur_PK1_stress = pPK1Stress_new[idx];
      Matrix3 cur_shape_inv = pShapeInv_new[idx];

      // Compute the internal force state
      Matrix3 cur_force_state = (cur_PK1_stress*cur_shape_inv)*pInfluence;

      // Compute the critical bond energy
      double hh = pHorizon[idx];
      double critical_bond_energy = 4*d_GIc/(M_PI*hh*hh*hh*hh);

      // Get the neighbor data
      NeighborList family = pNeighborList[idx];
      NeighborConnectivity connected = pNeighborConn[idx];
      NeighborBondEnergy bondEnergy_old = pNeighborBondEnergy[idx];

      // Loop through the neighbor list
      int neighborCount = pNeighborCount[idx];
      for (int ii = 0; ii < neighborCount; ii++) {
       
        // If the bond exists
        if (connected[ii]) {

          // Find the idx of the neighbor
          ParticleID familyID = family[ii];
          particleIndex family_idx;
          old_dw->getParticleIndex(familyIdMap, familyID, family_idx);

          // Find the position and displacement of the neighbor
          Point family_pos = pPosition_new_family[family_idx];
          Vector family_disp_old = pDisplacement_old_family[family_idx];
          Vector family_disp = pDisplacement_new_family[family_idx];
          Point family_ref_pos = family_pos - family_disp;

          // Compute relative reference position xi = x_family - x
          Vector xi = family_ref_pos - cur_ref_pos;

          // Compute relative displacement
          Vector eta_new = family_disp - cur_disp;

          // Compute increment of relative displacement
          Vector eta_old = family_disp_old - cur_disp_old;
          Vector eta_inc = eta_new - eta_old;

          // Find the first PK stress and shape tensor inverse
          Matrix3 family_PK1_stress = pPK1Stress_new_family[family_idx];
          Matrix3 family_shape_inv = pShapeInv_new_family[family_idx];

          // Compute the family internal force state
          Matrix3 family_force_state = (family_PK1_stress*family_shape_inv)*pInfluence;

          // Compute internal force difference
          Vector internal_force_diff = cur_force_state*xi - family_force_state*(-xi);

          // Compute the increment of energy
          double energy_inc = SCIRun::Dot(internal_force_diff, eta_inc);

          // Update the bond energy
          double bond_energy_new = bondEnergy_old[ii] + energy_inc;
          (pNeighborBondEnergy_new[idx])[ii] = bond_energy_new;

          // Compare the bond energy with the critical energy and break bonds
          if (bond_energy_new > critical_bond_energy) {
            (pNeighborConn_new[idx])[ii] = false;    // broken
          } else {
            (pNeighborConn_new[idx])[ii] = true;   // connected
          }

        } // Endif bond connected/broken
      }  // End neighbor particle loop

      // Now count the unbroken bonds
      int intact_bond_count = 0;
      for (int ii = 0; ii < neighborCount; ii++) {
        if ((pNeighborConn_new[idx])[ii]) {
          intact_bond_count++;
        }
      }
      double damaged_bond_ratio = 1.0 - (double)intact_bond_count/(double)neighborCount;
 
      // Update damage
      Matrix3 damageTensor(damaged_bond_ratio, 0.0, 0.0, 0.0, damaged_bond_ratio, 0.0, 
                           0.0, 0.0, damaged_bond_ratio); 
      pDamage_new[idx] = damageTensor;

    }  // End particle loop
  }  // End patch loop
}


void 
SphericalStrainEnergyDamageModel::addParticleState(std::vector<const Uintah::VarLabel*>& ,
                                                   std::vector<const Uintah::VarLabel*>& )
{
}
