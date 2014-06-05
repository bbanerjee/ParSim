#include <CCA/Components/Peridynamics/GradientComputer/PeridynamicsDefGradComputer.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/NeighborList.h>
#include <Core/Grid/Variables/NeighborConnectivity.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Util/DebugStream.h>

#include <iostream>

using namespace Vaango;

using Uintah::MaterialSubset;
using Uintah::Task;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::DataWarehouse;
using Uintah::ParticleSubset;
using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::Matrix3;
using Uintah::Ghost;
using Uintah::particleIndex;
using Uintah::long64;
using SCIRun::Vector;
using SCIRun::Point;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDDefGradDoing:+,PDDefGradDebug:+".....
//  bash     : export SCI_DEBUG="PDDefGradDoing:+,PDDefGradDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDDefGradDoing", false);
static DebugStream cout_dbg("PDDefGradDebug", false);


static const Matrix3 One(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

PeridynamicsDefGradComputer::PeridynamicsDefGradComputer(PeridynamicsFlags* flags,
                                                         PeridynamicsLabel* labels)
{
  d_flags = flags;
  d_labels = labels;
}

PeridynamicsDefGradComputer::~PeridynamicsDefGradComputer()
{
}

void
PeridynamicsDefGradComputer::addInitialComputesAndRequires(Task* task,
                                                           const PeridynamicsMaterial* matl,
                                                           const PatchSet*)
{
  cout_doing << "\t Scheduling initial compute variables in def grad computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  const MaterialSubset* matlset = matl->thisMaterial();
      
  // Computes (for explicit)
  task->computes(d_labels->pDefGradLabel,  matlset);
  task->computes(d_labels->pShapeTensorInvLabel, matlset);
}

void 
PeridynamicsDefGradComputer::initialize(const Patch* patch,
                                        PeridynamicsMaterial* matl,
                                        DataWarehouse* new_dw)
{
  cout_doing << "\t Initializing def grad computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  int matlIndex = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);
  ParticleVariable<Matrix3> pDefGrad;
  ParticleVariable<Matrix3> pShapeTensorInv;
  new_dw->allocateAndPut(pDefGrad,        d_labels->pDefGradLabel,        pset);
  new_dw->allocateAndPut(pShapeTensorInv, d_labels->pShapeTensorInvLabel, pset);
   
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++ ) {
    particleIndex idx = *iter;
    pDefGrad[idx] = One;
    pShapeTensorInv[idx] = One;
  }
}

void 
PeridynamicsDefGradComputer::addComputesAndRequires(Task* task,
                                                    const PeridynamicsMaterial* matl,
                                                    const PatchSet*)
{
  cout_doing << "\t Adding time stepping computes and requires in def grad computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  Ghost::GhostType gac = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  int numGhostCells = d_flags->d_numCellsInHorizon;

  // Get current body
  const MaterialSubset* matlset = matl->thisMaterial();

  // Requires 
  task->requires(Task::OldDW, d_labels->pPositionLabel,      matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_labels->pDisplacementLabel,  matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_labels->pVolumeLabel,        matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_labels->pParticleIDLabel,    matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_labels->pNeighborListLabel,  matlset, gnone);
  task->requires(Task::OldDW, d_labels->pNeighborConnLabel,  matlset, gnone);
  task->requires(Task::OldDW, d_labels->pNeighborCountLabel, matlset, gnone);

  // Computes 
  task->computes(d_labels->pDefGradLabel_preReloc,        matlset);
  task->computes(d_labels->pShapeTensorInvLabel_preReloc, matlset);
}


void 
PeridynamicsDefGradComputer::computeDeformationGradient(const Patch* patch,
                                                        const PeridynamicsMaterial* matl,
                                                        DataWarehouse* old_dw,
                                                        DataWarehouse* new_dw)
{
  cout_doing << "\t Computing deformation gradient in def grad computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  // F = sum[ w(xi) * OuterProduct(x,xi) * Vol ] * K_inv
  // K = sum[ w(xi) * OuterProduct(xi,xi) * Vol ]
  // w(xi) is influence function
  // xi = x'-x
  int matlIndex = matl->getDWIndex();

  // Get the set of particles completely contained in this patch
  ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

  // Get the set of particles contained in this patch + the ghost cells
  ParticleSubset* familySet = old_dw->getParticleSubset(matlIndex, patch, Ghost::AroundCells,
                                                        d_flags->d_numCellsInHorizon,
                                                        d_labels->pPositionLabel);

  // Create a map that takes particle IDs to the array index in the larger particle subset
  Uintah::ParticleIDMap familyIdMap;
  old_dw->createParticleIDMap(familySet, d_labels->pParticleIDLabel, familyIdMap);

  // TODO: Create factory for influence functions.  For now set to 1.
  const double pInfluence = 1.0;

  // Get the particle data for required variables
  constParticleVariable<Point> pPosition, pPosition_family;
  old_dw->get(pPosition,        d_labels->pPositionLabel, pset);
  old_dw->get(pPosition_family, d_labels->pPositionLabel, familySet);

  constParticleVariable<Vector> pDisp, pDisp_family;
  old_dw->get(pDisp,        d_labels->pDisplacementLabel, pset);
  old_dw->get(pDisp_family, d_labels->pDisplacementLabel, familySet);

  constParticleVariable<double> pVol_family;
  old_dw->get(pVol_family,  d_labels->pVolumeLabel, familySet);

  constParticleVariable<Uintah::NeighborList> pFamily;
  old_dw->get(pFamily,      d_labels->pNeighborListLabel,   pset);

  constParticleVariable<int> pFamilyCount;
  old_dw->get(pFamilyCount, d_labels->pNeighborCountLabel,  pset);

  constParticleVariable<Uintah::NeighborConnectivity> pFamilyConnected;
  old_dw->get(pFamilyConnected, d_labels->pNeighborConnLabel,  pset);

  // Allocate particle data for computed variables
  ParticleVariable<Matrix3> pDefGrad_new;
  new_dw->allocateAndPut(pDefGrad_new, d_labels->pDefGradLabel_preReloc, pset);

  ParticleVariable<Matrix3> pShapeTensorInv_new;
  new_dw->allocateAndPut(pShapeTensorInv_new, d_labels->pShapeTensorInvLabel_preReloc, pset);
   
  // Loop through the set of particles completely contained within this patch
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++ ) {

    particleIndex idx = *iter;

    Matrix3 defGrad_new(0.0);
    Matrix3 K(0.0);

    if (cout_dbg.active()) {
      cout_dbg << "\t\t Particle index = " << idx 
               << " current position = " << pPosition[idx] 
               << " displacement = " << pDisp[idx] 
               << std::endl;
    }

    for (int ii=0; ii < pFamilyCount[idx]; ii++) {

      // If the bond exists
      if (pFamilyConnected[idx][ii]) {

        // Get Particle Index from Neighbor List
        long64 pID_family = pFamily[idx][ii];
        particleIndex family_idx;
        old_dw->getParticleIndex(familyIdMap, pID_family, family_idx);

        if (cout_dbg.active()) {
          cout_dbg << "\t\t\t Family particle index = " << family_idx 
                   << " current position = " << pPosition_family[family_idx] 
                   << " displacement = " << pDisp_family[family_idx] 
                   << std::endl;
        }

        // Create Sums
        Vector x  = pPosition_family[family_idx] - pPosition[idx];
        Vector xi = x - (pDisp_family[family_idx] - pDisp[idx]);

        if (cout_dbg.active()) {
          cout_dbg << "\t\t\t Y = (yhat - y) = " << x
                   << " , xi = (xhat - x) =  " << xi << std::endl;
          cout_dbg << "\t\t\t xi otimes xi = " << Matrix3(xi, xi) << std::endl;
        }

        K += pInfluence * pVol_family[family_idx] * Matrix3(xi,xi);
        defGrad_new += pInfluence * pVol_family[family_idx] * Matrix3(x,xi); 
      }

    }

    pShapeTensorInv_new[idx] = K.Inverse();
    pDefGrad_new[idx] = defGrad_new * pShapeTensorInv_new[idx];

    if (cout_dbg.active()) {
      cout_dbg << "\t\t Particle index = " << idx 
               << " Def Grad = " << pDefGrad_new[idx] 
	       << " Shape tensor inverse = " << pShapeTensorInv_new[idx] << std::endl;
    }
  }
}

void
PeridynamicsDefGradComputer::addParticleState(std::vector<const Uintah::VarLabel*>& ,
                                              std::vector<const Uintah::VarLabel*>& )
{
  std::cout << "Deformation gradient computer::addParticleState called by mistake." << std::endl;
}

