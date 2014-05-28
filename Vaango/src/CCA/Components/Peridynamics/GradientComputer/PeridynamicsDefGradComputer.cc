#include <CCA/Components/Peridynamics/GradientComputer/PeridynamicsDefGradComputer.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/NeighborList.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <iostream>

using namespace Vaango;

static const Uintah::Matrix3 One(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

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
PeridynamicsDefGradComputer::addInitialComputesAndRequires(Uintah::Task* task,
                                                           const PeridynamicsMaterial* matl,
                                                           const Uintah::PatchSet*)
{
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();
      
  // Computes (for explicit)
  task->computes(d_labels->pDefGradLabel,  matlset);
  task->computes(d_labels->pShapeTensorInvLabel, matlset);
}

void 
PeridynamicsDefGradComputer::addComputesAndRequires(Uintah::Task* task,
                                                    const PeridynamicsMaterial* matl,
                                                    const Uintah::PatchSet*)
{
  Uintah::Ghost::GhostType gnone = Uintah::Ghost::None;
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();

  // Requires (for explicit)
  task->requires(Uintah::Task::OldDW, d_labels->pPositionLabel,      matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_labels->pDisplacementLabel,  matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_labels->pVolumeLabel,        matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_labels->pNeighborListLabel,  matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_labels->pNeighborCountLabel, matlset, gnone);
  task->requires(Uintah::Task::OldDW, d_labels->pParticleIDLabel,    matlset, gnone);

  // Computes (for explicit)
  task->computes(d_labels->pDefGradLabel_preReloc,      matlset);
  task->computes(d_labels->pShapeTensorInvLabel_preReloc, matlset);
}

void 
PeridynamicsDefGradComputer::initialize(const Uintah::Patch* patch,
                                        PeridynamicsMaterial* matl,
                                        Uintah::DataWarehouse* new_dw)
{
  int dwi = matl->getDWIndex();
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  Uintah::ParticleVariable<Uintah::Matrix3> pDefGrad;
  Uintah::ParticleVariable<Uintah::Matrix3> pShapeTensorInv;
  new_dw->allocateAndPut(pDefGrad,        d_labels->pDefGradLabel,      pset);
  new_dw->allocateAndPut(pShapeTensorInv, d_labels->pShapeTensorInvLabel, pset);
   
  Uintah::ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++ ) {
    Uintah::particleIndex idx = *iter;
    pDefGrad[idx] = One;
    pShapeTensorInv[idx] = One;
  }
}


void 
PeridynamicsDefGradComputer::computeDeformationGradient(const Uintah::Patch* patch,
                                                        const PeridynamicsMaterial* matl,
                                                        Uintah::DataWarehouse* old_dw,
                                                        Uintah::DataWarehouse* new_dw)
{
  // F = sum[ w(xi) * OuterProduct(x,xi) * Vol ] * K_inv
  // K = sum[ w(xi) * OuterProduct(xi,xi) * Vol ]
  // w(xi) is influence function
  // xi = x'-x
  int dwi = matl->getDWIndex();
  Uintah::ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

  const double pInfluence = 1.0;
  // TODO: Create factory for influence functions

  Uintah::constParticleVariable<SCIRun::Point> px;
  Uintah::constParticleVariable<SCIRun::Vector> pDisp;
  Uintah::constParticleVariable<double> pVol;
  Uintah::constParticleVariable<Uintah::NeighborList> pFamily;
  Uintah::constParticleVariable<int> pFamilyCount;
  Uintah::ParticleIDMap idMap;
  old_dw->get(px,           d_labels->pPositionLabel,       pset);
  old_dw->get(pDisp,        d_labels->pDisplacementLabel,   pset);
  old_dw->get(pVol,         d_labels->pVolumeLabel,         pset);
  old_dw->get(pFamily,      d_labels->pNeighborListLabel,   pset);
  old_dw->get(pFamilyCount, d_labels->pNeighborCountLabel,  pset);
  old_dw->createParticleIDMap(pset, d_labels->pParticleIDLabel, idMap);
   
  Uintah::ParticleVariable<Uintah::Matrix3> pDefGrad;
  Uintah::ParticleVariable<Uintah::Matrix3> pShapeTensorInv;
  new_dw->allocateAndPut(pDefGrad,      d_labels->pDefGradLabel_preReloc,      pset);
  new_dw->allocateAndPut(pShapeTensorInv, d_labels->pShapeTensorInvLabel_preReloc, pset);
   
  Uintah::ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++ ) {
    Uintah::particleIndex idx0 = *iter;
    Uintah::Matrix3 defGrad_new(0.0);
    Uintah::Matrix3 K(0.0);

    for (int ii=0; ii<pFamilyCount[idx0]; ii++) {
      // Get Particle Index from Neighbor List
      Uintah::long64 pID1 = pFamily[idx0][ii];
      Uintah::particleIndex idx1;
      old_dw->getParticleIndex(idMap, pID1, idx1);

      // Create Sums
      SCIRun::Vector x  = SCIRun::Vector(px[idx1]) - SCIRun::Vector(px[idx0]);
      SCIRun::Vector xi = x - (pDisp[idx1] - pDisp[idx0]);
      K += pInfluence * pVol[idx1] * Uintah::Matrix3(xi,xi);
      defGrad_new += pInfluence * pVol[idx1] * Uintah::Matrix3(x,xi); 
    }
    pShapeTensorInv[idx0] = K.Inverse();
    pDefGrad[idx0] = defGrad_new * pShapeTensorInv[idx0];
  }
}
