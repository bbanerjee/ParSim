/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardeningModel.h>

using namespace Uintah;
using namespace Vaango;

KinematicHardeningModel::KinematicHardeningModel()
{
  pBackStressLabel = VarLabel::create(
    "p.backStress", ParticleVariable<Matrix3>::getTypeDescription());
  pBackStressLabel_preReloc = VarLabel::create(
    "p.backStress+", ParticleVariable<Matrix3>::getTypeDescription());
}

KinematicHardeningModel::~KinematicHardeningModel()
{
  VarLabel::destroy(pBackStressLabel);
  VarLabel::destroy(pBackStressLabel_preReloc);
}

void
KinematicHardeningModel::addInitialComputesAndRequires(
  Task* task, const MPMMaterial* matl, const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pBackStressLabel, matlset);
}

void
KinematicHardeningModel::addComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pBackStressLabel, matlset, Ghost::None);
  task->computes(pBackStressLabel_preReloc, matlset);
}

void
KinematicHardeningModel::addComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patches,
                                                bool recurse) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::ParentOldDW, pBackStressLabel, matlset, Ghost::None);
}

void
KinematicHardeningModel::allocateCMDataAddRequires(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet* patch,
                                                   MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pBackStressLabel_preReloc, matlset, Ghost::None);
}

void
KinematicHardeningModel::allocateCMDataAdd(
  DataWarehouse* new_dw, ParticleSubset* addset,
  ParticleLabelVariableMap* newState, ParticleSubset* delset,
  DataWarehouse* old_dw)
{
  ParticleVariable<double> pBackStress_upd;
  constParticleVariable<double> pBackStress_old;

  new_dw->allocateTemporary(pBackStress_upd, addset);

  new_dw->get(pBackStress_old, pBackStressLabel_preReloc, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pBackStress_upd[*n] = pBackStress_old[*o];
  }

  (*newState)[pBackStressLabel] = pBackStress_upd.clone();
}

void
KinematicHardeningModel::addParticleState(std::vector<const VarLabel*>& from,
                                          std::vector<const VarLabel*>& to)
{
  from.push_back(pBackStressLabel);
  to.push_back(pBackStressLabel_preReloc);
}

void
KinematicHardeningModel::initializeBackStress(ParticleSubset* pset,
                                              DataWarehouse* new_dw)
{
  ParticleVariable<Matrix3> pBackStress_new;
  new_dw->allocateAndPut(pBackStress_new, pBackStressLabel, pset);
  for (int& iter : *pset) {
    pBackStress_new[iter] = 0.0;
  }
}

void
KinematicHardeningModel::getBackStress(
  ParticleSubset* pset, DataWarehouse* old_dw,
  constParticleVariable<Matrix3>& pBackStress)
{
  old_dw->get(pBackStress, pBackStressLabel, pset);
}

void
KinematicHardeningModel::allocateAndPutBackStress(
  ParticleSubset* pset, DataWarehouse* new_dw,
  ParticleVariable<Matrix3>& pBackStress_new)
{
  new_dw->allocateAndPut(pBackStress_new, pBackStressLabel_preReloc, pset);
}

void
KinematicHardeningModel::allocateAndPutRigid(ParticleSubset* pset,
                                             DataWarehouse* new_dw)
{
  ParticleVariable<Matrix3> pBackStress_new;
  new_dw->allocateAndPut(pBackStress_new, pBackStressLabel_preReloc, pset);
  for (int& iter : *pset) {
    pBackStress_new[iter] = 0.0;
  }
}
