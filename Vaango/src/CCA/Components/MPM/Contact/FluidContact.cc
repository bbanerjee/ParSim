/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/Contact/FluidContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>

#include <iostream>
#include <vector>

using namespace Uintah;
using std::string;
using std::vector;

FluidContact::FluidContact(const ProcessorGroup* myworld,
                           const MaterialManagerP& mat_manager,
                           const MPMLabel* labels,
                           const MPMFlags* flags,
                           ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  d_hydro_mpm_labels = std::make_unique<HydroMPMLabel>();
}

void
FluidContact::setContactMaterialAttributes()
{
}

void
FluidContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "fluid");
  d_matls.outputProblemSpec(contact_ps);
}

void
FluidContact::exchangeMomentum(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw,
                               const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  // Need to check whether we have null contact

  Ghost::GhostType gnone = Ghost::None;

  // Setting the fluid velocity equal to the rigid velocity,
  // to satisfy the no-slip boundary condition of impermeable wall-fluid contact
  // in computational fluid dynamics
  int numMatls = d_mat_manager->getNumMaterials();
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once, so store in
  // vectors of NCVariables
  constNCdoubleArray gMass(numMatls);
  constNCdoubleArray gVolume(numMatls);
  constNCVectorArray gVelocity(numMatls);
  constNCVectorArray gFluidVelocity(numMatls);
  constNCVectorArray gSurfNorm(numMatls);

  NCVectorArray gFluidVelocity_star(numMatls);
  NCVectorArray gFluidAcceleration(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    Vector dx       = patch->dCell();
    double cell_vol = dx.x() * dx.y() * dx.z();

    constNCVariable<double> NC_CCweight;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

    for (int m = 0; m < matls->size(); m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();

      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(
        gSurfNorm[m], d_mpm_labels->gSurfNormLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->get(gVelocity[m], gVelocity_label, dwi, patch, gnone, 0);
      new_dw->get(gFluidVelocity[m],
                  d_hydro_mpm_labels->gFluidVelocityLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      if (m != d_rigid_material) {
        new_dw->getModifiable(gFluidVelocity_star[m],
                              d_hydro_mpm_labels->gFluidVelocityStarLabel,
                              dwi,
                              patch);
        new_dw->getModifiable(gFluidAcceleration[m],
                              d_hydro_mpm_labels->gFluidAccelerationLabel,
                              dwi,
                              patch);
      }
    }

    // The normals are already computed if friction contact
    // It is stored in gSurfNorm

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector centerOfMassMom(0., 0., 0.);
      Point centerOfMassPos(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      for (int n = 0; n < numMatls; n++) {
        if (!d_matls.requested(n)) {
          continue;
        }
        centerOfMassMom += gVelocity[n][c] * gMass[n][c];
        centerOfMassMass += gMass[n][c];
        totalNodalVol += gVolume[n][c] * 8.0 * NC_CCweight[c];
      }
      centerOfMassPos /= centerOfMassMass;

      for (int n = 0; n < numMatls; n++) {
        // The rigid material does not change
        if (n == d_rigid_material) {
          continue;
        }

        // The rigid velocity
        Vector rigid_vel = gVelocity[d_rigid_material][c];

        if (!compare(gMass[d_rigid_material][c], 0.) &&
            (totalNodalVol / cell_vol) > d_vol_const) {
          gFluidVelocity_star[n][c] =
            gFluidVelocity[n][c] -
            Dot(gSurfNorm[n][c], gFluidVelocity[n][c] - rigid_vel) *
              gSurfNorm[n][c];
          gFluidAcceleration[n][c] =
            (gFluidVelocity_star[n][c] - gFluidVelocity[n][c]) / delT;
        }
      }
    } // nodeiterator
  }
}

void
FluidContact::addComputesAndRequires(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* ms,
                                     const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  Task* t = scinew Task("FluidContact::exchangeMomentum",
                        this,
                        &FluidContact::exchangeMomentum,
                        gVelocity_label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, gVelocity_label, Ghost::None);
  t->requires(
    Task::NewDW, d_hydro_mpm_labels->gFluidVelocityLabel, Ghost::None);
  t->modifies(d_hydro_mpm_labels->gFluidVelocityStarLabel, mss);
  t->modifies(d_hydro_mpm_labels->gFluidAccelerationLabel, mss);

  sched->addTask(t, patches, ms);

  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}
