/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

// SingleVel.cc
// One of the derived Contact classes.  This particular
// class contains methods for recapturing single velocity
// field behavior from objects belonging to multiple velocity
// fields.  The main purpose of this type of contact is to
// ensure that one can get the same answer using prescribed
// contact as can be gotten using "automatic" contact.

#include <CCA/Components/MPM/Contact/SingleVelContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace Uintah;
using std::vector;

SingleVelContact::SingleVelContact(const ProcessorGroup* myworld,
                                   const MaterialManagerP& mat_manager,
                                   const MPMLabel* labels,
                                   const MPMFlags* flags,
                                   ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  d_one_or_two_step = 2;

  ps->get("one_or_two_step", d_one_or_two_step);
  ps->getWithDefault("exclude_material", d_exclude_material, -999);
}

void
SingleVelContact::setContactMaterialAttributes()
{
}

void
SingleVelContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "single_velocity");
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  contact_ps->appendElement("exclude_material", d_exclude_material);
  d_matls.outputProblemSpec(contact_ps);
}

void
SingleVelContact::exchangeMomentum(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   const VarLabel* gVelocity_label)
{
  // If one_step, only do exchange for gVelocity_star
  if (d_one_or_two_step == 1 &&
      gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector zero(0.0, 0.0, 0.0);
    Vector centerOfMassVelocity(0.0, 0.0, 0.0);
    Vector centerOfMassMom(0.0, 0.0, 0.0);
    Vector Dvdt;
    double centerOfMassMass;

    // Retrieve necessary data from DataWarehouse
    std::vector<constNCVariable<double>> gMass(numMatls);
    std::vector<NCVariable<Vector>> gVelocity_star(numMatls);

    for (int m = 0; m < matls->size(); m++) {
      int matID = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, matID, patch, Ghost::None, 0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, matID, patch);
    }

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      centerOfMassMom  = zero;
      centerOfMassMass = 0.0;
      for (int n = 0; n < numMatls; n++) {
        if (d_matls.requested(n)) {
          centerOfMassMom += gVelocity_star[n][c] * gMass[n][c];
          centerOfMassMass += gMass[n][c];
        }
      }

      double excludeMass = 0.;
      if(d_exclude_material >=0){
        excludeMass = gMass[d_exclude_material][c];
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity = centerOfMassMom / centerOfMassMass;
      for (int n = 0; n < numMatls; n++) {
        if (d_matls.requested(n) && excludeMass < 1.0e-99) {
          Dvdt = (centerOfMassVelocity - gVelocity_star[n][c]) / delT;
          gVelocity_star[n][c] = centerOfMassVelocity;
        }
      }
    }
  }
}

void
SingleVelContact::addComputesAndRequires(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls,
                                         const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("SingleVelContact::exchangeMomentum",
                        this,
                        &SingleVelContact::exchangeMomentum,
                        gVelocity_label);

  const MaterialSubset* mss = matls->getUnion();
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);
  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);

  t->modifies(gVelocity_label, mss);

  sched->addTask(t, patches, matls);
}
