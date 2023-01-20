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

// SingleVel.cc
// One of the derived Contact classes.  This particular
// class contains methods for recapturing single velocity
// field behavior from objects belonging to multiple velocity
// fields.  The main purpose of this type of contact is to
// ensure that one can get the same answer using prescribed
// contact as can be gotten using "automatic" contact.
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/SingleVelContact.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <vector>

using namespace Uintah;
using std::vector;

SingleVelContact::SingleVelContact(const ProcessorGroup* myworld,
                                   ProblemSpecP& ps, MaterialManagerP& d_sS,
                                   MPMLabel* Mlb, MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  // Constructor
  d_mat_manager = d_sS;
  lb = Mlb;
  flag = MFlag;
}

SingleVelContact::~SingleVelContact()
{
}

void
SingleVelContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "single_velocity");
  d_matls.outputProblemSpec(contact_ps);
}

void
SingleVelContact::exchangeMomentum(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw, DataWarehouse* new_dw,
                                   const VarLabel* gVelocity_label)
{
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
    std::vector<constNCVariable<double>> gmass(numMatls);
    std::vector<NCVariable<Vector>> gvelocity_star(numMatls);

    for (int m = 0; m < matls->size(); m++) {
      int dwi = matls->get(m);
      new_dw->get(gmass[m], lb->gMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->getModifiable(gvelocity_star[m], gVelocity_label, dwi, patch);
    }

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      centerOfMassMom = zero;
      centerOfMassMass = 0.0;
      for (int n = 0; n < numMatls; n++) {
        if (d_matls.requested(n)) {
          centerOfMassMom += gvelocity_star[n][c] * gmass[n][c];
          centerOfMassMass += gmass[n][c];
        }
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity = centerOfMassMom / centerOfMassMass;
      for (int n = 0; n < numMatls; n++) {
        if (d_matls.requested(n)) {
          Dvdt = (centerOfMassVelocity - gvelocity_star[n][c]) / delT;
          gvelocity_star[n][c] = centerOfMassVelocity;
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
  Task* t = scinew Task("SingleVelContact::exchangeMomentum", this,
                        &SingleVelContact::exchangeMomentum,
                        gVelocity_label);

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);
  t->requires(Task::NewDW, lb->gMassLabel, Ghost::None);

  t->modifies(gVelocity_label, mss);

  sched->addTask(t, patches, matls);
}
