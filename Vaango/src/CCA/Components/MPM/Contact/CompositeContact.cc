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

#include <CCA/Components/MPM/Contact/CompositeContact.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

using namespace Uintah;

CompositeContact::CompositeContact(const ProcessorGroup* myworld,
                                   const MPMLabel* labels,
                                   const MPMFlags* flags,
                                   ProblemSpecP& ps)
  : Contact(myworld, nullptr, labels, flags, ps)
{
}

void
CompositeContact::setContactMaterialAttributes()
{
}

void
CompositeContact::outputProblemSpec(ProblemSpecP& ps)
{
  for (auto& contactModel : d_m) {
    contactModel->outputProblemSpec(ps);
  }
}

void
CompositeContact::add(std::unique_ptr<Contact> m)
{
  if (m->needNormals()) {
    d_need_normals = true;
  }
  if (m->useLogisticRegression()) {
    d_use_logistic_regression = true;
  }
  if (m->oneOrTwoStep() == 1) {
    d_one_or_two_step = 1;
  } else {
    d_one_or_two_step = 2;
  }

  d_m.push_back(std::move(m));
}

void
CompositeContact::exchangeMomentum(const ProcessorGroup* pg,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   const VarLabel* gVelocity_label)
{
  for (auto& contactModel : d_m) {
    contactModel->exchangeMomentum(
      pg, patches, matls, old_dw, new_dw, gVelocity_label);
  }
}

void
CompositeContact::addComputesAndRequires(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls,
                                         const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    Task* t = scinew Task(
      "Contact::initFriction", this, &CompositeContact::initFriction);
    t->computes(d_mpm_labels->frictionalWorkLabel);
    sched->addTask(t, patches, matls);
  }

  for (auto& contactModel : d_m) {
    contactModel->addComputesAndRequires(
      sched, patches, matls, gVelocity_label);
  }
}

void
CompositeContact::initFriction(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      NCVariable<double> frictionWork_m;
      new_dw->allocateAndPut(frictionWork_m,
                             d_mpm_labels->frictionalWorkLabel,
                             matls->get(m),
                             patch);
      frictionWork_m.initialize(0.);
    }
  }
}
