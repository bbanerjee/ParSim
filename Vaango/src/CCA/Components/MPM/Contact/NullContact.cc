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

// NullContact.cc
// One of the derived Contact classes.  This particular
// class is used when no contact is desired.  This would
// be used for example when a single velocity field is
// present in the problem, so doing contact wouldn't make
// sense.
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/NullContact.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Labels/MPMLabel.h>

using namespace Uintah;

NullContact::NullContact(const ProcessorGroup* myworld, SimulationStateP& d_sS,
                         MPMLabel* Mlb, MPMFlags* MFlags)
  : Contact(myworld, Mlb, MFlags, 0)
{
  // Constructor
  d_sharedState = d_sS;
  lb = Mlb;
  flag = MFlags;
}

NullContact::~NullContact()
{
}

void
NullContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "null");
  d_matls.outputProblemSpec(contact_ps);
}

void
NullContact::exchangeMomentum(const ProcessorGroup*, const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw, DataWarehouse* new_dw,
                              const VarLabel* gVelocity_label)
{
}

void
NullContact::addComputesAndRequires(SchedulerP& sched, const PatchSet* patches,
                                    const MaterialSet* matls,
                                    const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("NullContact::exchangeMomentum", this,
                        &NullContact::exchangeMomentum,
                        gVelocity_label);

  sched->addTask(t, patches, matls);
}
