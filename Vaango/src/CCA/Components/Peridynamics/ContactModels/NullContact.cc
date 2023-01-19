/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/ContactModels/NullContact.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>

#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>


using namespace Vaango;

using Uintah::ProcessorGroup;
using Uintah::ProblemSpecP;
using Uintah::SchedulerP;
using Uintah::PatchSubset;
using Uintah::MaterialSubset;
using Uintah::PatchSet;
using Uintah::MaterialSet;
using Uintah::Task;
using Uintah::DataWarehouse;
using Uintah::MaterialManagerP;

NullContact::NullContact(const ProcessorGroup* myworld,
                         MaterialManagerP& ss,
                         PeridynamicsLabel* labels,
                         PeridynamicsFlags* flags)
  : ContactModelBase(myworld, labels, flags, 0)
{
  d_mat_manager = ss;
}

NullContact::~NullContact()
{
}

void 
NullContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("ContactModel");
  contact_ps->appendElement("type","null");
  d_bodiesThatCanInteract.outputProblemSpec(contact_ps);
}

void 
NullContact::exchangeMomentumInterpolated(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset* matls,
                                          DataWarehouse* /*old_dw*/,
                                          DataWarehouse* new_dw)
{
}

void 
NullContact::exchangeMomentumIntegrated(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* matls,
                                        DataWarehouse* /*old_dw*/,
                                        DataWarehouse* new_dw)
{
}

void 
NullContact::addComputesAndRequiresInterpolated(SchedulerP & sched,
                                                const PatchSet* patches,
                                                const MaterialSet* ms)
{
  Task * t = scinew Task("NullContact::exchangeMomentumInterpolated", this, 
                          &NullContact::exchangeMomentumInterpolated);
  
  sched->addTask(t, patches, ms);
}

void 
NullContact::addComputesAndRequiresIntegrated(SchedulerP & sched,
                                              const PatchSet* patches,
                                              const MaterialSet* ms) 
{
  Task * t = scinew Task("NullContact::exchangeMomentumIntegrated", this, 
                         &NullContact::exchangeMomentumIntegrated);
  
  sched->addTask(t, patches, ms);
}
