/*
 * The MIT License
 *
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

#include <CCA/Components/Peridynamics/ContactModels/NullContact.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>

using namespace Vaango;

using Uintah::DataWarehouse;
using Uintah::MaterialManagerP;
using Uintah::MaterialSet;
using Uintah::MaterialSubset;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::ProblemSpecP;
using Uintah::ProcessorGroup;
using Uintah::SchedulerP;
using Uintah::Task;

NullContact::NullContact(const ProcessorGroup* myworld,
                         const MaterialManagerP& mat_manager,
                         PeridynamicsLabel* labels,
                         PeridynamicsFlags* flags)
  : ContactModelBase(myworld, mat_manager, labels, flags, nullptr)
{
}

NullContact::~NullContact() {}

void
NullContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("ContactModel");
  contact_ps->appendElement("type", "null");
  d_bodiesThatCanInteract.outputProblemSpec(contact_ps);
}

void
NullContact::exchangeMomentumInterpolated([[maybe_unused]] const ProcessorGroup*,
                                          [[maybe_unused]] const PatchSubset* patches,
                                          [[maybe_unused]] const MaterialSubset* matls,
                                          [[maybe_unused]] DataWarehouse* old_dw,
                                          [[maybe_unused]] DataWarehouse* new_dw)
{
}

void
NullContact::exchangeMomentumIntegrated([[maybe_unused]] const ProcessorGroup*,
                                        [[maybe_unused]] const PatchSubset* patches,
                                        [[maybe_unused]] const MaterialSubset* matls,
                                        [[maybe_unused]] DataWarehouse* old_dw,
                                        [[maybe_unused]] DataWarehouse* new_dw)
{
}

void
NullContact::addComputesAndRequiresInterpolated(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* ms)
{
  Task* t = scinew Task("NullContact::exchangeMomentumInterpolated",
                        this,
                        &NullContact::exchangeMomentumInterpolated);

  sched->addTask(t, patches, ms);
}

void
NullContact::addComputesAndRequiresIntegrated(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* ms)
{
  Task* t = scinew Task("NullContact::exchangeMomentumIntegrated",
                        this,
                        &NullContact::exchangeMomentumIntegrated);

  sched->addTask(t, patches, ms);
}
