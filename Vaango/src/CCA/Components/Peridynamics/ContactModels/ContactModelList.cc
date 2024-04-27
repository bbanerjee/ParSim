/*
 * The MIT License
 *
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

#include <CCA/Components/Peridynamics/ContactModels/ContactModelList.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>

using namespace Vaango;

using Uintah::DataWarehouse;
using Uintah::MaterialSet;
using Uintah::MaterialSubset;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::ProblemSpecP;
using Uintah::ProcessorGroup;
using Uintah::SchedulerP;

ContactModelList::ContactModelList(const Uintah::ProcessorGroup* myworld,
                                   PeridynamicsLabel* labels,
                                   PeridynamicsFlags* flags)
  : ContactModelBase(myworld, nullptr, labels, flags, nullptr)
{
}

ContactModelList::~ContactModelList()
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    delete *iter;
  }
}

void
ContactModelList::outputProblemSpec(ProblemSpecP& ps)
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->outputProblemSpec(ps);
  }
}

void
ContactModelList::add(ContactModelBase* model)
{
  d_modelList.push_back(model);
}

void
ContactModelList::exchangeMomentumInterpolated(const ProcessorGroup* pg,
                                               const PatchSubset* patches,
                                               const MaterialSubset* matls,
                                               DataWarehouse* old_dw,
                                               DataWarehouse* new_dw)
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->exchangeMomentumInterpolated(pg, patches, matls, old_dw, new_dw);
  }
}

void
ContactModelList::exchangeMomentumIntegrated(const ProcessorGroup* pg,
                                             const PatchSubset* patches,
                                             const MaterialSubset* matls,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->exchangeMomentumIntegrated(pg, patches, matls, old_dw, new_dw);
  }
}

void
ContactModelList::addComputesAndRequiresInterpolated(SchedulerP& sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls)
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->addComputesAndRequiresInterpolated(sched, patches, matls);
  }
}

void
ContactModelList::addComputesAndRequiresIntegrated(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls)
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->addComputesAndRequiresIntegrated(sched, patches, matls);
  }
}
