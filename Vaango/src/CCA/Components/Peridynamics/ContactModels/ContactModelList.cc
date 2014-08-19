#include <CCA/Components/Peridynamics/ContactModels/ContactModelList.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>

using namespace Vaango;

using Uintah::ProcessorGroup;
using Uintah::ProblemSpecP;
using Uintah::SchedulerP;
using Uintah::PatchSubset;
using Uintah::MaterialSubset;
using Uintah::PatchSet;
using Uintah::MaterialSet;
using Uintah::DataWarehouse;

ContactModelList::ContactModelList(const Uintah::ProcessorGroup* myworld, 
                                   PeridynamicsLabel* labels,
                                   PeridynamicsFlags* flags)
  : ContactModelBase(myworld, labels, flags, 0)
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
ContactModelList::addComputesAndRequiresInterpolated(SchedulerP & sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls) 
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->addComputesAndRequiresInterpolated(sched, patches, matls);
  }
}

void
ContactModelList::addComputesAndRequiresIntegrated(SchedulerP & sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls) 
{
  for (auto iter = d_modelList.begin(); iter != d_modelList.end(); iter++) {
    (*iter)->addComputesAndRequiresIntegrated(sched, patches, matls);
  }
}
