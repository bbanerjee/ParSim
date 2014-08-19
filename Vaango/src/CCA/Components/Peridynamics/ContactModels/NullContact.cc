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
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/SimulationStateP.h>
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
using Uintah::SimulationStateP;

NullContact::NullContact(const ProcessorGroup* myworld,
                         SimulationStateP& ss,
                         PeridynamicsLabel* labels,
                         PeridynamicsFlags* flags)
  : ContactModelBase(myworld, labels, flags, 0)
{
  d_sharedState = ss;
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
