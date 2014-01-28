#include <CCA/Components/Peridynamics/PeridynamicsSimulationState.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Vaango;

PeridynamicsSimulationState::PeridynamicsSimulationState(Uintah::ProblemSpecP &ps)
  : Uintah::SimulationState(ps)
{
  all_peridynamics_matls = 0;
}

void 
PeridynamicsSimulationState::registerPeridynamicsMaterial(PeridynamicsMaterial* matl)
{
  peridynamics_matls.push_back(matl);
  registerMaterial(matl);
}

void 
PeridynamicsSimulationState::registerPeridynamicsMaterial(PeridynamicsMaterial* matl,
                                                          unsigned int index)
{
  peridynamics_matls.push_back(matl);
  registerMaterial(matl, index);
}

void 
PeridynamicsSimulationState::finalizePeridynamicsMaterials()
{
  if (all_peridynamics_matls && all_peridynamics_matls->removeReference())
    delete all_peridynamics_matls;
  all_peridynamics_matls = scinew Uintah::MaterialSet();
  all_peridynamics_matls->addReference();
  vector<int> tmp_peridynamics_matls(peridynamics_matls.size());
  for( int i=0; i<(int)peridynamics_matls.size(); i++ ) {
    tmp_peridynamics_matls[i] = peridynamics_matls[i]->getDWIndex();
  }
  all_peridynamics_matls->addAll(tmp_peridynamics_matls);
}

void 
PeridynamicsSimulationState::clearPeridynamicsMaterials()
{
  if(all_peridynamics_matls && all_peridynamics_matls->removeReference())
    delete all_peridynamics_matls;
  }
  peridynamics_matls.clear();
  all_peridynamics_matls = 0;
}

PeridynamicsSimulationState::~PeridynamicsSimulationState()
{
  clearPeridynamicsMaterials();
}

const 
Uintah::MaterialSet* PeridynamicsSimulationState::allPeridynamicsMaterials() const
{
  ASSERT(all_peridynamics_matls != 0);
  return all_peridynamics_matls;
}
