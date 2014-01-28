#ifndef __VAANGO_PERIDYNAMICS_SIMULATION_STATE_H
#define __VAANGO_PERIDYNAMICS_SIMULATION_STATE_H

#include <Core/Grid/SimulationState.h>

namespace Vaango 
{
  class PeridynamicsMaterial;
   
  class PeridynamicsSimulationState : public Uintah::SimulationState 
  {

  public:
    PeridynamicsSimulationState(Uintah::ProblemSpecP &ps);
    ~PeridynamicsSimulationState();
    void finalizePeridynamicsMaterials();
    void clearPeridynamicsMaterials();

    void registerPeridynamicsMaterial(PeridynamicsMaterial* mat);
    void registerPeridynamicsMaterial(PeridynamicsMaterial* mat, unsigned int index);

    int getNumPeridynamicsMatls() const {
      return (int) peridynamics_matls.size();
    }

    PeridynamicsMaterial* getPeridynamicsMaterial(int idx) const {
      return peridynamics_matls[idx];
    }
  
    const Uintah::MaterialSet* allPeridynamicsMaterials() const;

  private:

    PeridynamicsSimulationState(const PeridynamicsSimulationState&);
    PeridynamicsSimulationState& operator=(const PeridynamicsSimulationState&);
      
    std::vector<PeridynamicsMaterial*>  peridynamics_matls;
    Uintah::MaterialSet* all_peridynamics_matls;

}; // end class PeridynamicsSimulationState

} // End namespace Vaango

#endif
