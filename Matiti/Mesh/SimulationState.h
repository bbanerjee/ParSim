#ifndef MATITI_SIMULATIONSTATE_H
#define MATITI_SIMULATIONSTATE_H

#include <Common/RefCounted.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/SimulationTime.h>
#include <Core/Grid/Ghost.h>
#include <Core/Geometry/Vector.h>
#include <Core/Math/MinMax.h>

#include <map>
#include <vector>
#include <iostream>

namespace Uintah {
  class VarLabel;
  class Material; 
}

namespace Matiti {

  class PeriMaterial;
   
  class SimulationState : public RefCounted {

    public:

      SimulationState(Uintah::ProblemSpecP &ps);
      ~SimulationState();

      void clearMaterials();

      const Uintah::VarLabel* get_delt_label() const {
        return delt_label;
      }

      void registerPeriMaterial(PeriMaterial*);
      void registerPeriMaterial(PeriMaterial*,unsigned int index);

      int getNumMatls() const {
        return (int)matls.size();
      }
      int getNumPeriMatls() const {
        return (int)peri_matls.size();
      }

      Uintah::Material* getMaterial(int idx) const {
        return matls[idx];
      }
      PeriMaterial* getPeriMaterial(int idx) const {
        return peri_matls[idx];
      }

      void finalizeMaterials();
      const MaterialSet* allPeriMaterials() const;
      const MaterialSet* allMaterials() const;

      double getElapsedTime() const { return d_elapsed_time; }
      void   setElapsedTime(double t) { d_elapsed_time = t; }

      // Returns the integer timestep index of the top level of the
      // simulation.  All simulations start with a top level time step
      // number of 0.  This value is incremented by one for each
      // simulation time step processed.  The 'set' function should only
      // be called by the SimulationController at the beginning of a
      // simulation.  The 'increment' function is called by the
      // SimulationController at the end of each timestep.
      int  getCurrentTimeStep() const { return d_timeStep; }
      void setCurrentTimeStep( int ts ) { d_timeStep = ts; }
      void incrementCurrentTimeStep() { d_timeStep++; }

      Uintah::Material* parseAndLookupMaterial(Uintah::ProblemSpecP& params,
                                       const std::string& name) const;
      Uintah::Material* getMaterialByName(const std::string& name) const;

      inline int getMaxMatlIndex() { return max_matl_index; }

      vector<vector<const Uintah::VarLabel* > > d_meshNodeState;
      vector<vector<const Uintah::VarLabel* > > d_meshNodeState_preReloc;

      double d_prev_delt;
      double d_current_delt;

      Uintah::SimulationTime* d_simTime;

    private:

      void registerMaterial(Material*);
      void registerMaterial(Material*,unsigned int index);

      SimulationState(const SimulationState&);
      SimulationState& operator=(const SimulationState&);
      
      const Uintah::VarLabel* delt_label;

      std::vector<Uintah::Material*>     matls;
      std::vector<PeriMaterial*> peri_matls;

      std::map<std::string, Uintah::Material*> named_matls;

      Uintah::MaterialSet    * all_peri_matls;
      Uintah::MaterialSet    * all_matls;

      // The time step that the top level (w.r.t. AMR) is at during a
      // simulation.  Usually corresponds to the Data Warehouse generation
      // number (it does for non-restarted, non-amr simulations).  I'm going to
      // attempt to make sure that it does also for restarts.
      int    d_timeStep;
      double d_elapsed_time;
  
  }; // end class SimulationState

} // End namespace Matiti

#endif
