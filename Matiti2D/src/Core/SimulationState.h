#ifndef MATITI_SIMULATIONSTATE_H
#define MATITI_SIMULATIONSTATE_H

#include <Vaango/Core/ProblemSpec/ProblemSpecP.h>
#include <Core/SimulationTime.h>
#include <Core/PeriMaterial.h>
#include <vector>

namespace Matiti {

  class SimulationState {

    public:
      double d_current_delt;
      double d_prev_delt;
      SimulationTime* d_simTime;

    public:
      SimulationState(Uintah::ProblemSpecP pspec);
      virtual ~SimulationState();

      void registerPeriMaterial(PeriMaterial* matl) {
        peri_matls.push_back(matl);
      }

      double getNumPeriMatls() const {
	return (int) peri_matls.size();
      }

      PeriMaterial* getPeriMaterial(int idx) const {
	return peri_matls[idx];
      }

    private:

      std::vector<PeriMaterial*> peri_matls;

      // Do not allow copying
      SimulationState(const SimulationState&);
      SimulationState& operator=(const SimulationState&);
   };

} // End namespace Matiti

#endif
