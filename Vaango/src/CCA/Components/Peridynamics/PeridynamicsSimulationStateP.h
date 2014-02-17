#ifndef __VAANGO_PERIDYNAMICS_SIMULATIONSTATEP_H
#define __VAANGO_PERIDYNAMICS_SIMULATIONSTATEP_H

namespace Uintah {
   template<class T> class Handle;
}

namespace Vaango {
   class PeridynamicsSimulationState;
   typedef Uintah::Handle<PeridynamicsSimulationState> PeridynamicsSimulationStateP;
}

#endif

