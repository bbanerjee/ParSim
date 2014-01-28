#ifndef __VAANGO_PERIDYNAMICS_COMMON_H__
#define __VAANGO_PERIDYNAMICS_COMMON_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsSimulationStateP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  class Uintah::ProcessorGroup;

  class PeridynamicsCommon {

  public:

    PeridynamicsCommon(const Uintah::ProcessorGroup* myworld);
    virtual ~PeridynamicsCommon();

    virtual void materialProblemSetup(const Uintah::ProblemSpecP& prob_spec,
                                      PeridynamicsSimulationStateP& sharedState,
                                      PeridynamicsFlags* flags);

   protected:
    const Uintah::ProcessorGroup* d_myworld;
  };
}

#endif
