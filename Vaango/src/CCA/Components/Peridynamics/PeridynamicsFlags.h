#ifndef __VAANGO_PERIDYNAMICS_FLAGS_H__
#define __VAANGO_PERIDYNAMICS_FLAGS_H__

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Ports/Output.h>
#include <Core/Geometry/Vector.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class PeridynamicsFlags
    \brief A structure that store the flags used for a Peridynamics simulation
  */
  /////////////////////////////////////////////////////////////////////////////


  class PeridynamicsFlags {

  public:

    enum IntegratorType {
      ForwardEuler,
      VelocityVerlet,
      BackwardEuler,
      None
    };

    const Uintah::ProcessorGroup* d_myworld;

    SCIRun::Vector d_gravity;
    std::string d_integrator_type; // Explicit or implicit time integration
    IntegratorType d_integrator;

    PeridynamicsFlags(const Uintah::ProcessorGroup* myworld);

    virtual ~PeridynamicsFlags();

    virtual void readPeridynamicsFlags(Uintah::ProblemSpecP& ps, Uintah::Output* dataArchive);
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

  private:

    PeridynamicsFlags(const PeridynamicsFlags& state);
    PeridynamicsFlags& operator=(const PeridynamicsFlags& state);
    
  };

} // End namespace Vaango

#endif  // __VAANGO_PERIDYNAMICS_FLAGS_H__ 

