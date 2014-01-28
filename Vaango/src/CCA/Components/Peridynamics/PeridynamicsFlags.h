#ifndef __VAANGO_PERIDYNAMICS_FLAGS_H__
#define __VAANGO_PERIDYNAMICS_FLAGS_H__

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Ports/Output.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class PeridynamicsFlags
    \brief A structure that store the flags used for a Peridynamics simulation
  */
  /////////////////////////////////////////////////////////////////////////////


  class PeridynamicsFlags {

  public:

    const Uintah::ProcessorGroup* d_myworld;
    Uintah::Vector d_gravity;

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

