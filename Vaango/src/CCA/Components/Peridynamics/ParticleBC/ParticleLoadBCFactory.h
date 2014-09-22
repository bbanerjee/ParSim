#ifndef __VAANGO_PARTICLE_LOAD_BC_FACTORY__
#define __VAANGO_PARTICLE_LOAD_BC_FACTORY__

#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <vector>


namespace Vaango {

  class ParticleLoadBCFactory
  {
  public:
    static void create(const Uintah::ProblemSpecP& ps);
    static void clean(); // delete all ParticleLoadBCs
    static std::vector<ParticleLoadBCBase*> particleLoadBCs;
  };

} // End namespace Vaango


#endif /* __VAANGO_PARTICLE_BC_FACTORY__ */

