#ifndef __VAANGO_PARTICLE_CREATORFACTORY_H_
#define __VAANGO_PARTICLE_CREATORFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>

namespace Vaango {

  class ParticleCreator;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;
  class PeridynamicsFlags;

  class ParticleCreatorFactory
  {
  public:
    static ParticleCreator* create(Uintah::ProblemSpecP& ps, PeridynamicsMaterial* mat,
                                   PeridynamicsFlags* flags);


  };
} // End namespace Vaango
      
#endif /* __VAANGO_PARTICLE_CREATORFACTORY_H_ */
