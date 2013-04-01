#ifndef MATITI_INTEGRATOR_H
#define MATITI_INTEGRATOR_H

namespace Matiti {

  class Integrator {

    public:

      Integrator();
      virtual ~Integrator();

      virtual void integrate() = 0;

  };
}

#endif
