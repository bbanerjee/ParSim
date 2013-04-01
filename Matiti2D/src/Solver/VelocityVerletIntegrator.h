#ifndef MATITI_VELOCITYVERLET_INTEGRATOR_H
#define MATITI_VELOCITYVERLET_INTEGRATOR_H

namespace Matiti {

  class VelocityVerletIntegrator {

    public:

      VelocityVerletIntegrator();
      virtual ~VelocityVerletIntegrator();

      void integrate();

  };
}

#endif
