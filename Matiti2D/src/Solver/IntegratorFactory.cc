#include <Solver/IntegratorFactory.h>
#include <Solver/Integrator.h>
#include <iostream>
#include <string>

using namespace Matiti;

static Integrator* IntegrationFactory::create()
{
  std::string d_integrator("explicit");
  std::string d_integration_scheme("velocity verlet");

  if (d_integrator == "explicit") {
    if (d_integration_scheme == "velocity verlet") {
      return (new VelocityVerletIntegrator());
    } else {
      std::cerr << "Integration scheme of type " << d_integration_scheme << " not implemented" 
                << std::endl;
      exit(0);
    }
  } else {
    std::cerr << "Integration of type " << d_integrator << " not implemented" << std::endl;
    exit(0);
  }
  return 0;
}
