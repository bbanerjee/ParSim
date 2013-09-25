#include <SimulationState.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <unistd.h>

using namespace Matiti;

SimulationState::SimulationState()
  : d_is_dynamic(true), d_modulus_type(ModulusType::Constant), d_horizon_factor(0.0)
{
}

SimulationState::~SimulationState()
{
}

void 
SimulationState::initialize(const Uintah::ProblemSpecP& ps)
{
  // Get peridynamics flags
  Uintah::ProblemSpecP peri_ps = ps->findBlock("Peridynamics");
  if (!peri_ps) {
    throw Exception("**ERROR** <Peridynamics> tag not found", __FILE__, __LINE__); 
  }

  std::string sim_type;
  peri_ps->require("simulation_type", sim_type);
  d_is_dynamic = true;
  if (sim_type == "static") d_is_dynamic = false;

  // Artificial viscosity coefficients for dynamic simulations 
  // Hard-coded for now 
  if (d_is_dynamic) {
    d_damping_coeff = {{0.2, 0.2, 0.2}};
  } else {
    d_damping_coeff = {{0.0, 0.0, 0.0}};
  }

  std::string mod_type;
  peri_ps->require("modulus_type", mod_type);
  d_modulus_type = ModulusType::Constant;
  if (mod_type == "conical") d_modulus_type = ModulusType::Conical;

  peri_ps->require("horizon_factor", d_horizon_factor);
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const SimulationState& state)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Dynamic = " << state.d_is_dynamic << " Modulus type = " << (int) (state.d_modulus_type) 
        << " Horizon factor = " << state.d_horizon_factor << std::endl;
    return out;
  }

}
