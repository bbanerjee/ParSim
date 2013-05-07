#include <SimulationState.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <unistd.h>

using namespace Emu2DC;

SimulationState::SimulationState()
  : d_max_time(0.0), d_delT(0.0), d_max_iter(0), d_output_folder_name(""), d_output_file_name(""),
    d_output_iter_interval(0), d_is_dynamic(true), d_modulus_type(ModulusType::Constant), 
    d_horizon_factor(0.0)
{
}

SimulationState::~SimulationState()
{
}

void 
SimulationState::initialize(const Uintah::ProblemSpecP& ps)
{
  // get simulation time information
  Uintah::ProblemSpecP time_ps = ps->findBlock("Time");
  if (!time_ps) {
    throw Exception("**ERROR** <Time> tag not found", __FILE__, __LINE__); 
  }

  time_ps->require("max_time", d_max_time);
  time_ps->require("max_iterations", d_max_iter);
  time_ps->require("delt", d_delT);
  
  // Get output files/interval
  Uintah::ProblemSpecP io_ps = ps->findBlock("Output");
  if (!io_ps) {
    throw Exception("**ERROR** <Output> tag not found", __FILE__, __LINE__); 
  } 

  io_ps->require("output_file", d_output_file_name);
  io_ps->require("output_iteration_interval", d_output_iter_interval);
  char buffer[2000];
  char * str = getcwd( buffer, 2000 );
  if (str == NULL) {
    throw Exception("**ERROR** Directory not returned by getcwd()", __FILE__, __LINE__); 
  } else {
    d_output_folder_name = std::string(buffer);
  }

  // Get peridynamics flags
  Uintah::ProblemSpecP peri_ps = ps->findBlock("Peridynamics");
  if (!io_ps) {
    throw Exception("**ERROR** <Peridynamics> tag not found", __FILE__, __LINE__); 
  }

  peri_ps->require("dimensions", d_dimensions);

  std::string sim_type;
  peri_ps->require("simulation_type", sim_type);
  d_is_dynamic = false;
  if (sim_type == "static") d_is_dynamic = false;

  std::string mod_type;
  peri_ps->require("modulus_type", mod_type);
  d_modulus_type = ModulusType::Constant;
  if (mod_type == "conical") d_modulus_type = ModulusType::Conical;

  peri_ps->require("horizon_factor", d_horizon_factor);
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const SimulationState& state)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Simulation state variables:" << std::endl;
    out << "Max T = " << state.d_max_time << " del T = " << state.d_delT
        << " Max iter = " << state.d_max_iter << std::endl;
    out << "Output dir = " << state.d_output_folder_name << " Output file = " << state.d_output_file_name
        << std::endl;
    out << "  Output iteration interval = " << state.d_output_iter_interval << std::endl;
    out << "Dynamic = " << state.d_is_dynamic << " Modulus type = " << (int) (state.d_modulus_type) 
        << " Horizon factor = " << state.d_horizon_factor << std::endl;
    return out;
  }

}
