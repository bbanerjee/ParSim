/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/SimulationState.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <unistd.h>

using namespace Matiti;

SimulationState::SimulationState()
  : d_is_dynamic(true), d_modulus_type(ModulusType::Constant), d_horizon_factor(1.001)
{
  d_damping_coeff = {{0.0, 0.0, 0.0}};
}

SimulationState::SimulationState(bool isDynamic, 
                                 Array3 dampingCoeff, 
                                 ModulusType modulusType, 
                                 double horizonFactor)
{
 d_is_dynamic = isDynamic;
 d_damping_coeff = dampingCoeff;
 d_modulus_type = modulusType;
 d_horizon_factor = horizonFactor;
}

SimulationState::~SimulationState()
{
}

void
SimulationState::clone(const SimulationState& state)
{
 d_is_dynamic = state.isDynamic();
 d_damping_coeff = state.dampingCoeff();
 d_modulus_type = state.modulusType();
 d_horizon_factor = state.horizonFactor();
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
