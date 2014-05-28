#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Vaango;

using Uintah::ProblemSpecP;

PeridynamicsFlags::PeridynamicsFlags(const Uintah::ProcessorGroup* myworld)
{
  d_myworld = myworld;
}

PeridynamicsFlags::~PeridynamicsFlags()
{
}

void
PeridynamicsFlags::readPeridynamicsFlags(ProblemSpecP& ps, Uintah::Output* dataArchive)
{
  ProblemSpecP root = ps->getRootNode();

  // Look for a <Peridynamics> block
  ProblemSpecP peridynamics_ps = root->findBlock("Peridynamics");
  if (!peridynamics_ps) {
    return;
  }

  // Find physical constants that are used by peridynamics
  d_gravity=SCIRun::Vector(0,0,0);
  peridynamics_ps->get("gravity",d_gravity);

  // Set the integrator type
  d_integrator_type = "forward_euler";
  d_integrator = ForwardEuler;
  peridynamics_ps->get("time_integrator", d_integrator_type);
  if (d_integrator_type == "forward_euler") {
    d_integrator = ForwardEuler;
  } else {
    d_integrator = ForwardEuler;  // default value
  }

  // Get the number of cells in the horizon
  d_numCellsInHorizon = 2.0;
  peridynamics_ps->get("num_cells_in_horizon", d_numCellsInHorizon);

}

void
PeridynamicsFlags::outputProblemSpec(ProblemSpecP& ps)
{
  ps->appendElement("gravity", d_gravity);
  ps->appendElement("time_integrator", d_integrator_type);
  ps->appendElement("num_cells_in_horizon", d_numCellsInHorizon);
}

