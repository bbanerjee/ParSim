#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Vaango;

PeridynamicsFlags::PeridynamicsFlags(const Uintah::ProcessorGroup* myworld)
{
  d_myworld = myworld;
}

PeridynamicsFlags::~PeridynamicsFlags()
{
}

void
PeridynamicsFlags::readPeridynamicsFlags(Uintah::ProblemSpecP& ps, Uintah::Output* dataArchive)
{
  Uintah::ProblemSpecP root = ps->getRootNode();
  Uintah::ProblemSpecP peridynamics_flag_ps = root->findBlock("Peridynamics");
  Uintah::ProblemSpecP phys_cons_ps = root->findBlock("PhysicalConstants");

  if(phys_cons_ps){
    phys_cons_ps->require("gravity",d_gravity);
  } else if (peridynamics_flag_ps) {
    peridynamics_flag_ps->require("gravity",d_gravity);
  } else{
    d_gravity=SCIRun::Vector(0,0,0);
  }
}

void
PeridynamicsFlags::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ps->appendElement("gravity", d_gravity);
}

