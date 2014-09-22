#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GeometryPiece/GeometryPiece.h>


using namespace Vaango;

ParticleForceBC::ParticleForceBC(Uintah::ProblemSpecP& ps)
  : ParticleLoadBCBase(ps)
{
  // Read and save the load curve information
  d_loadCurve = scinew ParticleLoadCurve<SCIRun::Vector>(ps);
}

ParticleForceBC::~ParticleForceBC()
{
  delete d_loadCurve;
}

void 
ParticleForceBC::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP press_ps = ps->appendChild("force");
  Uintah::ProblemSpecP geom_ps = press_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  press_ps->appendElement("numberOfParticlesOnLoadSurface",d_numParticlesOnLoadSurface);
  d_loadCurve->outputProblemSpec(press_ps);
}

std::string 
ParticleForceBC::getType() const
{
  return "Force";
}

