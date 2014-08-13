#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCFactory.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Vaango;

std::vector<ParticleLoadBCBase*> ParticleLoadBCFactory::particleLoadBCs;

void 
ParticleLoadBCFactory::create(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP test = ps->findBlock("ParticleBC");
  if (test){

    Uintah::ProblemSpecP current_ps = ps->findBlock("ParticleBC")->findBlock("Peridynamics");


    for(Uintah::ProblemSpecP child = current_ps->findBlock("force"); child != 0;
        child = child->findNextBlock("force") ) {
       particleLoadBCs.push_back(scinew ParticleForceBC(child));
    }

    for(Uintah::ProblemSpecP child = current_ps->findBlock("normal_force"); child != 0;
        child = child->findNextBlock("normal_force") ) {
       particleLoadBCs.push_back(scinew ParticleNormalForceBC(child));
    }

    for(Uintah::ProblemSpecP child = current_ps->findBlock("pressure"); child != 0;
        child = child->findNextBlock("pressure") ) {
       particleLoadBCs.push_back(scinew ParticlePressureBC(child));
    }

  }
}

void ParticleLoadBCFactory::clean()
{
  for (int i = 0; i < static_cast<int>(particleLoadBCs.size()); i++)
    delete particleLoadBCs[i];
  particleLoadBCs.clear();
}
