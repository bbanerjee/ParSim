#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreatorFactory.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsSimulationState.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Vaango;

ParticleCreator* ParticleCreatorFactory::create(Uintah::ProblemSpecP& ps, 
                                                PeridynamicsMaterial* mat,
                                                PeridynamicsFlags* flags)
{

  Uintah::ProblemSpecP cm_ps = ps->findBlock("material_model");
  string mat_type;
  cm_ps->getAttribute("type",mat_type);

  return scinew ParticleCreator(mat,flags);

}


