#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModelFactory.h>

#include <CCA/Components/Peridynamics/DamageModels/SphericalStrainEnergyDamageModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Vaango;

PeridynamicsDamageModel* 
PeridynamicsDamageModelFactory::create(Uintah::ProblemSpecP& ps,
                                       PeridynamicsLabel* labels,
                                       PeridynamicsFlags* flags)
{
  Uintah::ProblemSpecP child = ps->findBlock("damage_model");
  if(!child) {
    throw Uintah::ProblemSetupException("Cannot find damage_model tag in input file", __FILE__, __LINE__);
  }

  std::string damage_type;
  if(!child->getAttribute("type", damage_type)) {
    throw Uintah::ProblemSetupException("No type for damage_model in input file", __FILE__, __LINE__);
  }
   
  if (damage_type == "spherical_strain_energy") {

    return(scinew SphericalStrainEnergyDamageModel(child, labels, flags));

  } else {

    throw Uintah::ProblemSetupException("Unknown peridynamic damage model type ("+damage_type+") in input file", 
                                         __FILE__, __LINE__);
  }

  return 0;
}
