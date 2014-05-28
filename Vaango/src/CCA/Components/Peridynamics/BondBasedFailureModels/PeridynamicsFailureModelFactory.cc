#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModelFactory.h>

#include <CCA/Components/Peridynamics/FailureModels/BondDamageModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Vaango;

PeridynamicsFailureModel* 
PeridynamicsFailureModelFactory::create(Uintah::ProblemSpecP& ps,
                                        PeridynamicsFlags* flags)
{
  Uintah::ProblemSpecP child = ps->findBlock("failure_model");
  if(!child)
    throw Uintah::ProblemSetupException("Cannot find failure_model tag", __FILE__, __LINE__);
  std::string mat_type;
  if(!child->getAttribute("type", mat_type))
    throw Uintah::ProblemSetupException("No type for failure_model", __FILE__, __LINE__);
   
  if (mat_type == "bond_damage")
    return(scinew BondDamageModel(child, flags));
  else 
    throw Uintah::ProblemSetupException("Unknown peridynamic failure model type ("+mat_type+")", __FILE__, __LINE__);

  return 0;
}
