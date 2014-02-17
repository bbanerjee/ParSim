#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModelFactory.h>

#include <CCA/Components/Peridynamics/MaterialModels/LinearElasticBondModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Vaango;

PeridynamicsMaterialModel* 
PeridynamicsMaterialModelFactory::create(Uintah::ProblemSpecP& ps,
                                         PeridynamicsFlags* flags)
{
  Uintah::ProblemSpecP child = ps->findBlock("material_model");
  if(!child)
    throw Uintah::ProblemSetupException("Cannot find material_model tag", __FILE__, __LINE__);
  std::string mat_type;
  if(!child->getAttribute("type", mat_type))
    throw Uintah::ProblemSetupException("No type for material_model", __FILE__, __LINE__);
   
  if (mat_type == "linear_elastic_bond_based")
    return(scinew LinearElasticBondModel(child, flags));
  else 
    throw Uintah::ProblemSetupException("Unknown peridynamic material type ("+mat_type+")", __FILE__, __LINE__);

  return 0;
}
