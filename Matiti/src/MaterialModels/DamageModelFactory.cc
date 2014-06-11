#include <MaterialModels/DamageModelFactory.h>
#include <MaterialModels/DamageModelSimple.h>
#include <Pointers/DamageModelSP.h>

#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

DamageModelSP 
DamageModelFactory::create(Uintah::ProblemSpecP& damage_ps)
{
  // Get the attribute of the type
  std::string damage_type;
  if (!damage_ps->getAttribute("type", damage_type)) {
    std::ostringstream out;
    out << "**ERROR** Damage model does not have type attribute.";
    out << " Cannot determine which damage model to apply." ;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (damage_type == "isotropic") {
    return std::make_shared<DamageModelSimple>();
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown damage model type" << damage_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
