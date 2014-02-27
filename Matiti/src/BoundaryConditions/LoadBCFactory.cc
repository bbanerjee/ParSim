#include <BoundaryConditions/LoadBCFactory.h>
#include <BoundaryConditions/LoadBC.h>
#include <BoundaryConditions/ForceBC.h>
#include <BoundaryConditions/TractionBC.h>
#include <Pointers/LoadBCSP.h>

#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

LoadBCSP 
LoadBCFactory::create(Uintah::ProblemSpecP& load_ps)
{
  // Get the attribute of the type
  std::string load_type;
  if (!load_ps->getAttribute("type", load_type)) {
    std::ostringstream out;
    out << "**ERROR** Load BC does not have type attribute to determine if it is a force/traction/pressure etc." 
        << load_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (load_type == "force") {
    return (std::make_shared<ForceBC>());
  } else if (load_type == "traction") {
    return (std::make_shared<TractionBC>());
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown load boundary condition type" << load_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
