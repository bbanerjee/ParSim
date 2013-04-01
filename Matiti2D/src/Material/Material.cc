#include <Material/Material.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>
#include <iostream>
#include <sstream>

using namespace Uintah;


Material::Material()
{
}

Material::Material(ProblemSpecP& ps)
{
}

Material::~Material()
{
}

void Material::registerBondState(SimulationState* ss)
{
}
