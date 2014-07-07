#include <MaterialModels/DamageModelBase.h> 

#include <iostream>

using namespace Matiti;
  
DamageModelBase::DamageModelBase() 
{
}

DamageModelBase::~DamageModelBase() 
{
}

/* Make copies without having to read in stuff */
/*
DamageModelSP 
DamageModelBase::clone()
{
  std::cout << "**ERROR** Base class clone damage called" << std::endl;
}

void 
DamageModelBase::setVariation(double randomNumber, double coeffOfVar)
{
  std::cout << "**ERROR** Base set variation damage called" << std::endl;
}

void 
DamageModelBase::setAverage(const DamageModelBase* dam1, const DamageModelBase* dam2)
{
  std::cout << "**ERROR** Base class set average damage called" << std::endl;
}
*/

/*
void 
DamageModelBase::clone(const DamageModelBase* dam)
{
  std::cout << "**ERROR** Base class clone damage called" << std::endl;
}

void 
DamageModelBase::clone(const DamageModelBase* dam, double randomNumber, double coeffOfVar)
{
  std::cout << "**ERROR** Base class clone damage called" << std::endl;
}

void 
DamageModelBase::cloneAverage(const DamageModelBase* dam1, const DamageModelBase* dam2)
{
  std::cout << "**ERROR** Base class clone damage called" << std::endl;
}

*/

/**
 *  Compute the damage factor  given the damage index at a node.
 */
double 
DamageModelBase::computeDamageFactor(const double& damage_index) const
{
  std::cout << "**ERROR** Base class compute damage called" << std::endl;
}


