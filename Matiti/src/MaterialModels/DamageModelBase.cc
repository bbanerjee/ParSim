/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
  return 0;
}


