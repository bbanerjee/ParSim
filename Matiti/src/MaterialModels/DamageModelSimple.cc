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

#include <MaterialModels/DamageModelSimple.h> 
#include <Core/Node.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Matiti;

DamageModelSimple::DamageModelSimple() 
  : d_damage_viscosity(0.0, 0.0, 0.0), 
    d_damage_stretch(0.0, 0.0, 0.0), 
    d_damage_index_max(0.0)
{
}

DamageModelSimple::DamageModelSimple(const DamageModelSimple& dam)
  : d_damage_viscosity(dam.d_damage_viscosity),
    d_damage_stretch(dam.d_damage_stretch),
    d_damage_index_max(dam.d_damage_index_max)
{
}

DamageModelSimple::~DamageModelSimple()
{
}

DamageModelSP
DamageModelSimple::clone()
{
  return std::make_shared<DamageModelSimple>(*this);
}

void 
DamageModelSimple::setVariation(double randomNum, double coeffOfVar)
{
  d_damage_viscosity *= (std::abs(1.0+randomNum*coeffOfVar));
  d_damage_stretch *= (std::abs(1.0+randomNum*coeffOfVar));
  d_damage_index_max *= (std::abs(1.0+randomNum*coeffOfVar));
}

void 
//DamageModelSimple::setAverage(const DamageModelSimple* dam1, const DamageModelSimple* dam2)
DamageModelSimple::setAverage(const DamageModelBase* dam1i, const DamageModelBase* dam2i)
{
  const DamageModelSimple* dam1 = dynamic_cast<const DamageModelSimple*>(dam1i);
  const DamageModelSimple* dam2 = dynamic_cast<const DamageModelSimple*>(dam2i);
  d_damage_viscosity = (dam1->d_damage_viscosity+dam2->d_damage_viscosity)*0.5;
  d_damage_stretch = (dam1->d_damage_stretch+dam2->d_damage_stretch)*0.5;
  d_damage_index_max = 0.5*(dam1->d_damage_index_max+dam2->d_damage_index_max);
}

void
DamageModelSimple::clone(const DamageModelSimple* dam)
{
  d_damage_viscosity = dam->d_damage_viscosity;
  d_damage_stretch = dam->d_damage_stretch;
  d_damage_index_max = dam->d_damage_index_max;
}

void
DamageModelSimple::clone(const DamageModelSimple* dam,
                         double randomNum,
                         double coeffOfVar)
{
  d_damage_viscosity = dam->d_damage_viscosity*(std::abs(1.0+randomNum*coeffOfVar));
  d_damage_stretch = dam->d_damage_stretch*(std::abs(1.0+randomNum*coeffOfVar));
  d_damage_index_max = dam->d_damage_index_max*(std::abs(1.0+randomNum*coeffOfVar));
}

void
DamageModelSimple::cloneAverage(const DamageModelSimple* dam1, 
                                const DamageModelSimple* dam2)
{
  d_damage_viscosity = (dam1->d_damage_viscosity+dam2->d_damage_viscosity)*0.5;
  d_damage_stretch = (dam1->d_damage_stretch+dam2->d_damage_stretch)*0.5;
  d_damage_index_max = 0.5*(dam1->d_damage_index_max+dam2->d_damage_index_max);
}

void 
DamageModelSimple::initialize(Uintah::ProblemSpecP& dam_ps)
{
  // BB : This has now been moved up to Material.cc
  // Check for the <DamageModel> block
  // Uintah::ProblemSpecP dam_ps = ps->findBlock("DamageModel");
  // if (!dam_ps) return;

  // Get the material parameters
  Uintah::Vector viscosity(0.0, 0.0, 0.0);
  Uintah::Vector stretch(0.0, 0.0, 0.0);
  dam_ps->require("damage_viscosity", viscosity);
  dam_ps->require("damage_stretch", stretch);
  dam_ps->require("damage_index_max", d_damage_index_max);
  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_damage_viscosity[ii] = viscosity[ii];
    d_damage_stretch[ii] = stretch[ii];
  }
}

double 
DamageModelSimple::computeDamageFactor(const double& damage_index) const
{
  double coef3 = d_damage_stretch[2];
  if (damage_index > 0.9999) return coef3;

  double coef1 = d_damage_stretch[0];
  if (!(damage_index > coef1)) return 1.0;

  double coef2 = d_damage_stretch[1];
  if (!(coef2 > 0.0 && coef3 > 1.0)) return 1.0;
  
  double  damage_fac = 1.0 + coef2*(damage_index-coef1)/(1.0-damage_index);
  damage_fac = std::min(damage_fac, coef3);

  return damage_fac;
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const DamageModelSimple& dam)
  {
    //out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Damage model:" << std::endl;
    out << "  Viscosity = [" << dam.d_damage_viscosity[0] << ", " << dam.d_damage_viscosity[1] 
        << ", " << dam.d_damage_viscosity[2] << "]";
    out << "  Stretch = [" << dam.d_damage_stretch[0] << ", " << dam.d_damage_stretch[1] 
        << ", " << dam.d_damage_stretch[2] << "]";
    out << "  Damage index = " << dam.d_damage_index_max << std::endl;
    return out;
  }
}
