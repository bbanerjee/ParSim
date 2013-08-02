#include <DamageModel.h> 
#include <Node.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Matiti;

DamageModel::DamageModel()
  : d_damage_viscosity(0.0, 0.0, 0.0), d_damage_stretch(0.0, 0.0, 0.0), d_damage_index_max(0.0)
{
}

DamageModel::DamageModel(const DamageModel& dam)
  : d_damage_viscosity(dam.d_damage_viscosity),
    d_damage_stretch(dam.d_damage_stretch),
    d_damage_index_max(dam.d_damage_index_max)
{
}

DamageModel::~DamageModel()
{
}

void
DamageModel::clone(const DamageModelUP& dam)
{
  d_damage_viscosity = dam->d_damage_viscosity;
  d_damage_stretch = dam->d_damage_stretch;
  d_damage_index_max = dam->d_damage_index_max;
}

void
DamageModel::clone(const DamageModelUP& dam,
                   double randomNum,
                   double coeffOfVar)
{
  d_damage_viscosity = dam->d_damage_viscosity*(std::abs(1.0+randomNum*coeffOfVar));
  d_damage_stretch = dam->d_damage_stretch*(std::abs(1.0+randomNum*coeffOfVar));
  d_damage_index_max = dam->d_damage_index_max*(std::abs(1.0+randomNum*coeffOfVar));
}

void
DamageModel::cloneAverage(const DamageModelUP& dam1, const DamageModelUP& dam2)
{
  d_damage_viscosity = (dam1->d_damage_viscosity+dam2->d_damage_viscosity)*0.5;
  d_damage_stretch = (dam1->d_damage_stretch+dam2->d_damage_stretch)*0.5;
  d_damage_index_max = 0.5*(dam1->d_damage_index_max+dam2->d_damage_index_max);
}

void 
DamageModel::initialize(const Uintah::ProblemSpecP& ps)
{
  // Check for the <DamageModel> block
  Uintah::ProblemSpecP dam_ps = ps->findBlock("DamageModel");
  if (!dam_ps) return;

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
DamageModel::computeDamageFactor(const double& damage_index) const
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

  std::ostream& operator<<(std::ostream& out, const DamageModel& dam)
  {
    out.setf(std::ios::floatfield);
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
