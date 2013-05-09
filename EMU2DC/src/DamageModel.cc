#include <DamageModel.h> 
#include <Node.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Emu2DC;

DamageModel::DamageModel()
{
  d_damage_viscosity = {{0.0, 0.0, 0.0}};
  d_damage_index = 0.0;
  d_damage_stretch = {{0.0, 0.0, 0.0}};
}

DamageModel::DamageModel(const DamageModel& dam)
  : d_damage_viscosity(dam.d_damage_viscosity), d_damage_index(dam.d_damage_index),
    d_damage_stretch(dam.d_damage_stretch)
{
}

DamageModel::~DamageModel()
{
}

void
DamageModel::clone(const DamageModelUP& dam)
{
  d_damage_viscosity = dam->d_damage_viscosity;
  d_damage_index = dam->d_damage_index;
  d_damage_stretch = dam->d_damage_stretch;
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
  dam_ps->require("damage_index", d_damage_index);
  dam_ps->require("damage_stretch", stretch);
  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_damage_viscosity[ii] = viscosity[ii];
    d_damage_stretch[ii] = stretch[ii];
  }
}

void 
DamageModel::updateDamageIndex(const NodeP& node)
{
  int num_bonds_init = node->initialFamilySize();
  int num_bonds_cur = node->currentFamilySize();
  d_damage_index = (double) num_bonds_cur/(double) num_bonds_init;
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const DamageModel& dam)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Damage model:" << std::endl;
    out << "  Viscosity = [" << dam.d_damage_viscosity[0] << ", " << dam.d_damage_viscosity[1] 
        << ", " << dam.d_damage_viscosity[2] << "]";
    out << "  Stretch = [" << dam.d_damage_stretch[0] << ", " << dam.d_damage_stretch[1] 
        << ", " << dam.d_damage_stretch[2] << "]";
    out << "  Damage index = " << dam.d_damage_index << std::endl;
    return out;
  }
}
