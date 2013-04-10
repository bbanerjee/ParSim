#include <Material.h> 
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Emu2DC;

Material::Material()
  : d_id(0), d_have_name(false), d_name(""), d_density(0.0), d_young_modulus(0.0), d_fracture_energy(0.0),
    d_damage_model(new DamageModel())
{
}

Material::~Material()
{
}

// Assumes that the <Material> block has already been found 
// and a ps=mat_ps has been assigned
void 
Material::initialize(Uintah::ProblemSpecP& ps)
{
  // Find the material name
  if (ps->getAttribute("name", d_name)) {
    d_have_name = true;
  } else {
    d_have_name = false;
  }

  // Get the material parameters
  ps->require("density", d_density);
  ps->require("young_modulus", d_young_modulus);
  ps->require("fracture_energy", d_fracture_energy);

  // Get the damage model
  d_damage_model->initialize(ps);
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Material& mat)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Material model: id = " << mat.d_id << std::endl;
    out << "  Have name = " << mat.d_have_name << " Name = " << mat.d_name << std::endl;
    out << "  Density = " << mat.d_density << " Young = " << mat.d_young_modulus
        << " Fracture enrgy = " << mat.d_fracture_energy << std::endl;
    out << *(mat.damageModel());
    return out;
  }
}
