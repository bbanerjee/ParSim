#include <MaterialModels/Material.h> 
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

Material::Material()
  : d_id(0), d_have_name(false), d_name(""), d_density(0.0), d_young_modulus(0.0), d_fracture_energy(0.0),
    d_micro_modulus(0.0), d_strain(0.0), d_strain_energy(0.0),
    d_damage_model(new DamageModel())
{
  d_micro_modulus_model = MicroModulusModel::Constant;
}

Material::Material(const Material& mat)
  : d_id(mat.d_id), d_have_name(mat.d_have_name), d_name(mat.d_name), d_density(mat.d_density), 
    d_young_modulus(mat.d_young_modulus), d_fracture_energy(mat.d_fracture_energy),
    d_micro_modulus(mat.d_micro_modulus), d_strain(mat.d_strain), d_strain_energy(mat.d_strain_energy),
    d_damage_model(new DamageModel())
{
  d_micro_modulus_model = mat.d_micro_modulus_model;
  d_damage_model->clone(mat.d_damage_model);
}

Material::~Material()
{
}

Material& 
Material::operator=(const Material& mat)
{
  if (this != &mat) {
    d_id = mat.d_id;
    d_have_name = mat.d_have_name;
    d_name = mat.d_name; 
    d_density = mat.d_density; 
    d_young_modulus = mat.d_young_modulus;
    d_fracture_energy = mat.d_fracture_energy;
    d_micro_modulus_model = mat.d_micro_modulus_model;
    d_micro_modulus = mat.d_micro_modulus;
    d_strain = mat.d_strain;
    d_strain_energy = mat.d_strain_energy;
    d_damage_model->clone(mat.d_damage_model);
  }
  return *this;
}

void
Material::clone(const Material* mat)
{
  if (this != mat) {
    d_id = mat->d_id;
    d_have_name = mat->d_have_name;
    d_name = mat->d_name; 
    d_density = mat->d_density; 
    d_young_modulus = mat->d_young_modulus;
    d_fracture_energy = mat->d_fracture_energy;
    d_micro_modulus_model = mat->d_micro_modulus_model;
    d_micro_modulus = mat->d_micro_modulus;
    d_strain = mat->d_strain;
    d_strain_energy = mat->d_strain_energy;
    d_damage_model->clone(mat->d_damage_model);
  }
}

void
Material::clone(const Material* mat,
                double randomNum,
                double coeffOfVar)
{
  if (this != mat) {
    d_id = mat->d_id;
    d_have_name = mat->d_have_name;
    d_name = mat->d_name; 
    d_density = std::abs(mat->d_density*(1.0 + randomNum*coeffOfVar)); 
    d_young_modulus = std::abs(mat->d_young_modulus*(1.0 + randomNum*coeffOfVar));
    d_fracture_energy = std::abs(mat->d_fracture_energy*(1.0 + randomNum*coeffOfVar));
    d_micro_modulus_model = mat->d_micro_modulus_model;
    d_micro_modulus = mat->d_micro_modulus;
    d_strain = mat->d_strain;
    d_strain_energy = mat->d_strain_energy;
    d_damage_model->clone(mat->d_damage_model, randomNum, coeffOfVar);
  }
}

void
Material::cloneAverage(const Material* mat1, const Material* mat2)
{
  if (this != mat1) {
    d_id = mat1->d_id;
    d_have_name = mat1->d_have_name;
    d_name = mat1->d_name; 
    d_density = 0.5*(mat1->d_density+mat2->d_density); 
    d_young_modulus = 0.5*(mat1->d_young_modulus+mat2->d_young_modulus);
    d_fracture_energy = 0.5*(mat1->d_fracture_energy+mat2->d_fracture_energy);
    d_micro_modulus_model = mat1->d_micro_modulus_model;
    d_micro_modulus = 0.5*(mat1->d_micro_modulus+mat2->d_micro_modulus);
    d_strain = mat1->d_strain;
    d_strain_energy = mat1->d_strain_energy;
    d_damage_model->cloneAverage(mat1->d_damage_model,mat2->d_damage_model);
  }
}

//---------------------------------------------------------------------------------
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
  std::string micro_modulus_model;
  d_micro_modulus_model = MicroModulusModel::Constant;
  ps->require("micromodulus_type", micro_modulus_model);
  if (micro_modulus_model == "conical") d_micro_modulus_model = MicroModulusModel::Conical;

  ps->require("density", d_density);
  ps->require("young_modulus", d_young_modulus);
  ps->require("fracture_energy", d_fracture_energy);

  // Get the damage model
  d_damage_model->initialize(ps);
}

//---------------------------------------------------------------------------------
// This is an elastic model.  **TO DO** For more complex materials use factory concept.
void 
Material::computeForce(const Point3D& nodePos,
                       const Point3D& familyPos,
                       const Vector3D& nodeDisp,
                       const Vector3D& familyDisp,
                       const double& horizonSize,
                       Vector3D& force)
{
  // Compute new bond vector
  Vector3D xi = familyPos - nodePos;  
  Vector3D eta = familyDisp - nodeDisp;
  Vector3D pos_new = xi + eta;

  // Initial distance between the nodes
  double bond_length_init = xi.length();

  // New distance between the nodes
  double bond_length_new = pos_new.length();
  if ((bond_length_init <= 0.0) || (bond_length_new <= 0.0)) {
    std::ostringstream out;
    out << "**ERROR** Bond length is <= zero. " << std::endl
        << "    node 1 = " << nodePos << " node 2 = " << familyPos << std::endl
        << "    disp 1 = " << nodeDisp << " disp 2 = " << familyDisp << std::endl
        << "    xi = " << xi << " eta = " << eta << " pos = " << nodePos
        << "    pos_new = " << pos_new ;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
  if (bond_length_new > 1.0e16) {
    // diverges - break the bond and zero the force
    force.reset();
  }

  // Set the direction of the force
  force = pos_new/bond_length_new;

  // Compute bond micromodulus
  d_micro_modulus = computeMicroModulus(bond_length_init, horizonSize);
  //std::cout << "Micro modulus = " << d_micro_modulus << std::endl;
  if (std::isnan(d_micro_modulus)) {
    std::ostringstream out;
    out << "**ERROR** Micromodulus error.  Init bond length = " << bond_length_init
        << " Micromodulus = " << d_micro_modulus ;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // find bond strain and bond force 
  if (bond_length_new > 0.0 && bond_length_init < horizonSize) {

    // u is the displacement between the two nodes
    double disp = bond_length_new - bond_length_init; 
    if (bond_length_init > 0.0) {
      d_strain  = disp/bond_length_init;
    } else {
      d_strain = 0.0;
    }
 
    force *= (d_micro_modulus*d_strain);
    d_strain_energy = 0.5*d_micro_modulus*(disp*disp);
  } else {
    force *= 0.0;
    d_strain = 0.0;
    d_strain_energy = 0.0;
    d_micro_modulus = 0.0;
  }
  //std::cout << " xp = " << familyPos << " xi = " << nodePos 
  //          << " up = " << familyDisp << " ui = " << nodeDisp 
  //          << " force = " << force << " strain = " << d_strain << std::endl;
}

//---------------------------------------------------------------------------------
// Compute micromodulus
// See Hu et al. (2012), Bobaru et al (2009), and Silling and Askari (2005)
double
Material::computeMicroModulus(const double& bondLength, 
		              const double& horizonSize)
{
  double micromodulus = 0.0;
  double rad_cubed = horizonSize*horizonSize*horizonSize;

  // Assume Poisson's ratio = nu = 0.25
  if (d_micro_modulus_model == MicroModulusModel::Conical) {
    micromodulus = 32.0*d_young_modulus*(1.0-bondLength/horizonSize)/(M_PI*rad_cubed);
  } else {
    micromodulus = 12.0*d_young_modulus/(M_PI*rad_cubed*horizonSize);
  }
  return micromodulus;
}

//---------------------------------------------------------------------------------
double 
Material::computeCriticalStrain(const double& horizonSize) const
{
  double critical_strain = 0.0;
  if (d_micro_modulus_model == MicroModulusModel::Constant) {
    // For constant micro-modulus
    critical_strain = std::sqrt(5.0*d_fracture_energy/(6.0*horizonSize*d_young_modulus));
  } else if (d_micro_modulus_model == MicroModulusModel::Conical) {
    // For conical micro-modulus
    critical_strain = std::sqrt(5.0*M_PI*d_fracture_energy/(8.0*d_young_modulus*horizonSize));
  } 
  return critical_strain;
}

//---------------------------------------------------------------------------------
double 
Material::computeDamageFactor(const double& damageIndex) const
{
  return d_damage_model->computeDamageFactor(damageIndex);
}

namespace Matiti {

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
