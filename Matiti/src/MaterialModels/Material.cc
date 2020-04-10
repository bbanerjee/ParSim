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

#include <MaterialModels/Material.h> 
#include <MaterialModels/DamageModelFactory.h> 
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Matiti;

Material::Material()
  : d_id(0), d_have_name(false), d_name(""), d_density(0.0), d_young_modulus(0.0), d_fracture_energy(0.0),
    d_micro_modulus(0.0), d_strain(0.0), d_critical_strain(0.0), d_strain_energy(0.0), d_ring(0.0), d_earlywood_fraction(0.0),
    d_node_density(new Density()), d_wood(new Wood())
{
  d_damage_model = 0;
  d_micro_modulus_model = MicroModulusModel::Constant;
  d_coeffs.reserve(6);
}

Material::Material(const Material& mat)
  : d_id(mat.d_id), d_have_name(mat.d_have_name), d_name(mat.d_name), d_density(mat.d_density), 
    d_young_modulus(mat.d_young_modulus), d_fracture_energy(mat.d_fracture_energy), 
    d_micro_modulus(mat.d_micro_modulus), d_strain(mat.d_strain), d_critical_strain(mat.d_critical_strain), 
    d_strain_energy(mat.d_strain_energy), d_ring(mat.d_ring), d_node_density(new Density()), d_wood(new Wood())
{
  d_micro_modulus_model = mat.d_micro_modulus_model;
  d_damage_model = (mat.d_damage_model)->clone();
  d_node_density->clone(mat.d_node_density);
  d_wood->clone(mat.d_wood);
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
    d_critical_strain = mat.d_critical_strain;
    d_strain_energy = mat.d_strain_energy;
    d_ring = mat.d_ring;
    d_damage_model = (mat.d_damage_model)->clone();
    d_node_density->clone(mat.d_node_density);
    d_wood->clone(mat.d_wood);
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
    d_critical_strain = mat->d_critical_strain;
    d_strain_energy = mat->d_strain_energy;
    d_ring = mat->d_ring;
    d_damage_model = (mat->d_damage_model)->clone();
    d_node_density->clone(mat->d_node_density);
    d_wood->clone(mat->d_wood);
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
    d_critical_strain = mat->d_critical_strain;
    d_strain_energy = mat->d_strain_energy;
    d_damage_model = (mat->d_damage_model)->clone();
    d_damage_model->setVariation(randomNum, coeffOfVar);
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
    d_critical_strain = mat1->d_critical_strain;
    d_strain_energy = mat1->d_strain_energy;
    d_damage_model = (mat1->d_damage_model)->clone();
    d_damage_model->setAverage(mat1->d_damage_model.get(), mat2->d_damage_model.get());
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
   
    if ((d_have_name == true) && (d_name == "wood")) {
      WoodSP wood = std::make_shared<Wood>();
      wood->initialize(ps);
      setWood(wood);
      d_node_density = std::make_shared<Density>();
      d_node_density->initialize(ps);
    }
  
    else
   {
    ps->require("young_modulus", d_young_modulus);
    ps->require("fracture_energy", d_fracture_energy);

    // The mass density can be homogeneous or hetergeneous.
    Uintah::ProblemSpecP dens_ps = ps->findBlock("Density"); 
    if (!dens_ps) {
      std::ostringstream out;
      out << "**ERROR** No density information found for the material";
      throw Exception(out.str(), __FILE__, __LINE__);
    }
  
    if (!dens_ps->getAttribute("type", d_density_type)) {
       std::ostringstream out;
       out << "**ERROR** Density does not have type information" << d_density_type;
      throw Exception(out.str(), __FILE__, __LINE__);
    }
  
    if (d_density_type == "homogeneous") {
        dens_ps->require("density", d_density);
    } else if (d_density_type == "heterogeneous") {
        d_node_density = std::make_shared<Density>();
        d_node_density->initialize(dens_ps);
    } else {
        std::ostringstream out;
        out << "**ERROR** Unknown density type" << d_density_type;
        throw Exception(out.str(), __FILE__, __LINE__);
  }
 
 }

  // Get the damage model
  Uintah::ProblemSpecP damage_ps = ps->findBlock("DamageModel");
  d_damage_model = DamageModelFactory::create(damage_ps);
  d_damage_model->initialize(damage_ps);

}

void 
Material::initialize(int id, MicroModulusModel model, 
                     double density, double modulus, double fractureEnergy)
{
  d_id = id;
  d_micro_modulus_model = model;
  d_density = density;
  d_young_modulus = modulus;
  d_fracture_energy = fractureEnergy;
}

//---------------------------------------------------------------------------------
// This is an elastic model.  **TODO** For more complex materials use factory concept.
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

  // Compute initial and current bond lengths
  double bond_length_init = xi.length();
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
/**********************
 * for 2D simulation

  if (xi.z()!=0) {
    force.x(0);
    force.y(0);
    force.z(0);
  }
  else {
 **********************/
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
//  std::cout << " xp = " << familyPos << " xi = " << nodePos 
//            << " up = " << familyDisp << " ui = " << nodeDisp 
//            << " force = " << force << " strain = " << d_strain << std::endl;

/***********
 * for 2D simulation
 *}  //end of else
 ***********/
} 


void 
Material::computeForce(const Point3D& nodePos,
                       const Point3D& familyPos,
                       const Vector3D& nodeDisp,
                       const Vector3D& familyDisp,
                       const double& horizonSize,
                       const DensitySP density,
                       const WoodSP wood,
                       const Vector3D& gridSize,
                       Vector3D& force)
{
  // Compute new bond vector
  Vector3D xi = familyPos - nodePos;  
  Vector3D eta = familyDisp - nodeDisp;
  Vector3D pos_new = xi + eta;

  // Compute initial and current bond lengths
  double bond_length_init = xi.length();
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
      bool fiber_bond = ((xi.y() == 0));// && (xi.z() == 0));
      bool earlywood_node = earlywoodPoint(nodePos, density, wood);
      bool earlywood_family = earlywoodPoint(familyPos, density, wood);
//      std::cout << "fiber_bond= " << fiber_bond << "  earlywood_node= " << earlywood_node
//                                                << "  earlywood_family= " << earlywood_family << std::endl
//                                                << "  horizon size= " << horizonSize << " bond_length= "  
//                                                << bond_length_init << std::endl;  
      d_micro_modulus = wood->computeMicroModulus(bond_length_init, horizonSize, fiber_bond,
                                                    earlywood_node, earlywood_family, gridSize); 
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
//  double rad_cubed = horizonSize*horizonSize*horizonSize;
  double rad_cubed = pow(horizonSize, 3);

  // Assume Poisson's ratio = nu = 0.25
  if (d_micro_modulus_model == MicroModulusModel::Conical) {
    micromodulus = 32.0*d_young_modulus*(1.0-bondLength/horizonSize)/(M_PI*rad_cubed);
  } else {
    //micromodulus = 12.0*d_young_modulus/(M_PI*rad_cubed*horizonSize);
    //double num_bonds = 100;
    //micromodulus = 3.0*d_young_modulus/(horizonSize*num_bonds);
    micromodulus = 13.5*d_young_modulus/(M_PI*rad_cubed); // Test with 2d micromodulus
  } 
//  std::cout << std::setprecision(8) << "1.56^3= " << 1.56*1.56*1.56 << "= " << pow(1.56, 3) << std::endl;
//  std::cout << "horizon size= " << horizonSize << " size_cube= " << rad_cubed << " Pi= " << M_PI << std::endl;
//  std::cout << std::setprecision(12) << "Distance = " << bondLength << " Micro modulus = " << micromodulus << " modulus= " << d_young_modulus << std::endl;
  return micromodulus;
}

//---------------------------------------------------------------------------------
double 
Material::computeCriticalStrain(const double& horizonSize) 
{
  double critical_strain = 0.0;
  if (d_micro_modulus_model == MicroModulusModel::Constant) {
    // For constant micro-modulus
     critical_strain = std::sqrt(5.0*d_fracture_energy/(6.0*horizonSize*d_young_modulus));

   /*********** for 2D simulation
    *critical_strain = (2.0/3.0)*std::sqrt((2.0*M_PI*d_fracture_energy)/(3.0*horizonSize*d_young_modulus));
    ***********/

//     critical_strain = 0.001;
    //double num_bonds = 100;
    //critical_strain = std::sqrt(2.0*d_fracture_energy/d_young_modulus)*num_bonds;
  } else if (d_micro_modulus_model == MicroModulusModel::Conical) {
    // For conical micro-modulus
    critical_strain = std::sqrt(5.0*M_PI*d_fracture_energy/(8.0*d_young_modulus*horizonSize));
  }
  setCriticalStrain(critical_strain);
//  std::cout << "critical strain= " << d_critical_strain << "= " << criticalStrain() << "= " << critical_strain << std::endl; 
  return critical_strain;
}

//------------------------------------------------------------------------------------

bool 
Material::earlywoodPoint(const Point3D& xi, const DensitySP density, const WoodSP Wood)
{ 
 double ringWidth = density->ringWidth();
 double factor = Wood->earlywoodFraction();
//  std::cout << "ring width= " << ringWidth << " factor= " << factor << " y pos= " << xi.y() << std::endl;
//  std::cout << "name of material= " << d_name << std::endl;
  double pos = density->remind(ringWidth, xi.y()-2);
//  std::cout << "period pos= " << pos << std::endl;
  return ((pos >= 0.0) && (pos < ringWidth*factor));
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
    //out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Material model: id = " << mat.d_id << std::endl;
    out << "  Have name = " << mat.d_have_name << " Name = " << mat.d_name << std::endl;
    out << "  Density = " << mat.d_density << " Young = " << mat.d_young_modulus
        << " Fracture enrgy = " << mat.d_fracture_energy << std::endl;
    out << mat.damageModel();
    return out;
  }
}
