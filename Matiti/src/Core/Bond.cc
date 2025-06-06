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

#include <Core/Bond.h>
#include <Core/Node.h>
#include <MaterialModels/Material.h>
#include <Core/Exception.h>

#include <Containers/MaterialSPArray.h>

//#include <Pointers/WoodSP.h>
//#include <Containers/WoodSPArray.h> 
//#include <Woods/Wood.h>

#include <Pointers/DensitySP.h>
#include <Containers/DensitySPArray.h> 
#include <MaterialModels/Density.h>

#include <iostream>
#define _USE_MATH_DEFINE
#include <cmath>

using namespace Matiti;

Bond::Bond()
  :d_node1(0),d_node2(0),d_force(0.0,0.0,0.0),d_broken(false)
{
}

Bond::Bond(const NodeP node1, const NodeP node2)
  :d_node1(node1),d_node2(node2),d_mat(new Material()),d_force(0.0,0.0,0.0),d_broken(false)
{
  d_mat->clone(d_node1->material());
}

Bond::Bond(const NodeP node1, const NodeP node2, const Material* mat)
  :d_node1(node1),d_node2(node2),d_mat(new Material()),d_force(0.0,0.0,0.0),d_broken(false)
{
  d_mat->clone(mat);
}

Bond::Bond(const NodeP node1, const NodeP node2, const Material* mat1, const Material* mat2)
  :d_node1(node1),d_node2(node2),d_mat(new Material()),d_force(0.0,0.0,0.0),d_broken(false)
{
  d_mat->cloneAverage(mat1, mat2);
}

Bond::~Bond()
{
}

bool 
Bond::operator==(const Bond& bond) const
{
  return ((d_node1 == bond.d_node1 && d_node2 == bond.d_node2) ||
          (d_node1 == bond.d_node2 && d_node2 == bond.d_node1));
}

/**
 * Compute volume weighted internal force
 */

void
Bond::computeInternalForce()
{
//  std::cout << std::endl << "[" << d_node1->getID() << ", " << d_node2->getID() << "]   "
//                         << "Broken= " << d_broken << std::endl;
  if (d_broken) {
    d_force.reset();
  } else {

//    std::cout << "material name= " << d_mat->name() << std::endl;
    // Compute the bond internal force using the material model
    d_mat->computeForce(d_node1->position(), d_node2->position(),
                        d_node1->displacement(), d_node2->displacement(),
                        d_node1->horizonSize(), d_force);

    // Get the volumes associated with the two nodes
    double fam_volume = d_node2->volume();

    // Reduce volume if node2 is not fully within the horizon of current node.
    // **WARNING** Assuming a ball around node for calculating radius instead of a
    //             rectangular parallelepiped or a hex element as in EMUNE.
    //double node2_radius = d_node2->radius();
    double bond_length = (d_node2->position() - d_node1->position()).length();
    //double bond_length_plus_radius = bond_length + node2_radius;
    //double bond_length_minus_radius = bond_length - node2_radius;
    //double horizon_size = d_node1->horizonSize();
    double delta = d_node1->horizonSize();
    //double volume_frac = 1.0;
    //if (horizon_size < bond_length_minus_radius) {
      //volume_frac = 0.0;
    //} else if (horizon_size < bond_length_plus_radius) {

      // Compute the volume of the lens of intersection of the horizon and the ball around node 2 
      // (Source: http://mathworld.wolfram.com/Sphere-SphereIntersection.html)
      //double horizon_cap_height = (bond_length_plus_radius - horizon_size)*
                                  //(horizon_size - bond_length_minus_radius)/(2.0*bond_length);
      //double node2_cap_height = (horizon_size + bond_length_minus_radius)*
                                //(horizon_size - bond_length_minus_radius)/(2.0*bond_length);
      //double horizon_cap_volume = (M_PI/3.0)*horizon_cap_height*horizon_cap_height*(3.0*horizon_size - horizon_cap_height);
      //double node2_cap_volume = (M_PI/3.0)*node2_cap_height*node2_cap_height*(3.0*node2_radius - node2_cap_height);
      //double intersection_volume = horizon_cap_volume + node2_cap_volume;
      //volume_frac = intersection_volume/fam_volume;
    //} else {
      //volume_frac = 1.0;
    //}
    //fam_volume *= volume_frac;
//    std::cout << " endpoint position= " << d_node2->position() << " horizon size= " << d_node2->horizonSize() << std::endl;
//    std::cout << " volume fraction= " << volume_frac << " family volume= " << fam_volume << std::endl;
    
    // **TODO**  Compute family volume using element shapes surrounding a node
     double xi = bond_length;
     Array3 fam_interval = {{0.0, 0.0, 0.0}};
     //family_node->getInterval(fam_interval);
     fam_interval[0] = d_node2->getInterval()[0];
     fam_interval[1] = d_node2->getInterval()[1];
     fam_interval[2] = d_node2->getInterval()[2];  
  
     double fam_radij = 0.5*std::max(std::max(fam_interval[0], fam_interval[1]), fam_interval[2]);

     double volume_fac = 0.0;
     if (fam_radij > 0.0) {
       if (xi <= delta - fam_radij) {
         volume_fac = 1.0;
       } else if (xi <= delta + fam_radij) {
         volume_fac = (delta + fam_radij - xi)/(2.0*fam_radij);
       } else {
         volume_fac = 0.0;
       }
     } 
     fam_volume *= volume_fac;

     d_force *= fam_volume;

       
//    std::cout << std::endl  
//                           << " bond internal force 1= " << d_force << std::endl << std::endl;

    if (d_force.isnan()) {
      std::ostringstream out;
      out << "**ERROR**  Nan internal force" << d_force << " Bond = " << *this << " family vol = " << fam_volume;
      throw Exception(out.str(), __FILE__, __LINE__);
    }
  }
}



void
Bond::computeInternalForce(const MaterialSPArray& matList, const Vector3D& gridSize)
{
  if (d_broken) {
    d_force.reset();
  } else {

DensitySP density = matList.front()->getDensity();
WoodSP wood = matList.front()->getWood();    

//    std::cout << "material name= " << d_mat->name() << std::endl;
    // Compute the bond internal force using the material model
    d_mat->computeForce(d_node1->position(), d_node2->position(),
                        d_node1->displacement(), d_node2->displacement(),
                        d_node1->horizonSize(), density, wood, gridSize, d_force);

    // Get the volumes associated with the two nodes
    double fam_volume = d_node2->volume();

    // Reduce volume if node2 is not fully within the horizon of current node.
    // **WARNING** Assuming a ball around node for calculating radius instead of a
    //             rectangular parallelepiped or a hex element as in EMUNE.
    double node2_radius = d_node2->radius();
    double bond_length = (d_node2->position() - d_node1->position()).length();
    double bond_length_plus_radius = bond_length + node2_radius;
    double bond_length_minus_radius = bond_length - node2_radius;
    double horizon_size = d_node1->horizonSize();
    double volume_frac = 1.0;
    if (horizon_size < bond_length_minus_radius) {
      volume_frac = 0.0;
    } else if (horizon_size < bond_length_plus_radius) {

      // Compute the volume of the lens of intersection of the horizon and the ball around node 2 
      // (Source: http://mathworld.wolfram.com/Sphere-SphereIntersection.html)
      double horizon_cap_height = (bond_length_plus_radius - horizon_size)*
                                  (horizon_size - bond_length_minus_radius)/(2.0*bond_length);
      double node2_cap_height = (horizon_size + bond_length_minus_radius)*
                                (horizon_size - bond_length_minus_radius)/(2.0*bond_length);
      double horizon_cap_volume = (M_PI/3.0)*horizon_cap_height*horizon_cap_height*(3.0*horizon_size - horizon_cap_height);
      double node2_cap_volume = (M_PI/3.0)*node2_cap_height*node2_cap_height*(3.0*node2_radius - node2_cap_height);
      double intersection_volume = horizon_cap_volume + node2_cap_volume;
      volume_frac = intersection_volume/fam_volume;
    } else {
      volume_frac = 1.0;
    }
    fam_volume *= volume_frac;
    
    // **TODO**  Compute family volume using element shapes surrounding a node
    // double xi = bond_length_ref;
    // Array3 fam_interval = {{0.0, 0.0, 0.0}};
    // family_node->getInterval(fam_interval);
    // double fam_radij = 0.5*std::max(fam_interval[0], fam_interval[1]);

    // double volume_fac = 0.0;
    // if (fam_radij > 0.0) {
    //   if (xi <= delta - fam_radij) {
    //     volume_fac = 1.0;
    //   } else if (xi <= delta + fam_radij) {
    //     volume_fac = (delta + fam_radij - xi)/(2.0*fam_radij);
    //   } else {
    //     volume_fac = 0.0;
    //   }
    // } 
    // fam_volume *= volume_fac;

    d_force *= fam_volume;

    if (d_force.isnan()) {
      std::ostringstream out;
      out << "**ERROR**  Nan internal force" << d_force << " Bond = " << *this << " family vol = " << fam_volume;
      throw Exception(out.str(), __FILE__, __LINE__);
    }
  }
}

//-----------------------------------------------------------------------------
// Compute strain energy
double 
Bond::computeStrainEnergy() const
{
  return d_mat->strainEnergy()*(0.5*d_node2->volume());
}

//-----------------------------------------------------------------------------
// Compute volume weighted micro modulus
double 
Bond::computeMicroModulus() const
{
  return d_mat->microModulus()*d_node2->volume();
}

//-----------------------------------------------------------------------------
// Compute critical strain and flag bond as broken

bool 
Bond::checkAndFlagBrokenBond() 
{
//  std::cout << std::endl << " /////////////START OF BONDCHECK////////////////" << std::endl << std::endl;
  // Check if the bond is already broken
  if (d_node2->omit() || d_broken) return d_broken;
//  std::cout << std::endl << "omit= " << d_node2->omit() << "  Broken? " << d_broken << " Fly Node1? " << !d_node1->failureAllowed()
//                                                                                    << " Fly Node2? " << !d_node2->failureAllowed()
//                                                                                    << std::endl << std::endl;

  // Hack to prevent nodes flying off under the action of an external load
  if ((!d_node1->failureAllowed()) || (!d_node2->failureAllowed())) {
    d_broken = false;
    return d_broken;
  }
    
  // Compute critical stretch as a function of damage.
  double dmgij = std::max(d_node1->damageIndex(), d_node2->damageIndex());
  double damage_fac = d_mat->computeDamageFactor(dmgij);
  
/*  if (d_node1->getID() == 2) {
    std::cout << " node1 = " << d_node1->getID() << " damage index = " << d_node1->damageIndex()
               << " node2 = " << d_node2->getID() << " damage index = " << d_node2->damageIndex()
              << " damage fac = " << damage_fac << std::endl;
    Vector3D damage_stretch = d_mat->damageModel()->damageStretch();
    std::cout << "damage_stretch= (" << damage_stretch.x() << ", " << damage_stretch.y() << ", " << damage_stretch.z() << ") " << std::endl;
  }*/

  // Break bond if critical stretch exceeded.
  double critical_strain_cur;

  critical_strain_cur = d_mat->computeCriticalStrain(d_node1->horizonSize());

//  d_mat->computeCriticalStrain(d_node1->horizonSize());
  double ecr2 = critical_strain_cur*damage_fac;
  double str = d_mat->strain();
//  std::cout << "Critical Strain= " << critical_strain_cur << "  Damage Factor= " << damage_fac << std::endl;
//  std::cout << "Critical Strength= " << ecr2 << "  Strain= " << str << std::endl;
  if (str > ecr2) {
//    std::cout << "Breaking bond between nodes " << d_node1->getID() << " and " << d_node2->getID() 
//              << " critical strain = " << critical_strain_cur << " ecr2 = " << ecr2 << " strain = " << str << std::endl;
    d_broken = true;
  }

  //if (d_node1->getID() == 2) {
  // std::cout << "     crit_strain = " << critical_strain_cur << " scaled_crit_strain = " << ecr2
  //            << " strain = " << str << " broken = " << std::boolalpha << d_broken << std::endl;    
  //}
  return d_broken;
}


bool 
Bond::checkAndFlagBrokenBond(const MaterialSPArray& matList) 
{
  // Check if the bond is already broken
  if (d_node2->omit() || d_broken) return d_broken;

  // Hack to prevent nodes flying off under the action of an external load
  if ((!d_node1->failureAllowed()) || (!d_node2->failureAllowed())) {
    d_broken = false;
    return d_broken;
  }
    
  // Compute critical stretch as a function of damage.
  double dmgij = std::max(d_node1->damageIndex(), d_node2->damageIndex());
  double damage_fac = d_mat->computeDamageFactor(dmgij);
  
  //if (d_node1->getID() == 2) {
  //  std::cout << " node1 = " << d_node1->getID() << " damage index = " << d_node1->damageIndex()
  //             << " node2 = " << d_node2->getID() << " damage index = " << d_node2->damageIndex()
  //            << " damage fac = " << damage_fac << std::endl;
  //}

  // Break bond if critical stretch exceeded.
  double critical_strain_cur;

//  bool woodBond = ((d_mat->hasName() == true) && (d_mat->name() == "wood"));
  //std::cout << "woodBond= " << woodBond << std::endl;

//  if ((d_mat->hasName() == true) && (d_mat->name() == "wood"))
//   {
    DensitySP density = matList.front()->getDensity();
    WoodSP wood = matList.front()->getWood();
    bool earlywood_node_1 = d_mat->earlywoodPoint(d_node1->position(), density, wood);
    bool earlywood_node_2 = d_mat->earlywoodPoint(d_node2->position(), density, wood);
//    std::cout << "ring width= " << d_mat->getDensity()->ringWidth() << std::endl;
//   const WoodSPArray wood_array = d_mat->getWoodArray();
//   WoodSP wood = wood_array.front();
    critical_strain_cur = wood->computeCriticalStrain
                                        (d_node1, d_node2, earlywood_node_1, earlywood_node_2);
//   }
//  else
//   {
//    critical_strain_cur = d_mat->computeCriticalStrain(d_node1->horizonSize());
//   }
  d_mat->computeCriticalStrain(d_node1->horizonSize());
  double ecr2 = critical_strain_cur*damage_fac;
  double str = d_mat->strain();
  if (str > ecr2) {
//    std::cout << "Breaking bond between nodes " << d_node1->getID() << " and " << d_node2->getID() 
//              << " critical strain = " << critical_strain_cur << " ecr2 = " << ecr2 << " strain = " << str << std::endl;
    d_broken = true;
  }

  //if (d_node1->getID() == 2) {
  //  std::cout << "     crit_strain = " << critical_strain_cur << " scaled_crit_strain = " << ecr2
  //            << " strain = " << str << " broken = " << std::boolalpha << d_broken << std::endl;
  //}
  return d_broken;
}


namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Bond& bond)
  {
    //out.setf(std::ios::floatfield);
//    out.precision(3);
    out.precision(20);
    out << "Bond: [" << (bond.d_node1)->getID() << " - " << (bond.d_node2)->getID() 
        << "], broken = " << std::boolalpha << bond.d_broken 
        << ", force = " << bond.d_force
        << std::endl  
        << ", disp2 = " << (bond.d_node2)->displacement() 
        << ", disp1 = " <<  (bond.d_node1)->displacement()
        << std::endl 
        << ",material = " << (bond.d_mat)->id()
        << ", strain = " << (bond.d_mat)->strain() 
        << ", critical strain = " << (bond.d_mat)->criticalStrain() << std::endl
        << ",endpoit of Bond = (" << (bond.d_node2)->x() << ", " << (bond.d_node2)->y() << ", " << (bond.d_node2)->z() << ")"
        << std::endl;
    return out;
  }
}
