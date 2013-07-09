#include <Bond.h>
#include <Node.h>
#include <Material.h>
#include <Exception.h>

#include <iostream>

using namespace Matiti;

Bond::Bond()
  :d_node1(0),d_node2(0),d_force(0.0,0.0,0.0),d_broken(false)
{
}

Bond::Bond(const NodeP& node1, const NodeP& node2)
  :d_node1(node1),d_node2(node2),d_mat(new Material()),d_force(0.0,0.0,0.0),d_broken(false)
{
  d_mat->clone(d_node1->material());
}

Bond::Bond(const NodeP& node1, const NodeP& node2, const Material* mat)
  :d_node1(node1),d_node2(node2),d_mat(new Material()),d_force(0.0,0.0,0.0),d_broken(false)
{
  d_mat->clone(mat);
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
  if (d_broken) {
    d_force.reset();
  } else {
    d_mat->computeForce(d_node1->position(), d_node2->position(),
                        d_node1->displacement(), d_node2->displacement(),
                        d_node1->horizonSize(), d_force);
    double fam_volume = d_node2->volume();

    // **TODO**
    // Reduce volume if node2 is not fully within the horizon of current node.
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
  // Check if the bond is already broken
  if (d_node2->omit() || d_broken) {
    d_broken = true;
    return d_broken;
  }
    
  // Compute critical stretch as a function of damage.
  double dmgij = std::max(d_node1->damageIndex(), d_node2->damageIndex());
  double damage_fac = d_mat->computeDamageFactor(dmgij);
  //std::cout << " node1 = " << d_node1->getID() << " damage index = " << d_node1->damageIndex()
  //          << " node2 = " << d_node2->getID() << " damage index = " << d_node2->damageIndex()
  //          << " damage fac = " << damage_fac << std::endl;

  // Break bond if critical stretch exceeded.
  double critical_strain_cur = d_mat->computeCriticalStrain(d_node1->horizonSize());
  double ecr2 = critical_strain_cur*damage_fac;
  double str = d_mat->strain();
  if (str > ecr2) d_broken = true;
  //std::cout << "     crit_strain = " << critical_strain_cur << " scaled_crit_strain = " << ecr2
  //          << " strain = " << str << " broken = " << std::boolalpha << d_broken << std::endl;
  return d_broken;
}


namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Bond& bond)
  {
    out.setf(std::ios::floatfield);
    out.precision(3);
    out << "Bond: [" << (bond.d_node1)->getID() << " - " << (bond.d_node2)->getID() 
        << "], broken = " << std::boolalpha << bond.d_broken 
        << ", force = " << bond.d_force  
        << ", disp2 = " << (bond.d_node2)->displacement() 
        << " disp1 = " <<  (bond.d_node1)->displacement()
        << " material = " << (bond.d_mat)->id() << std::endl;
    return out;
  }
}
