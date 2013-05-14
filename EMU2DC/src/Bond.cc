#include <Bond.h>
#include <Node.h>
#include <Material.h>

using namespace Emu2DC;

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
// Compute micro modulus
double 
Bond::computeMicroModulus() const
{
  return d_mat->microModulus()*(d_node2->volume()/d_mat->density());
}
