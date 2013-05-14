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

