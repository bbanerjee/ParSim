#include <Node.h>
#include <NodeP.h>
#include <NodePArray.h>
#include <Material.h>

using namespace Emu2DC;

Node::Node()
  : d_dimension(3), d_id(0), d_mat_type(0), d_horizon_size(0.0), d_omit(false), d_surfaceNode(false),
    d_volume(0.0), d_material(new Material()),
    d_pos(0.0, 0.0, 0.0), d_disp(0.0, 0.0, 0.0), d_veloc(0.0, 0.0, 0.0), d_accel(0.0, 0.0, 0.0),
    d_new_veloc(0.0, 0.0, 0.0), d_new_disp(0.0, 0.0, 0.0), d_old_disp(0.0, 0.0, 0.0),
    d_force(0.0, 0.0, 0.0)
{
  d_adjacent_elements.reserve(10);
  d_neighbor_list.reserve(40);
  d_bond_materials.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const int id, const double xx, const double yy, const double zz, const int surfaceNode)
  : d_dimension(3), d_id(id), d_mat_type(0), d_horizon_size(0.0), d_omit(false),
    d_volume(0.0), d_material(new Material()),
    d_pos(xx, yy, zz), d_disp(0.0, 0.0, 0.0), d_veloc(0.0, 0.0, 0.0), d_accel(0.0, 0.0, 0.0),
    d_new_veloc(0.0, 0.0, 0.0), d_new_disp(0.0, 0.0, 0.0), d_old_disp(0.0, 0.0, 0.0),
    d_force(0.0, 0.0, 0.0)
{
  d_surfaceNode = false;
  if (surfaceNode) d_surfaceNode = true;

  d_adjacent_elements.reserve(10);
  d_neighbor_list.reserve(40);
  d_bond_materials.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const Node& node)
  : d_dimension(node.d_dimension), d_id(node.d_id), d_mat_type(node.d_mat_type), 
    d_horizon_size(node.d_horizon_size), d_omit(node.d_omit), d_surfaceNode(node.d_surfaceNode),
    d_volume(node.d_volume), d_adjacent_elements(node.d_adjacent_elements),
    d_neighbor_list(node.d_neighbor_list), 
    d_initial_family_size(node.d_initial_family_size),
    d_pos(node.d_pos),  d_disp(node.d_disp), d_veloc(node.d_veloc), d_accel(node.d_accel),
    d_new_veloc(node.d_new_veloc), d_new_disp(node.d_new_disp), d_old_disp(node.d_old_disp),
    d_force(node.d_force)
{
  d_material = MaterialUP(new Material()),
  d_material->clone(node.d_material.get());
  for (auto iter = node.d_bond_materials.begin(); iter != node.d_bond_materials.end(); ++iter) {
    MaterialUP mat_up(new Material());
    mat_up->clone((*iter).get());
    d_bond_materials.push_back(std::move(mat_up));
  }
  //std::cout << "created node " << d_id << std::endl;
}

Node::~Node()
{
  //std::cout << "deleted node " << d_id << std::endl;
}


void 
Node::setFamily(const NodePArray& fam) 
{
  d_neighbor_list = fam;

  // All the bonds have the same material at the beginning
  // **TO DO** Will have to think about what to do when the family changes
  for (auto iter = fam.begin(); iter != fam.end(); ++iter) {
    MaterialUP mat_up(new Material());
    mat_up->clone(d_material.get());
    d_bond_materials.push_back(std::move(mat_up));
  }
}

void
Node::assignMaterial(const Material* mat)
{
  d_material->clone(mat);
}

void 
Node::computeInitialDisplacement(const Vector3D& initVel, double delT)
{
  // Assume displacement is 0 at t = 0
  d_veloc = initVel;
  d_old_disp = d_veloc*delT;
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Node& node)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Node (" << node.d_id << ") = [" << node.d_pos.x() << " , " << node.d_pos.y() << " , " 
                    << node.d_pos.z() << "]" ;
    return out;
  }
}
