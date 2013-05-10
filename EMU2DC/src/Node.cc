#include <Node.h>
#include <NodeP.h>
#include <NodePArray.h>
#include <Material.h>

using namespace Emu2DC;

Node::Node()
  : d_dimension(3), d_id(0), d_mat_type(0), d_horizon_size(0.0), d_omit(false), d_surfaceNode(false),
    d_volume(0.0), d_material(new Material())
{
  d_pos = {{0.0, 0.0, 0.0}};
  d_disp = {{0.0, 0.0, 0.0}};
  d_veloc = {{0.0, 0.0, 0.0}};
  d_accel = {{0.0, 0.0, 0.0}};
  d_new_veloc = {{0.0, 0.0, 0.0}};
  d_new_disp = {{0.0, 0.0, 0.0}};
  d_old_disp = {{0.0, 0.0, 0.0}};
  d_force = {{0.0, 0.0, 0.0}};

  d_adjacent_elements.reserve(10);
  d_neighbor_list.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const int id, const double xx, const double yy, const double zz, const int surfaceNode)
  : d_dimension(3), d_id(id), d_mat_type(0), d_horizon_size(0.0), d_omit(false),
    d_volume(0.0), d_material(new Material()) 
{
  d_surfaceNode = false;
  if (surfaceNode) d_surfaceNode = true;

  d_pos = {{xx, yy, zz}};
  d_disp = {{0.0, 0.0, 0.0}};
  d_veloc = {{0.0, 0.0, 0.0}};
  d_accel = {{0.0, 0.0, 0.0}};
  d_new_veloc = {{0.0, 0.0, 0.0}};
  d_new_disp = {{0.0, 0.0, 0.0}};
  d_old_disp = {{0.0, 0.0, 0.0}};
  d_force = {{0.0, 0.0, 0.0}};

  d_adjacent_elements.reserve(10);
  d_neighbor_list.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const Node& node)
  : d_dimension(node.d_dimension), d_id(node.d_id), d_mat_type(node.d_mat_type), 
    d_horizon_size(node.d_horizon_size), d_omit(node.d_omit), d_surfaceNode(node.d_surfaceNode),
    d_volume(node.d_volume), d_material(new Material()), 
    d_adjacent_elements(node.d_adjacent_elements),
    d_neighbor_list(node.d_neighbor_list), d_initial_family_size(node.d_initial_family_size),
    d_pos(node.d_pos),  d_disp(node.d_disp), d_veloc(node.d_veloc), d_accel(node.d_accel),
    d_new_veloc(node.d_new_veloc), d_new_disp(node.d_new_disp), d_old_disp(node.d_old_disp),
    d_force(node.d_force)
{
  //std::cout << "created node " << d_id << std::endl;
}

Node::~Node()
{
  //std::cout << "deleted node " << d_id << std::endl;
}

void 
Node::assignMaterial(const Material* mat)
{
  d_material->clone(mat);
}

bool 
Node::operator<(const Node& node) const
{
  return (d_id < node.d_id);
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Node& node)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Node (" << node.d_id << ") = [" << node.d_pos[0] << " , " << node.d_pos[1] << " , " 
                    << node.d_pos[2] << "]" ;
    return out;
  }
}
