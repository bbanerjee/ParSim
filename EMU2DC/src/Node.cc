#include <Node.h>

using namespace Emu2DC;

Node::Node()
{
  d_id = 0;
  d_mat_type = 0;
  d_horizon_size = 0.0;
  d_iflag = false;

  d_volume = 0.0;
  d_density = 0.0;
  d_young = 0.0;

  d_strain_energy = 0;
  d_damage_index = 0;

  d_pos = {{0.0, 0.0, 0.0}};
  d_disp = {{0.0, 0.0, 0.0}};
  d_veloc = {{0.0, 0.0, 0.0}};
  d_accel = {{0.0, 0.0, 0.0}};
  d_new_veloc = {{0.0, 0.0, 0.0}};
  d_new_disp = {{0.0, 0.0, 0.0}};
  d_old_disp = {{0.0, 0.0, 0.0}};
  d_force = {{0.0, 0.0, 0.0}};

  d_bc = 0;
    
  d_nnodeelements = 0;
  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const int id, const double xx, const double yy, const double zz, const int hanging)
 : d_id(id)
{
  d_pos = {{xx, yy, zz}};
  d_iflag = false;
  if (hanging == 1) d_iflag = true;
}

Node::~Node()
{
  //std::cout << "deleted node " << d_id << std::endl;
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
    out << "Node (" << node.d_id << ") = [" << node.d_pos[0] << " , " << node.d_pos[1] << " , " << node.d_pos[2] << "]" << " hanging = " << node.d_iflag;
    return out;
  }
}
