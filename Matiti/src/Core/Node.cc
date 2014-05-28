#include <Core/Node.h>
#include <Core/Bond.h>
#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Containers/MaterialSPArray.h>
#include <MaterialModels/Material.h>
#include <algorithm>

#include <Core/Exception.h>

using namespace Matiti;

Node::Node()
  : d_id(0), d_mat_type(0), d_horizon_size(0.0), d_omit(false), d_surfaceNode(false),
    d_area(0.0), d_volume(0.0), d_radius(0.0), d_material(new Material()),
    d_pos(0.0, 0.0, 0.0), d_disp(0.0, 0.0, 0.0), d_disp_new(0.0, 0.0, 0.0), 
    d_vel(0.0, 0.0, 0.0), d_vel_new(0.0, 0.0, 0.0), d_accel(0.0, 0.0, 0.0),
    d_int_force(0.0, 0.0, 0.0), d_ext_force(0.0, 0.0, 0.0), d_damage_index(0.0), d_allow_failure(true)
{
  d_adjacent_elements.reserve(10);
  d_bonds.reserve(40);
  // d_neighbor_list.reserve(40);
  // d_bond_materials.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const int id, const double xx, const double yy, const double zz, const int surfaceNode)
  : d_id(id), d_mat_type(0), d_horizon_size(0.0), d_omit(false),
    d_area(0.0), d_volume(0.0), d_radius(0.0), d_material(new Material()),
    d_pos(xx, yy, zz), d_disp(0.0, 0.0, 0.0),  d_disp_new(0.0, 0.0, 0.0), 
    d_vel(0.0, 0.0, 0.0), d_vel_new(0.0, 0.0, 0.0), d_accel(0.0, 0.0, 0.0),
    d_int_force(0.0, 0.0, 0.0), d_ext_force(0.0, 0.0, 0.0), d_damage_index(0.0), d_allow_failure(true)
{
  d_surfaceNode = false;
  if (surfaceNode) d_surfaceNode = true;

  d_adjacent_elements.reserve(10);
  d_bonds.reserve(40);
  // d_neighbor_list.reserve(40);
  // d_bond_materials.reserve(40);

  //std::cout << "created node " << d_id << std::endl;
}

Node::Node(const Node& node)
  : d_id(node.d_id), d_mat_type(node.d_mat_type), 
    d_horizon_size(node.d_horizon_size), d_omit(node.d_omit), d_surfaceNode(node.d_surfaceNode),
    d_area(node.d_area), d_volume(node.d_volume), d_radius(node.d_radius), d_material(new Material()),
    d_adjacent_elements(node.d_adjacent_elements),
    // d_neighbor_list(node.d_neighbor_list), 
    d_bonds(node.d_bonds),
    d_initial_family_size(node.d_initial_family_size),
    d_pos(node.d_pos),  d_disp(node.d_disp), d_disp_new(node.d_disp_new), 
    d_vel(node.d_vel), d_vel_new(node.d_vel_new), d_accel(node.d_accel),
    d_int_force(node.d_int_force), d_ext_force(node.d_ext_force), d_damage_index(node.d_damage_index), 
    d_allow_failure(node.d_allow_failure)
{
  // d_material = MaterialUP(new Material()),
  d_material->clone(node.d_material.get());
  // for (auto iter = node.d_bond_materials.begin(); iter != node.d_bond_materials.end(); ++iter) {
  //   MaterialUP mat_up(new Material());
  //   mat_up->clone((*iter).get());
  //   d_bond_materials.push_back(std::move(mat_up));
  // }
  //std::cout << "created node " << d_id << std::endl;
}

Node::~Node()
{
  //std::cout << "deleted node " << d_id << std::endl;
}


// void 
// Node::setFamily(const NodePArray& fam) 
// {
//   d_neighbor_list = fam;
// 
//   // All the bonds have the same material at the beginning
//   // **TO DO** Will have to think about what to do when the family changes
//   for (auto iter = fam.begin(); iter != fam.end(); ++iter) {
//     MaterialUP mat_up(new Material());
//     mat_up->clone(d_material.get());
//     d_bond_materials.push_back(std::move(mat_up));
//   }
// }
// 

double
Node::computeStableTimestep(const double& factor) const
{
  // If the family is zero return large double 
  if (d_bonds.size() == 0) return 1.0e16;
  
  // Loop over the family of current node mi.
  double density = d_material->density();
  double denom = 0.0;
  for (auto family_iter = d_bonds.begin(); family_iter != d_bonds.end(); family_iter++) {
    BondP bond = *family_iter;
    double micromodulus = bond->computeMicroModulus();
    denom += micromodulus;
  }
  return factor*d_horizon_size*std::sqrt(2.0*density/denom);
}

void 
Node::computeInitialDisplacement(const Vector3D& initVel, double delT)
{
  // Assume displacement is 0 at t = 0
  d_vel = initVel;
  d_disp = d_vel*delT;
}

//-------------------------------------------------------------------------
// Find and erase broken bonds and update the damage index
//-------------------------------------------------------------------------
void Node::findAndDeleteBrokenBonds()
{
  // Check if bond strain exceeds critical strain and remove if true
  auto lambda_func =
    [&](const BondP& bond)
    {
      return bond->checkAndFlagBrokenBond();
    };
  d_bonds.erase(std::remove_if(d_bonds.begin(), d_bonds.end(), lambda_func), d_bonds.end());

  // Update the damage index
  //updateDamageIndex();
}

void Node::findAndDeleteBrokenBonds(const MaterialSPArray& matList)
{
  // Check if bond strain exceeds critical strain and remove if true
  auto lambda_func =
    [&](const BondP& bond)
    {
      return bond->checkAndFlagBrokenBond(matList);
    };
  d_bonds.erase(std::remove_if(d_bonds.begin(), d_bonds.end(), lambda_func), d_bonds.end());

  // Update the damage index
  //updateDamageIndex();
}

void 
Node::updateDamageIndex()
{
  int num_bonds_init = d_initial_family_size;
  int num_bonds_cur = currentFamilySize();
  //std::cout << "In " << __FILE__ << " line " << __LINE__ << " Node = " << d_id << std::endl;
  //std::cout << " initial bonds = " << num_bonds_init << std::endl;
  //std::cout << " current bonds = " << num_bonds_cur << std::endl;
  if (!(num_bonds_init > 0)) {
    std::ostringstream out;
    out << "**ERROR** Number of initial bonds is zero for node " << d_id;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
  d_damage_index = 1.0 - (double) num_bonds_cur/(double) num_bonds_init;
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Node& node)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);

    // Print node position
    out << "    Node (" << node.d_id << ") = [" << node.d_pos.x() << " , " << node.d_pos.y() << " , " 
                    << node.d_pos.z() << "]" 
        << std::boolalpha << " on surface = " << node.d_surfaceNode << std::endl;

    // Print adjacent elements
    out << "      Adjacent elements = " ;
    for (auto iter = (node.d_adjacent_elements).begin();
              iter != (node.d_adjacent_elements).end();
              ++iter) {
      out << "         " << *iter;
    }
    out << std::endl;
      
    // Print bonds
    out << "      Bonds = " << std::endl;
    for (auto iter = (node.d_bonds).begin();
              iter != (node.d_bonds).end();
              ++iter) {
      out << "         " << *(*iter);
    }
    out << std::endl;

    // Print volume and area
    out << "      Volume = " << node.d_volume << " radius = " << node.d_radius
        << " Area = " << node.d_area << std::endl;
    // Print family size
    out << "      Initial family size = " << node.d_initial_family_size 
        << " Damage index = " << node.d_damage_index << std::endl;

    // Print vectors
    out << "      Displacement: Old = " << node.d_disp << " New = " << node.d_disp_new << std::endl;
    out << "      Velocity: Old =" << node.d_vel << " New = " << node.d_vel_new << std::endl;
    out << "      Internal force = " << node.d_int_force << std::endl;
    out << "      External force = " << node.d_ext_force << std::endl;

    return out;
  }
}
