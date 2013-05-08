#include <VelocityBC.h> 
#include <Node.h>
#include <NodeP.h>
#include <Geometry/Point3D.h>
#include <Exception.h>
#include <ProblemSpecUtil.h>

#include <Core/ProblemSpec/ProblemSpec.h>


#include <vector>
#include <iostream>
#include <string>

using namespace Emu2DC;
  
VelocityBC::VelocityBC() {}
VelocityBC::~VelocityBC() {}

// Initialize velocity BC objects
void
VelocityBC::initialize(Uintah::ProblemSpecP& ps)
{
  // If ps is null return
  if (!(ps)) return;

  // Read the face on which this BC is to be applied
  std::string face_id;
  ps->require("face", face_id);
  if (face_id == "x-") d_face = FaceType::Xminus;
  else if (face_id == "x-") d_face = FaceType::Xplus;
  else if (face_id == "y-") d_face = FaceType::Yminus;
  else if (face_id == "y+") d_face = FaceType::Yplus;
  else if (face_id == "z-") d_face = FaceType::Zminus;
  else if (face_id == "z+") d_face = FaceType::Zplus;
  else {
    std::ostringstream out;
    out << "**ERROR** Unknown face id " << face_id << " for VelocityBC application." << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Read the boundary condition type
  std::string bc_type;
  ps->require("bc", bc_type);
  if (bc_type == "symmetry") d_bc = BCType::Symmetry;
  else if (bc_type == "periodic") d_bc = BCType::Periodic;
  else if (bc_type == "wall") d_bc = BCType::Wall;
  else if (bc_type == "outlet") d_bc = BCType::Outlet;
  else {
    std::ostringstream out;
    out << "**ERROR** Unknown BC type " << bc_type << " for VelocityBC." << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Read the points that describe the areas over which special velocity BCs are to be applied
  // and store in d_boundary
  for (Uintah::ProblemSpecP area_ps = ps->findBlock("Area"); area_ps != 0;
         area_ps = area_ps->findNextBlock("Area")) {

    // Read the BC type
    area_ps->require("bc", bc_type);
    if (bc_type == "symmetry") d_area_bc.push_back(BCType::Symmetry);
    else if (bc_type == "periodic") d_area_bc.push_back(BCType::Periodic);
    else if (bc_type == "wall") d_area_bc.push_back(BCType::Wall);
    else if (bc_type == "outlet") d_area_bc.push_back(BCType::Outlet);
    else {
      std::ostringstream out;
      out << "**ERROR** Unknown BC type " << bc_type << " for special VelocityBC over area." << std::endl;
      throw Exception(out.str(), __FILE__, __LINE__);
    }

    // Check for the <Boundary> block
    Uintah::ProblemSpecP boundary_ps = area_ps->findBlock("Boundary");
    if (boundary_ps) {

      // Get the boundary points and create linestring
      Uintah::ProblemSpecP node = boundary_ps->findBlock("point");
      if (node) { 
        Polygon3D boundary;
        Emu2DC_ProblemSpecUtil::readBoundary(node, boundary);
        d_area_boundary.push_back(boundary);
      } 
    }

    // Make sure that each area has a boundary and a bc
    if (d_area_bc.size() != d_area_boundary.size()) {
      std::ostringstream out;
      out << "**ERROR** Special VelocityBC area boundaries and bcs do not match.." << std::endl;
      throw Exception(out.str(), __FILE__, __LINE__);
    }
  }
}

// Apply boundary conditions to a set of nodes on the boundary of a body that have moved outside the 
// boundary of the computational domain
void
VelocityBC::applyVelocityBC(NodePArray& nodes, const Point3D& domain_min, const Point3D& domain_max) const
{
  // Loop through the nodes
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {

    // Apply bcs on to boundary nodes
    NodeP cur_node = *node_iter;

    // Look only at nodes that are on the surface of the body
    if (!(cur_node->onSurface())) continue; 

    // Get the old and current position of the node
    const Array3& pos = cur_node->position();
    const Array3& disp = cur_node->displacement();
    Point3D cur_pos(pos[0]+disp[0], pos[1]+disp[1], pos[2]+disp[2]);

    // If the node is inside the domain do nothing
    if (insideDomain(cur_pos, domain_min, domain_max)) continue;

    // Step 1: First look at the full face of the domain
    
    // The node is outside the domain push it back or let it go (apply BCs)
  }
}

bool 
VelocityBC::insideDomain(const Point3D& node_pos, const Point3D& domain_min, const Point3D& domain_max) const
{
  return node_pos.x() < domain_min.x() || node_pos.x() > domain_max.x() || 
      node_pos.y() < domain_min.y() || node_pos.y() > domain_max.y() || 
      node_pos.z() < domain_min.z() || node_pos.z() > domain_max.z();
}


