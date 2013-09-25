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

using namespace Matiti;
  
VelocityBC::VelocityBC()
  : d_restitution(0.0)
{
}

VelocityBC::~VelocityBC() 
{
}

// Initialize velocity BC objects
void
VelocityBC::initialize(Uintah::ProblemSpecP& ps)
{
  // If ps is null return
  if (!(ps)) return;

  // Get the coefficient of restitution
  d_restitution = 0.0;
  ps->get("coeff_of_restitution", d_restitution);

  // Read the face on which this BC is to be applied
  std::string face_id;
  ps->require("face", face_id);
  if (face_id == "x-") d_face = FaceType::Xminus;
  else if (face_id == "x+") d_face = FaceType::Xplus;
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
        Matiti_ProblemSpecUtil::readBoundary(node, boundary);
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
// **WARNING** None of the specialized boundary conditions are implemented yet. **TODO**
void
VelocityBC::apply(NodeP& node, 
                  const Point3D& hitPoint, 
                  const Point3D& domainMin,
                  const Point3D& domainMax) const
{
  // **WARNING** Just doing outlet and wall bcs for now
  if (d_bc == BCType::Symmetry) {
    std::ostringstream out;
    out << "**ERROR** Symmetry BCs not implemented yet.";
    throw Exception(out.str(), __FILE__, __LINE__);
  } else if (d_bc == BCType::Periodic) {
    std::ostringstream out;
    out << "**ERROR** Periodic BCs not implemented yet.";
    throw Exception(out.str(), __FILE__, __LINE__);
  } else if (d_bc == BCType::Outlet) {
    // Omit the current node  (?) **WARNING** maybe not
    node->omit(true);
  } 
  // Wall BCs below

  // Step 1: First which face is closest to intersection pt of ray through position
  // along the displacement direction

  // Set tolerance
  Vector3D domainRange = domainMax - domainMin;
  double tol = 0.001*std::abs(domainRange.min());

  // Find distance of hit point to domain face.  Values close to zero indicate that the
  // hit point is on a particular face.
  Vector3D dist_to_max = domainMax - hitPoint;
  Vector3D dist_to_min = hitPoint - domainMin;
    
  // The node is outside the domain push it back or let it go (apply BCs)
  // Otherwise apply reflective boundary condition
  // 1) Push back the node to the boundary along the boundary normal
  //    i.e. update displacement
  // 2) Reverse the normal velocity component 

  switch (d_face) {

    // Hit point is on x-
    case FaceType::Xminus:
      if (std::abs(dist_to_min.x()) < tol) {
        Vector3D normal(-1.0, 0.0, 0.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    // Hit point is on x+
    case FaceType::Xplus:
      if (std::abs(dist_to_max.x()) < tol) {
        Vector3D normal(1.0, 0.0, 0.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    // Hit point is on y-
    case FaceType::Yminus:
      if (std::abs(dist_to_min.y()) < tol) {
        Vector3D normal(0.0, -1.0, 0.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    // Hit point is on y+
    case FaceType::Yplus:
      if (std::abs(dist_to_max.y()) < tol) {
        Vector3D normal(0.0, 1.0, 0.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    // Hit point is on z-
    case FaceType::Zminus:
      if (std::abs(dist_to_min.z()) < tol) {
        Vector3D normal(0.0, 0.0, -1.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    // Hit point is on z+
    case FaceType::Zplus:
      if (std::abs(dist_to_max.z()) < tol) {
        Vector3D normal(0.0, 0.0, 1.0);
        updateVelocityAndPosition(node, hitPoint, normal); 
      }
      break;

    default:
      std::ostringstream out;
      out << "**ERROR** Facetype unknown. Aborting." ;
      throw Exception(out.str(), __FILE__, __LINE__); 
  } // end switch
}

// Update the velocity using reflection BCs
void
VelocityBC::updateVelocityAndPosition(NodeP& node, 
                                      const Point3D& hitPoint,
                                      const Vector3D& normal) const 
{
  Vector3D vel_new = node->newVelocity();
  double impluse = vel_new.dot(normal);
  std::cout << " normal = " << normal << " vel_new = " << vel_new - normal*((1.0+d_restitution)*impluse) << std::endl;
  node->newVelocity(vel_new - normal*((1.0+d_restitution)*impluse));
  node->newDisplacement(hitPoint - node->position());
} 
