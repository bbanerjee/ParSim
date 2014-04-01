#include <BoundaryConditions/DisplacementBC.h> 
#include <Geometry/Vector3D.h> 
#include <Core/Node.h>
#include <Pointers/NodeP.h>
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>
#include <iostream>

using namespace Matiti;
  
DisplacementBC::DisplacementBC() 
  : d_x_flag(false), d_y_flag(false), d_z_flag(false),
    d_x_value(0.0), d_y_value(0.0), d_z_value(0.0)
{
}

DisplacementBC::~DisplacementBC() 
{
}

// Initialize displacement boundary conditions on nodes
// TODO: Implement more general displacement BCs
void
DisplacementBC::initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes)
{
  // If ps is null return
  if (!(ps)) return;

  // Read the type of BC  (only symmetry bcs allowed for now)
  // symmetry x => u_1 = 0
  // symmetry y => u_2 = 0
  // symmetry z => u_3 = 0
  bool x_symmetry = false;
  bool y_symmetry = false;
  bool z_symmetry = false;
  ps->get("xsymmetry", x_symmetry);
  ps->get("ysymmetry", y_symmetry);
  ps->get("zsymmetry", z_symmetry);
  if (x_symmetry) {
    d_x_flag = true;
  } 
  if (y_symmetry) {
    d_y_flag = true;
  } 
  if (z_symmetry) {
    d_z_flag = true;
  } 
  if (!x_symmetry && !y_symmetry && !z_symmetry) {
    std::ostringstream out;
    out << "**ERROR** Displacement BC must be either xsymmetry, ysymmetry, or zsymmetry." 
        << " All symmetry BCs cannot be false." << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Read the region of application of the displaceemnt BC
  SCIRun::Vector box_min(0.0, 0.0, 0.0), box_max(0.0, 0.0, 0.0);
  ps->require("box_min", box_min);
  ps->require("box_max", box_max);

  // Make sure the box is OK
  for (int ii = 0; ii < 3; ii++) {
    if (box_max[ii] < box_min[ii]) {
      double temp = box_max[ii];
      box_max[ii] = box_min[ii];
      box_min[ii] = temp;
    }

  }
  if (box_min.x() == box_max.x() || box_min.y() == box_max.y() || box_min.z() == box_max.z()) {
    std::ostringstream out;
    out << "**ERROR** Box region for application of displacement BC is degenerate" << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Set the displacement BC flags at each node
  initializeDispBCSurfaceNodes(box_min, box_max, nodes);
}

void 
DisplacementBC::initializeDispBCSurfaceNodes(const SCIRun::Vector& boxMin, 
                                             const SCIRun::Vector& boxMax,
                                             NodePArray& nodes)
{
  std::cout << "Finding diplacement BC nodes" << std::endl;
  for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
    NodeP node = *iter;
    //if (!(node->onSurface())) continue;

    const Point3D& pos = node->position();
    if (pos.x() >= boxMin.x() && pos.y() >= boxMin.y() && pos.z() >= boxMin.z() && 
        pos.x() <= boxMax.x() && pos.y() <= boxMax.y() && pos.z() <= boxMax.z()) {
      d_surface_nodes.push_back(node);
    }
  }
} 

// This applies displacement BCs (and the related velocity BCs) 
// **WARNING** Make strong assumptions:
//   1) The coordinate system is the global one
//   2) Only symmetry BCs are allowed
void 
DisplacementBC::applyDisplacementBC()
{
  // Loop through surface nodes at which BC is to be applied
  for (auto iter = d_surface_nodes.begin(); iter != d_surface_nodes.end(); ++iter) {
    NodeP node = *iter;

    // Assign displacements, velocity, internal force
    if (d_x_flag) {
      node->xDisplacement(d_x_value);
      node->xVelocity(d_x_value);
      node->xNewDisplacement(d_x_value);
      node->xNewVelocity(d_x_value);
      node->xInternalForce(d_x_value);
    }
    if (d_y_flag) {
      node->yDisplacement(d_y_value);
      node->yVelocity(d_y_value);
      node->yNewDisplacement(d_y_value);
      node->yNewVelocity(d_y_value);
      node->yInternalForce(d_y_value);
    }
    if (d_z_flag) {
      node->zDisplacement(d_z_value);
      node->zVelocity(d_z_value);
      node->zNewDisplacement(d_z_value);
      node->zNewVelocity(d_z_value);
      node->zInternalForce(d_z_value);
      //std::cout << " Node = " << node->getID() << " Displacement = " << node->displacement() << std::endl;
      // Another way of doing things if needed.  This has been kept back so I can remember.
      //node->zExternalForce(-(node->internalForce())[2]);
    }
  }
}

