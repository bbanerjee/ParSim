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
{
}

DisplacementBC::~DisplacementBC() 
{
}

// Initialize displacement boundary conditions on nodes
void
DisplacementBC::initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes)
{
  // If ps is null return
  if (!(ps)) return;

  // Read the type of BC  (only symmetry bcs allowed for now)
  // symmetry x => u_1 = 0
  // symmetry y => u_2 = 0
  // symmetry z => u_3 = 0
  std::string bc_type;
  ps->require("symmetry", bc_type);
  if (bc_type == "x") {
    d_bc_type = BCType::XSymmetry;
    d_bc_val = 0.0;
  } else if (bc_type == "y") {
    d_bc_type = BCType::YSymmetry;
    d_bc_val = 0.0;
  } else if (bc_type == "z") {
    d_bc_type = BCType::ZSymmetry;
    d_bc_val = 0.0;
  } else {
    std::ostringstream out;
    out << "**ERROR** Symmetry BC type must be either x, y, or z." << std::endl;
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

  // Find the surface nodes
  NodePArray surface_nodes;
  findSurfaceNodesInBox(box_min, box_max, nodes, surface_nodes);
}


