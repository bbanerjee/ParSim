
#include <BoundaryConditions/TractionBC.h> 
#include <Geometry/Vector3D.h> 
#include <Core/Node.h>
#include <Pointers/NodeP.h>
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>
#include <iostream>

using namespace Matiti;
  
TractionBC::TractionBC() 
{
}

TractionBC::~TractionBC() 
{
}

// Initialize external forces on nodes
void
TractionBC::initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elems)
{
  // If ps is null return
  if (!(ps)) return;

  // Read the force vector
  SCIRun::Vector ext_traction(0.0, 0.0, 0.0);
  ps->require("traction", ext_traction);

  // Read the region of application of the force
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
    out << "**ERROR** Box region for application of external force is degenerate" << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Find the surface nodes
  NodePArray surface_nodes;
  findSurfaceNodesInBox(box_min, box_max, nodes, surface_nodes);

  //Find the max volume of surface nodes
  double maxVolume=0.0;
  findMaxVolume(surface_nodes, maxVolume);

  // Apply external force to the surface nodes
  computeExtForceDensity(ext_traction, surface_nodes, elems);
  
}

//********************************************************************
// subroutine ExtForceDenstiy
// Purpose : set the external force density array
//********************************************************************
void 
TractionBC::computeExtForceDensity(const SCIRun::Vector& extForce,
                                   NodePArray& surfaceNodes, 
                                   ElementPArray& elems)
{
  std::cerr << "ComputeExtForceDensity may not be implemented correctly yet. " << std::endl;
  for (auto node_iter = surfaceNodes.begin(); 
	    node_iter != surfaceNodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    double cur_node_vol = cur_node->volume();
    Vector3D ext_traction(-extForce[0]/cur_node_vol, -extForce[1]/cur_node_vol, -extForce[2]/cur_node_vol);
    cur_node->externalForce(ext_traction);
  }
}

