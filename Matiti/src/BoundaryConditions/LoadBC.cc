#include <BoundaryConditions/LoadBC.h> 
#include <Geometry/Vector3D.h> 
#include <Core/Node.h>
#include <Pointers/NodeP.h>
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>
#include <iostream>

using namespace Matiti;
  
LoadBC::LoadBC() 
{
}

LoadBC::~LoadBC() 
{
}


void
LoadBC::findSurfaceNodesInBox(const SCIRun::Vector& boxMin, 
                              const SCIRun::Vector& boxMax,
                              const NodePArray& nodes, 
                              NodePArray& surfaceNodes)
{
  for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
    NodeP node = *iter;
    if (!(node->onSurface())) continue;

    const Point3D& pos = node->position();
    if (pos.x() > boxMin.x() && pos.y() > boxMin.y() && pos.z() > boxMin.z() && 
        pos.x() < boxMax.x() && pos.y() < boxMax.y() && pos.z() < boxMax.z()) {
      surfaceNodes.push_back(node);
    }
  }
}

void
LoadBC::findMaxVolume(NodePArray& surfaceNodes, double& maxVol)
{
  for (auto node_iter = surfaceNodes.begin(); 
	    node_iter != surfaceNodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->numAdjacentElements()==4) {
        maxVol=cur_node->volume();
        break;
    }
  }
}



