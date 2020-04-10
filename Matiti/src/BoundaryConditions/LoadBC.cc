/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
LoadBC::findSurfaceNodesInBox(const Uintah::Vector& boxMin, 
                              const Uintah::Vector& boxMax,
                              const NodePArray& nodes, 
                              NodePArray& surfaceNodes)
{
  for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
    NodeP node = *iter;
    if (!(node->onSurface())) continue;
    //std::cout << "Node position = " << node->position() << " box min = " << boxMin.x() << "," << boxMin.y()
    //          << " box max = " << boxMax.x() << "," << boxMax.y() << std::endl;

    const Point3D& pos = node->position();
    if (pos.x() >= boxMin.x() && pos.y() >= boxMin.y() && pos.z() >= boxMin.z() && 
        pos.x() <= boxMax.x() && pos.y() <= boxMax.y() && pos.z() <= boxMax.z()) {
      surfaceNodes.push_back(node);
      node->allowFailure(false);
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



