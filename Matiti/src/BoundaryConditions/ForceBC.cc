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


#include <BoundaryConditions/ForceBC.h> 
#include <Geometry/Vector3D.h> 
#include <Core/Node.h>
#include <Pointers/NodeP.h>
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>
#include <iostream>

using namespace Matiti;
  
ForceBC::ForceBC() 
{
}

ForceBC::~ForceBC() 
{
}

// Initialize external forces on nodes
void
ForceBC::initialize(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray&)
{
  // If ps is null return
  if (!(ps)) return;

  // Read the force vector
  SCIRun::Vector ext_force(0.0, 0.0, 0.0);
  ps->require("force", ext_force);

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
  computeExtForceDensity(ext_force, surface_nodes, maxVolume);
  
}

//********************************************************************
// subroutine ExtForceDenstiy
// Purpose : set the external force density array
//********************************************************************
void 
ForceBC::computeExtForceDensity(const SCIRun::Vector& extForce,
                                NodePArray& surfaceNodes, double& maxVol)
{
  std::cerr << "ComputeExtForceDensity may not be implemented correctly yet. " << std::endl;
  for (auto node_iter = surfaceNodes.begin(); 
	    node_iter != surfaceNodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    double cur_node_vol = cur_node->volume();
    Vector3D ext_force(extForce[0]/cur_node_vol, extForce[1]/cur_node_vol, extForce[2]/cur_node_vol);
    cur_node->externalForce(ext_force);
  }
}

// Read the input data file to get
// 1) Surface node sets on which external forces are applied
// 2) external force vector
void 
ForceBC::initialize(std::string input_data_file)
{
  std::cerr << "initialize from input node file not implemented yet. " << std::endl;
  //for (int ii = 0; ii < d_num_ext_force_sets; ++ii) {
  //  input_stream << d_ext_force ;
  //  input_stream << d_ext_force_nodes;    
  //}
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void 
ForceBC::computeExtForceDensity(const NodePArray& nodes,
	  	                const Point3D& topLoc,
		                const Point3D& botLoc,
                                const Vector3D& extForce)
{
  // Set exyernal force magnitude
  double ext_force_mag = extForce[1];

  // Find the nodes at the top and bottom boundaries
  NodePArray top_nodes;
  NodePArray bot_nodes;
  findSurfaceNodes(nodes, topLoc, botLoc, top_nodes, bot_nodes);

  // Sort the nodes in increasing x-coordinate order and find span
  NodePArray sorted_top_nodes;
  std::vector<double> nodal_span_top;
  sortNodes(top_nodes, sorted_top_nodes, nodal_span_top);

  NodePArray sorted_bot_nodes;
  std::vector<double> nodal_span_bot;
  sortNodes(top_nodes, sorted_bot_nodes, nodal_span_bot);

  // Compute external force density
  int count = 0;
  Vector3D ext_force(0.0, 0.0, 0.0);
  for (auto node_iter=sorted_top_nodes.begin(); 
	    node_iter != sorted_top_nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    double cur_node_vol = cur_node->volume();
    ext_force.y(ext_force_mag*nodal_span_top[count]/(2.0*cur_node_vol));
    cur_node->externalForce(ext_force);
    count++;
  }

  count = 0;
  ext_force.reset();
  for (auto node_iter=sorted_bot_nodes.begin(); 
            node_iter != sorted_bot_nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    double cur_node_vol = cur_node->volume();
    ext_force.y(ext_force_mag*nodal_span_bot[count]/(2.0*cur_node_vol));
    cur_node->externalForce(ext_force);
    count++;
  }
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void
ForceBC::findSurfaceNodes(const NodePArray& nodes,
	  	           const Point3D& topLoc,
		           const Point3D& botLoc,
		           NodePArray& topNodes,
			   NodePArray& botNodes)
{
  double y_top = topLoc.y();
  double y_bot = botLoc.y();
  for (auto node_iter=nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    const Point3D& node_pos = cur_node->position();
    double y_cur = node_pos.y();
    if (std::abs(y_bot-y_cur) <= 1.0e-6) {
      botNodes.push_back(cur_node);
    } else {
      if (std::abs(y_top-y_cur) <= 1.0e-6) {
        topNodes.push_back(cur_node);
      } 
    }
  }
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void 
ForceBC::sortNodes(const NodePArray& surfaceNodes,
		   NodePArray& sortedSurfaceNodes,
		   std::vector<double>& nodalSpan)
{
  std::vector<double> xCoordsSurface;
  int count = 0;
  for (auto node_iter=surfaceNodes.begin(); node_iter != surfaceNodes.end(); 
            ++node_iter) {
    NodeP cur_node = *node_iter;
    const Point3D& node_pos = cur_node->position();
    xCoordsSurface.push_back(node_pos.x());
    sortedSurfaceNodes.push_back(cur_node);
    count++;
  }

  for (int ii=0; ii < count; ii++) {
    for (int jj=ii+1; jj < count; jj++) {
      if (xCoordsSurface[jj] > xCoordsSurface[ii]) {
        double temp = xCoordsSurface[ii];
        xCoordsSurface[ii] = xCoordsSurface[jj];
        xCoordsSurface[jj] = temp;
	NodeP tempNode = sortedSurfaceNodes[ii];
	sortedSurfaceNodes[ii] = sortedSurfaceNodes[jj];
	sortedSurfaceNodes[jj] = tempNode;
      }
    }
  }

  nodalSpan.push_back(xCoordsSurface[1] - xCoordsSurface[0]);
  for (int ii=1; ii < count-1; ii++) {
    nodalSpan.push_back(xCoordsSurface[ii+1] - xCoordsSurface[ii-1]);
  }
  nodalSpan.push_back(xCoordsSurface[count-1] - xCoordsSurface[count-2]);
}


