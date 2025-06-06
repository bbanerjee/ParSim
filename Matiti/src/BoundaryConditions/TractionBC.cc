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

#include <BoundaryConditions/TractionBC.h> 
#include <Geometry/Vector3D.h> 
#include <Core/Node.h>
#include <Core/Element.h>
#include <Core/Element2D.h>
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
  Uintah::Vector ext_traction(0.0, 0.0, 0.0);
  ps->require("traction", ext_traction);

  // Read the region of application of the force
  Uintah::Vector box_min(0.0, 0.0, 0.0), box_max(0.0, 0.0, 0.0);
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
//  std::cout << "Traction BCs: number of surface nodes = " << surface_nodes.size() << std::endl;

  if (surface_nodes.size() > 0) {

    // Find the area of support of each surface node
    findSurfaceNodeAreas(surface_nodes, elems);

    //Find the max volume of surface nodes
    double maxVolume=0.0;
    findMaxVolume(surface_nodes, maxVolume);

    // Apply external force to the surface nodes
    computeExtForceDensity(ext_traction, surface_nodes, elems);
  }
  
}

//-----------------------------------------------------------------
// Compute and set surface area covered by each surface node
// TODO: Need a better way of doig this
//-----------------------------------------------------------------
void
TractionBC::findSurfaceNodeAreas(NodePArray& surfaceNodes,
                                 ElementPArray& elems)
{
  // Loop thru the surface nodes
  for (auto iter = surfaceNodes.begin(); iter != surfaceNodes.end(); ++iter) {
    NodeP surf_node = *iter;

    // Find the elements adjacent to the nodes
    ElementPArray adj_elems = surf_node->getAdjacentElements();
//  std::cout << "Node ID: " << surf_node->getID() << " Number of adjacent elements: " << adj_elems.size() << std::endl;
    // Loop through elements
    std::vector<double> areas;
    for (auto e_iter = adj_elems.begin(); e_iter != adj_elems.end(); ++e_iter) {
      ElementP elem = *e_iter;
      NodePArray elem_nodes = elem->nodes();
//      std::cout << elem_nodes.size() << std::endl;

      // TODO: Create a list of element faces connected to the surface node
      switch(elem_nodes.size()) 
      {
        case 4: // Tetrahedron (nodes assumed to be in counter-clockwise order
                // base first, followed by vertex at peak)
        {
          Element2D face1(elem_nodes[0], elem_nodes[1], elem_nodes[2]); 
          Element2D face2(elem_nodes[0], elem_nodes[1], elem_nodes[3]); 
          Element2D face3(elem_nodes[0], elem_nodes[2], elem_nodes[3]); 
          Element2D face4(elem_nodes[1], elem_nodes[2], elem_nodes[3]); 
          if (face1.hasNode(surf_node) && face1.isSubset(surfaceNodes)) {
            areas.push_back(face1.area()/3.0);
          }
          if (face2.hasNode(surf_node) && face2.isSubset(surfaceNodes)) {
            areas.push_back(face2.area()/3.0);
          }
          if (face3.hasNode(surf_node) && face3.isSubset(surfaceNodes)) {
            areas.push_back(face3.area()/3.0);
          }
          if (face4.hasNode(surf_node) && face4.isSubset(surfaceNodes)) {
            areas.push_back(face4.area()/3.0);
          }
          break;
        }
        case 6: // Prism (nodes assumed to be in counter-clockwise order
                // base 1 followed by base 2
        {
          Element2D face1(elem_nodes[0], elem_nodes[1], elem_nodes[2]); 
          Element2D face2(elem_nodes[3], elem_nodes[4], elem_nodes[5]); 
          Element2D face3(elem_nodes[0], elem_nodes[1], elem_nodes[4], elem_nodes[3]); 
          Element2D face4(elem_nodes[1], elem_nodes[4], elem_nodes[5], elem_nodes[2]); 
          Element2D face5(elem_nodes[0], elem_nodes[3], elem_nodes[5], elem_nodes[2]); 
          if (face1.hasNode(surf_node) && face1.isSubset(surfaceNodes)) {
            areas.push_back(face1.area()/3.0);
          }
          if (face2.hasNode(surf_node) && face2.isSubset(surfaceNodes)) {
            areas.push_back(face2.area()/3.0);
          }
          if (face3.hasNode(surf_node) && face3.isSubset(surfaceNodes)) {
            areas.push_back(face3.area()/4.0);
          }
          if (face4.hasNode(surf_node) && face4.isSubset(surfaceNodes)) {
            areas.push_back(face4.area()/4.0);
          }
          if (face5.hasNode(surf_node) && face5.isSubset(surfaceNodes)) {
            areas.push_back(face5.area()/4.0);
          }
          break;
        }
        case 8: // Hex (nodes assumed to be in counter-clockwise order
                //  bottom first and then top)
        {
          Element2D face1(elem_nodes[0], elem_nodes[1], elem_nodes[2], elem_nodes[3]); 
          Element2D face2(elem_nodes[4], elem_nodes[5], elem_nodes[6], elem_nodes[7]); 
          Element2D face3(elem_nodes[0], elem_nodes[1], elem_nodes[5], elem_nodes[4]); 
          Element2D face4(elem_nodes[0], elem_nodes[4], elem_nodes[7], elem_nodes[3]); 
          Element2D face5(elem_nodes[1], elem_nodes[5], elem_nodes[6], elem_nodes[2]); 
          Element2D face6(elem_nodes[2], elem_nodes[6], elem_nodes[7], elem_nodes[3]); 
          if (face1.hasNode(surf_node) && face1.isSubset(surfaceNodes)) {
            areas.push_back(face1.area()/4.0);
          }
          if (face2.hasNode(surf_node) && face2.isSubset(surfaceNodes)) {
            areas.push_back(face2.area()/4.0);
          }
          if (face3.hasNode(surf_node) && face3.isSubset(surfaceNodes)) {
            areas.push_back(face3.area()/4.0);
          }
          if (face4.hasNode(surf_node) && face4.isSubset(surfaceNodes)) {
            areas.push_back(face4.area()/4.0);
          }
          if (face5.hasNode(surf_node) && face5.isSubset(surfaceNodes)) {
            areas.push_back(face5.area()/4.0);
          }
          if (face6.hasNode(surf_node) && face6.isSubset(surfaceNodes)) {
            areas.push_back(face6.area()/4.0);
          }
          break;
        }
        default:
        {
          std::ostringstream out;
          out << "**ERROR** Elements can only have 4, 6, or 8 nodes.";
          throw Exception(out.str(), __FILE__, __LINE__);
          break;
        }
      } // end of switch
    } // end e_iter

    if (areas.size() == 0) {
      std::ostringstream out;
      out << "**ERROR** Nodal area is zero.";
      throw Exception(out.str(), __FILE__, __LINE__);
    }

    // Add all areas and divide by the number of faces that contribute
    double nodeArea = 0.0;
    for (auto a_iter = areas.begin(); a_iter != areas.end(); a_iter++) {
      nodeArea += (*a_iter);
    }
    surf_node->area(nodeArea);


  } // end node loop
}

//********************************************************************
// subroutine ExtForceDenstiy
// Purpose : set the external force density array
//********************************************************************
void 
TractionBC::computeExtForceDensity(const Uintah::Vector& extTraction,
                                   NodePArray& surfaceNodes, 
                                   ElementPArray& elems)
{
  for (auto node_iter = surfaceNodes.begin(); 
	    node_iter != surfaceNodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    double cur_node_area = cur_node->area();
    double cur_node_vol = cur_node->volume();
    double fac = cur_node_area/cur_node_vol;
    Vector3D ext_force_den(extTraction[0]*fac, extTraction[1]*fac, extTraction[2]*fac);
    cur_node->externalForce(ext_force_den);

//    std::cout << "Surface node =" << *cur_node << std::endl;
//    std::cout << "area = " << cur_node->area() << "  ratio of area and volume = " << fac << std::endl;
  }
}

