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

#include <GeometryPiece/BoxGeometryPiece.h>
#include <Core/Exception.h>
#include <Core/Node.h>
#include <Core/Element.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

using namespace Matiti; 

BoxGeometryPiece::BoxGeometryPiece(Uintah::ProblemSpecP& ps,
                                   NodePArray& nodes,
                                   ElementPArray& elements, Vector3D& gridSize)
{
  d_name = "box";
  Uintah::Vector lower, upper;
  ps->require("lower", lower);
  ps->require("upper", upper);
  d_box = Box3D(Point3D(lower[0],lower[1],lower[2]), Point3D(upper[0],upper[1],upper[2]));
  if (d_box.isDegenerate()) {
    std::ostringstream out;
    out << "**ERROR** The box geometry piece is degenerate. Lower = " << d_box.lower()
        << " Upper = " << d_box.upper() << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
  Uintah::IntVector num_elem;
  ps->require("num_elements", num_elem);
  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_num_elements[ii] = num_elem(ii);
    if (d_num_elements[ii] < 1) d_num_elements[ii] = 1;
  }

  // Create nodes
  createNodes(nodes, gridSize);

  // Create elements
  createElements(elements);

  // Find adjacent elements for each node in the volume mesh
  findNodalAdjacentElements(elements);
}

BoxGeometryPiece::BoxGeometryPiece(const Point3D& lower,
                                   const Point3D& upper,
                                   const Uintah::IntVector& numElem,
                                   NodePArray& nodes,
                                   ElementPArray& elements, 
                                   Vector3D& gridSize)
{
  d_name = "box";
  d_box = Box3D(lower, upper);
  if (d_box.isDegenerate()) {
    std::ostringstream out;
    out << "**ERROR** The box geometry piece is degenerate. Lower = " << d_box.lower()
        << " Upper = " << d_box.upper() << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_num_elements[ii] = numElem(ii);
    if (d_num_elements[ii] < 1) d_num_elements[ii] = 1;
  }

  // Create nodes
  createNodes(nodes, gridSize);

  // Create elements
  createElements(elements);

  // Find adjacent elements for each node in the volume mesh
  findNodalAdjacentElements(elements);
}

BoxGeometryPiece::~BoxGeometryPiece()
{
}

Box3D 
BoxGeometryPiece::boundingBox() const
{
  return d_box;
}


bool 
BoxGeometryPiece::inside (const Point3D& pt) const
{
  return d_box.contains(pt); 
}


std::string 
BoxGeometryPiece::name() const
{
  return d_name;
}

void
BoxGeometryPiece::createNodes(NodePArray& nodes, Vector3D& gridSize)
{
  // Create nodes
  int nx = d_num_elements[0];
  int ny = d_num_elements[1];
  int nz = d_num_elements[2];
  double xmin = (d_box.lower())[0];
  double ymin = (d_box.lower())[1];
  double zmin = (d_box.lower())[2];
  Vector3D span = d_box.upper() - d_box.lower();
  double dx = span.x()/(double) nx;
  double dy = span.y()/(double) ny;
  double dz = span.z()/(double) nz;

  d_dx = dx;
  gridSize.x(dx);
  d_dy = dy;
  gridSize.y(dy);
  d_dz = dz;
  gridSize.z(dz);

  std::vector<double> xcoords, ycoords, zcoords;
  nx++;
  for (int ii=0; ii < nx; ++ii) {
    xcoords.emplace_back(xmin + dx*(double) ii);
  }
  ny++;
  for (int jj=0; jj < ny; ++jj) {
    ycoords.emplace_back(ymin + dy*(double) jj);
  }
  nz++;
  for (int kk=0; kk < nz; ++kk) {
    zcoords.emplace_back(zmin + dz*(double) kk);
  }
  int node_id = 0;
  int node_z = 0;
  for (auto ziter = zcoords.begin(); ziter != zcoords.end(); ++ziter) {
    ++node_z;
    bool on_z_surf = (node_z == 1) || (node_z == nz);
    int node_y = 0;
    for (auto yiter = ycoords.begin(); yiter != ycoords.end(); ++yiter) {
      ++node_y;
      bool on_y_surf = (node_y == 1) || (node_y == ny);
      int node_x = 0;
      for (auto xiter = xcoords.begin(); xiter != xcoords.end(); ++xiter) {
        ++node_x;
        bool on_x_surf = (node_x == 1) || (node_x == nx);
        ++node_id;
        NodeP node(new Node(node_id, *xiter, *yiter, *ziter, on_x_surf||on_y_surf||on_z_surf));
        nodes.emplace_back(node);
        //std::cout << "node = " << node_id << " " <<  node.position() ;

       // Add to the node ID -> node ptr map
       d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
      }
    }
  }
}

void
BoxGeometryPiece::createElements(ElementPArray& elements)
{
  // Create elements
  int nx = d_num_elements[0];
  int ny = d_num_elements[1];
  int nz = d_num_elements[2];

  int elem_id = 0;
  for (int kk=0; kk < nz; ++kk) {
    for (int jj=0; jj < ny; ++jj) {
      for (int ii=0; ii < nx; ++ii) {
        ++elem_id;

        std::vector<int> node_list;
        node_list.emplace_back(kk*(nx+1)*(ny+1)+jj*(nx+1)+ii+1);
        node_list.emplace_back(node_list[0]+1);
        node_list.emplace_back(node_list[1]+ (nx+1));
        node_list.emplace_back(node_list[2]-1);
        node_list.emplace_back(node_list[0] + (nx+1)*(ny+1));
        node_list.emplace_back(node_list[4]+1);
        node_list.emplace_back(node_list[5]+(nx+1));
        node_list.emplace_back(node_list[6]-1);

        // Find the node pointers
        NodePArray elem_nodes;
        if (d_id_ptr_map.empty()) {
          throw Exception("Could not find node id -> node ptr map", __FILE__, __LINE__);
        }
        for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
          int node_id = *iter;
          auto id_ptr_pair = d_id_ptr_map.find(node_id);
          if (id_ptr_pair == d_id_ptr_map.end()) {
            std::ostringstream out;
            out << "Could not find node id -> node ptr pair for node " << node_id;
            throw Exception(out.str(), __FILE__, __LINE__);
          }
          NodeP it = id_ptr_pair->second;
          elem_nodes.emplace_back(it); 
        }
     
        // Save the data
        ElementP elem(new Element(elem_id, elem_nodes));
        elem->computeVolume();
        elements.emplace_back(elem);
      }
    }
  }
}

void
BoxGeometryPiece::findNodalAdjacentElements(ElementPArray& elements)
{
  // Loop thru elements and find adjacent elements for each node
  for (auto elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter) {
    ElementP cur_elem = *elem_iter;

    // Loop thru nodes of each element
    NodePArray elem_nodes = cur_elem->nodes();
    for (auto elem_node_iter = elem_nodes.begin();
              elem_node_iter != elem_nodes.end(); ++elem_node_iter) {

      NodeP cur_elem_node = *elem_node_iter;
      cur_elem_node->addAdjacentElement(cur_elem);
    }
  }
}
