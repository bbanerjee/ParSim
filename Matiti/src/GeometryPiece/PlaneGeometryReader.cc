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

#include <GeometryPiece/PlaneGeometryReader.h>
#include <Core/Node.h>
#include <Core/Element.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace Matiti; 

PlaneGeometryReader::PlaneGeometryReader(Uintah::ProblemSpecP& ps,
                                         NodePArray& nodes,
                                         ElementPArray& elems)
{
  d_name = "plane_file";
  d_xmax = std::numeric_limits<double>::min();
  d_ymax = d_xmax;
  d_zmax = d_xmax;
  d_xmin = std::numeric_limits<double>::max();
  d_ymin = d_xmin;
  d_zmin = d_xmin;

  // Read input files
  readGeometryInputFiles(ps, nodes, elems);

  // Find adjacent elements for each node in the volume mesh
  findNodalAdjacentElements(elems);
}

PlaneGeometryReader::~PlaneGeometryReader()
{
}

std::string 
PlaneGeometryReader::name() const
{
  return d_name;
}

Box3D 
PlaneGeometryReader::boundingBox() const
{
  Point3D lower = Point3D(d_xmin, d_ymin, d_zmin);
  Point3D upper = Point3D(d_xmax, d_ymax, d_zmax);
  return Box3D(lower, upper);
}

bool 
PlaneGeometryReader::inside (const Point3D& pt) const
{
  Box3D box = boundingBox();
  return box.contains(pt);
}

void 
PlaneGeometryReader::readGeometryInputFiles(Uintah::ProblemSpecP& geom_ps,
                                            NodePArray& nodes,
                                            ElementPArray& elements)
{
  // Get the geometry (from input node and element files)
  std::string input_node_file;
  std::string input_element_file;
  geom_ps->require("input_node_file", input_node_file);
  geom_ps->require("input_element_file", input_element_file);
  std::cout << "Input geometry mesh files: " << input_node_file << ", " 
                                        << input_element_file << std::endl;

  // Read the input mesh file for nodes and elements
  readMeshNodesAndElements(input_node_file, input_element_file, nodes, elements);
}

//--------------------------------------------------------------------------------
// Read input two-dimensional mesh file in EMUNE format
//--------------------------------------------------------------------------------
void 
PlaneGeometryReader::readMeshNodesAndElements(const std::string& nodeFileName,
                                              const std::string& elementFileName,
                                              NodePArray& nodes,
                                              ElementPArray& elements)
{
  // Read the node file
  int num_nodes = readNodeFile(nodeFileName, nodes);

  // Read the element file
  readElementFile(elementFileName, elements, num_nodes);
}

int
PlaneGeometryReader::readNodeFile(const std::string& fileName,
                                  NodePArray& nodes)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open node input file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  int num_nodes = 0;
  std::string line;
  while (std::getline(file, line)) {

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
         std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    // Ignore empty lines
    if (line.empty()) continue;

    // Skip comment lines
    if (line[0] == '#') continue;

    // Read the data
    //std::cout << line << std::endl;
    std::istringstream data_stream(line);
    int node_id, surface_node_flag;
    double xcoord, ycoord, zcoord = 0.0;
    if (!(data_stream >> node_id >> xcoord >> ycoord >> surface_node_flag)) {
      throw Exception("Could not read node input data stream", __FILE__, __LINE__);
    } 

    // Save the data
    NodeP node(new Node(node_id, xcoord, ycoord, zcoord, surface_node_flag));
    nodes.emplace_back(node);
    ++num_nodes;

    // Find the bounding box
    d_xmax = (xcoord > d_xmax) ? xcoord : d_xmax;
    d_ymax = (ycoord > d_ymax) ? ycoord : d_ymax;
    d_zmax = (zcoord > d_zmax) ? zcoord : d_zmax;
    d_xmin = (xcoord < d_xmin) ? xcoord : d_xmin;
    d_ymin = (ycoord < d_ymin) ? ycoord : d_ymin;
    d_zmin = (zcoord < d_zmin) ? zcoord : d_zmin;

    // Add to the node ID -> node ptr map
    d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
  }

  // Create a duplicate layer a distance of min(x node spacing) away 
  // in the z-direction (assumes at least two nodes are available)
  auto iter = nodes.begin(); 
  double x_prev = (*iter)->x(); 
  iter++;
  double x_spacing = std::abs((*iter)->x() - x_prev);
  x_prev = (*iter)->x();
  iter++;
  for (; iter != nodes.end(); iter++) {
    x_spacing = std::min(x_spacing, std::abs((*iter)->x() - x_prev));
    x_prev = (*iter)->x();
  }
  //std::cout << "xspacing = " << x_spacing << std::endl;

  for (int ii = 0; ii < num_nodes; ii++) {

    // Save the data
    int node_id = num_nodes + ii + 1;
    double xcoord = nodes[ii]->x();
    double ycoord = nodes[ii]->y();
    double zcoord = x_spacing;
    bool surface_node_flag = nodes[ii]->onSurface();
    
    NodeP node(new Node(node_id, xcoord, ycoord, zcoord, surface_node_flag));
    nodes.emplace_back(node);

    // Add to the node ID -> node ptr map
    d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
  }

  return num_nodes;
}

void
PlaneGeometryReader::readElementFile(const std::string& fileName,
                                     ElementPArray& elements,
                                     int numNodes)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not element open input file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  std::string line;
  while (std::getline(file, line)) {

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
         std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    // Ignore empty lines
    if (line.empty()) continue;

    // Skip comment lines
    if (line[0] == '#') continue;

    // Read the element id
    std::istringstream data_stream(line);
    int element_id;
    if (!(data_stream >> element_id)) {
      throw Exception("Could not read element id from element input data stream", __FILE__, __LINE__);
    }

    // Read the node ids
    std::vector<int> node_list;
    int node;
    while (data_stream >> node) {
      node_list.emplace_back(node);
    }
    if (node_list.empty()) {
      throw Exception("Could not find nodes in element input data stream", __FILE__, __LINE__);
    }

    // Add the nodes in the parallel z-plane
    int nodes_in_elem = node_list.size();
    for (int ii =0; ii < nodes_in_elem; ii++) {
      node_list.emplace_back(node_list[ii]+numNodes);
    }

    // Find the node pointers
    if (d_id_ptr_map.empty()) {
      throw Exception("Could not find node id -> node ptr map", __FILE__, __LINE__);
    }
    NodePArray nodes;
    for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
      int node_id = *iter;
      auto id_ptr_pair = d_id_ptr_map.find(node_id);
      if (id_ptr_pair == d_id_ptr_map.end()) {
        std::ostringstream out;
        out << "Could not find node id -> node ptr pair for node " << node_id;
        throw Exception(out.str(), __FILE__, __LINE__);
      }
      NodeP it = id_ptr_pair->second;
      nodes.emplace_back(it); 
    }
     
    // Save the data
    ElementP elem(new Element(element_id, nodes));
    elem->computeVolume();
    elements.emplace_back(elem);
  }
}


void
PlaneGeometryReader::findNodalAdjacentElements(ElementPArray& elements)
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

