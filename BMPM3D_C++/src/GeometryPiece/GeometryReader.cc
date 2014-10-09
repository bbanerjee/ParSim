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

#include <GeometryPiece/GeometryReader.h>
#include <Node.h>
#include <Element.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace BrMPM; 

GeometryReader::GeometryReader(Uintah::ProblemSpecP& ps,
                               NodePArray& nodes,
                               ElementPArray& elems)
{
  d_xmax = std::numeric_limits<double>::min();
  d_ymax = d_xmax;
  d_zmax = d_xmax;
  d_xmin = std::numeric_limits<double>::max();
  d_ymin = d_xmin;
  d_zmin = d_xmin;
  d_num_buckets_x = 20;
  d_name = "file";

  // Read input files
  readGeometryInputFiles(ps, nodes, elems);

  // Find adjacent elements for each node in the volume mesh
  findNodalAdjacentElements(elems);

  // Find surface nodes
  findSurfaceNodes(nodes);
}

GeometryReader::~GeometryReader()
{
}

Box3D 
GeometryReader::boundingBox() const
{
  Point3D lower = Point3D(d_xmin, d_ymin, d_zmin);
  Point3D upper = Point3D(d_xmax, d_ymax, d_zmax);
  return Box3D(lower, upper);
}

bool 
GeometryReader::inside (const Point3D& pt) const
{
  Box3D box = boundingBox();
  return box.contains(pt);
}

std::string 
GeometryReader::name() const
{
  return d_name;
}

void 
GeometryReader::readGeometryInputFiles(Uintah::ProblemSpecP& geom_ps,
                                       NodePArray& nodes,
                                       ElementPArray& elements)
{
  // Get the geometry (from input node and element files)
  //Uintah::ProblemSpecP geom_ps = ps->findBlock("Geometry");
  std::string input_surface_mesh_file;
  std::string input_volume_mesh_file;
  geom_ps->require("input_surface_mesh_file", input_surface_mesh_file);
  geom_ps->require("input_volume_mesh_file", input_volume_mesh_file);
  std::cout << "Input geometry mesh files: " << input_surface_mesh_file << ", " 
                                        << input_volume_mesh_file << std::endl;

  // Read the input surface node file
  readSurfaceMeshNodes(input_surface_mesh_file);

  // Read the input volume mesh file for nodes and elements
  readVolumeMeshNodesAndElements(input_volume_mesh_file, nodes, elements);

}

//--------------------------------------------------------------------------------
// Read input triangulated surface mesh file in Abaqus format
//--------------------------------------------------------------------------------
void 
GeometryReader::readSurfaceMeshNodes(const std::string& fileName)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open node input surface mesh file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  std::string line;
  bool node_flag = false;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(),line.end(),' '),line.end());
    
    // Skip comment lines except *Node
    bool node_flag_old = node_flag;
    if (line[0] == '*') {
      node_flag = false;
      if (line.compare("*Node") == 0) {
        node_flag = true;
      }
      continue;
    }
    if (!node_flag && node_flag_old) {
      break;
    }

    // Read the nodal coordinates
    if (node_flag) {

      // Tokenize the string
      std::string data_piece;
      std::istringstream data_stream(line);
      std::vector<std::string> data;
      while (std::getline(data_stream, data_piece, ',')) {
        data.push_back(data_piece);
      }
      if (data.size() < 4) {
        std::ostringstream out;
        out << "Could not read node input surface mesh line " 
            << line << std::endl;
        throw Exception(out.str(), __FILE__, __LINE__);
      }
      
      auto iter = data.begin(); 
      //int node_id = std::stoi(*iter); 
      ++iter;
      double xcoord = std::stod(*iter); ++iter;
      double ycoord = std::stod(*iter); ++iter;
      double zcoord = std::stod(*iter);

      d_surf_pts.emplace_back(Point3D(xcoord, ycoord, zcoord)); 
      d_xmax = (xcoord > d_xmax) ? xcoord : d_xmax;
      d_ymax = (ycoord > d_ymax) ? ycoord : d_ymax;
      d_zmax = (zcoord > d_zmax) ? zcoord : d_zmax;
      d_xmin = (xcoord < d_xmin) ? xcoord : d_xmin;
      d_ymin = (ycoord < d_ymin) ? ycoord : d_ymin;
      d_zmin = (zcoord < d_zmin) ? zcoord : d_zmin;
      //std::cout << "id = " << node_id << "(" << xcoord << ", " 
      //          << ycoord << ", " << zcoord << ")" << std::endl;
    }
  }

  std::cout << "Body : min = [" << d_xmin << "," << d_ymin << "," << d_zmin
            << " max = [" << d_xmax << "," << d_ymax << "," << d_zmax
            << std::endl;

  // Create a box that surrounds the surface and
  // divide the box into a grid (hardcoded for now)
  double dx = (d_xmax - d_xmin)/(double) d_num_buckets_x;

  // Loop through nodes and add to buckets
  int node_id = 0;
  for (auto iter = d_surf_pts.begin(); iter != d_surf_pts.end(); ++iter) {
    double xx = (*iter).x();
    double yy = (*iter).y();
    double zz = (*iter).z();
    int x_cell_id = std::ceil((xx - d_xmin)/dx); 
    int y_cell_id = std::ceil((yy - d_ymin)/dx); 
    int z_cell_id = std::ceil((zz - d_zmin)/dx); 
    long64 cell_id = ((long64)x_cell_id << 16) | ((long64)y_cell_id << 32) | ((long64)z_cell_id << 48);
    d_bucket_to_node_map.insert(std::pair<long64, int>(cell_id, node_id));
    ++node_id;
  }

  // Print the map
  // for (auto  iter = d_bucket_to_node_map.begin(); iter != d_bucket_to_node_map.end(); iter++) {
  //   std::cout << "Cell = " << iter->first << " Node id = " << iter->second << std::endl;
  // }
}

//--------------------------------------------------------------------------------
// Read input tetrahedral volume mesh file in Abaqus format
//--------------------------------------------------------------------------------
void 
GeometryReader::readVolumeMeshNodesAndElements(const std::string& fileName,
                                               NodePArray& nodes,
                                               ElementPArray& elements)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open node input volume mesh file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  std::string line;
  bool node_flag = false;
  bool elem_flag = false;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(),line.end(),' '),line.end());
    
    // Skip comment lines except *Node
    if (line[0] == '*') {
      node_flag = false;
      elem_flag = false;
      if (line.compare("*Node") == 0) {
        node_flag = true;
      }
      std::string str("*Element");
      if (line.compare(0, str.length(), str) == 0) {
        elem_flag = true;
      }
      continue;
    }

    // Read the nodal coordinates
    if (node_flag) {
      readVolumeMeshNode(line, nodes);
    }

    // Read the element connectivity
    if (elem_flag) {
      readVolumeMeshElement(line, elements);
    }
  }
}

void
GeometryReader::readVolumeMeshNode(const std::string& inputLine,
                                   NodePArray& nodes)
{
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() < 4) {
    std::ostringstream out;
    out << "Could not read node input volume mesh line " 
        << inputLine << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  auto iter = data.begin(); 
  int node_id = std::stoi(*iter); ++iter;
  double xcoord = std::stod(*iter); ++iter;
  double ycoord = std::stod(*iter); ++iter;
  double zcoord = std::stod(*iter);

  // Save the data - surface_node_flag is 0
  bool on_surface = false;
  NodeP node(new Node(node_id, xcoord, ycoord, zcoord, on_surface));
  nodes.emplace_back(node);

  // Add to the node ID -> node ptr map
  d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
}

void
GeometryReader::readVolumeMeshElement(const std::string& inputLine,
                                      ElementPArray& elements)
{
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() < 5) {
    std::ostringstream out;
    out << "Could not read element input volume mesh line " 
        << inputLine << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Read the element id
  auto iter = data.begin(); 
  int element_id = std::stoi(*iter); ++iter;

  // Read the element node ids
  std::vector<int> node_list;
  for (; iter != data.end(); ++iter) {
    int node_id = std::stoi(*iter);
    node_list.emplace_back(node_id);
  }
  if (node_list.empty()) {
    throw Exception("Could not find nodes in element input data stream", __FILE__, __LINE__);
  }

  // Find the node pointers
  NodePArray elem_nodes;
  if (d_id_ptr_map.empty()) {
    throw Exception("Could not find node id -> node ptr map", __FILE__, __LINE__);
  }
  for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
    int node_id = *iter;
    auto id_ptr_pair = d_id_ptr_map.find(node_id);
    if (id_ptr_pair == d_id_ptr_map.end()) {
      std::string out = "Could not find node id -> node ptr pair for node " + node_id;
      throw Exception(out, __FILE__, __LINE__);
    }
    NodeP it = id_ptr_pair->second;
    elem_nodes.emplace_back(it); 
  }
     
  // Save the data
  ElementP elem(new Element(element_id, elem_nodes));
  elem->computeVolume();
  elements.emplace_back(elem);
}

void
GeometryReader::findNodalAdjacentElements(ElementPArray& elements)
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

void
GeometryReader::findSurfaceNodes(NodePArray& nodes)
{
  // Loop through nodes in volume mesh
  double dx = (d_xmax - d_xmin)/(double) d_num_buckets_x;
  int surf_node_count = 0;
  for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {

    // Get node position
    NodeP cur_node = *iter;
    Point3D pos = cur_node->position();

    // Compute bucket id
    double xx = pos.x();
    double yy = pos.y();
    double zz = pos.z();
    int x_cell_id = std::ceil((xx - d_xmin)/dx); 
    int y_cell_id = std::ceil((yy - d_ymin)/dx); 
    int z_cell_id = std::ceil((zz - d_zmin)/dx); 
    long64 cell_id = ((long64)x_cell_id << 16) | ((long64)y_cell_id << 32) | ((long64)z_cell_id << 48);

    // Find the cell id in the bucket-node map for the surface nodes
    auto range = d_bucket_to_node_map.equal_range(cell_id);
    // std::cout << "[" << x_cell_id << "," << y_cell_id << "," << z_cell_id <<"]"
    //           << cell_id <<  " =>";
    // for (auto range_iter = range.first; range_iter != range.second; ++range_iter) {
    //   std::cout << ' ' << range_iter->second;
    // }
    // std::cout << std::endl;

    // Loop over the range (surface nodes in the bucket)
    for (auto range_iter = range.first; range_iter != range.second; ++range_iter) {
      int node_id = (*range_iter).second;

      // Check if the point is on the surface (**WARNING** hardcoded for now)
      
      if (pos.distance(d_surf_pts[node_id]) < 1.0e-6) {

        // std::cout << "Node = " << pos << " Bucket = " << (*range_iter).first
        //           << " Node ID = " << (*range_iter).second 
        //           << " Pos = " << d_surf_pts[node_id] << std::endl;

        // Update the on-surface flags
        cur_node->onSurface(true);
        surf_node_count++;
        break;
      }
    }
    
  }
  std::cout << surf_node_count << std::endl;
}
