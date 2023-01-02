/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/AbaqusMeshGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>

namespace Uintah {

const string AbaqusMeshGeometryPiece::TYPE_NAME = "abaqus_mesh";

//---------------------------------------------------------------------------
// Read geometry data from file
//---------------------------------------------------------------------------
AbaqusMeshGeometryPiece::AbaqusMeshGeometryPiece(ProblemSpecP& ps) {
  // Set the default GeometryPiece type name
  d_name = "Unnamed " + TYPE_NAME + " from PS";

  // Get the input file name from .ups file
  ps->require("file_name", d_fileName);
  ps->getWithDefault("scaling_factor", d_scalefac, 1.0);
  ps->getWithDefault("translation_vector", d_translate, Vector(0.0, 0.0, 0.0));
  ps->getWithDefault("reflection_vector", d_reflect, Vector(1.0, 1.0, 1.0));
  ps->getWithDefault("axis_sequence", d_axis, IntVector(1, 2, 3));

  Vector rotate_row0, rotate_row1, rotate_row2;
  ps->getWithDefault("rotation_matrix_row0", rotate_row0, Vector(1, 0, 0));
  ps->getWithDefault("rotation_matrix_row1", rotate_row1, Vector(0, 1, 0));
  ps->getWithDefault("rotation_matrix_row2", rotate_row2, Vector(0, 0, 1));
  d_rotate = Matrix3(rotate_row0.x(),
                     rotate_row0.y(),
                     rotate_row0.z(),
                     rotate_row1.x(),
                     rotate_row1.y(),
                     rotate_row1.z(),
                     rotate_row2.x(),
                     rotate_row2.y(),
                     rotate_row2.z());

  // Decide whether to use the nodes or the gauss points at the MPM particle
  // positions (**warning** only for tets)
  d_use_gauss_pts = false;
  d_num_gauss_pts = 0;

  bool gauss_pts = false;
  ps->getWithDefault("use_gauss_points", d_use_gauss_pts, gauss_pts);
  if (d_use_gauss_pts) {
    int num_gauss_pts = 1;
    ps->getWithDefault("num_gauss_points", d_num_gauss_pts, num_gauss_pts);
    if (d_num_gauss_pts != 1 && d_num_gauss_pts != 4) {
      std::ostringstream out;
      out << "**ERROR** Gauss points can be used only for tetrahedral elements."
          << " Please input a value for num_gauss_pts of either 1 or 4.\n";
      throw ProblemSetupException(out.str(), __FILE__, __LINE__);
    }
  }

  proc0cout << "AbaqusMesh Geometry Piece: reading file " << d_fileName
            << std::endl;

  // Read the input file and save the data
  readMeshNodesAndElements(d_fileName);
}

//---------------------------------------------------------------------------
// Initialize geometry piece type name
//---------------------------------------------------------------------------
AbaqusMeshGeometryPiece::AbaqusMeshGeometryPiece(const string& /*file_name*/) {
  d_name = "Unnamed " + TYPE_NAME + " from file_name";
}

//---------------------------------------------------------------------------
// Write the output problem spec
//---------------------------------------------------------------------
void
AbaqusMeshGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("file_name", d_fileName);
  ps->appendElement("scaling_factor", d_scalefac);
  ps->appendElement("translation_vector", d_translate);
  ps->appendElement("reflection_vector", d_reflect);
  ps->appendElement("axis_sequence", d_axis);
  ps->appendElement("rotation_matrix_row0",
                    Vector(d_rotate(0, 0), d_rotate(0, 1), d_rotate(0, 2)));
  ps->appendElement("rotation_matrix_row1",
                    Vector(d_rotate(1, 0), d_rotate(1, 1), d_rotate(1, 2)));
  ps->appendElement("rotation_matrix_row2",
                    Vector(d_rotate(2, 0), d_rotate(2, 1), d_rotate(2, 2)));
  ps->appendElement("use_gauss_pts", d_use_gauss_pts);
  ps->appendElement("num_gauss_pts", d_num_gauss_pts);
}

//---------------------------------------------------------------------------
// Clone the geometry piece
//---------------------------------------------------------------------
GeometryPieceP
AbaqusMeshGeometryPiece::clone() const {
  return std::make_shared<AbaqusMeshGeometryPiece>(*this);
}

//--------------------------------------------------------------------------------
// Read input tetrahedral volume mesh file in Abaqus format
//--------------------------------------------------------------------------------
void
AbaqusMeshGeometryPiece::readMeshNodesAndElements(const std::string& fileName) {
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::ostringstream out;
    out << "Could not open node input volume mesh file " << fileName
        << " for reading \n";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Set up vectors to store nodes and elements
  std::vector<MeshNode> nodeArray;
  std::vector<SurfaceElement> surfElemArray;
  std::vector<VolumeElement> volElemArray;

  // Timing
  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  // Read file
  std::string line;
  bool node_flag                       = false;
  bool surf_elem_flag                  = false;
  [[maybe_unused]] bool line_elem_flag = false;
  bool vol_elem_flag                   = false;
  while (getline_safer(file, line)) {
    // std::cout << "line = " << line << std::endl;

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(), line.end(), ' '), line.end());

    // Skip comment lines except *Node and *Element
    if (line[0] == '*') {
      node_flag      = false;
      surf_elem_flag = false;
      vol_elem_flag  = false;
      if (line.compare("*Node") == 0 || line.compare("*NODE") == 0) {
        node_flag = true;
      }
      std::string str("*Element");
      std::string strCap("*ELEMENT");
      if (line.compare(0, str.length(), str) == 0 ||
          line.compare(0, str.length(), strCap) == 0) {
        std::regex line_token1(".*(T3D2).*");
        std::regex line_token2(".*(Line).*");
        std::regex surface_token1(".*(SURFACE).*");
        std::regex surface_token2(".*(Surface).*");
        // std::cout << "line = " << line ;
        if (std::regex_search(line.begin(), line.end(), line_token1) ||
            std::regex_search(line.begin(), line.end(), line_token2)) {
          line_elem_flag = true;
        } else if (std::regex_search(
                       line.begin(), line.end(), surface_token1) ||
                   std::regex_search(
                       line.begin(), line.end(), surface_token2)) {
          // std::cout << " contains the string SURFACE";
          surf_elem_flag = true;
        } else {
          // std::cout << " does not contain the string SURFACE: Volume
          // element";
          vol_elem_flag = true;
        }
        // std::cout << std::endl;
      }
      std::regex end_token(".*(End).*");
      if (std::regex_match(line.begin(), line.end(), end_token)) {
        // std::cout << "Line contains *End" << std::endl;
        node_flag      = false;
        surf_elem_flag = false;
        vol_elem_flag  = false;
        line_elem_flag = false;
      }
      continue;
    }

    // Read the nodal coordinates
    if (node_flag) {
      readMeshNode(line, nodeArray);
    }

    // Read the surface element connectivity
    if (surf_elem_flag) {
      readMeshSurfaceElement(line, surfElemArray);
    }

    // Read the volume element connectivity
    if (vol_elem_flag) {
      readMeshVolumeElement(line, volElemArray);
    }
  }

  // Timing
  endTime = std::chrono::system_clock::now();
  double time =
      std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done reading Abaqus mesh file in " << time << " millisecs"
            << std::endl;
  startTime = std::chrono::system_clock::now();

  // Renumber the volume elements starting from 1
  int lastSurfElemID = 0;
  if (surfElemArray.size() > 0) {
    lastSurfElemID = surfElemArray.back().id_;
  }
  // std::cout << "Last element id " << lastSurfElemID << std::endl;
  for (auto iter = volElemArray.begin(); iter != volElemArray.end(); ++iter) {
    (*iter).id_ -= lastSurfElemID;
    // std::cout << "Element id = " << (*iter).id_ << std::endl;
  }

  // Timing
  endTime = std::chrono::system_clock::now();
  time = std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done renumbering in " << time << " millisecs" << std::endl;
  startTime = std::chrono::system_clock::now();

  // Compute the volume of each volume element and sort in ascending order
  // of element id
  computeElementVolumes(nodeArray, volElemArray);
  std::sort(volElemArray.begin(),
            volElemArray.end(),
            [](const VolumeElement& a, const VolumeElement& b) {
              return b.id_ > a.id_;
            });

  // Create map of volume element index and element id
  std::map<int, int> volElemMap;
  int elemIndex = 0;
  for (auto elem : volElemArray) {
    volElemMap[elem.id_] = elemIndex;
    ++elemIndex;
  }

  // Timing
  endTime = std::chrono::system_clock::now();
  time = std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done computing element volumes in " << time << " millisecs"
            << std::endl;
  startTime = std::chrono::system_clock::now();

  // Print volume elements
  /*
  std::cout << "Volume elements : " << std::endl;
  for (auto iter = volElemArray.begin(); iter != volElemArray.end(); iter++) {
    VolumeElement volElem = *iter;
    std::cout << "\t ID = " << volElem.id_
              << ", nodes = (" << volElem.node1_ << ", " << volElem.node2_
              << ", " << volElem.node3_ << ", " << volElem.node4_ << ")"
              << ", vol = " << volElem.volume_ << std::endl;
  }

  std::cout << "Nodes : " << std::endl;
  for (auto iter = nodeArray.begin(); iter != nodeArray.end(); iter++) {
    MeshNode node = *iter;
    std::cout << "\t ID = " << node.id_
              << ", pos = (" << node.x_ << ", " << node.y_ << ", " << node.z_
              << "), vol = " << node.volume_ << std::endl;
  }
  */

  // Create particle locations and compute volumes
  std::vector<MeshNode> gaussPtArray;
  if (d_use_gauss_pts) {
    computeGaussPtVolumes(nodeArray, volElemArray, gaussPtArray);
  } else {
    computeNodalVolumes(nodeArray, volElemArray, volElemMap);
  }

  // Print nodes
  /*
  std::cout << "Nodes : " << std::endl;
  for (auto iter = nodeArray.begin(); iter != nodeArray.end(); iter++) {
    MeshNode node = *iter;
    std::cout << "\t ID = " << node.id_
              << ", pos = (" << node.x_ << ", " << node.y_ << ", " << node.z_
              << "), vol = " << node.volume_ << std::endl;
    std::cout << "\t Adj. Elems: ";
    for (auto e_iter = node.adjElements_.begin();
              e_iter != node.adjElements_.end();
              ++e_iter) {
      std::cout << *e_iter << ", ";
    }
    std::cout << std::endl;
  }
  */

  // Print surface elements
  /*
  std::cout << "Surface elements : " << std::endl;
  for (auto iter = surfElemArray.begin(); iter != surfElemArray.end(); iter++) {
    SurfaceElement surfElem = *iter;
    std::cout << "\t ID = " << surfElem.id_
              << ", nodes = (" << surfElem.node1_ << ", " << surfElem.node2_
              << ", " << surfElem.node3_ << ")" << std::endl;
  }
  */

  // Compute bounding box
  double xmax = std::numeric_limits<double>::min();
  double ymax = xmax;
  double zmax = xmax;
  double xmin = std::numeric_limits<double>::max();
  double ymin = xmin;
  double zmin = xmin;
  for (auto iter = nodeArray.begin(); iter != nodeArray.end(); iter++) {
    MeshNode node = *iter;
    xmax          = (xmax > node.x_) ? xmax : node.x_;
    ymax          = (ymax > node.y_) ? ymax : node.y_;
    zmax          = (zmax > node.z_) ? zmax : node.z_;
    xmin          = (xmin < node.x_) ? xmin : node.x_;
    ymin          = (ymin < node.y_) ? ymin : node.y_;
    zmin          = (zmin < node.z_) ? zmin : node.z_;
  }
  Point min(xmin, ymin, zmin);
  Point max(xmax, ymax, zmax);
  Vector fudge(1.e-5, 1.e-5, 1.e-5);
  min   = min - fudge;
  max   = max + fudge;
  d_box = Box(min, max);

  // Compute centroid (approx)
  Vector centroid(
      0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax));

  // Rotate points
  if (d_use_gauss_pts) {
    for (auto& node : gaussPtArray) {
      Vector point(node.x_, node.y_, node.z_);
      point -= centroid;
      point = d_rotate * point;
      point += centroid;
      node.x_ = point.x();
      node.y_ = point.y();
      node.z_ = point.z();
    }
  } else {
    for (auto& node : nodeArray) {
      Vector point(node.x_, node.y_, node.z_);
      point -= centroid;
      point = d_rotate * point;
      point += centroid;
      node.x_ = point.x();
      node.y_ = point.y();
      node.z_ = point.z();
    }
  }

  // Timing
  endTime = std::chrono::system_clock::now();
  time = std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done computing bounding box in " << time << " millisecs"
            << std::endl;

  std::cout << "Geometry object bounding box has min = " << min
            << " and max = " << max << std::endl;

  // Now push the coordinates and volumes of nodes into d_points and d_volumes
  if (d_use_gauss_pts) {
    for (auto node : gaussPtArray) {
      d_points.push_back(Point(node.x_, node.y_, node.z_));
      d_volume.push_back(node.volume_);
    }
  } else {
    for (auto iter = nodeArray.begin(); iter != nodeArray.end(); ++iter) {
      MeshNode node = *iter;
      d_points.push_back(Point(node.x_, node.y_, node.z_));
      d_volume.push_back(node.volume_);
    }
  }

  std::cout << "Number of MPM points = " << d_points.size()
            << " Number of MPM points/w volume = " << d_volume.size() << "\n";
}

void
AbaqusMeshGeometryPiece::readMeshNode(const std::string& inputLine,
                                      std::vector<MeshNode>& nodes) {
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() < 4) {
    std::ostringstream out;
    out << "Could not read nodal coordinates from " << inputLine << std::endl;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  auto iter   = data.begin();
  int node_id = std::stoi(*iter);
  ++iter;
  double xcoord = std::stod(*iter);
  ++iter;
  double ycoord = std::stod(*iter);
  ++iter;
  double zcoord = std::stod(*iter);

  int first        = d_axis.x() - 1;
  int second       = d_axis.y() - 1;
  int third        = d_axis.z() - 1;
  double coords[3] = {0.0, 0.0, 0.0};
  coords[first]    = xcoord;
  coords[second]   = ycoord;
  coords[third]    = zcoord;

  // Apply coordinate traanformations
  Vector point(d_translate.x() + (coords[0] * d_scalefac) * d_reflect.x(),
               d_translate.y() + (coords[1] * d_scalefac) * d_reflect.y(),
               d_translate.z() + (coords[2] * d_scalefac) * d_reflect.z());

  // Save the nodal coordinates
  nodes.emplace_back(MeshNode(node_id, point.x(), point.y(), point.z()));
}

void
AbaqusMeshGeometryPiece::readMeshVolumeElement(
    const std::string& inputLine, std::vector<VolumeElement>& elements) {
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() != 5) {
    std::ostringstream out;
    out << "**ERROR**Could not read volume element connectivity from input "
           "line: "
        << inputLine << std::endl;
    out << "\t\t Only four-noded tetrahedral elements are supported.\n";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Read the element id
  auto iter      = data.begin();
  int element_id = std::stoi(*iter);
  ++iter;

  // Read the element node ids
  std::vector<int> node_list;
  for (; iter != data.end(); ++iter) {
    int node_id = std::stoi(*iter);
    node_list.emplace_back(node_id);
  }
  if (node_list.empty()) {
    throw ProblemSetupException(
        "Could not find nodes in volume element input data stream",
        __FILE__,
        __LINE__);
  }

  // Save the data
  elements.emplace_back(VolumeElement(element_id, node_list));
}

void
AbaqusMeshGeometryPiece::readMeshSurfaceElement(
    const std::string& inputLine, std::vector<SurfaceElement>& elements) {
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() != 4) {
    std::ostringstream out;
    out << "Could not read surface element connectivity from input line: "
        << inputLine << std::endl;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Read the element id
  auto iter      = data.begin();
  int element_id = std::stoi(*iter);
  ++iter;

  // Read the element node ids
  std::vector<int> node_list;
  for (; iter != data.end(); ++iter) {
    int node_id = std::stoi(*iter);
    node_list.emplace_back(node_id);
  }
  if (node_list.empty()) {
    throw ProblemSetupException(
        "Could not find nodes in surface element input data stream",
        __FILE__,
        __LINE__);
  }

  // Save the data
  elements.emplace_back(SurfaceElement(element_id, node_list));
}

void
AbaqusMeshGeometryPiece::computeElementVolumes(
    std::vector<MeshNode>& nodes, std::vector<VolumeElement>& elements) {
  // Loop thru elements
  for (auto& elem : elements) {
    MeshNode node1 = nodes[elem.node1_ - 1];
    MeshNode node2 = nodes[elem.node2_ - 1];
    MeshNode node3 = nodes[elem.node3_ - 1];
    MeshNode node4 = nodes[elem.node4_ - 1];
    Point p0(node1.x_, node1.y_, node1.z_);
    Point p1(node2.x_, node2.y_, node2.z_);
    Point p2(node3.x_, node3.y_, node3.z_);
    Point p3(node4.x_, node4.y_, node4.z_);
    Vector vec01  = p1 - p0;
    Vector vec02  = p2 - p0;
    Vector vec03  = p3 - p0;
    Vector vec1x2 = Cross(vec01, vec02);
    double volume = 1.0 / 6.0 * std::abs(Dot(vec1x2, vec03));
    elem.volume_  = volume;
  }
}

void
AbaqusMeshGeometryPiece::computeGaussPtVolumes(
    std::vector<MeshNode>& nodes,
    std::vector<VolumeElement>& elements,
    std::vector<MeshNode>& gaussPts) {
  // Timing
  auto startTime = std::chrono::system_clock::now();

  // Set up parameter alpha
  double alpha = 0.25;
  if (d_num_gauss_pts == 4) {
    alpha = std::sqrt(5.0) / 5.0;
  }

  // Loop thru elements and quadrature point locations
  int pt_id = 0;
  for (auto elem : elements) {
    // Get the nodal coordinates
    Uintah::Vector e1(nodes[elem.node1_ - 1].x_,
                      nodes[elem.node1_ - 1].y_,
                      nodes[elem.node1_ - 1].z_);
    Uintah::Vector e2(nodes[elem.node2_ - 1].x_,
                      nodes[elem.node2_ - 1].y_,
                      nodes[elem.node2_ - 1].z_);
    Uintah::Vector e3(nodes[elem.node3_ - 1].x_,
                      nodes[elem.node3_ - 1].y_,
                      nodes[elem.node3_ - 1].z_);
    Uintah::Vector e4(nodes[elem.node4_ - 1].x_,
                      nodes[elem.node4_ - 1].y_,
                      nodes[elem.node4_ - 1].z_);

    // Compute locations of quadrature points
    if (d_num_gauss_pts == 1) {
      ++pt_id;
      Uintah::Vector q = (e1 + e2 + e3 + e4) * alpha;
      gaussPts.push_back(MeshNode(pt_id, q.x(), q.y(), q.z(), elem.volume_));

    } else {
      // Quadrature points (Journal of Computational and Applied Mathematics
      //                    Volume 236, Issue 17, November 2012, Pages
      //                    4348-4364)
      Uintah::FastMatrix A(4, 4);
      A(0, 0) = 0.5854101966249680;
      A(0, 1) = 0.1381966011250110;
      A(0, 2) = 0.1381966011250110;
      A(0, 3) = 0.1381966011250110;
      A(1, 0) = 0.1381966011250110;
      A(1, 1) = 0.5854101966249680;
      A(1, 2) = 0.1381966011250110;
      A(1, 3) = 0.1381966011250110;
      A(2, 0) = 0.1381966011250110;
      A(2, 1) = 0.1381966011250110;
      A(2, 2) = 0.5854101966249680;
      A(2, 3) = 0.1381966011250110;
      A(3, 0) = 0.1381966011250110;
      A(3, 1) = 0.1381966011250110;
      A(3, 2) = 0.1381966011250110;
      A(3, 3) = 0.5854101966249680;

      Uintah::Vector q =
          A(0, 0) * e1 + A(0, 1) * e2 + A(0, 2) * e3 + A(0, 3) * e4;
      ++pt_id;
      gaussPts.push_back(
          MeshNode(pt_id, q.x(), q.y(), q.z(), 0.25 * elem.volume_));

      q = A(1, 0) * e1 + A(1, 1) * e2 + A(1, 2) * e3 + A(1, 3) * e4;
      ++pt_id;
      gaussPts.push_back(
          MeshNode(pt_id, q.x(), q.y(), q.z(), 0.25 * elem.volume_));

      q = A(2, 0) * e1 + A(2, 1) * e2 + A(2, 2) * e3 + A(2, 3) * e4;
      ++pt_id;
      gaussPts.push_back(
          MeshNode(pt_id, q.x(), q.y(), q.z(), 0.25 * elem.volume_));

      q = A(3, 0) * e1 + A(3, 1) * e2 + A(3, 2) * e3 + A(3, 3) * e4;
      ++pt_id;
      gaussPts.push_back(
          MeshNode(pt_id, q.x(), q.y(), q.z(), 0.25 * elem.volume_));
    }
  }

  // Timing
  auto endTime = std::chrono::system_clock::now();
  auto time =
      std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done computing quadrature points in " << time << " millisecs"
            << std::endl;
  startTime = std::chrono::system_clock::now();
}

void
AbaqusMeshGeometryPiece::computeNodalVolumes(
    std::vector<MeshNode>& nodes,
    std::vector<VolumeElement>& elements,
    std::map<int, int>& elemMap) {
  // Timing
  auto startTime = std::chrono::system_clock::now();

  // Find the elements connected to each node
  findNodalAdjacentElements(nodes, elements);

  // Timing
  auto endTime = std::chrono::system_clock::now();
  auto time =
      std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done finding adjacent elements in " << time << " millisecs"
            << std::endl;
  startTime = std::chrono::system_clock::now();

  // Compute nodal volumes
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    MeshNode node   = *node_iter;
    double node_vol = 0.0;
    for (auto elem_id : node.adjElements_) {
      node_vol += 0.25 * elements[elemMap[elem_id]].volume_;
      /*
      std::cout << " Node = " << node.id_
                << " Elem = " << elem_id
                << " Elem Id = " << elem.id_
                << " vol += " <<  elem.volume_ << std::endl;
      */
    }
    (*node_iter).volume_ = node_vol;
  }

  // Timing
  endTime = std::chrono::system_clock::now();
  time = std::chrono::duration<double, std::milli>(endTime - startTime).count();
  std::cout << "Done computing nodal volumes in " << time << " millisecs"
            << std::endl;
}

void
AbaqusMeshGeometryPiece::findNodalAdjacentElements(
    std::vector<MeshNode>& nodes, std::vector<VolumeElement>& elements) {
  // Loop thru elements and find adjacent elements for each node
  for (VolumeElement cur_elem : elements) {
    // Loop thru nodes of each element
    /*
    std::cout << "Cur elem = " << cur_elem.id_ << " Nodes = "
              << cur_elem.node1_ << "," << cur_elem.node2_ << ","
              << cur_elem.node3_ << "," << cur_elem.node4_ << "\n";
    std::cout << "\t Ids: " << cur_elem.id_
              << " Node = " << cur_elem.node1_ << " id = " <<
    nodes[cur_elem.node1_-1].id_
              << " Node = " << cur_elem.node2_ << " id = " <<
    nodes[cur_elem.node2_-1].id_
              << " Node = " << cur_elem.node3_ << " id = " <<
    nodes[cur_elem.node3_-1].id_
              << " Node = " << cur_elem.node4_ << " id = " <<
    nodes[cur_elem.node4_-1].id_
              << "\n";
    */
    nodes[cur_elem.node1_ - 1].adjElements_.push_back(cur_elem.id_);
    nodes[cur_elem.node2_ - 1].adjElements_.push_back(cur_elem.id_);
    nodes[cur_elem.node3_ - 1].adjElements_.push_back(cur_elem.id_);
    nodes[cur_elem.node4_ - 1].adjElements_.push_back(cur_elem.id_);
    /*
    std::cout << "Node1: ";
    for (const auto& elem : nodes[cur_elem.node1_-1].adjElements_) {
      std::cout << elem << ",";
    }
    std::cout << "\n";
    std::cout << "Node4: ";
    for (const auto& elem : nodes[cur_elem.node4_-1].adjElements_) {
      std::cout << elem << ",";
    }
    std::cout << "\n";
    */
  }
}

//______________________________________________________________________
// *TODO* Implement this
bool
AbaqusMeshGeometryPiece::inside(const Point& p) const {
  // proc0cout << "**WARNING: 'inside' for Abaqus Mesh Geometry not implemented
  // yet."
  //           << std::endl;
  // Check p with the lower coordinates
  if (p == Max(p, d_box.lower()) && p == Min(p, d_box.upper()))
    return true;
  else
    return false;
}
//______________________________________________________________________
//
Box
AbaqusMeshGeometryPiece::getBoundingBox() const {
  return d_box;
}

//______________________________________________________________________
//
unsigned int
AbaqusMeshGeometryPiece::createPoints() {
  std::cout << "AbaqusMeshGeometryPiece: d_points.size() = " << d_points.size()
            << "\n";
  return d_points.size();
}

// Read a line without end of line characters messing things up
// From:
// https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf/6089413#6089413
std::istream&
AbaqusMeshGeometryPiece::getline_safer(std::istream& is, std::string& t) {
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if (sb->sgetc() == '\n') sb->sbumpc();
        return is;
      case std::streambuf::traits_type::eof():
        // Also handle the case when the last line has no line ending
        if (t.empty()) is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}

} // end namespace Uintah