#include <Core/GeometryPiece/AbaqusMeshGeometryPiece.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <regex>

using namespace Uintah;
using namespace SCIRun;
using std::cerr;
using std::cout;
using std::endl;

const string AbaqusMeshGeometryPiece::TYPE_NAME = "abaqus_mesh";

//---------------------------------------------------------------------------
// Read geometry data from file
//---------------------------------------------------------------------------
AbaqusMeshGeometryPiece::AbaqusMeshGeometryPiece(ProblemSpecP & ps)
{
  // Set the default GeometryPiece type name
  name_ = "Unnamed " + TYPE_NAME + " from PS";

  // Get the input file name from .ups file
  ps->require("file_name",d_fileName);
  
  proc0cout << "AbaqusMesh Geometry Piece: reading file " << d_fileName << std::endl;
  
  // Read the input file and dave the data
  readMeshNodesAndElements(d_fileName);
}

//---------------------------------------------------------------------------
// Initialize geometry piece type name
//---------------------------------------------------------------------------
AbaqusMeshGeometryPiece::AbaqusMeshGeometryPiece(const string& /*file_name*/)
{
  name_ = "Unnamed " + TYPE_NAME + " from file_name";
}

AbaqusMeshGeometryPiece::~AbaqusMeshGeometryPiece()
{
}

//---------------------------------------------------------------------------
// Write the output problem spec
//---------------------------------------------------------------------
void
AbaqusMeshGeometryPiece::outputHelper( ProblemSpecP & ps ) const
{
  ps->appendElement("file_name",d_fileName);
}

//---------------------------------------------------------------------------
// Clone the geometry piece
//---------------------------------------------------------------------
GeometryPieceP
AbaqusMeshGeometryPiece::clone() const
{
  return scinew AbaqusMeshGeometryPiece(*this);
}

//--------------------------------------------------------------------------------
// Read input tetrahedral volume mesh file in Abaqus format
//--------------------------------------------------------------------------------
void 
AbaqusMeshGeometryPiece::readMeshNodesAndElements(const std::string& fileName)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::ostringstream out;
    out << "Could not open node input volume mesh file " 
        << fileName << " for reading \n";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Set up vectors to store nodes and elements
  std::vector<MeshNode> nodeArray;
  std::vector<SurfaceElement> surfElemArray;
  std::vector<VolumeElement> volElemArray;

  // Read file
  std::string line;
  bool node_flag = false;
  bool surf_elem_flag = false;
  bool vol_elem_flag = false;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(),line.end(),' '),line.end());
    
    // Skip comment lines except *Node and *Element
    if (line[0] == '*') {
      node_flag = false;
      surf_elem_flag = false;
      vol_elem_flag = false;
      if (line.compare("*Node") == 0) {
        node_flag = true;
      }
      std::string str("*Element");
      if (line.compare(0, str.length(), str) == 0) {
        std::regex surface_token("SURFACE");
        if (std::regex_match(line.begin(), line.end(), surface_token)) {
          surf_elem_flag = true;
        } else {
          vol_elem_flag = true;
        }
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

  // Compute the volume of each volume element
  computeElementVolumes(nodeArray, volElemArray);

  // Compute nodal volumes
  computeNodalVolumes(nodeArray, volElemArray);

  // Compute bounding box
  double xmax = std::numeric_limits<double>::min();
  double ymax = xmax;
  double zmax = xmax;
  double xmin = std::numeric_limits<double>::max();
  double ymin = xmin;
  double zmin = xmin;
  for (auto iter = nodeArray.begin(); iter != nodeArray.end(); iter++) {
    MeshNode node = *iter;
    xmax = (xmax > node.x_) ? xmax : node.x_;
    ymax = (ymax > node.y_) ? ymax : node.y_;
    zmax = (zmax > node.z_) ? zmax : node.z_;
    xmin = (xmin < node.x_) ? xmin : node.x_;
    ymin = (ymin < node.y_) ? ymin : node.y_;
    zmin = (zmin < node.z_) ? zmin : node.z_;
  }
  Point min(xmin, ymin, zmin);
  Point max(xmax, ymax, zmax);
  Vector fudge(1.e-5,1.e-5,1.e-5);
  min = min - fudge;
  max = max + fudge;
  d_box = Box(min,max);  
}

void
AbaqusMeshGeometryPiece::readMeshNode(const std::string& inputLine,
                                      std::vector<MeshNode>& nodes)
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
    out << "Could not read nodal coordinates from " 
        << inputLine << std::endl;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  auto iter = data.begin(); 
  int node_id = std::stoi(*iter); ++iter;
  double xcoord = std::stod(*iter); ++iter;
  double ycoord = std::stod(*iter); ++iter;
  double zcoord = std::stod(*iter);

  // Save the nodal coordinates
  nodes.emplace_back(MeshNode(node_id, xcoord, ycoord, zcoord));
}

void
AbaqusMeshGeometryPiece::readMeshVolumeElement(const std::string& inputLine,
                                               std::vector<VolumeElement>& elements)
{
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() != 5) {
    std::ostringstream out;
    out << "Could not read element connectivity from input line: " 
        << inputLine << std::endl;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
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
    throw ProblemSetupException("Could not find nodes in volume element input data stream", __FILE__, __LINE__);
  }

  // Save the data
  elements.emplace_back(VolumeElement(node_list));
}

void
AbaqusMeshGeometryPiece::readMeshSurfaceElement(const std::string& inputLine,
                                                std::vector<SurfaceElement>& elements)
{
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
  auto iter = data.begin(); 
  int element_id = std::stoi(*iter); ++iter;

  // Read the element node ids
  std::vector<int> node_list;
  for (; iter != data.end(); ++iter) {
    int node_id = std::stoi(*iter);
    node_list.emplace_back(node_id);
  }
  if (node_list.empty()) {
    throw ProblemSetupException("Could not find nodes in surface element input data stream", __FILE__, __LINE__);
  }

  // Save the data
  elements.emplace_back(SurfaceElement(node_list));
}

void
AbaqusMeshGeometryPiece::computeElementVolumes(std::vector<MeshNode>& nodes,
                                               std::vector<VolumeElement>& elements)
{

  // Loop thru elements
  for (auto iter = elements.begin(); iter != elements.end(); iter++) {
    VolumeElement elem = *iter;
    MeshNode node1 = nodes[elem.node1_-1];
    MeshNode node2 = nodes[elem.node2_-1];
    MeshNode node3 = nodes[elem.node3_-1];
    MeshNode node4 = nodes[elem.node4_-1];
    Point p0(node1.x_, node1.y_, node1.z_);
    Point p1(node2.x_, node2.y_, node2.z_);
    Point p2(node3.x_, node3.y_, node3.z_);
    Point p3(node4.x_, node4.y_, node4.z_);
    Vector vec01 = p1 - p0;
    Vector vec02 = p2 - p0;
    Vector vec03 = p3 - p0;
    Vector vec1x2 = Cross(vec01, vec02);
    double volume = 1.0/6.0*std::abs(Dot(vec1x2,vec03));
    (*iter).volume_ = volume; 
  }
}

void
AbaqusMeshGeometryPiece::computeNodalVolumes(std::vector<MeshNode>& nodes,
                                             std::vector<VolumeElement>& elements)
{

  // Find the elements connected to each node
  findNodalAdjacentElements(nodes, elements);

  // Compute nodal volumes
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    MeshNode node = *node_iter;
    double node_vol = 0.0;
    for (auto elem_iter = node.adjElements_.begin(); elem_iter != node.adjElements_.end();
          ++ elem_iter) {
      node_vol += 0.25*(*elem_iter)->volume_ ;
    }
    node.volume_ = node_vol;
  }
}

void
AbaqusMeshGeometryPiece::findNodalAdjacentElements(std::vector<MeshNode>& nodes,
                                                   std::vector<VolumeElement>& elements)
{
  // Loop thru elements and find adjacent elements for each node
  for (auto elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter) {
    VolumeElement cur_elem = *elem_iter;

    // Loop thru nodes of each element
    nodes[cur_elem.node1_-1].adjElements_.push_back(&cur_elem);
    nodes[cur_elem.node2_-1].adjElements_.push_back(&cur_elem);
    nodes[cur_elem.node3_-1].adjElements_.push_back(&cur_elem);
    nodes[cur_elem.node4_-1].adjElements_.push_back(&cur_elem);
  }
}


//______________________________________________________________________
//
bool
AbaqusMeshGeometryPiece::inside(const Point& p) const
{
  //Check p with the lower coordinates
  if (p == Max(p,d_box.lower()) && p == Min(p,d_box.upper()) )
    return true;
  else
    return false;
}
//______________________________________________________________________
//
Box
AbaqusMeshGeometryPiece::getBoundingBox() const
{
  return d_box;
}


//______________________________________________________________________
//
unsigned int
AbaqusMeshGeometryPiece::createPoints()
{
  cerr << "You should be reading points .. not creating them" << endl;  
  return 0;
}
