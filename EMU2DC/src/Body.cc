#include <Body.h>
#include <Node.h>
#include <Element.h>
#include <Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Emu2DC;

   
Body::Body()
  : d_id(0), d_mat_id(0), d_nodes(0), d_elements(0)
{
}

Body::~Body()
{
}

void 
Body::initialize(Uintah::ProblemSpecP& ps,
                 const MaterialSPArray& matList)
{
  if (!ps) return;
  
  // Find the material name
  Uintah::ProblemSpecP mat_ps = ps->findBlock("material");
  if (!mat_ps) {
    throw Exception("A body must have a material assigned to it", __FILE__, __LINE__);
  }

  std::string mat_name;
  if (!(mat_ps->getAttribute("name", mat_name))) {
    throw Exception("A material assigned to a body does not have a name", __FILE__, __LINE__);
  }

  // Check if the name corresponds to an actual material
  bool found_mat_name = false;
  for (auto iter = matList.begin(); iter != matList.end(); iter++) {
    Material* mat = (*iter).get();
    if (mat_name == mat->name()) {
      found_mat_name = true;
      d_mat_id = mat->id();
      break;
    }
  }
  if (!found_mat_name) {
    throw Exception("Material assigned to the body does not have a known name", __FILE__, __LINE__);
  }

  // Get the geometry (from input node and element files)
  Uintah::ProblemSpecP geom_ps = ps->findBlock("Geometry");
  std::string input_node_file;
  std::string input_element_file;
  geom_ps->require("input_node_file", input_node_file);
  geom_ps->require("input_element_file", input_element_file);
  std::cout << "Input geometry files: " << input_node_file << ", " << input_element_file << std::endl;

  // Read the input node file
  readNodeFile(input_node_file);

  // Read the input element file
  readElementFile(input_element_file);

}

void
Body::readNodeFile(const std::string& fileName)
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open node input file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  std::string line;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
         std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    // Skip comment lines
    if (line[0] == '#') continue;

    // Read the data
    std::istringstream data_stream(line);
    int node_id, hanging;
    double xcoord, ycoord, zcoord;
    if (!(data_stream >> node_id >> xcoord >> ycoord >> zcoord >> hanging)) {
      throw Exception("Could not read node input data stream", __FILE__, __LINE__);
    }

    // Save the data
    NodeP node(new Node(node_id, xcoord, ycoord, zcoord, hanging));
    d_nodes.emplace_back(node);

    // Add to the node ID -> node ptr map
    d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
  }
}

void
Body::readElementFile(const std::string& fileName)
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

    // Ignore empty lines
    if (line.empty()) continue;

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
         std::not1(std::ptr_fun<int, int>(std::isspace))));
    
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

    // Find the node pointers
    if (d_id_ptr_map.empty()) {
      throw Exception("Could not find node id -> node ptr map", __FILE__, __LINE__);
    }
    NodePArray nodes;
    for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
      int node_id = *iter;
      auto id_ptr_pair = d_id_ptr_map.find(node_id);
      if (id_ptr_pair == d_id_ptr_map.end()) {
        std::string out = "Could not find node id -> node ptr pair for node " + node_id;
        throw Exception(out, __FILE__, __LINE__);
      }
      NodeP it = id_ptr_pair->second;
      nodes.emplace_back(it); 
    }
     
    // Save the data
    ElementP elem(new Element(element_id, nodes));
    d_elements.emplace_back(elem);
  }
}


namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Body& body)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Body:" << body.d_id << std::endl;
    out << "  Material = " << body.d_mat_id << std::endl;
    for (auto iter = (body.d_nodes).begin(); iter != (body.d_nodes).end(); ++iter) {
      out << *(*iter) << std::endl ;
    }
    for (auto iter = (body.d_elements).begin(); iter != (body.d_elements).end(); ++iter) {
      out << *(*iter) << std::endl ;
    }
    return out;
  }
}
