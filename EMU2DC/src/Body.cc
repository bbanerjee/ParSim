#include <Body.h>
#include <Node.h>
#include <NodePArray.h>
#include <Element.h>
#include <CrackSP.h>
#include <Crack.h>
#include <ForceBC.h>
#include <Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Emu2DC;

   
Body::Body()
  : d_id(0), d_mat_id(0), d_initial_velocity({{0.0,0.0,0.0}})
{
  d_nodes.reserve(1000);
  d_elements.reserve(1000);
}

Body::~Body()
{
}

void 
Body::initialize(Uintah::ProblemSpecP& ps,
                 const Domain& domain,
                 const SimulationState& state,
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
  readNodeFile(input_node_file, state.dimensions());
  setInitialNodeHorizon(domain.horizon());

  // Read the input element file
  readElementFile(input_element_file);

  // Compute nodal volumes
  computeNodalVolumes();

  // Create the cell-node map for the family computer
  initializeFamilyComputer(domain);

  // Read the initial conditions for each body
  Uintah::ProblemSpecP ic_ps = ps->findBlock("InitialConditions");
  Uintah::Vector initial_velocity;
  ic_ps->require("velocity", initial_velocity);
  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_initial_velocity[ii] = initial_velocity[ii];
  }

  // Get the initial crack information
  for (Uintah::ProblemSpecP crack_ps = ic_ps->findBlock("Crack"); crack_ps != 0;
       crack_ps = crack_ps->findNextBlock("Crack")) {
    CrackSP crack = std::make_shared<Crack>();
    crack->initialize(crack_ps); 
    d_cracks.emplace_back(crack);
  }

  // Read the external force boundary conditions (displacement/velocity bcs may be added later)
  // The BC information is used to update the external force on particles/nodes.
  Uintah::ProblemSpecP bc_ps = ps->findBlock("BoundaryConditions");
  for (Uintah::ProblemSpecP force_ps = bc_ps->findBlock("ForceBC"); force_ps != 0;
       force_ps = force_ps->findNextBlock("ForceBC")) {
    ForceBC ext_force;
    ext_force.initialize(force_ps, d_nodes); 
  }
    
}

void
Body::readNodeFile(const std::string& fileName, const int dim)
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
    int node_id, boundary_node;
    double xcoord, ycoord, zcoord;
    if (!(data_stream >> node_id >> xcoord >> ycoord >> zcoord >> boundary_node)) {
      throw Exception("Could not read node input data stream", __FILE__, __LINE__);
    }

    // Save the data
    NodeP node(new Node(node_id, xcoord, ycoord, zcoord, boundary_node));
    node->dimension(dim);
    d_nodes.emplace_back(node);

    // Add to the node ID -> node ptr map
    d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
  }
}

void
Body::setInitialNodeHorizon(const double horizon)
{
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    (*iter)->horizonSize(horizon);
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
    elem->computeVolume();
    d_elements.emplace_back(elem);
  }

  // Loop thru elements and find adjacent elements for each node
  for (auto elem_iter = d_elements.begin(); elem_iter != d_elements.end(); ++elem_iter) { 
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


// general method to calculate the volume of each node in any type of grid (uniform and non-uniform)
void 
Body::computeNodalVolumes()
{
  for (auto node_iter = d_nodes.begin(); node_iter != d_nodes.end(); ++node_iter) {

    const ElementPArray& adj_elems = (*node_iter)->getAdjacentElements();
    double vol = 0.0;

    for (auto elem_iter = adj_elems.begin(); elem_iter != adj_elems.end(); ++elem_iter) {

      double elem_vol = (*elem_iter)->volume();
      int num_nodes = (*elem_iter)->numNodes();
      vol += elem_vol/(double) num_nodes;
    }

    (*node_iter)->volume(vol);
  }
}

void 
Body::initializeFamilyComputer(const Domain& domain)
{
  if (d_nodes.size() == 0) {
    std::ostringstream out;
    out << "The body has to be discretized into nodes before a cell-node map can be created" << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }  
  d_family_computer.createCellNodeMap(domain, d_nodes);
  //d_family_computer.printCellNodeMap();
}

void 
Body::createInitialFamily(const Domain& domain)
{
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    NodeP cur_node = *iter;
    NodePArray neighbor_list;
    d_family_computer.getInitialFamily(cur_node, domain, neighbor_list);
    cur_node->setFamily(neighbor_list);
    cur_node->initialFamilySize((int) neighbor_list.size());
  }
}

void 
Body::updateFamily(const Domain& domain)
{
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    NodeP cur_node = *iter;
    NodePArray neighbor_list;
    d_family_computer.getCurrentFamily(cur_node, domain, neighbor_list);
    cur_node->setFamily(neighbor_list);
  }
}

void 
Body::printFamily()
{
  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    NodeP cur_node = *iter;
    NodePArray neighbor_list = cur_node->getFamily();
    int count = 0;
    std::cout << "Current node = " << *cur_node << std::endl;
    for (constNodePIterator fam_iter = neighbor_list.begin(); fam_iter != neighbor_list.end(); fam_iter++) {
      std::cout << " Neigbor node (" << ++count << ") = "<< *(*fam_iter) << std::endl;
    } 
  }
}

void
Body::removeBondsIntersectedByCracks()
{
  // Loop through body nodes
  for (auto node_iter = d_nodes.begin(); node_iter != d_nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    NodePArray neighbor_list = cur_node->getFamily();

    // Loop through cracks in the body
    for (auto crack_iter = d_cracks.begin(); crack_iter != d_cracks.end(); ++crack_iter) {
      (*crack_iter)->breakBonds(cur_node, neighbor_list);
    }
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
    for (auto iter = (body.d_cracks).begin(); iter != (body.d_cracks).end(); ++iter) {
      out << *(*iter) << std::endl ;
    }
    return out;
  }
}
