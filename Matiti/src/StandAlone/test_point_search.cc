// Test point search with hashed unordered_multimap containers
#include <Core/Domain.h>
#include <Core/Node.h>
#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Types/CellNodePMap.h>
#include <Core/FamilyComputer.h>
#include <map>
#include <tr1/unordered_map>
#include <vector>
#include <cstdint>
#include <ctime>
#include <random>
#include <memory>
#include <chrono>

using namespace Matiti;

typedef int64_t long64;
typedef std::tr1::unordered_multimap<long64, Node*, Hash64> CellNodeMap;
// typedef std::tr1::unordered_multimap<long64, Node*, HashCell> CellNodeMap;
typedef CellNodeMap::iterator CellNodeMapIterator;
typedef CellNodeMap::const_iterator constCellNodeMapIterator;
typedef std::pair<long64, Node*> CellNodePair; 
typedef std::pair<CellNodeMapIterator, CellNodeMapIterator> CellNodePairIterator; 
typedef std::pair<constCellNodeMapIterator, constCellNodeMapIterator> constCellNodePairIterator; 

typedef std::vector<Node*> NodeArray;
typedef std::vector<Node*>::iterator NodeIterator;
typedef std::vector<Node*>::const_iterator constNodeIterator;

//typedef std::shared_ptr<Node> NodeP;
typedef std::tr1::unordered_multimap<long64, NodeP, Hash64> CellNodePMap;
// typedef std::tr1::unordered_multimap<long64, NodeP, HashCell> CellNodePMap;
typedef CellNodePMap::iterator CellNodePMapIterator;
typedef CellNodePMap::const_iterator constCellNodePMapIterator;
typedef std::pair<long64, NodeP> CellNodePPair; 
typedef std::pair<CellNodePMapIterator, CellNodePMapIterator> CellNodePPairIterator; 
typedef std::pair<constCellNodePMapIterator, constCellNodePMapIterator> constCellNodePPairIterator; 

//typedef std::vector<NodeP> NodePArray;
//typedef std::vector<NodeP>::iterator NodePIterator;
//typedef std::vector<NodeP>::const_iterator constNodePIterator;

void test_point_search_raw_pointer(); 
void test_point_search_shared_pointer(); 
void test_point_search_FamilyComputer();
void generate_random_points(const Domain& domain, const int& numPoints, NodeArray& nodes);
void generate_random_points(const Domain& domain, const int& numPoints, NodePArray& nodes);

int main()
{
  auto t1 = std::chrono::high_resolution_clock::now();
  test_point_search_raw_pointer();
  auto t2 = std::chrono::high_resolution_clock::now();
  test_point_search_shared_pointer();
  auto t3 = std::chrono::high_resolution_clock::now();
  test_point_search_FamilyComputer(); 
  auto t4 = std::chrono::high_resolution_clock::now();
  std::cout << "Node*: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " NodeP: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() 
            << " Family: " << std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count() << std::endl;
  return 0;
}

void test_point_search_raw_pointer() 
{
  // inputs
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 3.0, ymax = 2.0, zmax = 1.0;
  double horizon = 0.5;
  double num_points = 100;

  // set up a box domain with horizon
  Point3D lower(xmin, ymin, zmin);
  Point3D upper(xmax, ymax, zmax);

  Domain domain(lower, upper, horizon);
  std::cout << domain;

  // Create a set of randomly distributed points in the box
  NodeArray node_list;
  generate_random_points(domain, num_points, node_list);

  // Find which cells the points sit in and create a unordered map that maps nodes to cells
  CellNodeMap cell_node_map;
  std::cout << cell_node_map.bucket_count() << std::endl;
  int count = 0;
  for (constNodeIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    //std::cout << " Node (" << ++count << ") = "<< *(*iter);
    IntArray3 cell;
    Node* node = *iter;
    const Point3D& pos = node->position();
    domain.findCellIndex(pos, cell);
    //std::cout << " Cell = [" << cell[0] << ", " << cell[1] << ", " << cell[2] << "]" << std::endl;
    long64 cellID = ((long64)cell[0] << 16) | ((long64)cell[1] << 32) | ((long64)cell[2] << 48);
    cell_node_map.insert(CellNodePair(cellID, node));
  }
  //std::cout << cell_node_map.bucket_count() << std::endl;
  //cell_node_map.rehash((size_t)domain.totalCells());
  //std::cout << cell_node_map.bucket_count() << std::endl;

  // Read the data for a particular cell, say, [2, 1, 1];
  //long64 cell211 = ((long64)2 << 16) | ((long64)1 << 32) | ((long64)1 << 48);
  //CellNodePairIterator nodes = cell_node_map.equal_range(cell211);
  //for (auto it = nodes.first; it != nodes.second; ++it) {
  //  std::cout << "Cell (2,1,1): key = " << it->first << " value = " << *(it->second) << std::endl;
  //}

  // Print out all the data
  //for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
  //  std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  //}

  // Check the buckets
  //std::cout << "cell_node_map's buckets contain:\n";
  //for ( unsigned i = 0; i < cell_node_map.bucket_count(); ++i) {
  //  std::cout << "bucket #" << i << " contains:";
  //  for ( auto local_it = cell_node_map.begin(i); local_it!= cell_node_map.end(i); ++local_it ) {
  //    std::cout << " " << local_it->first << ":" << *(local_it->second);
  //  }
  //  std::cout << std::endl;
  //}

  // Loop through nodes and find cell range in horizon
  IntArray3 num_cells = domain.numCells();
  for (NodeIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    Node* cur_node = *iter;
    IntArray3 cur_cell;
    domain.findCellIndex(cur_node->position(), cur_cell);
    int iimin = std::max(1, cur_cell[0]-1);
    int jjmin = std::max(1, cur_cell[1]-1);
    int kkmin = std::max(1, cur_cell[2]-1);
    int iimax = std::min(cur_cell[0]+1, num_cells[0]);
    int jjmax = std::min(cur_cell[1]+1, num_cells[1]);
    int kkmax = std::min(cur_cell[2]+1, num_cells[2]);
    IntArray3 cell_min({{iimin, jjmin, kkmin}});
    IntArray3 cell_max({{iimax, jjmax, kkmax}});
    std::cout << "Node = " << *cur_node << " Min cell = [" << iimin << ", " << jjmin << ", " << kkmin << "]";
    std::cout << " Max cell = [" << iimax << ", " << jjmax << ", " << kkmax << "]" << std::endl;

    // Find the nodes inside the cells within the range
    NodeArray neighbor_list;
    for (int ii=iimin; ii <= iimax; ++ii) {
      for (int jj=jjmin; jj <= jjmax; ++jj) {
        for (int kk=kkmin; kk <= kkmax; ++kk) {
          long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
          CellNodePairIterator nodes = cell_node_map.equal_range(cellID);
          for (auto it = nodes.first; it != nodes.second; ++it) {
            Node* near_node = it->second;
            if (cur_node == near_node) continue;
            //double dist = cur_node->distance(*near_node);
            //std::cout << "Cell (" << ii <<"," << jj << "," << kk 
            //          << ") : key = " << it->first << " value = " << *(it->second) 
            //          << " distance = " << dist ;
            if (cur_node->distance(*near_node) < horizon) {
              neighbor_list.push_back(near_node);
	      //std::cout << " in." << std::endl;
            } //else {
	      //std::cout << " out." << std::endl;
            //}
          }
        }
      }
    }
    count = 0;
    for (constNodeIterator iter = neighbor_list.begin(); iter != neighbor_list.end(); iter++) {
      std::cout << " Neigbor node (" << ++count << ") = "<< *(*iter) << std::endl;
    }
  }

  // Delete the nodes
  for (NodeIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    delete *iter;
  }
}

void test_point_search_shared_pointer() 
{
  // inputs
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 3.0, ymax = 2.0, zmax = 1.0;
  double horizon = 0.5;
  double num_points = 100;

  // set up a box domain with horizon
  Point3D lower(xmin, ymin, zmin);
  Point3D upper(xmax, ymax, zmax);

  Domain domain(lower, upper, horizon);
  std::cout << domain;

  // Create a set of randomly distributed points in the box
  NodePArray node_list;
  generate_random_points(domain, num_points, node_list);

  // Find which cells the points sit in and create a unordered map that maps nodes to cells
  CellNodePMap cell_node_map;
  std::cout << cell_node_map.bucket_count() << std::endl;
  int count = 0;
  for (constNodePIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    //std::cout << " Node (" << ++count << ") = "<< *(*iter);
    IntArray3 cell;
    NodeP node = *iter;
    const Point3D& pos = node->position();
    domain.findCellIndex(pos, cell);
    //std::cout << " Cell = [" << cell[0] << ", " << cell[1] << ", " << cell[2] << "]" << std::endl;
    long64 cellID = ((long64)cell[0] << 16) | ((long64)cell[1] << 32) | ((long64)cell[2] << 48);
    cell_node_map.insert(CellNodePPair(cellID, node));
  }
  //std::cout << cell_node_map.bucket_count() << std::endl;
  //cell_node_map.rehash((size_t)domain.totalCells());
  //std::cout << cell_node_map.bucket_count() << std::endl;

  // Read the data for a particular cell, say, [2, 1, 1];
  //long64 cell211 = ((long64)2 << 16) | ((long64)1 << 32) | ((long64)1 << 48);
  //CellNodePPairIterator nodes = cell_node_map.equal_range(cell211);
  //for (auto it = nodes.first; it != nodes.second; ++it) {
  //  std::cout << "Cell (2,1,1): key = " << it->first << " value = " << *(it->second) << std::endl;
  //}

  // Print out all the data
  //for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
  //  std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  //}

  // Check the buckets
  //std::cout << "cell_node_map's buckets contain:\n";
  //for ( unsigned i = 0; i < cell_node_map.bucket_count(); ++i) {
  //  std::cout << "bucket #" << i << " contains:";
  //  for ( auto local_it = cell_node_map.begin(i); local_it!= cell_node_map.end(i); ++local_it ) {
  //    std::cout << " " << local_it->first << ":" << *(local_it->second);
  //  }
  //  std::cout << std::endl;
  //}

  // Loop through nodes and find cell range in horizon
  IntArray3 num_cells = domain.numCells();
  for (NodePIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    NodeP cur_node = *iter;
    IntArray3 cur_cell;
    domain.findCellIndex(cur_node->position(), cur_cell);
    int iimin = std::max(1, cur_cell[0]-1);
    int jjmin = std::max(1, cur_cell[1]-1);
    int kkmin = std::max(1, cur_cell[2]-1);
    int iimax = std::min(cur_cell[0]+1, num_cells[0]);
    int jjmax = std::min(cur_cell[1]+1, num_cells[1]);
    int kkmax = std::min(cur_cell[2]+1, num_cells[2]);
    IntArray3 cell_min({{iimin, jjmin, kkmin}});
    IntArray3 cell_max({{iimax, jjmax, kkmax}});
    std::cout << "Node = " << *cur_node << " Min cell = [" << iimin << ", " << jjmin << ", " << kkmin << "]";
    std::cout << " Max cell = [" << iimax << ", " << jjmax << ", " << kkmax << "]" << std::endl;

    // Find the nodes inside the cells within the range
    NodePArray neighbor_list;
    for (int ii=iimin; ii <= iimax; ++ii) {
      for (int jj=jjmin; jj <= jjmax; ++jj) {
        for (int kk=kkmin; kk <= kkmax; ++kk) {
          long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
          CellNodePPairIterator nodes = cell_node_map.equal_range(cellID);
          for (auto it = nodes.first; it != nodes.second; ++it) {
            NodeP near_node = it->second;
            if (cur_node == near_node) continue;
            //double dist = cur_node->distance(*near_node);
            //std::cout << "Cell (" << ii <<"," << jj << "," << kk 
            //          << ") : key = " << it->first << " value = " << *(it->second) 
            //          << " distance = " << dist ;
            if (cur_node->distance(*near_node) < horizon) {
              neighbor_list.push_back(near_node);
	      //std::cout << " in." << std::endl;
            } //else {
	      //std::cout << " out." << std::endl;
            //}
          }
        }
      }
    }
    count = 0;
    for (constNodePIterator iter = neighbor_list.begin(); iter != neighbor_list.end(); iter++) {
      std::cout << " Neigbor node (" << ++count << ") = "<< *(*iter) << std::endl;
    }
  }
}

// Test the FamilyComputer classss
void test_point_search_FamilyComputer() 
{
  // inputs
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 3.0, ymax = 2.0, zmax = 1.0;
  double horizon = 0.5;
  double num_points = 100;

  // set up a box domain with horizon
  Point3D lower(xmin, ymin, zmin);
  Point3D upper(xmax, ymax, zmax);

  Domain domain(lower, upper, horizon);
  std::cout << domain;

  // Create a set of randomly distributed points in the box
  NodePArray node_list;
  generate_random_points(domain, num_points, node_list);

  // Create a family computer
  FamilyComputer fc;
  fc.createCellNodeMap(domain, node_list);
  fc.printCellNodeMap();

  // Loop through nodes and find cell range in horizon
  for (NodePIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    NodeP cur_node = *iter;
    NodePArray neighbor_list;
    fc.getInitialFamily(cur_node, domain, neighbor_list);
    int count = 0;
    for (constNodePIterator iter = neighbor_list.begin(); iter != neighbor_list.end(); iter++) {
      std::cout << " Neigbor node (" << ++count << ") = "<< *(*iter) << std::endl;
    }
  }
}


void generate_random_points(const Domain& domain, const int& numPoints, NodeArray& nodeList)
{
  double xmin = domain.lower().x();
  double ymin = domain.lower().y();
  double zmin = domain.lower().z();
  double xmax = domain.upper().x();
  double ymax = domain.upper().y();
  double zmax = domain.upper().z();

  // set up uniformly distributed random numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<double> xrand(xmin, xmax);
  std::uniform_real_distribution<double> yrand(ymin, ymax);
  std::uniform_real_distribution<double> zrand(zmin, zmax);

  // Generate points
  for (int pt = 0;  pt < numPoints; ++pt) {
    Point3D pos(xrand(rand_gen), yrand(rand_gen), zrand(rand_gen));
    //std::cout << "pt(" << pt <<") = [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << std::endl;
    Node* node = new Node();
    long64 partID = (((long64) pos.x()<<16)|((long64) pos.y()<<32)|(long64)pos.z()<<48)|(long64)(pt+1);
    node->setID(partID);
    node->position(pos);
    nodeList.push_back(node); 
  }
}

void generate_random_points(const Domain& domain, const int& numPoints, NodePArray& nodeList)
{
  double xmin = domain.lower().x();
  double ymin = domain.lower().y();
  double zmin = domain.lower().z();
  double xmax = domain.upper().x();
  double ymax = domain.upper().y();
  double zmax = domain.upper().z();

  // set up uniformly distributed random numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<double> xrand(xmin, xmax);
  std::uniform_real_distribution<double> yrand(ymin, ymax);
  std::uniform_real_distribution<double> zrand(zmin, zmax);

  // Generate points
  for (int pt = 0;  pt < numPoints; ++pt) {
    Point3D pos(xrand(rand_gen), yrand(rand_gen), zrand(rand_gen));
    //std::cout << "pt(" << pt <<") = [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << std::endl;
    NodeP node = std::make_shared<Node>();
    long64 partID = (((long64) pos.x()<<16)|((long64) pos.y()<<32)|(long64)pos.z()<<48)|(long64)(pt+1);
    node->setID(partID);
    node->position(pos);
    //nodeList.push_back(node); 
    nodeList.emplace_back(node); 
  }
}

