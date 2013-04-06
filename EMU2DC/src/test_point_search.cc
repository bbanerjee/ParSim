// Test point search with hashed unordered_multimap containers
#include <Domain.h>
#include <Node.h>
#include <NodeP.h>
#include <NodePArray.h>
#include <CellNodePMap.h>
#include <map>
#include <tr1/unordered_map>
#include <vector>
#include <cstdint>
#include <ctime>
#include <random>
#include <memory>
#include <chrono>

using namespace Emu2DC;

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
void generate_random_points(const Domain& domain, const int& numPoints, NodeArray& nodes);
void generate_random_points(const Domain& domain, const int& numPoints, NodePArray& nodes);
void test_map_node();
void test_map_node_pointer();
void test_multimap_node_pointer();
void test_key_generation();
bool bit(int ii, int jj);
void test_random_number_generator();
void test_random_point_generator();
void test_unorderedmap_node_pointer();
void test_unorderedmultimap_node_pointer();

int main()
{
  auto t1 = std::chrono::high_resolution_clock::now();
  test_point_search_raw_pointer();
  auto t2 = std::chrono::high_resolution_clock::now();
  test_point_search_shared_pointer();
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "Node*: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " NodeP: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() << std::endl;
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
  Array3 lower = {{xmin, ymin, zmin}};
  Array3 upper = {{xmax, ymax, zmax}};

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
    Array3 pos = node->position();
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
  Array3 lower = {{xmin, ymin, zmin}};
  Array3 upper = {{xmax, ymax, zmax}};

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
    Array3 pos = node->position();
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

void generate_random_points(const Domain& domain, const int& numPoints, NodeArray& nodeList)
{
  double xmin = (domain.lower())[0];
  double ymin = (domain.lower())[1];
  double zmin = (domain.lower())[2];
  double xmax = (domain.upper())[0];
  double ymax = (domain.upper())[1];
  double zmax = (domain.upper())[2];

  // set up uniformly distributed random numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<double> xrand(xmin, xmax);
  std::uniform_real_distribution<double> yrand(ymin, ymax);
  std::uniform_real_distribution<double> zrand(zmin, zmax);

  // Generate points
  for (int pt = 0;  pt < numPoints; ++pt) {
    Array3 pos = {{xrand(rand_gen), yrand(rand_gen), zrand(rand_gen)}};
    //std::cout << "pt(" << pt <<") = [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << std::endl;
    Node* node = new Node();
    long64 partID = (((long64) pos[0]<<16)|((long64) pos[1]<<32)|(long64)pos[2]<<48)|(long64)(pt+1);
    node->setID(partID);
    node->position(pos);
    nodeList.push_back(node); 
  }
}

void generate_random_points(const Domain& domain, const int& numPoints, NodePArray& nodeList)
{
  double xmin = (domain.lower())[0];
  double ymin = (domain.lower())[1];
  double zmin = (domain.lower())[2];
  double xmax = (domain.upper())[0];
  double ymax = (domain.upper())[1];
  double zmax = (domain.upper())[2];

  // set up uniformly distributed random numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<double> xrand(xmin, xmax);
  std::uniform_real_distribution<double> yrand(ymin, ymax);
  std::uniform_real_distribution<double> zrand(zmin, zmax);

  // Generate points
  for (int pt = 0;  pt < numPoints; ++pt) {
    Array3 pos = {{xrand(rand_gen), yrand(rand_gen), zrand(rand_gen)}};
    //std::cout << "pt(" << pt <<") = [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << std::endl;
    NodeP node = std::make_shared<Node>();
    long64 partID = (((long64) pos[0]<<16)|((long64) pos[1]<<32)|(long64)pos[2]<<48)|(long64)(pt+1);
    node->setID(partID);
    node->position(pos);
    //nodeList.push_back(node); 
    nodeList.emplace_back(node); 
  }
}

void test_map_node()
{
  typedef std::map<Node, Node> Bond;
  typedef std::map<Node, Node>::iterator BondIterator;
  typedef std::map<Node, Node>::const_iterator constBondIterator;

  typedef std::vector<Node> NodeArray;
  typedef std::vector<Node>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  for (int ii = 0; ii < 10; ii++) {
    Node node;
    Array3 pos = {{(double)ii, 0.0, 0.0}};
    node.setID(ii);
    node.position(pos);
    nodelist.push_back(node); 
  }

  // Create two bonds
  Bond bondList;
  for (int ii = 0; ii < 9; ii++) {
    bondList.insert(std::pair<Node, Node>(nodelist[ii], nodelist[ii+1]));
  }

  for (constBondIterator iter = bondList.begin(); iter != bondList.end(); iter++) {
    std::cout << "Master node = " << iter->first << " Slave Node = " << iter->second << std::endl;
  }

}

void test_map_node_pointer()
{
  typedef std::map<Node*, Node*> Bond;
  typedef std::map<Node*, Node*>::iterator BondIterator;
  typedef std::map<Node*, Node*>::const_iterator constBondIterator;

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  for (int ii = 0; ii < 10; ii++) {
    Node* node = new Node();
    Array3 pos = {{(double)ii, 0.0, 0.0}};
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // Create two bonds
  Bond bondList;
  for (int ii = 0; ii < 9; ii++) {
    bondList.insert(std::pair<Node*,Node*>(nodelist[ii], nodelist[ii+1]));
  }

  for (constBondIterator iter = bondList.begin(); iter != bondList.end(); iter++) {
    std::cout << "Master node = " << *(iter->first) << " Slave Node = " << *(iter->second) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_multimap_node_pointer()
{
  typedef std::multimap<Node*, Node*> Bond;
  typedef std::multimap<Node*, Node*>::iterator BondIterator;
  typedef std::multimap<Node*, Node*>::const_iterator constBondIterator;

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  for (int ii = 0; ii < 10; ii++) {
    Node* node = new Node();
    Array3 pos = {{(double)ii, 0.0, 0.0}};
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // Create two bonds
  Bond bondList;
  for (int ii = 1; ii < 10; ii++) {
    bondList.insert(std::pair<Node*,Node*>(nodelist[0], nodelist[ii]));
  }

  for (constBondIterator iter = bondList.begin(); iter != bondList.end(); iter++) {
    std::cout << "Master node = " << *(iter->first) << " Slave Node = " << *(iter->second) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_key_generation()
{
  typedef int64_t long64;

  // Divide the domain into m x n x p cells
  int mm = 3;
  int nn = 2;
  int pp = 1;
  int nbits = 31;
  for (int ii = 0; ii < mm; ii++) {
    for (int jj = 0; jj < nn; jj++) {
      for (int kk = 0; kk < pp; kk++) {
        std::cout << "Cell = ["<< ii << "," << jj << "," << kk << "];" ;
        long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
        short int particleID = 10;
        std::cout << "cellID = " << cellID << " + Particle = " << (cellID | (long64) particleID) << std::endl;

        long64 placebit = (long64)1<<nbits;
        //std::cout << "placebit = " << placebit << std::endl;
        long64 bitsum = 0;
        for (int bb=0; bb < nbits; bb++) {
          long64 i_bit = ((long64)ii & ((long64)1 << bb));
          long64 j_bit = ((long64)jj & ((long64)1 << bb));
          long64 k_bit = ((long64)kk & ((long64)1 << bb));
          //std::cout << "bit(ix,j) = " << i_bit << " bit(iy,j) = " << j_bit << " bit(iz,j) = " << k_bit 
          //          << std::endl;
          //bitsum += ((long64)1<<(3*bb))*((i_bit<<2)+(j_bit<<1)+k_bit);
          bitsum += ((long64)1<<(3*bb))*((i_bit<<2)|(j_bit<<1)|k_bit);
        } 
        //std::cout << "bitsum = " << bitsum << std::endl;
        std::cout << "key = " << placebit+bitsum << std::endl;
      }
    }
  }
}

// get the j-th bit from i
bool bit(int ii, int jj)
{
  return ii & (1 << jj);
}

void test_random_number_generator()
{
  // set up a box domain
  //Array3 lower = {{0.0, 0.0, 0.0}};
  //Array3 upper = {{3.0, 2.0, 1.0}};

  // seed the random number generator
  std::srand(0);
  int num_points = 10;
  for (int ii = 0; ii < num_points; ii++) {
    std::cout << std::rand() << ", " ;
  }
  std::cout << std::endl;

  // Alternative uniform distribution from <random>
  unsigned int seed = 1;
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(0.0,3.0);

  std::cout << "some random numbers between 1 and 3: ";
  for (int i=0; i<10; ++i)
    std::cout << distribution(generator) << " ";

  std::cout << std::endl;
}

void test_random_point_generator()
{
  // set up a box domain
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 3.0, ymax = 2.0, zmax = 1.0;
  //Array3 lower = {{0.0, 0.0, 0.0}};
  //Array3 upper = {{3.0, 2.0, 1.0}};

  // set up uniformly distributed random numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_real_distribution<double> xrand(xmin, xmax);
  std::uniform_real_distribution<double> yrand(ymin, ymax);
  std::uniform_real_distribution<double> zrand(zmin, zmax);

  // Generate 20 points
  int num_points = 20;
  for (int pt = 0;  pt < num_points; ++pt) {
    Array3 node = {{xrand(rand_gen), yrand(rand_gen), zrand(rand_gen)}};
    std::cout << "pt(" << pt <<") = [" << node[0] << ", " << node[1] << ", " << node[2] << "]" << std::endl;
  }
}

void test_unorderedmap_node_pointer()
{
  typedef int64_t long64;
  typedef std::tr1::unordered_map<long64, Node*> CellNodeMap;
  typedef CellNodeMap::iterator CellNodeMapIterator;
  typedef CellNodeMap::const_iterator constCellNodeMapIterator;

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  int num_nodes = 10;
  for (int ii = 0; ii < num_nodes; ii++) {
    Node* node = new Node();
    Array3 pos = {{(double)ii, 0.0, 0.0}};
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // set up a list to get uniformly distributed random node numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_int_distribution<int> dist(0, num_nodes-1);

  // Create a unordered map
  CellNodeMap cell_node_map;
  int nx = 3, ny = 2, nz = 1;
  for (int ii = 0; ii < nx; ii++) {
    for (int jj = 0; jj < ny; jj++) {
      for (int kk = 0; kk < nz; kk++) {
        long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
        int nodeIndex = dist(rand_gen);
        cell_node_map[cellID] = nodelist[nodeIndex];
        std::cout << "cell = " << cellID << " node = " << *(nodelist[nodeIndex]) ;
        nodeIndex = dist(rand_gen);
        cell_node_map[cellID] = nodelist[nodeIndex];
        std::cout << " node = " << *(nodelist[nodeIndex]) << std::endl ;
      }
    }
  }

  // Read the data for a particular cell, say, [2, 1, 0];
  long64 cell210 = ((long64)2 << 16) | ((long64)1 << 32) | ((long64)0 << 48);
  CellNodeMapIterator it = cell_node_map.find(cell210);
  if (it != cell_node_map.end()) {
    std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  } else {
    std::cout << "key not found " << std::endl;
  }

  // Print out all the data
  for (constCellNodeMapIterator it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_unorderedmultimap_node_pointer()
{

  // Create ten nodes
  NodeArray nodelist;
  int num_nodes = 10;
  for (int ii = 0; ii < num_nodes; ii++) {
    Node* node = new Node();
    Array3 pos = {{(double)ii, 0.0, 0.0}};
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // set up a list to get uniformly distributed random node numbers
  unsigned int seed = 1;
  std::default_random_engine rand_gen(seed);
  std::uniform_int_distribution<int> dist(0, num_nodes-1);

  // Create a unordered map
  CellNodeMap cell_node_map;
  int nx = 3, ny = 2, nz = 1;
  for (int ii = 0; ii < nx; ii++) {
    for (int jj = 0; jj < ny; jj++) {
      for (int kk = 0; kk < nz; kk++) {
        long64 cellID = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
        int nodeIndex = dist(rand_gen);
        cell_node_map.insert(CellNodePair(cellID, nodelist[nodeIndex]));
        std::cout << "cell = " << cellID << " node = " << *(nodelist[nodeIndex]) ;
        nodeIndex = dist(rand_gen);
        cell_node_map.insert(CellNodePair(cellID, nodelist[nodeIndex]));
        std::cout << " node = " << *(nodelist[nodeIndex]) << std::endl ;
      }
    }
  }

  // Read the data for a particular cell, say, [2, 1, 0];
  long64 cell210 = ((long64)2 << 16) | ((long64)1 << 32) | ((long64)0 << 48);
  CellNodePairIterator nodes = cell_node_map.equal_range(cell210);
  for (auto it = nodes.first; it != nodes.second; ++it) {
    std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  }

  // Print out all the data
  for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  }

  // Check the buckets
  std::cout << "cell_node_map's buckets contain:\n";
  for ( unsigned i = 0; i < cell_node_map.bucket_count(); ++i) {
    std::cout << "bucket #" << i << " contains:";
    for ( auto local_it = cell_node_map.begin(i); local_it!= cell_node_map.end(i); ++local_it )
      std::cout << " " << local_it->first << ":" << *(local_it->second);
    std::cout << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}
