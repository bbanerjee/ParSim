// Test point search with hashed unordered_multimap containers
#include <Domain.h>
#include <Node.h>
#include <map>
#include <tr1/unordered_map>
#include <vector>
#include <cstdint>
#include <ctime>
#include <random>

using namespace Emu2DC;

// lookup3 hash function
typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef uint64_t u128;

inline u32 lookup3( const u8 *key, u32 len, u32 seed ) {
  #if defined(_MSC_VER)
    #define rot(x,k) _rotl(x,k)
  #else
    #define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
  #endif

  #define mix(a,b,c) \
  { \
    a -= c;  a ^= rot(c, 4);  c += b; \
    b -= a;  b ^= rot(a, 6);  a += c; \
    c -= b;  c ^= rot(b, 8);  b += a; \
    a -= c;  a ^= rot(c,16);  c += b; \
    b -= a;  b ^= rot(a,19);  a += c; \
    c -= b;  c ^= rot(b, 4);  b += a; \
  }

  #define final(a,b,c) \
  { \
    c ^= b; c -= rot(b,14); \
    a ^= c; a -= rot(c,11); \
    b ^= a; b -= rot(a,25); \
    c ^= b; c -= rot(b,16); \
    a ^= c; a -= rot(c,4);  \
    b ^= a; b -= rot(a,14); \
    c ^= b; c -= rot(b,24); \
  }

  u32 a, b, c;
  a = b = c = 0xdeadbeef + len + seed;

  const u32 *k = (const u32 *)key;

  while ( len > 12 ) {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    len -= 12;
    k += 3;
  }

  switch( len ) {
    case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
    case 11: c+=k[2]&0xffffff; b+=k[1]; a+=k[0]; break;
    case 10: c+=k[2]&0xffff; b+=k[1]; a+=k[0]; break;
    case 9 : c+=k[2]&0xff; b+=k[1]; a+=k[0]; break;
    case 8 : b+=k[1]; a+=k[0]; break;
    case 7 : b+=k[1]&0xffffff; a+=k[0]; break;
    case 6 : b+=k[1]&0xffff; a+=k[0]; break;
    case 5 : b+=k[1]&0xff; a+=k[0]; break;
    case 4 : a+=k[0]; break;
    case 3 : a+=k[0]&0xffffff; break;
    case 2 : a+=k[0]&0xffff; break;
    case 1 : a+=k[0]&0xff; break;
    case 0 : return c;
  }

  final(a,b,c);
  return c;
}

// function object class for Hashing with lookup3
struct Hash64 {
  std::size_t operator() (const long64& key) const {
    return lookup3((const u8*) &key, sizeof(key), 13 );
  }
};

// function object class for Hashing with cell id
struct HashCell {
  std::size_t operator() (const long64& cellID) const {
    return std::hash<int>()(cellID >> 16) ^ std::hash<int>()(cellID >> 32) ^ std::hash<int>()(cellID >> 48);
  }
};

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

void generate_random_points(const Domain& domain, const int& numPoints, NodeArray& nodes);
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

  // inputs
  double xmin = 0.0, ymin = 0.0, zmin = 0.0;
  double xmax = 3.0, ymax = 2.0, zmax = 1.0;
  double horizon = 0.1;
  double num_points = 10000;

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
  std::cout << cell_node_map.bucket_count() << std::endl;
  //cell_node_map.rehash((size_t)domain.totalCells());
  //std::cout << cell_node_map.bucket_count() << std::endl;

  // Read the data for a particular cell, say, [2, 1, 1];
  long64 cell211 = ((long64)2 << 16) | ((long64)1 << 32) | ((long64)1 << 48);
  CellNodePairIterator nodes = cell_node_map.equal_range(cell211);
  for (auto it = nodes.first; it != nodes.second; ++it) {
    std::cout << "Cell (2,1,1): key = " << it->first << " value = " << *(it->second) << std::endl;
  }

  // Print out all the data
  //for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
  //  std::cout << "key = " << it->first << " value = " << *(it->second) << std::endl;
  //}

  // Check the buckets
  std::cout << "cell_node_map's buckets contain:\n";
  for ( unsigned i = 0; i < cell_node_map.bucket_count(); ++i) {
    std::cout << "bucket #" << i << " contains:";
    for ( auto local_it = cell_node_map.begin(i); local_it!= cell_node_map.end(i); ++local_it )
      std::cout << " " << local_it->first << ":" << *(local_it->second);
    std::cout << std::endl;
  }

  // Delete the nodes
  for (NodeIterator iter = node_list.begin(); iter != node_list.end(); iter++) {
    delete *iter;
  }
  return 0;
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
    node->setPosition(pos);
    nodeList.push_back(node); 
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
    node.setPosition(pos);
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
    node->setPosition(pos);
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
    node->setPosition(pos);
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
    node->setPosition(pos);
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
    node->setPosition(pos);
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
