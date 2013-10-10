// Test the map, mutimap, unorderd_map, unordered_multimap containers
#include <Core/Node.h>
#include <map>
#include <tr1/unordered_map>
#include <vector>
#include <cstdint>
#include <ctime>
#include <random>
#include <memory>
#include <chrono>

using namespace Matiti;


void test_map_node();
void test_map_node_pointer(const int numNodes);
void test_map_node_handle(const int numNodes);
void test_multimap_node_pointer();
void test_key_generation();
bool bit(int ii, int jj);
void test_random_number_generator();
void test_random_point_generator();
void test_unorderedmap_node_pointer();
void test_unorderedmultimap_node_pointer();
void test_unorderedmultimap_with_hash();
void test_pack_unpack();

// CrapWow64 hash function
typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef uint64_t u128;
#define NBITS 32
inline u64 CrapWow64(const u8 *key, u64 len, u64 seed ) {
  const u64 m = 0x95b47aa3355ba1a1, n = 0x8a970be7488fda55;
  #define cwfold( a, b, lo, hi ) { p = (u64)(a) * (u128)(b); lo ^= (u64)p; hi ^= (u64)(p >> NBITS); }
  #define cwmixa( in ) { cwfold( in, m, k, h ); }
  #define cwmixb( in ) { cwfold( in, n, h, k ); }

  const u64 *key8 = (const u64 *)key;
  u64 h = len, k = len + seed + n;
  u128 p;

  while ( len >= 16 ) { cwmixb(key8[0]) cwmixa(key8[1]) key8 += 2; len -= 16; }
  if ( len >= 8 ) { cwmixb(key8[0]) key8 += 1; len -= 8; }
  if ( len ) { cwmixa( key8[0] & ( ( (u64)1 << ( len * 8 ) ) - 1 ) ) }
  cwmixb( h ^ (k + n) )
  return k ^ h;
}

// lookup3 hash function
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
    //return std::hash<int>()((cellID >> 16) & 0xffff) ^ std::hash<int>()((cellID >> 32) & 0xffff) ^ std::hash<int>()(cellID >> 48);
    return std::hash<int>()(cellID >> 16) ^ std::hash<int>()(cellID >> 32) ^ std::hash<int>()(cellID >> 48);
  }
};

int main()
{
  //test_map_node();
  int num_nodes = 1000000;
  auto t1 = std::chrono::high_resolution_clock::now();
  test_map_node_pointer(num_nodes);
  auto t2 = std::chrono::high_resolution_clock::now();
  test_map_node_handle(num_nodes);
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "Node*: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << "NodeP: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() << std::endl;
  //test_multimap_node_pointer();
  //test_key_generation();
  //test_random_number_generator();
  //test_random_point_generator();
  //test_unorderedmap_node_pointer();
  //test_unorderedmultimap_node_pointer();
  //test_unorderedmultimap_with_hash();
  //test_pack_unpack();
  return 0;
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
    Point3D pos((double)ii, 0.0, 0.0);
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

void test_map_node_pointer(const int numNodes)
{
  typedef std::map<Node*, Node*> Bond;
  typedef std::map<Node*, Node*>::iterator BondIterator;
  typedef std::map<Node*, Node*>::const_iterator constBondIterator;

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create nodes
  NodeArray nodelist;
  for (int ii = 0; ii < numNodes; ii++) {
    Node* node = new Node();
    Point3D pos((double)ii, 0.0, 0.0);
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // Create two bonds
  Bond bondList;
  for (int ii = 0; ii < numNodes-1; ii++) {
    bondList.insert(std::pair<Node*,Node*>(nodelist[ii], nodelist[ii+1]));
  }

  for (constBondIterator iter = bondList.begin(); iter != bondList.end(); iter++) {
    // std::cout << "Master node = " << *(iter->first) << " Slave Node = " << *(iter->second) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_map_node_handle(const int numNodes)
{
  typedef std::shared_ptr<Node> NodeP;
  typedef std::map<NodeP, NodeP> Bond;
  typedef std::map<NodeP, NodeP>::iterator BondIterator;
  typedef std::map<NodeP, NodeP>::const_iterator constBondIterator;

  typedef std::vector<NodeP> NodeArray;
  typedef std::vector<NodeP>::iterator NodeIterator;
  typedef std::vector<NodeP>::const_iterator constNodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  for (int ii = 0; ii < numNodes; ii++) {
    //NodeP node(new Node());
    NodeP node = std::make_shared<Node>();
    Point3D pos((double)ii, 0.0, 0.0);
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
    //nodelist.emplace_back(node); 
  }

  //int count = 0;
  //for (constNodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
  //  std::cout << "Node (" << (++count) << ") = " << *(*iter) << std::endl;
  //}

  //for (int ii = 0; ii < count; ++ii) {
  //  std::cout << "Node (" << ii+1 << ") = " << nodelist[ii] << std::endl;
  //}

  // Create two bonds
  Bond bondList;
  for (int ii = 0; ii < numNodes-1; ii++) {
    bondList.insert(std::pair<NodeP,NodeP>(nodelist[ii], nodelist[ii+1]));
  }

  for (constBondIterator iter = bondList.begin(); iter != bondList.end(); iter++) {
   //  std::cout << "Master node = " << *(iter->first) << " Slave Node = " << *(iter->second) << std::endl;
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
    Point3D pos((double)ii, 0.0, 0.0);
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
    Point3D pos((double)ii, 0.0, 0.0);
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
  typedef int64_t long64;
  typedef std::tr1::unordered_multimap<long64, Node*> CellNodeMap;
  typedef CellNodeMap::iterator CellNodeMapIterator;
  typedef CellNodeMap::const_iterator constCellNodeMapIterator;
  typedef std::pair<long64, Node*> CellNodePair; 
  typedef std::pair<CellNodeMapIterator, CellNodeMapIterator> CellNodePairIterator; 
  typedef std::pair<constCellNodeMapIterator, constCellNodeMapIterator> constCellNodePairIterator; 

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  int num_nodes = 10;
  for (int ii = 0; ii < num_nodes; ii++) {
    Node* node = new Node();
    Point3D pos((double)ii, 0.0, 0.0);
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

  // Check the hash function
  CellNodeMap::hasher hash_func = cell_node_map.hash_function();
  for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "key = " << it->first << " hashes as = " << hash_func(it->first) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_unorderedmultimap_with_hash()
{
  typedef int64_t long64;
  //typedef std::tr1::unordered_multimap<long64, Node*, Hash64> CellNodeMap;
  typedef std::tr1::unordered_multimap<long64, Node*, HashCell> CellNodeMap;
  typedef CellNodeMap::iterator CellNodeMapIterator;
  typedef CellNodeMap::const_iterator constCellNodeMapIterator;
  typedef std::pair<long64, Node*> CellNodePair; 
  typedef std::pair<CellNodeMapIterator, CellNodeMapIterator> CellNodePairIterator; 
  typedef std::pair<constCellNodeMapIterator, constCellNodeMapIterator> constCellNodePairIterator; 

  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  // Create ten nodes
  NodeArray nodelist;
  int num_nodes = 10000;
  for (int ii = 0; ii < num_nodes; ii++) {
    Node* node = new Node();
    Point3D pos((double)ii, 0.0, 0.0);
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

  // Erase a bit of data and Print out all the data
  for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "Before: key = " << it->first << " value = " << *(it->second) << std::endl;
  }
  cell_node_map.erase(65536);
  for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "After: key = " << it->first << " value = " << *(it->second) << std::endl;
  }

  // Check the hash function
  CellNodeMap::hasher hash_func = cell_node_map.hash_function();
  for (auto it = cell_node_map.begin(); it != cell_node_map.end(); ++it) {
    std::cout << "key = " << it->first << " hashes as = " << hash_func(it->first) << std::endl;
  }

  // Delete the stuff
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

void test_pack_unpack()
{
  // pack
  //int ii = 15, jj = 22, kk = 27;
  int ii = 2, jj = 1, kk = 0;
  long64 cell = ((long64)ii << 16) | ((long64)jj << 32) | ((long64)kk << 48);
  std::cout << "[" << ii << ", " << jj << ", " << kk << "] packed into " << cell << std::endl;

  // unpack one by one
  int kk_up = (cell >> 48);
  std::cout << " kk_up = " << kk_up << std::endl;
  int jj_up = (cell >> 32) & 0xffff ;
  std::cout << " jj_up = " << jj_up << std::endl;
  int ii_up = (cell >> 16) & 0xffff;
  std::cout << " ii_up = " << ii_up << std::endl;
 
  // partial unpack
  kk_up = (cell >> 48);
  std::cout << " kk_up = " << kk_up << std::endl;
  jj_up = (cell >> 32);
  std::cout << " jj_up = " << jj_up << std::endl;
  ii_up = (cell >> 16);
  std::cout << " ii_up = " << ii_up << std::endl;
}
