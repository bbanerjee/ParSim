// Test the map, mutimap, unorderd_map, unordered_multimap containers
#include <Node.h>
#include <map>
#include <tr1/unordered_map>
#include <vector>

using namespace Emu2DC;


void test_map_node();
void test_map_node_pointer();
void test_multimap_node_pointer();
void test_key_generation();
bool bit(int ii, int jj);

int main()
{
  //test_map_node();
  //test_map_node_pointer();
  //test_multimap_node_pointer();
  test_key_generation();
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
  // Divide the domain into m x n x p cells
  int mm = 100;
  int nn = 100;
  int pp = 10;
  int nbits = 31;
  for (int ii = 0; ii < mm; ii++) {
    for (int jj = 0; jj < nn; jj++) {
      for (int kk = 0; kk < pp; kk++) {
        std::cout << "Cell = ["<< ii << "," << jj << "," << kk << "];" ;
        long placebit = 1<<nbits;
        //std::cout << "placebit = " << placebit << std::endl;
        long bitsum = 0;
        for (int bb=0; bb < nbits; bb++) {
          int i_bit = (int) (ii & (1 << bb));
          int j_bit = (int) (jj & (1 << bb));
          int k_bit = (int) (kk & (1 << bb));
          //std::cout << "bit(ix,j) = " << i_bit << " bit(iy,j) = " << j_bit << " bit(iz,j) = " << k_bit 
          //          << std::endl;
          bitsum += (1<<(3*bb))*((i_bit<<2)+(j_bit<<1)+k_bit);
        } 
        //std::cout << "bitsum = " << bitsum << std::endl;
        std::cout << "key = " << placebit+bitsum << std::endl;
      }
    }
  }
}

bool bit(int ii, int jj)
{
  return ii & (1 << jj);
}

/*
void test_unorderedmap_node_pointer()
{
  typedef std::tr1::unordered_map<Node*, Node*> Bond;
  typedef std::tr1::unordered_map<Node*, Node*>::iterator BondIterator;
  typedef std::tr1::unordered_map<Node*, Node*>::const_iterator constBondIterator;

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
*/
