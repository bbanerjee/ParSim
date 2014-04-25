// Test boost adjacency_list graph
#include <Core/Node.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <vector>

using namespace Matiti;

void test_graph_node(const int numNodes);

int main()
{
  int num_nodes = 10;
  test_graph_node(num_nodes);
  return 0;
}

void test_graph_node(const int numNodes)
{
  typedef std::vector<Node*> NodeArray;
  typedef std::vector<Node*>::iterator NodeIterator;

  typedef boost::property<boost::edge_weight_t, double> Weight;
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                         boost::no_property, Weight> UndirectedGraph;

  typename boost::graph_traits<UndirectedGraph>::vertex_descriptor node1, node2;
  typename boost::graph_traits<UndirectedGraph>::out_edge_iterator out, out_end;
  typename boost::graph_traits<UndirectedGraph>::in_edge_iterator in, in_end;

  // Initialize undirected graph
  UndirectedGraph graph(numNodes);

  // Create nodes
  NodeArray nodelist;
  for (int ii = 0; ii < numNodes; ii++) {
    Node* node = new Node();
    Point3D pos((double)ii, 0.0, 0.0);
    node->setID(ii);
    node->position(pos);
    nodelist.push_back(node); 
  }

  // Create bonds (edges of the graph)
  for (int ii = 0; ii < numNodes-1; ii++) {
    node1 = boost::vertex(ii, graph);
    node2 = boost::vertex(ii+1, graph);
    boost::add_edge(node1, node2, graph);
  }

  // Output graph
  for (int ii = 0; ii < numNodes; ii++) {
    std::cout << "out_edges: " << ii << std::endl;
    node1 = boost::vertex(ii, graph);
    for (boost::tie(out, out_end) = out_edges(node1, graph);
         out != out_end; ++out) {
      std::cout << ' ' << *out;
    }
    std::cout << std::endl;
    std::cout << "in_edges: " << ii << std::endl;
    for (boost::tie(in, in_end) = in_edges(node1, graph);
         in != in_end; ++in) {
      std::cout << ' ' << *in;
    }
    std::cout << std::endl;
  }

  // Delete the nodes
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

