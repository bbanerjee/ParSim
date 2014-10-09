/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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
  //UndirectedGraph graph(numNodes);
  UndirectedGraph graph;

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

  // Remove bonds (edges of the graph)
  for (int ii = 0; ii < numNodes-1; ii++) {
    if (ii == 4) {
      node1 = boost::vertex(ii, graph);
      node2 = boost::vertex(ii+1, graph);
      boost::remove_edge(node1, node2, graph);
    }
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
  }

  // Delete the nodes
  for (NodeIterator iter = nodelist.begin(); iter != nodelist.end(); iter++) {
    delete *iter;
  }
}

