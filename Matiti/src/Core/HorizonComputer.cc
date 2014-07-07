#include <Core/HorizonComputer.h>
#include <Core/Body.h>
#include <Core/Node.h>
#include <Core/Element.h>
#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>

using namespace Matiti;

void 
HorizonComputer::operator()(const BodySP body,
                            SimulationState& state)
{
  NodePArray nodes = body->nodes();
  ElementPArray elements = body->elements();

    double COEF=std::sqrt(3);

  // Loop thru all nodes and find maximum 
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {

    NodeP cur_node = *node_iter;
    ElementPArray adj_elems = cur_node->getAdjacentElements();

    double max_length = 0.0;

    // Loop thru elements attached to this node
    for (auto elem_iter = adj_elems.begin(); elem_iter != adj_elems.end(); ++elem_iter) {
        
      // Loop thru nodes of the current element
      ElementP cur_elem = *elem_iter;
      NodePArray cur_elem_nodes = cur_elem->nodes();
      for (auto elem_node_iter = cur_elem_nodes.begin(); 
                elem_node_iter != cur_elem_nodes.end(); ++elem_node_iter) {
          
        NodeP cur_elem_node = *elem_node_iter;
        double dist = cur_node->distance(*cur_elem_node);
        max_length = std::max(max_length, dist);
      }
    }
    cur_node->horizonSize(state.horizonFactor()*max_length/COEF);
  }
}



