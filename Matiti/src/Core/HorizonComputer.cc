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

  // double COEF=std::sqrt(3);

  // Loop thru all nodes and find maximum 
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {

    NodeP cur_node = *node_iter;
    ElementPArray adj_elems = cur_node->getAdjacentElements();

    double max_length = 0.0;
    double max_length_xyz = 0.0;
    double max_length_x = 0.0;
    double max_length_y = 0.0;
    double max_length_z = 0.0;

    // Loop thru elements attached to this node
    for (auto elem_iter = adj_elems.begin(); elem_iter != adj_elems.end(); ++elem_iter) {
        
      // Loop thru nodes of the current element
      ElementP cur_elem = *elem_iter;
      NodePArray cur_elem_nodes = cur_elem->nodes();
      for (auto elem_node_iter = cur_elem_nodes.begin(); 
                elem_node_iter != cur_elem_nodes.end(); ++elem_node_iter) {
          
        NodeP cur_elem_node = *elem_node_iter;

        double dist_x = std::abs(cur_node->x()-cur_elem_node->x());
        double dist_y = std::abs(cur_node->y()-cur_elem_node->y());
        double dist_z = std::abs(cur_node->z()-cur_elem_node->z());

        double dist = cur_node->distance(*cur_elem_node);
        max_length = std::max(max_length, dist);

        max_length_x = std::max(max_length_x, dist_x);
        max_length_y = std::max(max_length_y, dist_y);
        max_length_z = std::max(max_length_z, dist_z);

        max_length_xyz = std::max(std::max(max_length_x, max_length_y), max_length_z);
      }
    }
   
//   std::cout << std::endl << std::endl << "Max_length= " << max_length << "   Max_length_xyz= " << max_length_xyz 
//             << "  horizon_factor= " << state.horizonFactor() << std::endl << std::endl;
   
    //cur_node->horizonSize(state.horizonFactor()*max_length/COEF);
    
//    cur_node->horizonSize(state.horizonFactor()*max_length);
     
     cur_node->horizonSize(state.horizonFactor()*max_length_xyz);

     Array3 interval = {{0.0, 0.0, 0.0}};
     interval[0] = max_length_x;  
     interval[1] = max_length_y;
     interval[2] = max_length_z;
     cur_node->setInterval(interval);  
//   std::cout << "node's_horizon_size= " << cur_node->horizonSize() << std::endl;
//   std::cout << "node's interval vector= (" << cur_node->getInterval()[0] << ", "
//                                            << cur_node->getInterval()[1] << ", " 
//                                            << cur_node->getInterval()[2] << ")" << std::endl;                           
  }














}



