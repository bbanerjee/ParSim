#include <Horizon.h>
#include <Node.h>

using namespace Emu2DC;

Horizon::Horizon() 
{
}

Horizon::~Horizon() 
{
}

// ****************************************************************************************************
// This subroutine will calculate the horizon)size for each node of any uniform and non-uniform grid. 
// This value is the largest edge connceted to this node times an scalar defined by the user
//*****************************************************************************************************
void 
Horizon::calculateHorizon(const GlobalState* state,
                          const DataWarehouse* old_dw)
{
  NodeArray nodes;
  old_dw->getNodeList(nodes);
  
  ElemArray elements;
  old_dw->getElementList(elements);

  // Loop thru all nodes
  double radnod_max = -1.0e-6;
  for (NodeArrayIter node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {

    double max_length = 0.0;
    double max_interval_x = 0.0;
    double max_interval_y = 0.0;
    int count = 0;

    Node* cur_node = *node_iter;
    int cur_node_id = cur_node->getID();

    cur_node->setHorizonSize(0.0);

    // Loop thru elements 
    for (ElemArrayIter elem_iter = elements.begin(); elem_iter != elements.end(); ++elem_iter) {
        
      Element* cur_elem = *elem_iter;

      // Compute area, xlength, ylength of the element
      double area = 0.0;
      double xlength = 0.0;
      double ylength = 0.0;
      cur_elem->computeGeometry(area, xlength, ylength);

      // Loop thru nodes of each element
      NodeArrayIter elem_node_iter = cur_elem->d_elementnodes.begin();
      for (; elem_node_iter != cur_elem->d_elementnodes.end(); ++elem_node_iter) {
          
        Node* cur_elem_node = *elem_node_iter;
        int cur_elem_node_id = cur_elem_node->getId();

        if (cur_elem_node_id == cur_node_id) {

          cur_node->addAdjacentElement(cur_elem->d_id);

          // If this node is a hanging node
          NodeArrayIter inner_elem_node_iter = cur_elem->d_elementnodes.begin();
          for (; inner_elem_node_iter != cur_elem->d_elementnodes.end(); ++inner_elem_node_iter) {
            if ((*inner_elem_node_iter)->isHangingNode()) {
              cur_node->updateAdjacentElementNumHangingNodes(count);
            }
          }

          cur_node->updateAdjacentElementDepth(cur_elem->d_depth);
          cur_node->updateAdjacentElementSize(area);
           
	  max_interval_x = std::max(max_interval_x, xlength);            
          max_interval_y = std::max(max_interval_y, ylength);            

          count++;
	}
      }
    }

    Array3 interval = {max_interval_x, max_interval_y, 0.0};
    (state->interval).push_back(interval);

    max_length = std:max(max_interval_x, max_interval_y);

    cur_node->setHorizonSize(state->horizon_factor*max_length);
    double radnod = 0.5*max_length;
    radnod_max =  (radnod > radnod_max) ? radnod : radnod_max;
    (state->radnod).push_back(radnod);
    state->radnod_max = radnod_max;
  }
}



