#include <InitialConditions.h>
#include <Node.h>

using namespace Emu2DC;
  
InitialConditions::InitialConditions() {};
InitialConditions::~InitialConditions() {};

void 
InitialConditions::apply_initial_conditions(const Domain& domain,
                                            NodeArray& nodes)
{
  Array3 zero = {{0.0, 0.0, 0.0}};

  for (NodeArrayIter node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    Node* cur_node = *node_iter;
    cur_node->setDisplacement(zero);
    cur_node->setVelocity(zero);
    cur_node->setNewDisplacement(zero);
    cur_node->setNewVelocity(zero);
  }
}
