#include <ForceBC.h> 
#define EMU2DC_FORCEBC_H

using namespace Emu2DC;
  
ForceBC::ForceBC() {}
ForceBC::~ForceBC() {}

// Read the input data file to get
// 1) Boundary node sets on which external forces are applied
// 2) external force vector
void 
ForceBC:initialize(std::string input_data_file)
{
  std::cerr << "initialize not implemented yet. " << endl;
  //for (int ii = 0; ii < d_num_ext_force_sets; ++ii) {
  //  input_stream << d_ext_force ;
  //  input_stream << d_ext_force_nodes;    
  //}
}

//********************************************************************
// subroutine ExtForceDenstiy
// Purpose : set the external force density array
//********************************************************************
void 
ForceBC::computeExtForceDensity(NodeArray& nodes)
{
  std::cerr << "ComputeExtForceDensity not implemented yet. " << endl;
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void 
ForceBC::computeExtForceDensity(const NodeArray& nodes,
	  	                const Array3& topLoc,
		                const Array3& botLoc)
{
  // Find the nodes at the top and bottom boundaries
  NodeArray top_nodes;
  NodeArray bot_nodes;
  findBoundaryNodes(nodes, topLoc, botLoc, top_nodes, bot_nodes);

  // Sort the nodes in increasing x-coordinate order and find span
  NodeArray sorted_top_nodes;
  vector<double> nodal_span_top;
  sortNodes(top_nodes, sorted_top_nodes, nodal_span_top);

  NodeArray sorted_bot_nodes;
  vector<double> nodal_span_bot;
  sortNodes(top_nodes, sorted_bot_nodes nodal_span_bot);

  // Compute external force density
  int count = 0;
  Array3 ext_force = {{0.0, 0.0, 0.0}};
  for (NodeArrayIter node_iter=sorted_top_nodes.begin(); 
		  node_iter != sorted_top_nodes.end(); ++node_iter) {
    Node* cur_node = *node_iter;
    double cur_node_vol = cur_node->getVolume();
    ext_force[1] = -ext_force_mag*nodal_span_top[count]/(2.0*cur_node_vol);
    cur_node->setExternalForce(ext_force);
    count++;
  }

  count = 0;
  ext_force = {{0.0, 0.0, 0.0}};
  for (NodeArrayIter node_iter=sorted_bot_nodes.begin(); 
		  node_iter != sorted_bot_nodes.end(); ++node_iter) {
    Node* cur_node = *node_iter;
    double cur_node_vol = cur_node->getVolume();
    ext_force[1] = -ext_force_mag*nodal_span_bot[count]/(2.0*cur_node_vol);
    cur_node->setExternalForce(ext_force);
    count++;
  }
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void
ForceBC::findBoundaryNodes(const NodeArray& nodes,
	  	           const Array3& topLoc,
		           const Array3& botLoc,
		           NodeArray& topNodes,
			   NodeArray& botNodes)
{
  double y_top = topLoc[1];
  double y_bot = botLoc[1];
  Array3 node_pos = {{0.0, 0.0, 0.0}};
  for (NodeArrayIter node_iter=nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    Node* cur_node = *node_iter;
    cur_node->getPosition(node_pos);
    double y_cur = node_pos[1];
    if (std::abs(y_bot-y_pos) <= 1.0e-6) {
      botNodes->push_back(cur_node);
    } else {
      if (std::abs(y_top-y_pos) <= 1.0e-6) {
        topNodes->push_back(cur_node);
      } 
    }
  }
}

// Special treatment for rectangular axis-aligned domain with forces applied at the top
// and bottom boundaries (for testing compatibility with EMUNE)
void 
ForceBC::sortNodes(const NodeArray& boundaryNodes,
		   NodeArray& sortedBoundaryNodes,
		   std::vector<double>& nodalSpan)
{
  vector<double> xCoordsBoundary;
  int count = 0;
  Array3 node_pos = {{0.0, 0.0, 0.0}};
  for (NodeArrayIter node_iter=boundaryNodes.begin(); node_iter != boundaryNodes.end(); 
		  ++node_iter) {
    Node* cur_node = *node_iter;
    cur_node->getPosition(node_pos);
    xCoordsBoundary->push_back(node_pos[0]);
    sortedBoundaryNodes->push_back(cur_node);
    count++;
  }

  for (int ii=0; ii < count; ii++) {
    for (int jj=ii+1; jj < count; jj++) {
      if (xCoordsBoundary[jj] > xCoordsBoundary[ii]) {
        double temp = xCoordsBoundary[ii];
        xCoordsBoundary[ii] = xCoordsBoundary[jj];
        xCoordsBoundary[jj] = temp;
	Node* tempNode = sortedBoundaryNodes[ii];
	sortedBoundaryNodes[ii] = sortedBoundaryNodes[jj];
	sortedBoundaryNodes[jj] = tempNode;
      }
    }
  }

  nodalSpan->push_back(xCoordsBoundary[1] - xCoordsBoundary[0]);
  for (int ii=1; ii < count-1; ii++) {
    nodalSpan->push_back(xCoordsBoundary[ii+1] - xCoordsBoundary[ii-1]);
  }
  nodalSpan->push_back(xCoordsBoundary[count-1] - xCoordsBoundary[count-2]);
}


