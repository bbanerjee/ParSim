#include <DispBC.h>
#include <Node.h>

using namespace Emu2DC;

DispBC::DispBC() {};
DispBC::~DispBC() {};

// Apply displacement or velocity boundary conditions.
// (for now only when a node touches the domain boundary)
void
DispBC::apply_boundary_cond(const Domain& domain,
                            NodeArray& nodes)
{
  Array3 pos_old = {{0.0, 0.0, 0.0}};
  Array3 pos_new = {{0.0, 0.0, 0.0}};
  Array3 disp_old = {{0.0, 0.0, 0.0}};
  Array3 disp_new = {{0.0, 0.0, 0.0}};
  Array3 vel_old = {{0.0, 0.0, 0.0}};
  Array3 vel_new = {{0.0, 0.0, 0.0}};

  for (NodeArrayIter node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    Node* cur_node = *node_iter;

    // Get the data
    cur_node->getPosition(pos_old);
    cur_node->getDisplacement(disp_old);
    cur_node->getNewDisplacement(disp_new);
    cur_node->getVelocity(vel_old);
    cur_node->getNewVelocity(vel_new);

    Array3 pos_new = {{0.0, 0.0, 0.0}};
    pos_new[0] = pos_old[0] + disp_old[0];
    pos_new[1] = pos_old[1] + disp_old[1];
    pos_new[2] = pos_old[2] + disp_old[2];

    // Find if the node is inside domain (then do nothing)
    if (domain.inside(pos_new)) continue;

    // Otherwise apply reflective boundary condition
    // 1) Push back the node to the boundary along the boundary normal
    //    i.e. update displacement
    // 2) Reverse the normal velocity component 
    

  }

  // This should be a separate task
  /*
    ! Replace the values of displacement, velocity and strain energy with new ones
    do i=1,nnodes
      nodes(i)%disp(1)=nodes(i)%new_disp(1)
      nodes(i)%disp(2)=nodes(i)%new_disp(2)
      nodes(i)%veloc(1)=nodes(i)%new_veloc(1)
      nodes(i)%veloc(2)=nodes(i)%new_veloc(2)
      nodes(i)%strain_energy=wt(i)
    enddo
  */
    
}

