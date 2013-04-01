#include <Dynamics.h>
#include <Node.h>

Dynamics::Dynamics() {}
Dynamics::~Dynamics() {}

void
Dynamics::do_dynamics(const GlobalFlags* flags,
                      const GlobalState& old_state,
                      GlobalState& new_state,
                      DamageModel& damage_model,
                      BondFamily& bond_family,
                      NodeArray& node_list)
{
  
  std::cerr << "....Begin Solver...." << std::endl;

  // Get the time step and iteration information
  int nt = old_state.numTimesteps;
  int dt = old_state.delT;
  int output_freq = old_state.snapshot_frequency;
  
  // Local variables
  double time1 = 0.0, time2 = 0.0;
  // call cpu_time(time1)
  damageModel.compute_critical_strain();
  bond_family.sort_ref();
  if (flags->introduce_precrack) {
    initial_broken_bond();
  }
  // call cpu_time(time2)
  //std::cerr << "damage, computational time(seconds):" << time2-time1 << std::endl;

  // Output nodal information every snapshots_frequency iteration   
  int output_file_count = 0
  write_output(output_file_count);

  // Set up initial external force
  // call cpu_time(time1)
  ForceBC* extForceBC = new ForceBC();
  for (NodeArrayIter nodeiter = node_list.begin(); nodeiter != node_list.end(); ++nodeiter) {
    Node* cur_node = node_list[*nodeiter];
    Array3 extForce = {0.0, 0.0, 0.0};
    extForceBC->computeExtForceDensity(); 
  }
  // call cpu_time(time2)
  // std::cerr << "ext, computational time(seconds):" << time2-time1 << std::endl;

  // Start time iteration
  VelocityVerlet* integrator = new VelocityVerlet();
  DispBC* dispBC = new DispBC();
  for (int iter = 0; iter < nt; ++iter) {

    // Velocity-Verlet scheme
    // call cpu_time(time1)
    integrator->peri_motion_velocity_verlet2();
    // call cpu_time(time2)
    // std::cerr << "peri_motion, computational time(seconds):" << time2-time1 << std::endl;

    // Apply boundary conditions
    disp_bc->apply_boundary_cond();
  
    // Output nodal information every snapshots_frequency iteration   
    if (std::mod(iter, output_freq) == 0) {
      output_file_count++;
      write_output(output_file_count);
    }
  }
  
  // clean up
  delete extForceBC;
  delete dispBC;
  delete integrator;
}

void
Dynamics::write_output(const int& output_file_count)
{
  // count valid nodes
  int valid_node_count = 0;
  for (NodeArrayIter nodeiter = node_list.begin(); nodeiter != node_list.end(); ++nodeiter) {
    if (node_list[*nodeiter]->omit) continue;  // skip this node
    valid_node_count++;
  }

  // Write the output to individual files
  string output_file_name = old_state.output_file_name;
  string current_output_file_name = output_file_name + toint(output_file_count) + ".tec";
  //write(current_output_file_name, fmt='(A,I5.5,A)') trim(output_file_name), output_file_count, '.tec'

  ofstream output_file(current_output_file_name);
  output_file << "TITLE=\"simulation results\" " << endl;
  output_file << "VARIABLES=\"X\",\"Y\",\"DX\",\"DY\",\"VX\",\"VY\",\"DAM\",\"W\" " << endl;
  output_file << "ZONE I=" << valid_node_count << " SOLUTIONTIME=" << iter*dt << " F=POINT" << endl;
  for (NodeArrayIter nodeiter = node_list.begin(); nodeiter != node_list.end(); ++nodeiter) {
    Node* cur_node = node_list[*nodeiter];
    if (cur_node->omit) continue;  // skip this node
    double xdisp = cur_node->displacement[0];
    double ydisp = cur_node->displacement[1];
    double cur_x_pos = cur_node->position[0] + xdisp;
    double cur_y_pos = cur_node->position[1] + ydisp;
    output_file << cur_x_pos << " " << cur_y_pos << " " << xdisp << " " << ydisp 
                << cur_node->velocity[0] << cur_node->velocity[1]
                << cur_node->damage_index << cur_node->wt << endl;
  }
  output_file.close();
}
