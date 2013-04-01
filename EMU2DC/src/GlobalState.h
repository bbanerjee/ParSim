#ifndef EMU2DC_GLOBALSTATE_H
#define EMU2DC_GLOBALSTATE_H

#include <Node.h>
#include <Element.h>
#include <Line.h>
#include <string>
#include <vector>
#include <array>

namespace Emu2DC {

  class GlobalState {

  public:
 
     GlobalState();
     ~GlobalState();

     //  Output file folder and name --
     std::string output_folder_name;
     std::string output_file_name;

     //_____________________ Refinement Global variables _____________________
     std::vector<int> temp;

     int nbc, nforce;                    // number of nodes for prescribing boundary conditions and tractions  
    
     int nroot_elements;
     Element aux_elements[8];

     bool grid_ok = true;  // In the case of two neighbor elements with difference in the level of refinement 
                           // larger than 2, this variable recevies .false. 
  
     int global_depth;
     int max_depth;        // maximum depth of an element on the grid
     int min_depth;        // minimum depth of an element on the grid
     int max_size;         // maximum size of an element on the grid
     int min_size;         // minimum size of an element on the grid
     int dim;              // dimension of the problem: if dim=2, then is a 2D problem
     int iteration=0;
     int iterations;
     int mnode;
     double x_min_lim, y_min_lim, z_min_lim, x_max_lim, y_max_lim, z_max_lim,cont;

     int nbox_x,nbox_y,nboxes;
     double max_strain_energy;
     double damage_index_cri;
  
     bool modified_mesh= false;  // This variable tells us if the mesh have been modified( refined or coarsed) 

     //_____________________ Solver Global variables _____________________
     double critical_stretch;
     double fracture_energy;
     double dc1, dc2, dc3;
     int myn, mbd, nbd1, nbd2;
     double delta, dt;
     bool brokij, brokij_old, broke_first;

     int nbroke, ncyc;

     // Material properties
     double young;
     double  denst;

     // Octree parameters (make vector of vector)
     static const int max_global_nodes=2000000;   // Maximum number of nodes 
     static const int maxfam=2000; 
     static const int ncell1 =10;
     static const int ncell2 =10;
     static const int ncell_def1=10;
     static const int ncell_def2=10;

     int cell_location_in_bin[ncell1][ncell2], number_of_nodes_in_cell[ncell1][ncell2];
     int cell_location_in_bin_def[ncell_def1][ncell_def2], number_of_nodes_in_cell_def[ncell_def1][ncell_def2];

     int mfam, m_def_fam;

     double rad_search, radnod_max;

     // parameters
     double x1_lo, x1_hi;
     double x2_lo, x2_hi;
     double y1_lo, y1_hi;
     double y2_lo, y2_hi;
     std::vector<double> xt1, xt2, yt1, yt2;
 
     std::vector<double> a,b;
     std::vector<int> jsnr,irnr;
     int mctr;
     int  n_elements;
     std::vector<int> info_element;

     // Number of deleted regions
     static const int ndel=1;

     int mdtype[ndel];
     double x1dlo[ndel], x1dhi[ndel], x2dlo[ndel], x2dhi[ndel];
     double  x1dcen[ndel], x2dcen[ndel], raddl[ndel];

     std::vector<Line> Crack_Line;
     int ncracks=0;
     bool newiteration;

     // Time steps
     int iter;
     int nt;

     int snapshots_frequence;

     // external forces
     std::vector<double> xtempbottom, ytempbottom, xtemptop, ytemptop;
     std::vector<int> ncorrespond_bottom, ncorrespond_top;

     std::vector<double> ext1, ext2, RadHorizon;
     std::vector<std::array<double, 3> > interval;
     double shortrange_dist_fac_nom[10],  shortrange_dist_fac_init[10],  shortrange_force_fac[10]; 

     double speed_iter=0;
     std::vector<double> top_crack_speed, bottom_crack_speed; 
     double bc_velocity;
     double force_mag;     // external force magnitude
     double visco, visdk;  // damage viscosity parameters
     double horizon_factor;

     // Decision variables
     int do_adap_refinement;
     int introduce_precrack;
     int dynamic_or_static;
     int sw_micro;

     // Initial Crack
     double initial_crack_x0,initial_crack_y0,initial_crack_x1,initial_crack_y1;
   
} // end namespace
#endif
