#include <MPMPatch.h>
#include <Node.h>
#include <NodePArray.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

#include <VelocityBC.h>
#include <Body.h>
#include <cmath>

using namespace BrMPM;

MPMPatch::MPMPatch() 
  : d_lower(0.0, 0.0, 0.0), d_upper(1.0, 1.0, 1.0),
    d_node_counts(1.0, 1.0, 1.0),  d_cellsize(0.0, 0.0, 0.0),
    d_nof_ghost(0), d_tol(1e-4),
    d_nof_particles_per_cell(0), xcoord(0.0), ycoord(0), zcoord(0),
    d_num_grids(0.0, 0.0, 0.0), d_grids(0.0, 0.0, 0.0)
{
  //d_num_cells = {{1, 1}};
}

MPMPatch::~MPMPatch() {}

MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper):d_lower(lower), d_upper(upper)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());   
  d_num_cells = {{1, 1, 1}};
 // d_horizon = std::max(d_xrange, d_yrange);
}

MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper, const IntArray3& numCells):d_lower(lower),
                                                                             d_upper(upper),
                                                                             d_num_cells(numCells)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
 // d_horizon =std::max(d_xrange/(double)d_num_cells[0], d_yrange/(double)d_num_cells[1]));                           
}
    
MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper, const double& horizon):d_lower(lower), 
                                                                          d_upper(upper),
                                                                          d_horizon(horizon)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_num_cells[0] = (int) (d_xrange/horizon);
  d_num_cells[1] = (int) (d_yrange/horizon);
  d_num_cells[2] = (int) (d_zrange/horizon);
}

void 
MPMPatch::initialize(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP dom_ps = ps->findBlock("Domain");
  if (!dom_ps) return;

  Uintah::Vector lower(0.0, 0.0, 0.0);
  Uintah::Vector upper(0.0, 0.0, 0.0);
  Uintah::IntVector num_cells(1, 1, 1);
//  Uintah::IntVector num_grids(1, 1, 1);
  dom_ps->require("min", lower);
  dom_ps->require("max", upper);
  dom_ps->require("num_cells", num_cells);
  dom_ps->require("num_ghost", d_nof_ghost);
 // dom_ps->require("domain_thick", d_thick);
  dom_ps->require("particlesPerElement", d_nof_particles_per_cell);
  
  d_lower.x(lower[0]);
  d_lower.y(lower[1]);
  d_lower.z(lower[2]);
  d_upper.x(upper[0]);
  d_upper.y(upper[1]);
  d_upper.z(upper[2]);

  d_num_cells[0] = num_cells[0];
  d_num_cells[1] = num_cells[1];
  d_num_cells[2] = num_cells[2];
  d_num_grids[0] = 1 + num_cells[0] + 2*d_nof_ghost;
  d_num_grids[1] = 1 + num_cells[1] + 2*d_nof_ghost;
  d_num_grids[2] = 1 + num_cells[2] + 2*d_nof_ghost;

  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());

  double dx = d_xrange/(d_num_grids[0] + 1);
  double dy = d_yrange/(d_num_grids[1] + 1);
  double dz = d_zrange/(d_num_grids[2] + 1);

  for (int x_id = 0; x_id < d_num_grids[0] + 2; x_id++) {
       xcoord = d_lower.x() + x_id * dx;
       for (int y_id = 0; y_id < d_num_grids[1] + 2; y_id++) {
            ycoord = d_lower.y() + y_id * dy;
            for (int z_id = 0; z_id < d_num_grids[2] + 2; z_id++) {
                 zcoord = d_lower.z() + z_id * dz;
                 d_grids.x(xcoord);
                 d_grids.y(ycoord);
                 d_grids.z(zcoord);
                 d_gridsPosition.emplace_back(d_grids);
            }
       }
  }
 
  d_cellsize.x(d_xrange/d_num_grids[0]);
  d_cellsize.y(d_yrange/d_num_grids[1]);
  d_cellsize.z(d_zrange/d_num_grids[2]);

}



bool
MPMPatch::insidePatch(const Point3D& point) const
{
  if ((point.x() < d_lower.x()) || (point.y() < d_lower.y()) || (point.z() < d_lower.z())) {
     return false;
  }
  else if ((point.x() > d_upper.x()) || (point.x() > d_upper.x()) || (point.x() > d_upper.x())) {
     return false;
  }   
  else {
     return true;
  }
}



bool
MPMPatch::allInsidePatch(const std::vector<Point3D> points) const
{
  for (auto iter = points.begin(); iter != points.end(); iter++)
      {
        Point3D point = *iter;
        if (!insidePatch(point)) {
            return false;
            break;
           }
       }
        return true;
      }
}




 
 









 



