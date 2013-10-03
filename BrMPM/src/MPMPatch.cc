#include <MPMPatch.h>
#include <Node.h>
#include <NodePArray.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

#include <VelocityBC.h>
#include <Body.h>
#include <cmath>

using namespace Matiti;

MPMPatch::MPMPatch() 
  : d_lower(0.0, 0.0), d_upper(1.0, 1.0), d_xrange(1.0), 
    d_yrange(1.0), d_horizon(1.0)
{
  d_num_cells = {{1, 1}};
}

MPMPatch::~MPMPatch() {}

MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper):d_lower(lower), d_upper(upper)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());  
  d_num_cells = {{1, 1}};
  d_horizon = std::max(d_xrange, d_yrange);
}

MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper, const IntArray2& numCells):d_lower(lower), 
                                                                             d_upper(upper),
                                                                             d_num_cells(numCells)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_horizon =std::max(d_xrange/(double)d_num_cells[0], d_yrange/(double)d_num_cells[1]));                           
}
    
MPMPatch::MPMPatch(const Point3D& lower, const Point3D& upper, const double& horizon):d_lower(lower), 
                                                                          d_upper(upper),
                                                                          d_horizon(horizon)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_num_cells[0] = (int) (d_xrange/horizon);
  d_num_cells[1] = (int) (d_yrange/horizon);
}

void 
MPMPatch::initialize(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP dom_ps=ps->findBlock("Domain");
  if (!dom_ps) return;

  Uintah::Vector lower(0.0, 0.0);
  Uintah::Vector upper(0.0, 0.0);
  Uintah::IntVector num_cells(1,1);
//  Uintah::IntVector num_grids(1,1);
  dom_ps->require("min", lower);
  dom_ps->require("max", upper);
  dom_ps->require("num_cells", num_cells);
  dom_ps->require("num_ghost", d_ghost);
  dom_ps->require("domain_thick", d_thick);
  dom_ps->require("particlesPerElement", d_particlesperelement); 
  
  d_lower.x(lower[0]);
  d_lower.y(lower[1]);
  d_upper.x(upper[0]);
  d_upper.y(upper[1]);
  d_num_cells[0] = num_cells[0];
  d_num_cells[1] = num_cells[1];
  d_num_grids[0] = 1+num_cells[0]+2*d_ghost;
  d_num_grids[1] = 1+num_cells[1]+2*d_ghost;

  d_xrange=std::abs(d_upper.x() - d_lower.x());
  d_yrange=std::abs(d_upper.y() - d_lower.y());





