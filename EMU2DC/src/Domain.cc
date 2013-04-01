#include <Domain.h>
#include <cmath>

using namespace Emu2DC;

Domain::Domain() 
{
  d_lower = {{0.0, 0.0, 0.0}};
  d_upper = {{1.0, 1.0, 1.0}};
  d_numcells = {{1, 1, 1}};
  d_xrange = 1.0;
  d_yrange = 1.0;
  d_zrange = 1.0;
  d_horizon = 1.0;
}

Domain::~Domain() {}

Domain::Domain(const Array3& lower, const Array3& upper):d_lower(lower), d_upper(upper)
{
  d_xrange = std::abs(d_upper[0] - d_lower[0]);
  d_yrange = std::abs(d_upper[1] - d_lower[1]);
  d_zrange = std::abs(d_upper[2] - d_lower[2]);
  d_numcells = {{1, 1, 1}};
  d_horizon = std::max(std::max(d_xrange, d_yrange), d_zrange);
}

Domain::Domain(const Array3& lower, const Array3& upper, const IntArray3& numcells):d_lower(lower), 
                                                                             d_upper(upper),
                                                                             d_numcells(numcells)
{
  d_xrange = std::abs(d_upper[0] - d_lower[0]);
  d_yrange = std::abs(d_upper[1] - d_lower[1]);
  d_zrange = std::abs(d_upper[2] - d_lower[2]);
  d_horizon = std::max(std::max(d_xrange/(double)d_numcells[0], d_yrange/(double)d_numcells[1]),
                           d_zrange/(double)d_numcells[2]);
}
    
Domain::Domain(const Array3& lower, const Array3& upper, const double& horizon):d_lower(lower), 
                                                                          d_upper(upper),
                                                                          d_horizon(horizon)
{
  d_xrange = std::abs(d_upper[0] - d_lower[0]);
  d_yrange = std::abs(d_upper[1] - d_lower[1]);
  d_zrange = std::abs(d_upper[2] - d_lower[2]);
  d_numcells[0] = (int) (d_xrange/horizon);
  d_numcells[1] = (int) (d_yrange/horizon);
  d_numcells[2] = (int) (d_zrange/horizon);
}

void Domain::findCellIndex(const Array3& point,
                       IntArray3& cell) const
{
  cell[0] = 1 + (int)((point[0] - d_lower[0])/d_horizon);
  cell[1] = 1 + (int)((point[1] - d_lower[1])/d_horizon);
  cell[2] = 1 + (int)((point[2] - d_lower[2])/d_horizon);
}

bool Domain::inside(const Array3& point) const
{
  return point[0] > d_lower[0] && point[1] > d_lower[1] && point[2] > d_lower[2]
             && point[0] < d_upper[0] && point[1] < d_upper[1] && point[2] < d_upper[2];
}

