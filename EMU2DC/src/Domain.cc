#include <Domain.h>
#include <VelocityBC.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>

using namespace Emu2DC;

Domain::Domain() 
  : d_lower(0.0, 0.0, 0.0), d_upper(1.0, 1.0, 1.0), d_xrange(1.0), 
    d_yrange(1.0), d_zrange(1.0), d_horizon(1.0)
{
  d_num_cells = {{1, 1, 1}};
}

Domain::~Domain() {}

Domain::Domain(const Point3D& lower, const Point3D& upper):d_lower(lower), d_upper(upper)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_num_cells = {{1, 1, 1}};
  d_horizon = std::max(std::max(d_xrange, d_yrange), d_zrange);
}

Domain::Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells):d_lower(lower), 
                                                                             d_upper(upper),
                                                                             d_num_cells(numCells)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_horizon = std::max(std::max(d_xrange/(double)d_num_cells[0], d_yrange/(double)d_num_cells[1]),
                           d_zrange/(double)d_num_cells[2]);
}
    
Domain::Domain(const Point3D& lower, const Point3D& upper, const double& horizon):d_lower(lower), 
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
Domain::initialize(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP dom_ps = ps->findBlock("Domain");
  if (!dom_ps) return;

  Uintah::Vector lower(0.0, 0.0, 0.0);
  Uintah::Vector upper(0.0, 0.0, 0.0);
  Uintah::IntVector num_cells(1,1,1);
  dom_ps->require("min", lower);
  dom_ps->require("max", upper);
  dom_ps->get("num_cells", num_cells);

  d_lower.x(lower[0]);
  d_lower.y(lower[1]);
  d_lower.z(lower[2]);
  d_upper.x(upper[0]);
  d_upper.y(upper[1]);
  d_upper.z(upper[2]);
  d_num_cells[0] = num_cells(0);
  d_num_cells[1] = num_cells(1);
  d_num_cells[2] = num_cells(2);

  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_horizon = std::min(std::max(d_xrange/(double)d_num_cells[0], d_yrange/(double)d_num_cells[1]),
                           d_zrange/(double)d_num_cells[2]);

  // Read the velocity boundary conditions for this domain
  Uintah::ProblemSpecP bc_ps = dom_ps->findBlock("BoundaryConditions");
  for (Uintah::ProblemSpecP vel_ps = bc_ps->findBlock("VelocityBC"); vel_ps != 0;
       vel_ps = vel_ps->findNextBlock("VelocityBC")) {
    VelocityBCSP vel_BC = std::make_shared<VelocityBC>();
    vel_BC->initialize(vel_ps); 
    d_vel_BC.emplace_back(vel_BC);
  }
}

const Point3D& Domain::lower() const
{
  return d_lower;
}

const Point3D& Domain::upper() const
{
  return d_upper;
}

const double& Domain::horizon() const
{
  return d_horizon;
}

const double& Domain::xrange() const
{
  return d_xrange;
}

const double& Domain::yrange() const
{
  return d_yrange;
}

const double& Domain::zrange() const
{
  return d_zrange;
}

const IntArray3& Domain::numCells() const
{
  return d_num_cells;
}

const double Domain::totalCells() const
{
  return (d_num_cells[0]*d_num_cells[1]*d_num_cells[2]);
}

void Domain::findCellIndex(const Point3D& point,
                           IntArray3& cell) const
{
  cell[0] = 1 + (int)((point.x() - d_lower.x())/d_horizon);
  cell[1] = 1 + (int)((point.y() - d_lower.y())/d_horizon);
  cell[2] = 1 + (int)((point.z() - d_lower.z())/d_horizon);
}

void Domain::findCellIndex(const long64& cell_key,
                           IntArray3& cell) const
{
  int kk_up = (cell_key >> 48);
  int jj_up = (cell_key >> 32) & 0xffff;
  int ii_up = (cell_key >> 16) & 0xffff;
  cell[0] = ii_up;
  cell[1] = jj_up;
  cell[2] = kk_up;
}

bool Domain::inside(const Point3D& point) const
{
  return point.x() > d_lower.x() && point.y() > d_lower.y() && point.z() > d_lower.z()
             && point.x() < d_upper.x() && point.y() < d_upper.y() && point.z() < d_upper.z();
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Domain& domain)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Computational domain:" << std::endl;
    out << "  Lower = [" << domain.d_lower.x() << ", " << domain.d_lower.y() << ", " << domain.d_lower.z() << "]";
    out << "  Upper = [" << domain.d_upper.x() << ", " << domain.d_upper.y() << ", " << domain.d_upper.z() << "]" << std::endl;
    out << "  Horizon size = " << domain.d_horizon << std::endl;
    out << "  Range = [" << domain.d_xrange << ", " << domain.d_yrange << ", " << domain.d_zrange << "]" << std::endl;
    out << "  Cells = [" << domain.d_num_cells[0] << ", " << domain.d_num_cells[1] << ", " << domain.d_num_cells[2] << "]" << std::endl;
    return out;
  }
}
