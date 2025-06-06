/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Domain.h>
#include <BoundaryConditions/VelocityBC.h>
#include <Core/Body.h>
#include <Containers/NodePArray.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>

using namespace Matiti;

Domain::Domain() 
  : d_lower(0.0, 0.0, 0.0), d_upper(1.0, 1.0, 1.0), d_xrange(1.0), 
    d_yrange(1.0), d_zrange(1.0), d_cell_size(1.0, 1.0, 1.0)
{
  d_num_cells = {{1, 1, 1}};
}

Domain::~Domain() {}

Domain::Domain(const Point3D& lower, const Point3D& upper)
  : d_lower(lower), d_upper(upper)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_num_cells = {{1, 1, 1}};
  d_cell_size.x(d_xrange);
  d_cell_size.y(d_yrange); 
  d_cell_size.z(d_zrange);
}

Domain::Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells)
  : d_lower(lower), d_upper(upper), d_num_cells(numCells)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_cell_size.x(d_xrange/(double)d_num_cells[0]);
  d_cell_size.y(d_yrange/(double)d_num_cells[1]);
  d_cell_size.z(d_zrange/(double)d_num_cells[2]);
}
    
Domain::Domain(const Point3D& lower, const Point3D& upper, const Uintah::Vector& cellSize)
  : d_lower(lower), d_upper(upper), d_cell_size(cellSize)
{
  d_xrange = std::abs(d_upper.x() - d_lower.x());
  d_yrange = std::abs(d_upper.y() - d_lower.y());
  d_zrange = std::abs(d_upper.z() - d_lower.z());
  d_num_cells[0] = (int) (d_xrange/cellSize[0]);
  d_num_cells[1] = (int) (d_yrange/cellSize[1]);
  d_num_cells[2] = (int) (d_zrange/cellSize[2]);
}

void
Domain::clone(const Domain& domain)
{
  d_lower = domain.d_lower;
  d_upper = domain.d_upper;
  d_xrange = domain.d_xrange;
  d_yrange = domain.d_yrange;
  d_zrange = domain.d_zrange;
  d_cell_size = domain.d_cell_size;
  d_num_cells = domain.d_num_cells;

  // *TODO* Add velocity BC
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
  d_cell_size.x(d_xrange/(double)d_num_cells[0]);
  d_cell_size.y(d_yrange/(double)d_num_cells[1]);
  d_cell_size.z(d_zrange/(double)d_num_cells[2]);

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

const Uintah::Vector& Domain::cellSize() const
{
  return d_cell_size;
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
  cell[0] = 1 + (int)((point.x() - d_lower.x())/d_cell_size[0]);
  cell[1] = 1 + (int)((point.y() - d_lower.y())/d_cell_size[1]);
  cell[2] = 1 + (int)((point.z() - d_lower.z())/d_cell_size[2]);
  //std::cout << " point = " << point << " d_lower = " << d_lower << " cell size = " << d_cell_size 
  //          << " cell id = " << cell[0] << "," << cell[1] << "," << cell[2] << std::endl;
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

void 
Domain::applyVelocityBC(BodySP body) const
{
  // Loop through nodes in the body
  const NodePArray& nodes = body->nodes();
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
 
    // Apply bcs on to boundary nodes
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;

    // Get the old and current position of the node
    const Point3D& pos = cur_node->position();
    const Vector3D& disp_new = cur_node->newDisplacement();
    Point3D cur_pos(pos+disp_new);

    // If the node is inside the domain do nothing 
    if (inside(cur_pos)) continue;

    // If the node has not moved, do nothing
    if (disp_new.length() < std::numeric_limits<double>::epsilon()) continue;

    // Find the point of intersection of the domain with a ray through the node along
    // the velocity direction
    Point3D hit_point;
    std::cout << "Node = " << cur_node->getID() << " old_pos = " << pos << " new_pos = " << cur_pos
              << " Lower = " << d_lower << " Upper = " << d_upper << " Disp_new = " << disp_new;
    intersection(pos, disp_new, hit_point);
    std::cout << " Hit point = " << hit_point << std::endl;

    // Apply appropriate velocity boundary conditions
    //std::cout << "Before apply BC: Node = " << cur_node->getID() << " vel = " << cur_node->velocity() << " disp = " << cur_node->displacement() << std::endl;
    for (auto iter = d_vel_BC.begin(); iter != d_vel_BC.end(); ++iter) {
      (*iter)->apply(cur_node, hit_point, d_lower, d_upper);
    }
    //std::cout << "After apply BC:  Node = " << cur_node->getID() << " vel = " << cur_node->newVelocity() << " disp = " << cur_node->newDisplacement() << std::endl;
  }
}

// Find the intersection point of particle position segment (t_n+1-t_n) with domain box
// For algorithm see: http://www.cs.utah.edu/~awilliam/box/
bool
Domain::intersection(const Point3D& point, const Vector3D& ray,
                     Point3D& hitPoint) const
{
  //std::cout << "    Domain intersection with boundary:" << " ray = " << ray
  //          << " upper = " << d_upper << " lower = " << d_lower ;
  Vector3D inv_ray = ray.invDirection();
  Vector3D t1 = (d_lower - point)*inv_ray;
  Vector3D t2 = (d_upper - point)*inv_ray;
  Vector3D tn = Matiti::min(t1, t2);
  Vector3D tf = Matiti::max(t1, t2);
  double tnear = tn.max();
  double tfar = tf.min();
  double tt = (tnear < 0.0 || tnear > 1.0) ? tfar : tnear;
  //std::cout << " tnear = " << tnear << " tfar = " << tfar << std::endl;
  hitPoint = point + ray*tt;
  return !(tt < 0.0 || tt > 1.0);
  //if(tnear <= tfar){
  //  hitPoint = point + ray*tnear;
  //  return true;
  //} else {
  //  return false;
  //}
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Domain& domain)
  {
    //out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Computational domain:" << std::endl;
    out << "  Lower = [" << domain.d_lower.x() << ", " << domain.d_lower.y() << ", " << domain.d_lower.z() << "]";
    out << "  Upper = [" << domain.d_upper.x() << ", " << domain.d_upper.y() << ", " << domain.d_upper.z() << "]" << std::endl;
    out << "  Range = [" << domain.d_xrange << ", " << domain.d_yrange << ", " << domain.d_zrange << "]" << std::endl;
    out << "  Cells = [" << domain.d_num_cells[0] << ", " << domain.d_num_cells[1] << ", " << domain.d_num_cells[2] << "]" << std::endl;
    out << "  Cell size = [ " << domain.d_cell_size[0] << "," << domain.d_cell_size[1] << "," << domain.d_cell_size[2] << "]" << std::endl;
    return out;
  }
}
