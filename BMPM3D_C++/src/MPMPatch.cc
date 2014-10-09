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

#include <MPMPatch.h>
#include <ShapeFunctions/MPMShapeFunctionFactory.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <cmath>

using namespace BrMPM;

MPMPatch::MPMPatch() 
: d_num_particles_per_cell(2, 2, 2),
  d_lower(0.0, 0.0, 0.0),
  d_upper(1.0, 1.0, 1.0),
  d_node_counts(2, 2, 2),
  d_num_ghost(0.0, 0.0, 0.0),  // TODO: should be intvector (?)
  d_cell_size(1.0, 1.0, 1.0),
  d_time(0.0),
  d_time_final(1.0),
  d_delT(0.1),
  d_iteration(0),
  d_tol(1.0e-6)
{
}

MPMPatch::~MPMPatch() {}

void 
MPMPatch::initialize(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP patch_ps = ps->findBlock("Patch");
  if (!patch_ps) return;

  // Get the basic grid data
  Uintah::Vector lower(0.0, 0.0, 0.0);
  Uintah::Vector upper(0.0, 0.0, 0.0);
  Uintah::IntVector num_cells(1, 1, 1);
  Uintah::IntVector num_ghost(1, 1, 1);
  Uintah::IntVector num_ppc(1, 1, 1);
  patch_ps->require("min", lower);
  patch_ps->require("max", upper);
  patch_ps->require("num_cells", num_cells);
  patch_ps->require("num_ghost_cells", num_ghost);
  patch_ps->require("particle_per_grid_cell", num_ppc);

  // Copy the data to local variables (we are using BrMPM::Vector3D instead of Uintah::Vector)
  for (int ii = 0; ii < 3; ++ii) {
    d_lower[ii] = lower[ii];
    d_upper[ii] = upper[ii];
    d_num_particles_per_cell[ii] = num_ppc[ii];
    d_node_counts[ii] = num_cells[ii]+1+2*num_ghost[ii];
    d_num_ghost[ii] = num_ghost[ii];
  }

  // Read the shape function information from the input file and create
  // the right type of shape function
  d_shape = MPMShapeFunctionFactory::create(ps);
}

void
MPMPatch::initGrid(DoubleNodeData& gx)
{
  // TODO:  Hooman to complete this.
}

bool
MPMPatch::insidePatch(const Point3D& point) const
{
  return (point.x() >= d_lower.x() && point.y() >= d_lower.y() && point.z() >= d_lower.z() &&
          point.x() <= d_upper.x() && point.y() <= d_upper.y() && point.z() <= d_upper.z());
}

bool
MPMPatch::allInsidePatch(const Point3DParticleData& points) const
{
  for (auto iter = points.begin(); iter != points.end(); iter++)
  {
    if (!insidePatch(*iter)) return false;
  }
  return true;
}
















