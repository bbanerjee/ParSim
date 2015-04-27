/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#include <Core/RigidBody.h>
#include <Core/Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

//#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Matiti;

RigidBody::RigidBody()
  : d_id(0), d_density(0.0), d_com(0.0, 0.0, 0.0), d_radius(1.0),
    d_init_vel(0.0, 0.0, 0.0), d_init_acc(0.0, 0.0, 0.0), 
    d_init_ang_vel(0.0, 0.0, 0.0), d_init_ang_acc(0.0, 0.0, 0.0),
    d_body_force(0.0, 0.0, 0.0), 
    d_ext_force(0.0, 0.0, 0.0), d_ext_torque(0.0, 0.0, 0.0),
    d_rot_center(0.0, 0.0, 0.0), d_rot_vel(0.0, 0.0, 0.0)
{
  computeInertialProperties();
}

RigidBody::~RigidBody()
{
}

void
RigidBody::computeInertialProperties()
{
  // Compute volume and mass
  d_volume = 4.0/3.0*M_Pi*d_radius*d_radius*d_radius;
  d_mass = d_density*d_volume;

  // Compute the square of the radius of gyration
  rog_sq = 0.6*d_radius*d_radius;
  d_rog_sq = Vector3D(rog_sq, rog_sq, rog_sq);
}

void 
RigidBody::initialize(Uintah::ProblemSpecP& ps)
{
  if (!ps) return;
  
  // Read the density
  ps->require("density", d_density);

  // Read the position and geometry 
  // (**TODO** More general geometries later)
  // Assume sphere
  ps->require("center_of_mass", d_com);
  ps->require("radius", d_radius);

  // Compute inertial properties
  computeInertialProperties();

  // Read the initial conditions 
  // **TODO** Assumes initial spin is zero. Make more general.
  ps->require("initial_velocity", d_init_vel);
  ps->require("initial_acceleration", d_init_acc);
  ps->require("initial_angular_velocity", d_init_ang_vel);
  ps->require("initial_angular_acceleration", d_init_ang_acc);

  // Read the external force boundary conditions 
  ps->require("external_force", d_ext_force);
  ps->require("external_torque", d_ext_torque);

  // Read the initial body forces
  ps->require("body_force", d_body_force);

  // For a rotating coordinate system
  ps->require("center_of_rotation", d_rot_center);
  ps->require("angular_velocity_of_rotation", d_rot_vel);
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const RigidBody& body)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "RigidBody: ID" << body.d_id << std::endl;
    return out;
  }
}
