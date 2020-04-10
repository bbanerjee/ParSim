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

#include <Core/ConvexHullRigidBody.h>
#include <Core/Exception.h>

//#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Matiti;

ConvexHullRigidBody::ConvexHullRigidBody()
  : d_id(0), d_com(0.0, 0.0, 0.0), d_vel(0.0, 0.0, 0.0), 
    d_body_force(0.0, 0.0, 0.0), 
    d_rot_center(0.0, 0.0, 0.0), d_rot_vel(0.0, 0.0, 0.0)
{
}

ConvexHullRigidBody::~ConvexHullRigidBody()
{
}

void 
ConvexHullRigidBody::initialize(const int&                         id,
                                const std::vector<double>&         masses,
                                const std::vector<double>&         volumes,
                                const std::vector<Uintah::Vector>& positions,
                                const std::vector<Uintah::Vector>& velocities,
                                const double&                      velScaleFactor, 
                                const Uintah::Vector&              bodyForce,
                                const Uintah::Vector&              centerOfRotation,
                                const Uintah::Vector&              angularVelOfRotation)
{
  // Save the id
  d_id = id;

  // Save mass, volume, center of mass, and average velocity
  d_volume = 0.0;
  for (auto vol : volumes) {
    d_volume += vol;
  }
  int count = 0;
  d_mass = 0.0;
  for (auto mass : masses) {
    d_mass += mass;
    d_com  += positions[count]*mass;
    d_vel  += velocities[count];
    ++count;
  }
  d_com /= d_mass;
  d_vel /= (count/velScaleFactor);

  // Copy positions to local Uintah::Vector
  d_positions = positions;

  // Save body force quantities
  d_body_force = bodyForce;
  d_rot_center = centerOfRotation;
  d_rot_vel    = angularVelOfRotation;
}


namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const ConvexHullRigidBody& body)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "ConvexHullRigidBody: ID" << body.d_id << std::endl;
    return out;
  }
}
