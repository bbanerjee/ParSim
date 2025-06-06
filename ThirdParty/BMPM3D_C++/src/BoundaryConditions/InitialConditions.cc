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

#include <InitialConditions.h>
#include <Node.h>
#include <Crack.h>
#include <CrackSP.h>
#include <BondPArray.h>

#include <Core/ProblemSpec/ProblemSpec.h>

using namespace BrMPM;
  
InitialConditions::InitialConditions()
  : d_initial_velocity(0.0,0.0,0.0),
    d_body_force(0.0, 0.0, 0.0)
{
}

InitialConditions::~InitialConditions()
{
}

void
InitialConditions::readInitialConditions(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP ic_ps = ps->findBlock("InitialConditions");

  // Get the initial velocity and body force
  Uintah::Vector initial_velocity(0.0, 0.0, 0.0), gravity(0.0, 0.0, 0.0);
  ic_ps->require("velocity", initial_velocity);
  ic_ps->require("gravity", gravity);
  for (unsigned int ii = 0; ii < 3; ++ii) {
    d_initial_velocity[ii] = initial_velocity[ii];
    d_body_force[ii] = gravity[ii];
  }

  // Get the initial crack information
  for (Uintah::ProblemSpecP crack_ps = ic_ps->findBlock("Crack"); crack_ps != 0;
       crack_ps = crack_ps->findNextBlock("Crack")) {
    CrackSP crack = std::make_shared<Crack>();
    crack->initialize(crack_ps); 
    d_cracks.emplace_back(crack);
  }

}

void 
InitialConditions::applyInitialVelocity(NodePArray& nodes)
{
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;

    cur_node->velocity(d_initial_velocity);
  }
}

void 
InitialConditions::applyBodyForce(NodePArray& nodes)
{
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;

    double density = cur_node->density();
    double volume = cur_node->volume();
    Vector3D ext_force = cur_node->externalForce();
    ext_force += d_body_force*(density*volume);
    cur_node->externalForce(ext_force);
  }
}

void
InitialConditions::removeBondsIntersectedByCracks(NodePArray& nodes)
{
  // Loop through body nodes
  for (auto node_iter = nodes.begin(); node_iter != nodes.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;

    // Loop through cracks in the body
    BondPArray& bonds = cur_node->getBonds();
    for (auto crack_iter = d_cracks.begin(); crack_iter != d_cracks.end(); ++crack_iter) {
      (*crack_iter)->breakBonds(cur_node, bonds);
    }
  }
}

namespace BrMPM {

  std::ostream& operator<<(std::ostream& out, const InitialConditions& ic)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Initial Conditions: " ;
    out << "Initial velocity = " << ic.d_initial_velocity
        << " Gravity = " << ic.d_body_force << std::endl;
    for (auto iter = (ic.d_cracks).begin(); iter != (ic.d_cracks).end(); ++iter) {
      out << *(*iter) << std::endl ;
    }
    return out;
  }
}
