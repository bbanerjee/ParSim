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

#ifndef MATITI_RIGID_BODY_DYNAMICS_H
#define MATITI_RIGID_BODY_DYNAMICS_H

#include <Core/SimulationState.h>
#include <Core/Time.h>
#include <InputOutput/OutputVTK.h>
#include <Core/Domain.h>
#include <Containers/BodySPArray.h>
#include <BoundaryConditions/VelocityBC.h>

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {

  class RigidBodyDynamics {

  public:
    RigidBodyDynamics();
    ~RigidBodyDynamics();

    void problemSetup(Uintah::ProblemSpecP& ps);
    
    void problemSetup(Time& time,
                      OutputVTK& output,
                      SimulationState& state,
                      Domain& domain,
                      MaterialSPArray& matList,
                      BodySPArray& bodyList);

    void run();

  protected:

    void applyInitialConditions();

    void computeInternalForceDensity(const NodeP node,
                                     Vector3D& internalForce, const Vector3D& gridSize);

    void computeInternalForceDensity(const NodeP cur_node,
                                     Vector3D& internalForce, const Vector3D& gridSize, const int& iter,
                                     std::vector <int>& posBrokenBond);

    void integrateNodalAcceleration(const Vector3D& velocityOld,
                                    const Vector3D& accelerationOld,
                                    double delT,
                                    Vector3D& velocityNew);

    void integrateNodalVelocity(const Vector3D& displacementOld,
                                const Vector3D& velocityOld,
                                double delT,
                                Vector3D& displacementNew);

    void checkMemoryUsage(double& resident_mem, double& shared_mem);

  private:

    Time d_time;
    OutputVTK d_output;
    SimulationState d_state;
    Domain d_domain;
    MaterialSPArray d_mat_list;
    WoodSPArray d_wood_list;
    BodySPArray d_body_list;
    VelocityBC d_velocitybc;

  }; // end class
} // end namespace

#endif

