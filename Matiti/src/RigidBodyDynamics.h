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

#include <Core/Time.h>
#include <InputOutput/OutputVTK.h>
#include <Core/Domain.h>
#include <Containers/RigidBodySPArray.h>

#include <btBulletDynamicsCommon.h>
#include <BulletDynamics/MLCPSolvers/btDantzigSolver.h>
#include <BulletDynamics/MLCPSolvers/btSolveProjectedGaussSeidel.h>
#include <BulletDynamics/MLCPSolvers/btMLCPSolver.h>
#include <BulletCollision/CollisionShapes/btShapeHull.h>


#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {

  class RigidBodyDynamics {

  public:
    RigidBodyDynamics();
    ~RigidBodyDynamics();

    void problemSetup(Uintah::ProblemSpecP& ps);
    
    void problemSetup(Time& time,
                      OutputVTK& output,
                      Domain& domain,
                      RigidBodySPArray& bodyList,
                      ConvexHullRigidBodySPArray& convexBodyList);

    void run();

    void checkMemoryUsage(double& resident_mem, double& shared_mem);

    void createWalls();

    void createGround();

    void createRigidBodies(const double& radius);

    void createRigidBodies(const std::vector<Point3D>& points);

    void createConvexHullRigidBodies();

  private:

    Time d_time;
    OutputVTK d_output;
    Domain d_domain;
    RigidBodySPArray d_body_list;
    ConvexHullRigidBodySPArray d_convex_body_list;

    Uintah::Vector d_ground_min;
    Uintah::Vector d_ground_max;

    struct Wall {
      Uintah::Vector box_min;
      Uintah::Vector box_max;
    }; 
    std::vector<Wall> d_walls;

    // Bullet setup 
    btDefaultCollisionConfiguration* d_config; // Collision configuration
    btCollisionDispatcher* d_dispatch; // Collision dispatcher
    btBroadphaseInterface* d_broadphase; // Broad phase (for the interface)
    btSequentialImpulseConstraintSolver* d_solver; // Constraint solver
    //btMLCPSolver* d_solver; // Constraint solver
    btDiscreteDynamicsWorld* d_world; // Dynamics world
  
    // Create empty array
    btAlignedObjectArray<btCollisionShape*> d_collisionShapes;

  private:

    void readPointsFromFile(const std::string& fileName,
                            std::vector<Uintah::Vector>& positions,
                            std::vector<Uintah::Vector>& velocities,
                            std::vector<double>& masses,
                            std::vector<double>& volumes);
    void initializeBullet();
    void setupBulletRigidBodies();
    void deleteBulletRigidBodies();
    void deleteBulletShapes();

  }; // end class
} // end namespace

#endif

