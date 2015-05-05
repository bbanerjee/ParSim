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

#include <RigidBodyDynamics.h>
#include <Core/RigidBody.h>
#include <Core/Exception.h>

#include <Pointers/RigidBodySP.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>
#include <chrono>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>

#include <errno.h>
#include <fenv.h>

using namespace Matiti;

RigidBodyDynamics::RigidBodyDynamics() 
{
  // Set up bullet physics
  initializeBullet();
}

RigidBodyDynamics::~RigidBodyDynamics() 
{
  deleteBulletRigidBodies();
  deleteBulletShapes();
  delete d_world;
  delete d_solver;
  delete d_interface;
  delete d_dispatch;
  delete d_config;
}

void 
RigidBodyDynamics::initializeBullet()
{
  // Set up default collision configuration
  d_config = new btDefaultCollisionConfiguration();

  // Set up collision dispatcher
  d_dispatch = new btCollisionDispatcher(d_config);

  // Set up broad phase (for the interface)
  d_interface = new btDbvtBroadphase();

  // Set up constraint solver
  d_solver = new btSequentialImpulseConstraintSolver();

  // Create the world
  d_world = new btDiscreteDynamicsWorld(d_dispatch, d_interface, d_solver, d_config);
}

void
RigidBodyDynamics::problemSetup(Uintah::ProblemSpecP& ps)
{
  // Set up the time information
  d_time.initialize(ps);
  // std::cout << d_time ;

  // Set up the output information
  d_output.initialize(ps);
  // std::cout << d_output ;

  // Set up the domain
  d_domain.initialize(ps);
  // std::cout << d_domain ;

  // Set up the body information
  unsigned long count = 0;
  for (Uintah::ProblemSpecP body_ps = ps->findBlock("RigidBody"); body_ps != 0;
       body_ps = body_ps->findNextBlock("RigidBody")) {

    // std::cout << count << endl;

    // Initialize the body (nodes, elements, cracks)
    RigidBodySP body = std::make_shared<RigidBody>();
    body->initialize(body_ps); 
    body->id(count);
    d_body_list.emplace_back(body);
    
    ++count;
    // std::cout << *body;
  }

  // Set up bullet
  setupBulletRigidBodies();
}

void 
RigidBodyDynamics::problemSetup(Time& time,
                                OutputVTK& output,
                                Domain& domain,
                                RigidBodySPArray& bodyList)
{
  d_time.clone(time);
  d_output.clone(output);
  d_domain.clone(domain);
  d_body_list = bodyList;

  // Set up bullet
  setupBulletRigidBodies();
}

void
RigidBodyDynamics::setupBulletRigidBodies()
{
  // Create the ground for bullet
  Vector3D groundMin(d_domain.lower().x(), d_domain.lower().y(), d_domain.lower().z()); 
  Vector3D groundMax(d_domain.upper().x(), d_domain.upper().y(), 
                     d_domain.lower().z()+0.001*d_domain.zrange());
  createGround(groundMin, groundMax);

  // Create the rigid body list for bullet
  auto iter = d_body_list.begin();
  double radius = (*iter)->radius();
  createRigidBodies(radius);
}

// The shape is a box **NOTE** Generalize later
void
RigidBodyDynamics::createGround(const Vector3D& boxMin, const Vector3D& boxMax)
{
  // Get box dimensions
  double xLen = boxMax[0] - boxMin[0];
  double yLen = boxMax[1] - boxMin[1];
  double zLen = boxMax[2] - boxMin[2];

  // Create shape
  btCollisionShape* shape = 
    new btBoxShape(btVector3(btScalar(xLen), btScalar(yLen), btScalar(zLen)));
  d_collisionShapes.push_back(shape);

  // Move the shape to the right location
  btTransform groundTransform;
  groundTransform.setIdentity();
  groundTransform.setOrigin(btVector3(boxMin[0], boxMin[1], boxMin[2]));

  // Create motion state
  btDefaultMotionState* motionState = new btDefaultMotionState(groundTransform);

  // The ground does not move
  btScalar mass(0.0);
  btVector3 localInertia(0.0, 0.0, 0.0);

  // Create the body
  btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
  btRigidBody* body = new btRigidBody(rbInfo);

  // Add the body to the dynamics world
  d_world->addRigidBody(body);
}

// The shape is a sphere **NOTE** Generalize later
void
RigidBodyDynamics::createRigidBodies(const double& radius)
                                     
{
  // Create shape
  btCollisionShape* shape = new btSphereShape(radius);
  d_collisionShapes.push_back(shape); 

  // Create transforms
  btTransform bodyTransform;

  // Initialize local inertia
  btVector3 localInertia(0.0, 0.0, 0.0);

  // Loop thru the list of bodies
  for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {    
    
    // Set the position
    Vector3D pos =  (*body_iter)->centerOfMass();
    bodyTransform.setIdentity();
    bodyTransform.setOrigin(btVector3(pos[0], pos[1], pos[2]));

    // Create motion state
    btDefaultMotionState* motionState = new btDefaultMotionState(bodyTransform);

    // Get the mass
    btScalar mass((*body_iter)->mass());

    // Compute local inertia
    shape->calculateLocalInertia(mass, localInertia);

    // Create the body
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);

    // Add the body to the dynamics world
    d_world->addRigidBody(body);
  }
}

void
RigidBodyDynamics::run()
{
  std::cerr << "....Begin Solver...." << std::endl;

  // Constants
  Vector3D Zero(0.0, 0.0, 0.0);

  // Write the output at the beginning of the simulation
  d_output.write(d_time, d_domain, d_body_list);

  // Do time incrementation
  double cur_time = 0.0;
  int cur_iter = 1;

  while (cur_time < d_time.maxTime() && cur_iter < d_time.maxIter()) {
   
    auto t1 = std::chrono::high_resolution_clock::now();

    // Get the current delT
    double delT = d_time.delT();
    d_world->stepSimulation(delT, 10);

    // Set the body force
    Vector3D body_force = d_body_list[0]->bodyForce();
    d_world->setGravity(btVector3(body_force.x(), body_force.y(), body_force.z()));

    // Loop through the rigid bodies
    for (int jj = d_world->getNumCollisionObjects()-1; jj >= 0; jj--) {
      btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
      btRigidBody* body = btRigidBody::upcast(obj);

      if (body && body->getMotionState()) {
        btTransform trans;
        body->getMotionState()->getWorldTransform(trans);

        if (jj > 0) {
          d_body_list[jj-1]->position(Vector3D(trans.getOrigin().getX(),
                                               trans.getOrigin().getY(),
                                               trans.getOrigin().getZ()));
        }
      }
    }


    //  Get memory usage
    double res_mem = 0.0, shar_mem = 0.0;
    checkMemoryUsage(res_mem, shar_mem);

    //  Print out current time and iteration
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time = " << cur_time << " Iteration = " << cur_iter 
              << " Compute time = " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
              << " ms"
              << " Memory = " << (res_mem-shar_mem)/1024 << " MB" << std::endl;

    // Update the current time and current iteration
    cur_time = d_time.incrementTime(delT);
    ++cur_iter;
        
    // Output nodal information every snapshots_frequency iteration   
    int output_freq = d_output.outputIteratonInterval();
    if (cur_iter%output_freq == 0) {
      d_output.write(d_time, d_domain, d_body_list);
       //std::cout << "Wrote out data at time " << d_time << std::endl;
    }
  }
}

void
RigidBodyDynamics::deleteBulletRigidBodies()
{
  // Loop through the rigid bodies
  for (int jj = d_world->getNumCollisionObjects()-1; jj >= 0; jj--) {
    btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
    btRigidBody* body = btRigidBody::upcast(obj);

    if (body && body->getMotionState()) {
      delete body->getMotionState();
    }
    d_world->removeCollisionObject(obj);
    delete obj;
  }
}

void
RigidBodyDynamics::deleteBulletShapes()
{
  // Loop through the shapes
  for (int jj = 0; jj < d_collisionShapes.size(); jj++) {
    btCollisionShape* shape = d_collisionShapes[jj];
    d_collisionShapes[jj] = 0;
    delete shape;
  }
  d_collisionShapes.clear();
}

void
RigidBodyDynamics::checkMemoryUsage(double& resident_mem, double& shared_mem)
{
  int tSize = 0, resident = 0, share = 0;
  std::ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  resident_mem = resident * page_size_kb;
  shared_mem = share * page_size_kb;
}

