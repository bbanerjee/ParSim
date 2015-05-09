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
  delete d_broadphase;
  delete d_dispatch;
  delete d_config;
}

void 
RigidBodyDynamics::initializeBullet()
{
  // Set up default collision configuration
  d_config = new btDefaultCollisionConfiguration();
  d_config->setConvexConvexMultipointIterations(7); // Needed for contact detection
                                                    // for small objects

  // Set up collision dispatcher
  d_dispatch = new btCollisionDispatcher(d_config);

  // Set up broad phase (for the interface)
  d_broadphase = new btDbvtBroadphase();
  //btVector3 worldAabbMin(-1000,-1000,-1000);
  //btVector3 worldAabbMax(1000,1000,1000);

  //btHashedOverlappingPairCache* pairCache = new btHashedOverlappingPairCache();
  //d_broadphase = new btAxisSweep3(worldAabbMin,worldAabbMax,3500,pairCache);

  // Set up constraint solver
  d_solver = new btSequentialImpulseConstraintSolver();
  //btDantzigSolver* mlcp = new btDantzigSolver();
  //d_solver = new btMLCPSolver(mlcp);

  // Create the world
  d_world = new btDiscreteDynamicsWorld(d_dispatch, d_broadphase, d_solver, d_config);
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

  // Get the walls (**TODO** Use the domain wall BCs in the next version)
  for (Uintah::ProblemSpecP wall_ps = ps->findBlock("Wall"); wall_ps != 0;
         wall_ps = wall_ps->findNextBlock("Wall")) {

    // Get the wall box
    SCIRun::Vector wall_min, wall_max;
    wall_ps->require("wall_min", wall_min);
    wall_ps->require("wall_max", wall_max);
    Wall wall;
    wall.box_min = wall_min;
    wall.box_max = wall_max;
    d_walls.push_back(wall);
  }

  // Set up the body information
  Uintah::ProblemSpecP file_ps = ps->findBlock("RigidBodyFile");
  if (file_ps) {
    // Get the file name
    std::string fileName;
    file_ps->require("rigid_body_file", fileName);

    // Get the particle stride (number to skip in file)
    int stride = 1;
    file_ps->require("particle_stride", stride);

    // Get the ground box
    file_ps->require("ground_min", d_ground_min);
    file_ps->require("ground_max", d_ground_max);

    // Get the velocity scale factor
    double vel_scale_fac = 1.0;
    file_ps->require("velocity_scale_factor", vel_scale_fac);

    // Get the body force and rotation information
    // Read the initial body forces
    Uintah::Vector body_force;
    file_ps->require("body_force", body_force);

    // For a rotating coordinate system
    Uintah::Vector rot_center, rot_vel;
    file_ps->require("center_of_rotation", rot_center);
    file_ps->require("angular_velocity_of_rotation", rot_vel);

    // Try to open file
    std::ifstream file(fileName);
    if (!file.is_open()) {
      std::string out = "Could not open node input rigid body file " + fileName + " for reading \n";
      throw Exception(out, __FILE__, __LINE__);
    }

    // Read the file
    unsigned long count = 0;
    std::string line;
    while (std::getline(file, line)) { 
 
      // Ignore empty lines
      if (line.empty()) continue;

      // Tokenize the line
      std::string data_piece;
      std::istringstream data_stream(line);
      std::vector<std::string> data;
      while (std::getline(data_stream, data_piece, ' ')) {
        data.push_back(data_piece);
      }
      if (data.size() < 9) {
        std::ostringstream out;
        out << "Could not read node input rigid body file line "
            << line << std::endl;
        throw Exception(out.str(), __FILE__, __LINE__);
      }

      if (count%stride == 0) {
        // Put the read data into variables
        auto iter = data.begin();
        //double time = std::stod(*iter); 
        ++iter;
        double com_x = std::stod(*iter); ++iter;
        double com_y = std::stod(*iter); ++iter;
        double com_z = std::stod(*iter); ++iter;
        double vel_x = std::stod(*iter)*vel_scale_fac; ++iter;
        double vel_y = std::stod(*iter)*vel_scale_fac; ++iter;
        double vel_z = std::stod(*iter)*vel_scale_fac; ++iter;
        double mass = std::stod(*iter); ++iter;
        double vol = std::stod(*iter); 

        // Create a Rigid body
        RigidBodySP body = std::make_shared<RigidBody>();
        body->initialize(mass, vol, SCIRun::Vector(com_x, com_y, com_z),
                         SCIRun::Vector(vel_x, vel_y, vel_z),
                         body_force, rot_center, rot_vel); 
        body->id(count);
        d_body_list.emplace_back(body);
      }
      ++count;
    }

  } else {
    // Set up the ground for bullet
    d_ground_min = SCIRun::Vector(d_domain.lower().x(), d_domain.lower().y(), 
                                 d_domain.lower().z()); 
    d_ground_max = SCIRun::Vector(d_domain.upper().x(), d_domain.upper().y(), 
                       d_domain.lower().z()+0.01*d_domain.zrange());

    // Read the rigid bodies
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

  // Set up the ground for bullet
  d_ground_min = SCIRun::Vector(d_domain.lower().x(), d_domain.lower().y(), 
                                 d_domain.lower().z()); 
  d_ground_max = SCIRun::Vector(d_domain.upper().x(), d_domain.upper().y(), 
                       d_domain.lower().z()+0.01*d_domain.zrange());

  // Set up bullet
  setupBulletRigidBodies();
}

void
RigidBodyDynamics::setupBulletRigidBodies()
{
  // Create walls
  createWalls();

  // Create the ground
  createGround();

  // Create the rigid body list for bullet
  auto iter = d_body_list.begin();
  double radius = (*iter)->radius();
  createRigidBodies(radius);
}

void 
RigidBodyDynamics::createWalls()
{
  for (auto iter = d_walls.begin(); iter != d_walls.end(); iter++) {
    SCIRun::Vector wall_min = (*iter).box_min;  
    SCIRun::Vector wall_max = (*iter).box_max;  

    // Get box dimensions
    double xLen = std::abs(wall_max[0] - wall_min[0]);
    double yLen = std::abs(wall_max[1] - wall_min[1]);
    double zLen = std::abs(wall_max[2] - wall_min[2]);

    // Create shape
    btCollisionShape* shape = 
      new btBoxShape(btVector3(btScalar(xLen), btScalar(yLen), btScalar(zLen)));
    //shape->setMargin(0.05);
    d_collisionShapes.push_back(shape);

    // Move the shape to the right location
    btTransform wallTransform;
    wallTransform.setIdentity();
    wallTransform.setOrigin(btVector3(0.5*(wall_min[0]+wall_max[0]), 
                                      0.5*(wall_min[1]+wall_max[1]), 
                                      0.5*(wall_min[2]+wall_max[2])));

    // Create motion state
    btDefaultMotionState* motionState = new btDefaultMotionState(wallTransform);

    // The wall does not move
    btScalar mass(0.0);
    btVector3 localInertia(0.0, 0.0, 0.0);

    // Create the body
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);
    body->setFriction(0.1);
    body->setRollingFriction(0.1);

    // Add the body to the dynamics world
    d_world->addRigidBody(body);
  }
}

// The shape is a box **NOTE** Generalize later
void
RigidBodyDynamics::createGround()
{
  // Get box dimensions
  double xLen = d_ground_max[0] - d_ground_min[0];
  double yLen = d_ground_max[1] - d_ground_min[1];
  double zLen = d_ground_max[2] - d_ground_min[2];

  // Create shape
  btCollisionShape* shape = 
    new btBoxShape(btVector3(btScalar(xLen), btScalar(yLen), btScalar(zLen)));
  //shape->setMargin(0.05);
  d_collisionShapes.push_back(shape);

  // Move the shape to the right location
  btTransform groundTransform;
  groundTransform.setIdentity();
  groundTransform.setOrigin(btVector3(d_ground_min[0], d_ground_min[1], d_ground_min[2]));

  // Create motion state
  btDefaultMotionState* motionState = new btDefaultMotionState(groundTransform);

  // The ground does not move
  btScalar mass(0.0);
  btVector3 localInertia(0.0, 0.0, 0.0);

  // Create the body
  btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
  btRigidBody* body = new btRigidBody(rbInfo);
  body->setFriction(1);
  body->setRollingFriction(1);

  // Add the body to the dynamics world
  d_world->addRigidBody(body);
}

// The shape is a sphere **NOTE** Generalize later
void
RigidBodyDynamics::createRigidBodies(const double& radius)
                                     
{
  // Create shape
  btCollisionShape* shape = new btSphereShape(radius);
  std::cout << shape->getMargin() << std::endl;
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

    // Set the friction coeffs
    body->setFriction(1);
    body->setRollingFriction(1);

    // Set the initial velocity
    Vector3D init_vel = (*body_iter)->velocity();
    body->setLinearVelocity(btVector3(init_vel.x(), init_vel.y(), init_vel.z()));

    // Set the initial body force acceleration
    Vector3D body_force = (*body_iter)->bodyForce();
    body->setGravity(btVector3(body_force.x(), body_force.y(), body_force.z()));

    // Add the body to the dynamics world
    d_world->addRigidBody(body);
    d_world->setGravity(btVector3(0.0, 0.0, 0.0));
  }
}

void
RigidBodyDynamics::run()
{
  std::cout << "....Begin Solver...." << std::endl;

  // Write the output at the beginning of the simulation
  d_output.write(d_time, d_domain, d_body_list);

  // Get the ground + walls count
  int static_bodies = (int) d_walls.size() + 1;

  // Do time incrementation
  double cur_time = 0.0;
  int cur_iter = 1;
  while (cur_time < d_time.maxTime() && cur_iter < d_time.maxIter()) {
   
    auto t1 = std::chrono::high_resolution_clock::now();

    // Get the current delT
    double delT = d_time.delT();
    d_world->stepSimulation(delT, 10, 1.0/480.0); // Third argument is needed for
                                                  // contact detetion of small objects

    // Loop through the rigid bodies
    std::cout << "Num collision objects = " << d_world->getNumCollisionObjects() << std::endl;
    std::cout << "Num rigid bodies = " << d_body_list.size() << std::endl;
    for (int jj = d_world->getNumCollisionObjects()-1; jj >= 0; jj--) {
      btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
      btRigidBody* body = btRigidBody::upcast(obj);

      if (body && body->getMotionState()) {
        btTransform trans;
        body->getMotionState()->getWorldTransform(trans);

        if (jj > static_bodies) {
          // Update position and velocity
          Vector3D pos(trans.getOrigin().getX(), trans.getOrigin().getY(),
                       trans.getOrigin().getZ());
          Vector3D vel(body->getLinearVelocity().getX(), body->getLinearVelocity().getY(),
                       body->getLinearVelocity().getZ());
          d_body_list[jj-static_bodies]->position(pos);
          d_body_list[jj-static_bodies]->velocity(vel);
      
          // Get the rotation origin and velocity
          Vector3D oo = d_body_list[jj-static_bodies]->rotatingCoordCenter();
          Vector3D omega = d_body_list[jj-static_bodies]->rotatingCoordAngularVelocity();
          Vector3D rr = pos - oo;
          Vector3D omegaxr = omega.cross(rr);
          Vector3D omegaxomegaxr = omega.cross(omegaxr);
          Vector3D omegaxv = omega.cross(vel);

          // Update the body force
          Vector3D body_force = d_body_list[jj-static_bodies]->bodyForce();
          body_force = body_force - omegaxomegaxr - omegaxv*2.0;

          body->setGravity(btVector3(body_force.x(), body_force.y(), body_force.z()));

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

  // Print position, mass, volume
  for (int jj = d_world->getNumCollisionObjects()-1; jj >= 0; jj--) {
    btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
    btRigidBody* body = btRigidBody::upcast(obj);

    if (body && body->getMotionState()) {
      btTransform trans;
      body->getMotionState()->getWorldTransform(trans);

      if (jj > static_bodies) {
        std::cerr << trans.getOrigin().getX() << "," << trans.getOrigin().getY() << "," 
                  << trans.getOrigin().getZ() << "," 
                  << d_body_list[jj-static_bodies]->mass() << ","
                  << d_body_list[jj-static_bodies]->volume() << std::endl;
      }
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

