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
#include <Core/SphereRigidBody.h>
#include <Core/ConvexHullRigidBody.h>
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

  // Set up the ground for bullet
  d_ground_min = SCIRun::Vector(d_domain.lower().x(), d_domain.lower().y(), 
                                d_domain.lower().z()); 
  d_ground_max = SCIRun::Vector(d_domain.upper().x(), d_domain.upper().y(), 
                                d_domain.lower().z()+0.01*d_domain.zrange());

  // Set up the body information
  unsigned long count = 0;

  // A file containing a list of speherical rigid bodies
  Uintah::ProblemSpecP file_ps = ps->findBlock("RigidBodyFile");
  if (file_ps) {

    std::cout << "Reading rigid body file " << std::endl;
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

    // Read the points
    std::vector<Uintah::Vector> positions;
    std::vector<Uintah::Vector> velocities;
    std::vector<double> masses;
    std::vector<double> volumes;
    readPointsFromFile(fileName, positions, velocities, masses, volumes);
   
    // Create spherical rigid bodies
    int count = 0;
    for (auto position : positions) {

      if (count%stride == 0) {
        // Create a Rigid body
        RigidBodySP body = std::make_shared<SphereRigidBody>();
        body->initialize(masses[count], volumes[count], position,
                         velocities[count]*vel_scale_fac,
                         body_force, rot_center, rot_vel); 
        body->id(count);
        d_body_list.emplace_back(body);
      }
      ++count;
    }
    std::cout << "Done reading rigid body file " << std::endl;

  } else {

    // Read the rigid bodies
    for (Uintah::ProblemSpecP body_ps = ps->findBlock("RigidBody"); body_ps != 0;
         body_ps = body_ps->findNextBlock("RigidBody")) {

      // std::cout << count << endl;

      // Initialize the body 
      RigidBodySP body = std::make_shared<SphereRigidBody>();
      body->initialize(body_ps); 
      body->id(count);
      d_body_list.emplace_back(body);
    
      ++count;
      // std::cout << *body;
    }
  }

  // Read convex rigid bodies
  // **WARNING** Assumes the ground location has already been defined
  int convex_body_count = 0;
  for (Uintah::ProblemSpecP body_ps = ps->findBlock("ConvexHullRigidBody"); body_ps != 0;
       body_ps = body_ps->findNextBlock("ConvexHullRigidBody")) {

    std::cout << "Reading convex hull file " << std::endl;
    // Get the file name
    std::string fileName;
    body_ps->require("rigid_body_file", fileName);

    // Get the velocity scale factor
    double vel_scale_fac = 1.0;
    body_ps->require("velocity_scale_factor", vel_scale_fac);

    // Get the body force and rotation information
    // Read the initial body forces
    Uintah::Vector body_force;
    body_ps->require("body_force", body_force);

    // For a rotating coordinate system
    Uintah::Vector rot_center, rot_vel;
    body_ps->require("center_of_rotation", rot_center);
    body_ps->require("angular_velocity_of_rotation", rot_vel);

    // Read the points
    std::vector<Uintah::Vector> positions;
    std::vector<Uintah::Vector> velocities;
    std::vector<double> masses;
    std::vector<double> volumes;
    readPointsFromFile(fileName, positions, velocities, masses, volumes);

    // Create convex hull rigid bodies
    ConvexHullRigidBodySP body = std::make_shared<ConvexHullRigidBody>();
    body->initialize(convex_body_count,
                     masses, volumes, positions, velocities, vel_scale_fac,
                     body_force, rot_center, rot_vel); 
    d_convex_body_list.emplace_back(body);
    ++convex_body_count;
    std::cout << "Done reading convex hull file " << std::endl;
  }

  // Set up bullet
  setupBulletRigidBodies();
}

void 
RigidBodyDynamics::readPointsFromFile(const std::string& fileName,
                                      std::vector<SCIRun::Vector>& positions,
                                      std::vector<SCIRun::Vector>& velocities,
                                      std::vector<double>& masses,
                                      std::vector<double>& volumes)
{

  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open node input rigid body file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read the file
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

    //std::cout << " line = " << line << std::endl;
    // Put the read data into variables
    auto iter = data.begin();
    //double time = std::stod(*iter); 
    ++iter;
    double com_x = std::stod(*iter); ++iter;
    double com_y = std::stod(*iter); ++iter;
    double com_z = std::stod(*iter); ++iter;
    double vel_x = std::stod(*iter); ++iter;
    double vel_y = std::stod(*iter); ++iter;
    double vel_z = std::stod(*iter); ++iter;
    double mass = std::stod(*iter); ++iter;
    double vol = std::stod(*iter); 

    positions.push_back(Uintah::Vector(com_x, com_y, com_z));
    velocities.push_back(Uintah::Vector(vel_x, vel_y, vel_z));
    masses.push_back(mass);
    volumes.push_back(vol);
  }
}
                        

void 
RigidBodyDynamics::problemSetup(Time& time,
                                OutputVTK& output,
                                Domain& domain,
                                RigidBodySPArray& bodyList,
                                ConvexHullRigidBodySPArray& convexBodyList)
{
  d_time.clone(time);
  d_output.clone(output);
  d_domain.clone(domain);
  d_body_list = bodyList;
  d_convex_body_list = convexBodyList;

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
  if (d_body_list.size() > 0) {
    auto iter = d_body_list.begin();
    double radius = (*iter)->radius();
    createRigidBodies(radius);
  }

  // Create the convex hull rigid body list for bullet
  createConvexHullRigidBodies();
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

// The shape is a convex hull
void
RigidBodyDynamics::createConvexHullRigidBodies()
{
  for (auto& convex_body : d_convex_body_list) {

    // Create shape
    btConvexHullShape* originalConvexShape = new btConvexHullShape();

    // Add the points to the shape
    for (auto position : convex_body->getPositions()) {
      originalConvexShape->addPoint(btVector3(position.x(), position.y(), position.z()));
    }

    // Create a reduced convex hull shape
    btShapeHull* hull = new btShapeHull(originalConvexShape);
    btScalar margin = originalConvexShape->getMargin();
    hull->buildHull(margin);
    btConvexHullShape* shape = 
      new btConvexHullShape((const btScalar*) hull->getVertexPointer(), hull->numVertices());
    
    std::cout << shape->getMargin() << std::endl;
    d_collisionShapes.push_back(shape); 

    // Create transforms (no transformation for now)
    btTransform bodyTransform;
    bodyTransform.setIdentity();

    // Create motion state
    btDefaultMotionState* motionState = new btDefaultMotionState(bodyTransform);

    // Get the mass
    btScalar mass(convex_body->mass());

    // Initialize local inertia
    btVector3 localInertia(0.0, 0.0, 0.0);

    // Compute local inertia
    shape->calculateLocalInertia(mass, localInertia);

    // Create the body
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);

    // Set the friction coeffs
    body->setFriction(1);
    body->setRollingFriction(1);

    // Set the initial velocity
    Uintah::Vector init_vel = convex_body->velocity();
    body->setLinearVelocity(btVector3(init_vel.x(), init_vel.y(), init_vel.z()));

    // Set the initial body force acceleration
    Uintah::Vector body_force = convex_body->bodyForce();
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
  d_output.write(d_time, d_domain, d_body_list, d_convex_body_list);

  // Get the ground + walls count
  int static_bodies = (int) d_walls.size() + 1;

  // Get the sphere and convex rigid body count
  int sphere_bodies = (int) d_body_list.size();
  int convex_bodies = (int) d_convex_body_list.size();
  //std::cout << "Static objects = " << static_bodies
  //          << " Sphere objects = " << sphere_bodies
  //           << " Convex objects = " << convex_bodies << std::endl;
  std::cout << "Num collision objects = " << d_world->getNumCollisionObjects() << std::endl;
  std::cout << "Num rigid bodies = " << d_body_list.size() << std::endl;
  std::cout << "Num convex hull rigid bodies = " << d_convex_body_list.size() << std::endl;

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
    for (int jj = d_world->getNumCollisionObjects()-1; jj > 0; jj--) {

      // For sphere rigid bodies
      if (jj >= static_bodies && jj < static_bodies + sphere_bodies) {
        //std::cout << " jj = " << jj;
        //std::cout << " Processing sphere" << std::endl;

        btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
        btRigidBody* body = btRigidBody::upcast(obj);

        if (body && body->getMotionState()) {
          btTransform trans;
          body->getMotionState()->getWorldTransform(trans);

          // Update position and velocity
          Vector3D pos(trans.getOrigin().getX(), trans.getOrigin().getY(),
                       trans.getOrigin().getZ());
          Vector3D vel(body->getLinearVelocity().getX(), body->getLinearVelocity().getY(),
                       body->getLinearVelocity().getZ());

          int index = jj - static_bodies;

          d_body_list[index]->position(pos);
          d_body_list[index]->velocity(vel);
      
          // Get the rotation origin and velocity
          Vector3D oo = d_body_list[index]->rotatingCoordCenter();
          Vector3D omega = d_body_list[index]->rotatingCoordAngularVelocity();
          Vector3D rr = pos - oo;
          Vector3D omegaxr = omega.cross(rr);
          Vector3D omegaxomegaxr = omega.cross(omegaxr);
          Vector3D omegaxv = omega.cross(vel);

          // Update the body force
          Vector3D body_force = d_body_list[index]->bodyForce();
          body_force = body_force - omegaxomegaxr - omegaxv*2.0;

          body->setGravity(btVector3(body_force.x(), body_force.y(), body_force.z()));

        } // end if body && body->getMotionState

      } else { // For convex hull rigid bodies
        if (jj >= static_bodies) {
          //std::cout << " jj = " << jj;
          //std::cout << " Processing convex hull" << std::endl;

          btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];

          btRigidBody* body = btRigidBody::upcast(obj);

          if (body && body->getMotionState()) {
            btTransform trans;
            body->getMotionState()->getWorldTransform(trans);

            btCollisionShape* shape = obj->getCollisionShape();

            int index = jj - static_bodies - sphere_bodies;

            if (shape->isConvex()) {

              // Get the translation of the origin
              //btVector3 translation = trans.getOrigin();

              // Get the rotation 
              //btMatrix3x3 rotation = trans.getBasis();

              // Update positions of vertices
              std::vector<Uintah::Vector> positions;
              std::vector<Uintah::Vector> old_positions = d_convex_body_list[index]->getPositions();
              Uintah::Vector old_com = d_convex_body_list[index]->centerOfMass();
              for (auto position : old_positions) {

                // Rotate and translate the points
                btVector3 old_pos(position.x(), position.y(), position.z());
                btVector3 pos = trans*old_pos;
                positions.push_back(Uintah::Vector(pos.x(), pos.y(), pos.z()));

                /*
                // Get relative positions
                Uintah::Vector old_rel_pos = position - old_com;
                btVector3 old_rel_pos_vec(old_rel_pos.x(), old_rel_pos.y(), old_rel_pos.z());

                // Rotate the points
                btVector3 rot_mat_row1 = rotation[0];
                btVector3 rot_mat_row2 = rotation[1];
                btVector3 rot_mat_row3 = rotation[2];
                btScalar new_rel_pos1 = rot_mat_row1.dot(old_rel_pos_vec);
                btScalar new_rel_pos2 = rot_mat_row2.dot(old_rel_pos_vec);
                btScalar new_rel_pos3 = rot_mat_row3.dot(old_rel_pos_vec);
                btVector3 new_rel_pos_vec(new_rel_pos1, new_rel_pos2, new_rel_pos3);

                // Translate the points
                btVector3 pos = new_rel_pos_vec + translation;
                positions.push_back(Uintah::Vector(pos.x(), pos.y(), pos.z()));
                */
              }
              d_convex_body_list[index]->setPositions(positions);

              // Update velocity
              Uintah::Vector com_pos(trans.getOrigin().getX(), trans.getOrigin().getY(),
                                     trans.getOrigin().getZ());
              d_convex_body_list[index]->setCenterOfMass(com_pos);

              Uintah::Vector vel(body->getLinearVelocity().getX(), body->getLinearVelocity().getY(),
                                 body->getLinearVelocity().getZ());
              d_convex_body_list[index]->setVelocity(vel);
      
              // Get the rotation origin and velocity
              Uintah::Vector oo = d_convex_body_list[index]->rotatingCoordCenter();
              Uintah::Vector omega = d_convex_body_list[index]->rotatingCoordAngularVelocity();
              Uintah::Vector rr = com_pos - oo;
              Uintah::Vector omegaxr = Cross(omega, rr);
              Uintah::Vector omegaxomegaxr = Cross(omega, omegaxr);
              Uintah::Vector omegaxv = Cross(omega, vel);

              // Update the body force
              Uintah::Vector body_force = d_convex_body_list[index]->bodyForce();
              body_force = body_force - omegaxomegaxr - omegaxv*2.0;

              body->setGravity(btVector3(body_force.x(), body_force.y(), body_force.z()));

            } // end if shape isConvex
          } // end if body && body->getMotionState
          //std::cout << " Done processing convex hull" << std::endl;
        } // ceck if jj > static_bodies
      } // end if jj is a convex body index
    } // End loop over collision objects

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
      d_output.write(d_time, d_domain, d_body_list, d_convex_body_list);
      //std::cout << "Wrote out data at time " << d_time << std::endl;
    }

    // Print position, mass, volume
    /*
    for (int jj = d_world->getNumCollisionObjects()-1; jj > 0; jj--) {
      btCollisionObject* obj = d_world->getCollisionObjectArray()[jj];
      btRigidBody* body = btRigidBody::upcast(obj);

      if (body && body->getMotionState()) {
        btTransform trans;
        body->getMotionState()->getWorldTransform(trans);

        //std::cout << " jj = " << jj << " static = " << static_bodies
        //          << " sphere = " << sphere_bodies << std::endl;
        if (jj >= static_bodies) {
          if (jj < static_bodies + sphere_bodies) {
            int index = jj - static_bodies;
            std::cerr << trans.getOrigin().getX() << "," << trans.getOrigin().getY() << "," 
                      << trans.getOrigin().getZ() << "," 
                      << d_body_list[index]->mass() << ","
                      << d_body_list[index]->volume() << std::endl;
          } else {
            //std::cout << " Printing convex hull data " << std::endl;
            int index = jj - static_bodies - sphere_bodies;
            std::cerr << trans.getOrigin().getX() << "," << trans.getOrigin().getY() << "," 
                      << trans.getOrigin().getZ() << "," 
                      << d_convex_body_list[index]->mass() << ","
                      << d_convex_body_list[index]->volume() << std::endl;
          }
        }
      }
    }
    */
  } // End time loop
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

