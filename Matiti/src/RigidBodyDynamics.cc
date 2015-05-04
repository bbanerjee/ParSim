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
  delete d_config;
  delete d_dispatch;
  delete d_interface;
  delete d_solver;
  delete d_world;
}

void 
RigidBodyDynamics::initializeBullet()
{
  // Set up default collision configuration
  btDefaultCollisionConfiguration* d_config = new btDefaultCollisionConfiguration();

  // Set up collision dispatcher
  btCollisionDispatcher* d_dispatch = new btCollisionDispatcher(d_config);

  // Set up broad phase (for the interface)
  btBroadphaseInterface* d_interface = new btDbvtBroadphase();

  // Set up constraint solver
  btSequentialImpulseConstraintSolver* d_solver = new btSequentialImpulseConstraintSolver();

  // Create the world
  btDiscreteDynamicsWorld* d_world = new btDiscreteDynamicsWorld(d_dispatch, 
                                                                 d_interface, 
                                                                 d_solver, 
                                                                 d_config);
}

void
RigidBodyDynamics::problemSetup(Uintah::ProblemSpecP& ps)
{
  // Set up the simulation state data
  d_state.initialize(ps);
  // std::cout << d_state ;

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

}

void 
RigidBodyDynamics::problemSetup(Time& time,
                                OutputVTK& output,
                                SimulationState& state,
                                Domain& domain,
                                RigidBodySPArray& bodyList)
{
  d_time.clone(time);
  d_output.clone(output);
  d_state.clone(state);
  d_domain.clone(domain);
  d_body_list = bodyList;

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
       
    // Do the computations separately for each body
    for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {    


    } // end of for (body)
    
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

