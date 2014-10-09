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

// Unit test for 8 particle body
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <Peridynamics.h>
#include <MaterialModels/Material.h>
#include <Core/Body.h>
#include <Core/Exception.h>

#include <iostream>
#include <string>
#include <chrono>


using namespace boost::unit_test;
using namespace Matiti;

bool init_unit_test();
bool init_ut8Particle();
void ut8Particle();

int main(int argc, char* argv[])
{
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

bool init_unit_test()
{
  init_ut8Particle();
  return true;
}

bool init_ut8Particle()
{  
  bool success = true;
  framework::master_test_suite().
        add(BOOST_TEST_CASE(&ut8Particle));
  return success;
}

void ut8Particle()
{
  // Create a SimulationState object
  bool isDynamic = true;
  Array3 dampingCoeff = {{0.0, 0.0, 0.0}};
  SimulationState::ModulusType modulusType = SimulationState::ModulusType::Constant;
  double horizonFactor = 1.001;
  SimulationState state(isDynamic, dampingCoeff, modulusType, horizonFactor);  

  // Create a Time object
  double maxTime = 1.0;
  double delT = 1.0e-2;
  int maxIter = 4;
  double factor = 1.0;
  double curTime = 0.0;
  Time time(maxTime, delT, maxIter, factor, curTime);

  // Create a Output object
  std::string outFile("ut8ParticleTest");
  int outInterval = 1000;
  OutputVTK output(outFile, outInterval); 

  // Create a Domain object
  Point3D lower(0.0, 0.0, 0.0);
  Point3D upper(1.0, 1.0, 1.0);
  IntArray3 numCells = {{1,1,1}};
  Domain domain(lower, upper, numCells);

  // Create a Material array
  MaterialSP mat = std::make_shared<Material>();
  mat->initialize(0, Material::MicroModulusModel::Constant, 1.0e3, 1.0e7, 1.0);
  MaterialSPArray matList;
  matList.emplace_back(mat);

  // Create the geometry (nodes, elements)
  Point3D boxLower(0.0, 0.0, 0.0);
  Point3D boxUpper(1.0, 1.0, 0.001);
  Uintah::IntVector numElem(1,1,1);
  NodePArray nodes;
  ElementPArray elements;
  Vector3D gridSize;
  BoxGeometryPiece geom(boxLower, boxUpper, numElem, nodes, elements, gridSize);

  // Create the initial conditions
  Vector3D initVel(0.0, 0.0, 0.0);
  Vector3D gravity(0.0, 0.0, 0.0);
  InitialConditions ic(initVel, gravity);

  // Create a Body array
  BodySP body = std::make_shared<Body>();
  body->initialize(0, nodes, elements, gridSize, domain, state, matList, ic);
  BodySPArray bodyList;
  bodyList.emplace_back(body);

  // Compute the horizons of nodes in the body
  HorizonComputer compute_horizon;
  compute_horizon(body, state);

  // Create the initial family of each node and remove any bonds
  // intersected by initial cracks
  for (auto iter = bodyList.begin(); iter != bodyList.end(); ++iter) {

    // Create family
    (*iter)->createInitialFamily(domain);

    // After all bond deletions have been completed update nodal damage indices
    (*iter)->updateDamageIndex();
  }

  // Create a Peridynamics object
  Peridynamics peri;

  // Initialize
  peri.problemSetup(time, output, state, domain, matList, bodyList);

  // Run
  peri.run();

}
