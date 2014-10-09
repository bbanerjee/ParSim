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

// Test the XML input file reader

#include <Core/SimulationState.h>
#include <Core/Time.h>
#include <InputOutput/Output.h>
#include <Core/Domain.h>
#include <InputOutput/ProblemSpecReader.h>
#include <Core/Exception.h>
#include <MaterialModels/Material.h>
#include <Containers/MaterialSPArray.h>
#include <Core/Body.h>
#include <Containers/BodySPArray.h>
#include <Core/HorizonComputer.h>
#include <iostream>
#include <string>
#include <chrono>

using namespace Matiti;

void test_reader(const std::string& filename);

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "Usage: test_reader <filename>" << std::endl;
    exit(0);
  }

  std::string filename = argv[1];
  try {
    test_reader(filename);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }

  return 0;
}

void test_reader(const std::string& filename)
{
  Uintah::ProblemSpecP ps = ProblemSpecReader().readInputFile(filename);

  if (!ps) {
    throw Exception("Cannot read input file", __FILE__, __LINE__);
  }
  if (ps->getNodeName() != "Vaango") {
    throw Exception("Input file is not a Vaango file specification.", __FILE__, __LINE__);
  }

  // Get the basic simulation information
  SimulationState ss;
  ss.initialize(ps);
  std::cout << ss ;

  // Get the time information
  Time time(ps);
  std::cout << time;
  
  // Get the output information
  Output output(ps);
  std::cout << output;

  // Get the domain information
  Domain domain;
  domain.initialize(ps);
  std::cout << domain;

  // Get the material information
  MaterialSPArray mat_list;
  int count = 0;
  for (Uintah::ProblemSpecP mat_ps = ps->findBlock("Material"); mat_ps != 0;
       mat_ps = mat_ps->findNextBlock("Material")) {

    MaterialSP mat = std::make_shared<Material>();
    mat->initialize(mat_ps); 
    mat->id(count);
    mat_list.emplace_back(mat);
    ++count;
    // std::cout << *mat << std::endl;
  }
  
  // Get the body information
  BodySPArray body_list;
  count = 0;
  for (Uintah::ProblemSpecP body_ps = ps->findBlock("Body"); body_ps != 0;
       body_ps = body_ps->findNextBlock("Body")) {

    // Initialize the body (nodes, elements, cracks)
    auto t1 = std::chrono::high_resolution_clock::now();
    BodySP body = std::make_shared<Body>();
    body->initialize(body_ps, domain, ss, mat_list); 
    body->id(count);
    body_list.emplace_back(body);
    ++count;
    auto t2 = std::chrono::high_resolution_clock::now();
 
    // Compute the horizons of nodes in the body
    HorizonComputer compute_horizon;
    compute_horizon(body, ss);
    auto t3 = std::chrono::high_resolution_clock::now();

    // std::cout << *body << std::endl;
    std::cout << "Body " << count << " : init = " << (t2-t1).count() 
                                  << " horizon = " << (t3-t2).count() << std::endl;

  }

  // Create the family of each node and store
  // Also remove any bonds intersected by cracks
  for (auto iter = body_list.begin(); iter != body_list.end(); ++iter) {

    // Create family
    auto t1 = std::chrono::high_resolution_clock::now();
    (*iter)->createInitialFamily(domain);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Body " << *iter << " : family = " << (t2-t1).count();

    // Remove bonds
    t2 = std::chrono::high_resolution_clock::now();
    (*iter)->removeBondsIntersectedByCracks();
    auto t3 = std::chrono::high_resolution_clock::now();

    std::cout << "Body " << *iter << " : bond removal = " << (t3-t2).count() << std::endl;
  }

  ps = 0;  // give up memory held by ps
}
