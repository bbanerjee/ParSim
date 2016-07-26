/*
 * The MIT License
 *
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

// Test the reader, run, and output

#include <RigidBodyDynamics.h>
#include <InputOutput/ProblemSpecReader.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>
#include <chrono>

using namespace Matiti;

void test_bullet(const std::string& filename);

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "Usage: test_bullet <filename>" << std::endl;
    exit(0);
  }

  std::string filename = argv[1];
  try {
    test_bullet(filename);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }

  return 0;
}

void test_bullet(const std::string& filename)
{
  Uintah::ProblemSpecP ps = ProblemSpecReader().readInputFile(filename);

  if (!ps) {
    throw Exception("Cannot read input file", __FILE__, __LINE__);
  }
  if (ps->getNodeName() != "Vaango") {
    throw Exception("Input file is not a Vaango file specification.", __FILE__, __LINE__);
  }

  // Create a Rigid Body dynamics object
  RigidBodyDynamics rigid;

  // Initialize
  auto t1 = std::chrono::high_resolution_clock::now();
  rigid.problemSetup(ps);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Bullet rigid body dynamics: Problem setup time = " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
            << " msec" << std::endl;

  // Run
  rigid.run();
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "Bullet rigid body dynamics: Run time = " 
            << std::chrono::duration_cast<std::chrono::seconds>(t3-t2).count() 
            << " sec" << std::endl;

  ps = 0;  // give up memory held by ps
}
