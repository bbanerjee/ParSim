// Test the reader, run, and output

#include <Peridynamics.h>
#include <InputOutput/ProblemSpecReader.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>
#include <chrono>

using namespace Matiti;

void test_peri(const std::string& filename);

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "Usage: test_peri <filename>" << std::endl;
    exit(0);
  }

  std::string filename = argv[1];
  try {
    test_peri(filename);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return -1;
  }

  return 0;
}

void test_peri(const std::string& filename)
{
  Uintah::ProblemSpecP ps = ProblemSpecReader().readInputFile(filename);

  if (!ps) {
    throw Exception("Cannot read input file", __FILE__, __LINE__);
  }
  if (ps->getNodeName() != "Vaango") {
    throw Exception("Input file is not a Vaango file specification.", __FILE__, __LINE__);
  }

  // Create a Peridynamics object
  Peridynamics peri;

  // Initialize
  auto t1 = std::chrono::high_resolution_clock::now();
  peri.problemSetup(ps);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Peridynamics: Problem setup time = " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
            << " msec" << std::endl;

  // Run
  peri.run();
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "Peridynamics: Run time = " 
            << std::chrono::duration_cast<std::chrono::seconds>(t3-t2).count() 
            << " sec" << std::endl;

  ps = 0;  // give up memory held by ps
}
