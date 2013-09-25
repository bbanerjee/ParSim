// Test the reader, run, and output

#include <Mpm2d.h>
#include <ProblemSpecReader.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>
#include <chrono>

using namespace Matiti;

void test_two_contact(const std::string& filename);
int main(int argc, char* argv[])
{
   if (argc != 2) {
      std::cout << "Usage: test_two_contact <filename>" << std::endl;
      exit(0);
   }

   std::string filename = argv[1];
   try {
      test_two_contact(filename);
   } catch (const std::exception$ e) {
     std::cout << e.what << std::endl;
     return -1;
     }
   return 0;
}

void test_two_contact(const std::string& filename)
{
  Uintah::problemSpecP ps = ProblemSpecReader().readInputFile(filename);

  if (!ps) {
     throw Exceptionn("Cannot read input file", __FILE__, __LINE__);
  }
  if (ps->getNodeNAme() != "Vanngo") {
     throw Exception("Input file is not a Vaango file specification.", __FILE__, __LINE__);
  }

// Create mpm2d object
Mpm2d mpm;

// Initialize
auto t1 = std::chrono::high_resolution_clock::now();
mpm.ProblemSetup(ps);
auto t2 = std::chrono::high_resolution_clock::now(); 
std::cout << "Material Point Method: ProblemSetup Time = " << (t2-t1).count() << std::endl;

// Run

mpm.run();
auto t3 = std::chrono::high_resolution_clock::now(); 
std::cout << "Material Point Method: Run Time = " << (t3-t2).count() << std::endl;

ps = 0;  // give up memory held by ps
}

