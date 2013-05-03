// Test the XML input file reader

#include <SimulationState.h>
#include <Domain.h>
#include <ProblemSpecReader.h>
#include <Exception.h>
#include <Material.h>
#include <MaterialSPArray.h>
#include <Body.h>
#include <BodySPArray.h>
#include <HorizonComputer.h>
#include <iostream>
#include <string>

using namespace Emu2DC;

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

  // Get the domain information
  Domain domain;
  domain.initialize(ps);
  std::cout << domain ;

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
    BodySP body = std::make_shared<Body>();
    body->initialize(body_ps, domain, mat_list); 
    body->id(count);
    body_list.emplace_back(body);
    ++count;
 
    // Compute the horizons of nodes in the body
    HorizonComputer compute_horizon;
    compute_horizon(body, ss);

    // std::cout << *body << std::endl;
  }

  // Create the family of each node and store
  for (auto iter = body_list.begin(); iter != body_list.end(); ++iter) {
    (*iter)->createInitialFamily(domain);
  }

  // Do some processing to remove bonds intersected by cracks
  ps = 0;  // give up memory held by ps
}
