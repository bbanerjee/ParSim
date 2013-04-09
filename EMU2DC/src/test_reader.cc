// Test the XML input file reader

#include <Domain.h>
#include <ProblemSpecReader.h>
#include <Exception.h>
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

  // Get input and output files/interval
  Uintah::ProblemSpecP io_ps = ps->findBlock("InputOutput");
  std::string input_node_file;
  std::string input_element_file;
  std::string output_file;
  int output_iter_interval = 0;
  io_ps->get("input_node_file", input_node_file);
  io_ps->get("input_element_file", input_element_file);
  io_ps->get("output_file", output_file);
  io_ps->get("output_iteration_interval", output_iter_interval);
  std::cout << "Input files: " << input_node_file << ", " << input_element_file << std::endl;
  std::cout << "Output file: " << output_file << std::endl;
  std::cout << "  Interval: " << output_iter_interval << std::endl;

  // Get the domain information
  Domain domain;
  domain.initialize(ps);
  std::cout << domain ;
  
  ps = 0;  // give up memory held by ps
}
