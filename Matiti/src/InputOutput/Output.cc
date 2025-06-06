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

#include <InputOutput/Output.h>
#include <Core/Exception.h>
#include <Containers/NodePArray.h>
#include <Core/Body.h>
#include <Core/Node.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <unistd.h>
#include <fstream>

using namespace Matiti;

Output::Output()
{
}

Output::Output(const std::string& fileName,
               int iterInterval)
{
  d_output_file_name = fileName;
  d_output_iter_interval = iterInterval;
  char buffer[2000];
  char * str = getcwd( buffer, 2000 );
  if (str == NULL) {
    throw Exception("**ERROR** Directory not returned by getcwd()", __FILE__, __LINE__); 
  } else {
    d_output_folder_name = std::string(buffer) + ".vtk";
  }
  d_output_file_count = 0;
}

Output::Output(const Uintah::ProblemSpecP& ps)
{
  initialize(ps);
}

Output::~Output()
{
}

void
Output::clone(const Output& output)
{
 d_output_file_name = output.d_output_file_name; 
 d_output_iter_interval = output.d_output_iter_interval; 
 d_output_folder_name = output.d_output_folder_name;
 d_output_file_count = output.d_output_file_count;
}

void
Output::initialize(const Uintah::ProblemSpecP& ps)
{
  // Get output files/interval
  Uintah::ProblemSpecP io_ps = ps->findBlock("Output");
  if (!io_ps) {
    throw Exception("**ERROR** <Output> tag not found", __FILE__, __LINE__); 
  } 

  io_ps->require("output_file", d_output_file_name);
  io_ps->require("output_iteration_interval", d_output_iter_interval);
  char buffer[2000];
  char * str = getcwd( buffer, 2000 );
  if (str == NULL) {
    throw Exception("**ERROR** Directory not returned by getcwd()", __FILE__, __LINE__); 
  } else {
    d_output_folder_name = std::string(buffer);
  }
  std::cout << "Output folder = " << d_output_folder_name << std::endl;
  d_output_file_count = 0;
}

void
Output::write(const Time& time, const Domain&, 
              const RigidBodySPArray& bodyList,
              const ConvexHullRigidBodySPArray& convexBodyList) 
{
  std::cout << "**WARNING** Not implemented yet." << std::endl;
}

void
Output::write(const Time& time, const Domain& , const BodySPArray& bodyList) 
{
  // Write the output to individual files
  // std::string output_file_name = d_output_file_name;
  // std::string current_output_file_name = d_output_file_name + std::to_string(d_output_file_count) + ".tec";
  //write(current_output_file_name, fmt='(A,I5.5,A)') trim(output_file_name), output_file_count, '.tec'

  // std::ofstream output_file(current_output_file_name);

  // Write the output to individual files
  std::ostringstream of_name;
  //of_name.setf(std::ios::basefield);
  of_name.precision(5);
  of_name << outputFile() << outputFileCount() << ".vtu"; 
  std::ofstream output_file(of_name.str());

  // Loop through bodies
  for (auto body_iter = bodyList.begin(); body_iter != bodyList.end(); ++body_iter) {

    const NodePArray& node_list = (*body_iter)->nodes();

    // count valid nodes
    int valid_node_count = 0;
    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
      if ((*node_iter)->omit()) continue;  // skip this node
      valid_node_count++;
    }

    // write the headers
    output_file << "TITLE=\"simulation results\" " << std::endl;
    output_file << "VARIABLES=\"X\",\"Y\",\"Z\",\"DX\",\"DY\",\"DZ\",\"VX\",\"VY\",\"VZ\",\"DAM\"" << std::endl;
    output_file << "ZONE I=" << valid_node_count << " SOLUTIONTIME=" << time.currentTime() 
                << " F=POINT" << std::endl;

    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
      NodeP cur_node = *node_iter;
      if (cur_node->omit()) continue;  // skip this node
      double xdisp = cur_node->displacement()[0];
      double ydisp = cur_node->displacement()[1];
      double zdisp = cur_node->displacement()[2];
      double cur_x_pos = cur_node->position().x() + xdisp;
      double cur_y_pos = cur_node->position().y() + ydisp;
      double cur_z_pos = cur_node->position().z() + zdisp;
      output_file << cur_x_pos << " " << cur_y_pos << " " << cur_z_pos << " " 
                  << xdisp << " " << ydisp << " " << zdisp
                  << cur_node->velocity()[0] << cur_node->velocity()[1] 
                  << cur_node->velocity()[2] 
                  << cur_node->damageIndex() << std::endl;
    }
  }

  output_file.close();
 
  // Increment the output file count
  incrementOutputFileCount();

  //of_name.setf(std::ios::floatfield, std::ios::basefield);
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Output& output)
  {
    //out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Output dir = " << output.d_output_folder_name << " Output file = " << output.d_output_file_name
        << std::endl;
    out << "  Output iteration interval = " << output.d_output_iter_interval << std::endl;
    return out;
  }

}
