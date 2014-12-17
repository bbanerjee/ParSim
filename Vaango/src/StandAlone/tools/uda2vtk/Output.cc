/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2014-2015 Parresia Research Limited, New Zealand
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

#include <StandAlone/tools/uda2vtk/Output.h>
#include <unistd.h>
#include <fstream>

using namespace Vaango;

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
    d_output_folder_name = std::string(buffer);
  }
  d_output_file_count = 0;
}

Output::~Output()
{
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
  d_output_file_count = 0;
}

namespace Vaango {

  std::ostream& operator<<(std::ostream& out, const Output& output)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Output dir = " << output.d_output_folder_name << " Output file = " << output.d_output_file_name
        << std::endl;
    out << "  Output iteration interval = " << output.d_output_iter_interval << std::endl;
    return out;
  }

}
