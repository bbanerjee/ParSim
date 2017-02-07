/*
 * The MIT License
 *
 * Copyright (c) 2017- Parresia Research Limited, New Zealand
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

#include <unistd.h>
#include <fstream>
#include <exception>

using namespace dem;

Output::Output(const std::string& fileName, int iterInterval) {
  d_output_file_name = fileName;
  d_output_iter_interval = iterInterval;
  char buffer[2000];
  char* str = getcwd(buffer, 2000);
  if (str == NULL) {
    std::cout << "**ERROR** Directory not returned by getcwd()" << __FILE__ <<
                    __LINE__ << "\n";
  } else {
    d_output_folder_name = std::string(buffer) + ".vtk";
  }
  d_output_file_count = 0;

  d_cartComm = nullptr;
  d_mpiProcX = 1;
  d_mpiProcY = 1;
  d_mpiProcZ = 1;
}

Output::~Output() {}

void Output::clone(const Output& output) {
  d_output_file_name = output.d_output_file_name;
  d_output_iter_interval = output.d_output_iter_interval;
  d_output_folder_name = output.d_output_folder_name;
  d_output_file_count = output.d_output_file_count;
  d_mpiProcX = output.d_mpiProcX;
  d_mpiProcY = output.d_mpiProcY;
  d_mpiProcZ = output.d_mpiProcZ;
}

void Output::write() {
  std::cout << "Please call the writer routine from the correct derived class.\n";
}

namespace dem {

std::ostream& operator<<(std::ostream& out, const Output& output) {
  out.setf(std::ios::floatfield);
  out.precision(6);
  out << "Output dir = " << output.d_output_folder_name
      << " Output file = " << output.d_output_file_name << std::endl;
  out << "  Output iteration interval = " << output.d_output_iter_interval
      << std::endl;
  return out;
}
}
