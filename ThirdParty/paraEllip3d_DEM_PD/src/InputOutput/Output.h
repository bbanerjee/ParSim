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

#ifndef __DEM_OUTPUT_H__
#define __DEM_OUTPUT_H__

#include <mpi.h>

#include <iostream>
#include <string>

namespace dem {

class Output
{
public:
  friend std::ostream& operator<<(std::ostream& out, const dem::Output& output);

public:
  Output() = delete;
  Output(const std::string& fileName, int iterInterval);
  virtual ~Output();

  void clone(const Output& output);

  virtual void write() = 0;

  inline void outputFolder(const std::string& folder)
  {
    d_output_folder_name = folder;
  }
  inline std::string outputFolder() const { return d_output_folder_name; }
  inline std::string outputFile() const { return d_output_file_name; }
  inline int outputIteratonInterval() const { return d_output_iter_interval; }

  int outputFileCount() const { return d_output_file_count; }

  // Set the processor distribution
  void setMPIComm(const MPI_Comm& cartComm) { d_cartComm = cartComm; }

  void setMPIProc(size_t x, size_t y, size_t z)
  {
    d_mpiProcX = x;
    d_mpiProcY = y;
    d_mpiProcZ = z;
  }

protected:
  void incrementOutputFileCount() { d_output_file_count++; }

  // Processor distribution
  MPI_Comm d_cartComm;
  std::size_t d_mpiProcX;
  std::size_t d_mpiProcY;
  std::size_t d_mpiProcZ;

private:
  //  Output file folder and name
  std::string d_output_folder_name;
  std::string d_output_file_name;
  int d_output_iter_interval;

  int d_output_file_count;

  Output(const Output&) = delete;
}; // end class

} // end namespace

#endif
