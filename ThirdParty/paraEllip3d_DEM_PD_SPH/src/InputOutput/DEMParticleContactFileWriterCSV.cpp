/*
 * The MIT License
 *
 * Copyright (c) 2017-2018 Parresia Research Limited, New Zealand
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

#include <InputOutput/DEMParticleContactFileWriterCSV.h>
#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMContact.h>
#include <iomanip>

using namespace dem;

DEMParticleContactFileWriterCSV::DEMParticleContactFileWriterCSV(MPI_Comm world,
                                                                 const std::string& outputFileName) 
{
  std::string fileName = outputFileName + ".csv";
  int err = MPI_File_open(world, fileName.c_str(), 
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &d_contactFile);
  MPI_Comm_rank(world, &d_rank);
  if (err) {
    char err_string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(err, err_string, &resultlen);
    std::cerr << "Process: " << d_rank 
              << " Could not open file " << fileName
              << " for writing particle contact information." << std::endl
              << " : MPI Error: = " << err_string << std::endl;
    throw false;            
  }
}

DEMParticleContactFileWriterCSV::~DEMParticleContactFileWriterCSV() 
{
  MPI_File_close(&d_contactFile);
}

/**
 * Write the particle contacts in CSV format (contacts are per patch)
 */
void
DEMParticleContactFileWriterCSV::write(const DEMContactArray& contacts)
{
  std::stringstream dataStream;
  dataStream.setf(std::ios::scientific, std::ios::floatfield);
  for (const auto& contact : contacts) {
    contact.write(dataStream);
    std::cout << "Rank = " << d_rank 
              << " data = " << dataStream.str() << "\n";
  }

  //int length = (OWID * 28 + 1) * d_contacts.size();
  int length = dataStream.str().length();

  // write a file at a location specified by a shared file pointer (blocking,
  // collective) // note MPI_File_write_shared is non-collective
  MPI_Status status;
  int err = MPI_File_write_ordered(d_contactFile, 
                                   const_cast<char*>(dataStream.str().c_str()),
                                   length, MPI_CHAR, &status);
  if (err) {
    char err_string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(err, err_string, &resultlen);
    std::cerr << "Rank: " << d_rank << " err = " << err_string << "\n";
  }
}
