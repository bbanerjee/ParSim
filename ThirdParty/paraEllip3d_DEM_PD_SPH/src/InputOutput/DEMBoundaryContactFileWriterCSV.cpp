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

#include <InputOutput/DEMBoundaryContactFileWriterCSV.h>
#include <Boundary/Boundary.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Const/Constants.h>
#include <iomanip>

using namespace dem;

DEMBoundaryContactFileWriterCSV::DEMBoundaryContactFileWriterCSV(const std::string& outputFileName) 
    : d_outputStream(outputFileName+".csv")
{
  if (!d_outputStream) {
    debugInf << "Could not create output stream for " 
              << outputFileName << std::endl;
    exit(-1);
  }
  d_outputStream.setf(std::ios::scientific, std::ios::floatfield);
  d_outputStream.precision(OPREC);
}

DEMBoundaryContactFileWriterCSV::~DEMBoundaryContactFileWriterCSV() 
{
  d_outputStream.close();
}

/**
 * Write the boundary contacts in CSV format
 */
void
DEMBoundaryContactFileWriterCSV::write(const BoundaryPArray& boundaries)
{
  for (const auto& boundary : boundaries) {

    int boundaryID = static_cast<int>(boundary->getId());
    const BoundaryContactArray& contacts = boundary->getBoundaryContacts();

    d_outputStream << std::setw(OWID) << boundaryID << std::endl;
    d_outputStream << std::setw(OWID) << contacts.size() << std::endl;
    d_outputStream << std::setw(OWID) << "pos_x" << std::setw(OWID) << "pos_y"
      << std::setw(OWID) << "pos_z" << std::setw(OWID) << "normal_x"
      << std::setw(OWID) << "normal_y" << std::setw(OWID) << "normal_z"
      << std::setw(OWID) << "tangt_x" << std::setw(OWID) << "tangt_y"
      << std::setw(OWID) << "tangt_z" << std::setw(OWID) << "pentr" << std::endl;

    for (const auto& contact : contacts) {
      contact.write(d_outputStream);
    }
  }
}
