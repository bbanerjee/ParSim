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

#include <InputOutput/DEMParticleContactFileWriterXML.h>
#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMContact.h>
#include <Core/Const/Constants.h>
#include <iomanip>

using namespace dem;

DEMParticleContactFileWriterXML::DEMParticleContactFileWriterXML(MPI_Comm world,
                                                                 const std::string& outputFileName) 
  : d_doc("Ellip3D_contact"),
    d_outputFileName(outputFileName+".xml")
{
  int err = MPI_File_open(world, d_outputFileName.c_str(), 
                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &d_contactFile);
  MPI_Comm_rank(world, &d_rank);
  if (err) {
    char err_string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(err, err_string, &resultlen);
    std::cerr << "Process: " << d_rank 
              << " Could not open file " << d_outputFileName
              << " for writing particle contact information." << std::endl
              << " : MPI Error: = " << err_string << std::endl;
    throw false;            
  }
}

DEMParticleContactFileWriterXML::~DEMParticleContactFileWriterXML() 
{
  MPI_File_close(&d_contactFile);
}

/**
 * Write the particle contacts in XML format
 */
void
DEMParticleContactFileWriterXML::write(const DEMContactArray& contacts)
{
  // Create a proxy output
  zen::XmlOut xml(d_doc);

  // Write the title
  //std::string title = "Ellip3D particle contact XML file";
  //xml["Meta"]["title"](title);

  // Write particle contacts
  for (const auto& contact : contacts) {

    contact.write(xml);
  }

  // Serialize the xml document
  std::string dataString; 
  zen::serialize(d_doc.root(), dataString, "\r\n", "  ", 0);

  //int length = (OWID * 28 + 1) * d_contacts.size();
  int length = dataString.length();

  // write a file at a location specified by a shared file pointer (blocking,
  // collective) // note MPI_File_write_shared is non-collective
  MPI_Status status;
  MPI_File_write_ordered(d_contactFile, 
                         const_cast<char*>(dataString.c_str()),
                         length, MPI_CHAR, &status);

}