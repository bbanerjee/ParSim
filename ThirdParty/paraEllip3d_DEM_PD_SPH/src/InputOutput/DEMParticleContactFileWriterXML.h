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

#ifndef ELLIP3D_DEM_PARTICLE_CONTACT_FILE_WRITER_XML_H
#define ELLIP3D_DEM_PARTICLE_CONTACT_FILE_WRITER_XML_H

#include <DiscreteElements/DEMContainers.h>
#include <Core/Const/Constants.h>
#include <InputOutput/zenxml/xml.h>

namespace dem {

class DEMParticleContactFileWriterXML
{

public:

  DEMParticleContactFileWriterXML(MPI_Comm world,
                                  const std::string& outputFileName);
  ~DEMParticleContactFileWriterXML(); 

  void write(const DEMContactArray& contacts);

private:

  MPI_File d_contactFile;

  zen::XmlDoc d_doc;
  std::string d_outputFileName;

  DEMParticleContactFileWriterXML() = delete;
  DEMParticleContactFileWriterXML(DEMParticleContactFileWriterXML const&) = delete;
  void operator=(DEMParticleContactFileWriterXML const&) = delete;
};
}

#endif
