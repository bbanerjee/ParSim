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

#ifndef ELLIP3D_DEM_PARTICLE_FILE_WRITER_H
#define ELLIP3D_DEM_PARTICLE_FILE_WRITER_H

#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/Gradation.h>

#include <Core/Types/RealTypes.h>

#include <InputOutput/zenxml/xml.h>

namespace dem {

class DEMParticleFileWriter
{

public:
  DEMParticleFileWriter() = default;
  ~DEMParticleFileWriter() = default;

  void writeCSV(const DEMParticlePArray& particles,
                const Gradation& gradation,
                const std::string& outputFileName) const;

  void writeXML(const DEMParticlePArray& particles,
                const Gradation& gradation,
                const std::string& outputFileName) const;

private:

  template <typename T, int numComponents>
  void writeParticleValues(zen::XmlElement& element, 
                           const std::string& name,
                           const std::string& unit,
                           std::vector<T>& particleData) const;

  DEMParticleFileWriter(DEMParticleFileWriter const&) = delete;
  void operator=(DEMParticleFileWriter const&) = delete;
};
}
#endif
