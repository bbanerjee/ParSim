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

#ifndef __DEM_OUTPUT_TECPLOT_H__
#define __DEM_OUTPUT_TECPLOT_H__

#include <Core/Geometry/Box.h>
#include <DiscreteElements/DEMContainers.h>
#include <Peridynamics/PeriContainers.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/Output.h>

#include <iostream>
#include <sstream>

namespace dem {

template <typename TArray>
class OutputTecplot : public Output
{
public:
  OutputTecplot() = delete;
  OutputTecplot(const OutputTecplot&) = delete;

  OutputTecplot(const std::string& folderName, int iterInterval);
  virtual ~OutputTecplot();

  void write(int frame, REAL time = 0.0);

  void setDomain(const Box* domain) { d_domain = domain; }
  void setPatchBox(const Box* patchBox) { d_patchBox = patchBox; }
  void setParticles(const TArray* particles)
  {
    d_particles = particles;
  }

  void writeDomain(const Box* domain, REAL time = 0.0);
  void writePatchBoxGrid(const Box* patchBox, REAL time = 0.0);
  void writeParticles(const TArray* particles, int frame, REAL time = 0.0);
  void writeSieves(const Gradation* gradation);

private:

  const Box* d_domain;
  const Box* d_patchBox;
  const TArray* d_particles;

}; // end class

} // end namespace

#endif
