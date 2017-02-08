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
#include <DiscreteElements/Containers.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/Output.h>

#include <iostream>
#include <sstream>

namespace dem {

class OutputTecplot : public Output
{
public:
  OutputTecplot() = delete;
  OutputTecplot(const OutputTecplot&) = delete;

  OutputTecplot(const std::string& fileName, int iterInterval);
  virtual ~OutputTecplot();

  void write();

  void setDomain(const Box* domain) { d_domain = domain; }
  void setGrid(const Box* grid) { d_grid = grid; }
  void setParticles(const ParticlePArray* particles)
  {
    d_particles = particles;
  }

  void getFileNames(std::ostringstream& domainFileName,
                    std::ostringstream& nodeFileName,
                    std::ostringstream& particleFileName);

  void writeDomain(const Box* domain, const std::string& fileName);

  void writeGrid(const Box* grid, const std::string& fileName);

  void writeParticles(const ParticlePArray* particles,
                      const std::string& fileName);

  void writeSieves(const Gradation* gradation, const std::string& fileName);

private:
  std::ostringstream d_output_dir;

  const Box* d_domain;
  const Box* d_grid;
  const ParticlePArray* d_particles;

}; // end class

} // end namespace

#endif
