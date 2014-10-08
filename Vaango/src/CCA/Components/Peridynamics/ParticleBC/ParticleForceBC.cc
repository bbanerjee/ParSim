/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GeometryPiece/GeometryPiece.h>


using namespace Vaango;

ParticleForceBC::ParticleForceBC(Uintah::ProblemSpecP& ps)
  : ParticleLoadBCBase(ps)
{
  // Read and save the load curve information
  d_loadCurve = scinew ParticleLoadCurve<SCIRun::Vector>(ps);
}

ParticleForceBC::~ParticleForceBC()
{
  delete d_loadCurve;
}

void 
ParticleForceBC::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP press_ps = ps->appendChild("force");
  Uintah::ProblemSpecP geom_ps = press_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  press_ps->appendElement("numberOfParticlesOnLoadSurface",d_numParticlesOnLoadSurface);
  d_loadCurve->outputProblemSpec(press_ps);
}

std::string 
ParticleForceBC::getType() const
{
  return "Force";
}

