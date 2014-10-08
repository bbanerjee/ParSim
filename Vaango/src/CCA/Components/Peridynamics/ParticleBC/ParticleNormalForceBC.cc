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

#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>

using namespace Vaango;

ParticleNormalForceBC::ParticleNormalForceBC(Uintah::ProblemSpecP& ps)
  : ParticleLoadBCBase(ps)
{
  // Read and save the load curve information
  d_loadCurve = scinew ParticleLoadCurve<double>(ps);
}

ParticleNormalForceBC::~ParticleNormalForceBC()
{
  delete d_loadCurve;
}


void 
ParticleNormalForceBC::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP press_ps = ps->appendChild("normal_force");
  Uintah::ProblemSpecP geom_ps = press_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  press_ps->appendElement("numberOfParticlesOnLoadSurface",d_numParticlesOnLoadSurface);
  d_loadCurve->outputProblemSpec(press_ps);
}

// Get the type of this object for BC application
std::string 
ParticleNormalForceBC::getType() const
{
  return "NormalForce";
}

// Calculate the force vector to be applied to a particular
// particle location
SCIRun::Vector
ParticleNormalForceBC::getForceVector(const SCIRun::Point& px, double force,
                                      const double time) const
{
  SCIRun::Vector normal_force(0.0,0.0,0.0);
  if (d_surfaceType == "box") {
    Uintah::BoxGeometryPiece* gp = dynamic_cast<Uintah::BoxGeometryPiece*>(d_surface);
    SCIRun::Vector normal(0.0, 0.0, 0.0);
    normal[gp->thicknessDirection()] = 1.0;
    normal_force = normal*force;
  } else if (d_surfaceType == "cylinder") {
    Uintah::CylinderGeometryPiece* gp = dynamic_cast<Uintah::CylinderGeometryPiece*>(d_surface);
    SCIRun::Vector normal = gp->radialDirection(px);
    normal_force = normal*force;
  } else if (d_surfaceType == "sphere") {
    Uintah::SphereGeometryPiece* gp = dynamic_cast<Uintah::SphereGeometryPiece*>(d_surface);
    SCIRun::Vector normal = gp->radialDirection(px);
    normal_force = normal*force;
  } else {
    throw Uintah::ParameterNotFound("ERROR: Unknown surface specified for pressure BC",
                            __FILE__, __LINE__);
  }
  return normal_force;
}

