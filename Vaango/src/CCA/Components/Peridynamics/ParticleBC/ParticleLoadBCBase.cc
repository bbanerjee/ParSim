/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCBase.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Vaango;

// Store the geometry object and the load curve
ParticleLoadBCBase::ParticleLoadBCBase(Uintah::ProblemSpecP& ps)
{
  // First read the geometry information
  // d_surface is the geometry object containing the surface to be loaded.
  // **WARNING** Currently allows only for box, cylinder or sphere.
  Uintah::ProblemSpecP parent = ps->findBlock("geom_object");
  Uintah::ProblemSpecP child  = parent->findBlock();
  std::string go_type         = child->getNodeName();
  // std::cerr << "ParticleLoadBCBase::go_type = " << go_type << std::endl;

  if (go_type == "box") {
    d_surface     = scinew Uintah::BoxGeometryPiece(child);
    d_surfaceType = "box";
  } else if (go_type == "sphere") {
    d_surface     = scinew Uintah::SphereGeometryPiece(child);
    d_surfaceType = "sphere";
  } else if (go_type == "cylinder") {
    d_surface     = scinew Uintah::CylinderGeometryPiece(child);
    d_surfaceType = "cylinder";
  } else {
    throw Uintah::ParameterNotFound(
      "** ERROR ** No surface specified for pressure BC.", __FILE__, __LINE__);
  }

  d_numParticlesOnLoadSurface = 0; // this value is read in on a restart
  ps->get("numberOfParticlesOnLoadSurface", d_numParticlesOnLoadSurface);
}

// Destroy the pressure BCs
ParticleLoadBCBase::~ParticleLoadBCBase()
{
  delete d_surface;
}

// Locate and flag the particles to which this pressure BC is
// to be applied. Assumes that the "checkForSurface" function in
// ParticleCreator.cc has been used to identify this particle as being on the
// surface of the body. WARNING : For this logic to work, the surface object
// should be a box (zero volume), cylinder, sphere geometry piece that touches
// or contains the surface on which the pressure is to be applied.
bool
ParticleLoadBCBase::flagSurfaceParticle(const Uintah::Point& p,
                                        const Uintah::Vector& dxpp)
{
  bool flag = false;
  if (d_surfaceType == "box") {

    // Create box that is min-dxpp, max+dxpp;
    Uintah::Box box = d_surface->getBoundingBox();
    Uintah::GeometryPiece* volume =
      scinew Uintah::BoxGeometryPiece(box.lower() - dxpp, box.upper() + dxpp);

    if (volume->inside(p)) {
      flag = true;
    }
    delete volume;

  } else if (d_surfaceType == "cylinder") {

    double tol = 0.9 * dxpp.minComponent();
    Uintah::CylinderGeometryPiece* cgp =
      dynamic_cast<Uintah::CylinderGeometryPiece*>(d_surface);

    // Create a cylindrical annulus with radius-|dxpp|, radius+|dxpp|
    std::shared_ptr<Uintah::GeometryPiece> outer =
      std::make_shared<Uintah::CylinderGeometryPiece>(
        cgp->top(), cgp->bottom(), cgp->radius() + tol);
    std::shared_ptr<Uintah::GeometryPiece> inner =
      std::make_shared<Uintah::CylinderGeometryPiece>(
        cgp->top(), cgp->bottom(), cgp->radius() - tol);

    Uintah::GeometryPiece* volume =
      scinew Uintah::DifferenceGeometryPiece(outer, inner);
    if (volume->inside(p)) {
      flag = true;
    }
    delete volume;

  } else if (d_surfaceType == "sphere") {
    // Create a spherical shell with radius-|dxpp|, radius+|dxpp|
    double tol = dxpp.length();
    Uintah::SphereGeometryPiece* sgp =
      dynamic_cast<Uintah::SphereGeometryPiece*>(d_surface);

    std::shared_ptr<Uintah::GeometryPiece> outer =
      std::make_shared<Uintah::SphereGeometryPiece>(sgp->origin(), sgp->radius() + tol);
    std::shared_ptr<Uintah::GeometryPiece> inner =
      std::make_shared<Uintah::SphereGeometryPiece>(sgp->origin(), sgp->radius() - tol);
    Uintah::GeometryPiece* volume =
      scinew Uintah::DifferenceGeometryPiece(outer, inner);
    if (volume->inside(p)) {
      flag = true;
    }

    delete volume;

  } else {
    throw Uintah::ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }

  return flag;
}

// Calculate the area of the surface on which the pressure or traction BC
// is applied
double
ParticleLoadBCBase::getSurfaceArea() const
{
  double area = 0.0;
  if (d_surfaceType == "box") {
    Uintah::BoxGeometryPiece* gp =
      dynamic_cast<Uintah::BoxGeometryPiece*>(d_surface);
    area = gp->volume() / gp->smallestSide();
  } else if (d_surfaceType == "cylinder") {
    Uintah::CylinderGeometryPiece* gp =
      dynamic_cast<Uintah::CylinderGeometryPiece*>(d_surface);
    area = gp->surfaceArea();
  } else if (d_surfaceType == "sphere") {
    Uintah::SphereGeometryPiece* gp =
      dynamic_cast<Uintah::SphereGeometryPiece*>(d_surface);
    area = gp->surfaceArea();
  } else {
    throw Uintah::ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }
  return area;
}
