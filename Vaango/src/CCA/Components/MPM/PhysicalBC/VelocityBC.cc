/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/BBox.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <iostream>

using namespace Uintah;


// Store the geometry object and the load curve
VelocityBC::VelocityBC(ProblemSpecP& ps,
                       const GridP& grid,
                       [[maybe_unused]] const MPMFlags* flags)
{
  // First read the geometry information
  // d_surface is the geometry object containing the surface to be loaded.
  // **WARNING** Currently allows only for box, cylinder or sphere.
  d_dxpp = Vector(1.0, 1.0, 1.0); // Only needed for axisymmetric end, see below
  ProblemSpecP parent = ps->findBlock("geom_object");
  ProblemSpecP child  = parent->findBlock();
  std::string go_type = child->getNodeName();

  // std::cerr << "VelocityBC::go_type = " << go_type << std::endl;

  if (go_type == "box") {
    d_surface     = scinew BoxGeometryPiece(child);
    d_surfaceType = "box";
  } else if (go_type == "sphere") {
    d_surface     = scinew SphereGeometryPiece(child);
    d_surfaceType = "sphere";
  } else if (go_type == "cylinder") {
    d_surface     = scinew CylinderGeometryPiece(child);
    d_surfaceType = "cylinder";
    CylinderGeometryPiece* cgp =
      dynamic_cast<CylinderGeometryPiece*>(d_surface);
    d_cylinder_end      = cgp->cylinder_end();
    d_axisymmetric_end  = cgp->axisymmetric_end();
    d_axisymmetric_side = cgp->axisymmetric_side();
    if (d_axisymmetric_end) {
      ps->require("res", d_res);
      Vector dx = grid->getLevel(0)->dCell();
      d_dxpp    = Vector(dx.x() / ((double)d_res.x()),
                      dx.y() / ((double)d_res.y()),
                      dx.z() / ((double)d_res.z()));
    }
  } else {
    throw ParameterNotFound(
      "** ERROR ** No surface specified for velocity BC.", __FILE__, __LINE__);
  }

  // Read the scaling function for the load curve
  ps->getWithDefault(
    "load_curve_scaling_function", d_scaling_function_expr, "1.0");

  // Parse the expression
  d_symbol_table.add_variable("t", d_time);
  d_symbol_table.add_variable("X", d_pos_x);
  d_symbol_table.add_variable("Y", d_pos_y);
  d_symbol_table.add_variable("Z", d_pos_z);
  d_symbol_table.add_constants();
  d_expression.register_symbol_table(d_symbol_table);
  if (!d_parser.compile(d_scaling_function_expr, d_expression)) {
    std::ostringstream out;
    out << "** ERROR ** Failed to parse load_curve_scaling_function"
        << d_scaling_function_expr << ".  Parser error was " << d_parser.error()
        << "." << std::endl;
    for (std::size_t i = 0; i < d_parser.error_count(); ++i) {
      exprtk::parser_error::type error = d_parser.get_error(i);

      out << "\t Error: " << i << " Position: " << error.token.position
          << " Type: " << exprtk::parser_error::to_str(error.mode)
          << " Msg: " << error.diagnostic << std::endl;
    }
    out << "Please check your input file." << std::endl;
    throw ParameterNotFound(out.str(), __FILE__, __LINE__);
  }

  // Read and save the load curve information
  d_loadCurve = scinew LoadCurve<Vector>(ps);
}

// Destroy the velocity BCs
VelocityBC::~VelocityBC()
{
  delete d_surface;
  delete d_loadCurve;
}

void
VelocityBC::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP vel_ps  = ps->appendChild("velocity");
  ProblemSpecP geom_ps = vel_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  d_loadCurve->outputProblemSpec(vel_ps);
  vel_ps->appendElement("res", d_res);
  vel_ps->appendElement("load_curve_scaling_function", d_scaling_function_expr);
}

// Get the type of this object for BC application
std::string
VelocityBC::getType() const
{
  return "Velocity";
}

// Locate and flag the material points to which this velocity BC is
// to be applied. Assumes that the "checkForSurface" function in
// ParticleCreator.cc has been used to identify this material point as being on
// the surface of the body. WARNING : For this logic to work, the surface object
// should be a box (zero volume), cylinder, sphere geometry piece that touches
// or contains the surface on which the velocity is to be applied.
bool
VelocityBC::flagMaterialPoint(const Point& p, const Vector& dxpp)
{
  bool flag = false;
  if (d_surfaceType == "box") {
    // Create box that is min-dxpp, max+dxpp;
    Box box               = d_surface->getBoundingBox();
    GeometryPieceP volume = std::make_shared<BoxGeometryPiece>(
      box.lower() - dxpp, box.upper() + dxpp);

    if (volume->inside(p))
      flag = true;

  } else if (d_surfaceType == "cylinder") {
    double tol = 0.9 * dxpp.minComponent();
    CylinderGeometryPiece* cgp =
      dynamic_cast<CylinderGeometryPiece*>(d_surface);

    if (!d_cylinder_end && !d_axisymmetric_end) { // Not a cylinder end
      // Create a cylindrical annulus with radius-|dxpp|, radius+|dxpp|
      GeometryPieceP outer = std::make_shared<CylinderGeometryPiece>(
        cgp->top(), cgp->bottom(), cgp->radius() + tol);
      GeometryPieceP inner = std::make_shared<CylinderGeometryPiece>(
        cgp->top(), cgp->bottom(), cgp->radius() - tol);

      GeometryPieceP volume =
        std::make_shared<DifferenceGeometryPiece>(outer, inner);
      if (volume->inside(p)) {
        flag = true;
      }
    } else if (d_cylinder_end || d_axisymmetric_end) {
      Vector add_ends = tol * (cgp->top() - cgp->bottom()) /
                        (cgp->top() - cgp->bottom()).length();

      GeometryPieceP end = std::make_shared<CylinderGeometryPiece>(
        cgp->top() + add_ends, cgp->bottom() - add_ends, cgp->radius());
      if (end->inside(p)) {
        flag = true;
      }
    }
  } else if (d_surfaceType == "sphere") {
    // Create a spherical shell with radius-|dxpp|, radius+|dxpp|
    double tol               = dxpp.length();
    SphereGeometryPiece* sgp = dynamic_cast<SphereGeometryPiece*>(d_surface);
    GeometryPieceP outer =
      std::make_shared<SphereGeometryPiece>(sgp->origin(), sgp->radius() + tol);
    GeometryPieceP inner =
      std::make_shared<SphereGeometryPiece>(sgp->origin(), sgp->radius() - tol);
    GeometryPieceP volume =
      std::make_shared<DifferenceGeometryPiece>(outer, inner);
    if (volume->inside(p))
      flag = true;

  } else {
    throw ParameterNotFound(
      "ERROR: Unknown surface specified for velocity BC", __FILE__, __LINE__);
  }

  return flag;
}

Vector
VelocityBC::velocity(double t, const Point& pX)
{
  Vector velocity = d_loadCurve->getLoad(t);

  d_time              = t;
  d_pos_x             = pX.x();
  d_pos_y             = pX.y();
  d_pos_z             = pX.z();
  double scale_factor = d_expression.value();
  // std::cout << " t = " << t << " x = " << d_pos_x << " scale_factor = " <<
  // scale_factor << " V = "
  //           << velocity*scale_factor << std::endl;

  return (velocity * scale_factor);
}

// Calculate the velocity vector to be applied to a particular
// material point location
Vector
VelocityBC::getVelocityVector(const Point& px,
                              const Vector& pDisp,
                              const double time)
{
  Point pX = px - pDisp;
  // std::cout << "t = " << time << " x = " << px << " X = " << pX << " U = " <<
  // pDisp << std::endl;
  Vector vel = velocity(time, pX);
  return vel;
}

// Update the load curve
void
VelocityBC::updateLoadCurve(const std::vector<double>& time,
                            const std::vector<Vector>& velocity)
{
  d_loadCurve->setTimeLoad(time, velocity);
}

namespace Uintah {
// A method to print out the velocity bcs
ostream&
operator<<(std::ostream& out, const VelocityBC& bc)
{
  out << "Begin MPM Velocity BC # = " << bc.loadCurveID() << std::endl;
  std::string surfType = bc.getSurfaceType();
  out << "    Surface of application = " << surfType << std::endl;
  if (surfType == "box") {
    Box box = (bc.getSurface())->getBoundingBox();
    out << "        " << box << std::endl;
  } else if (surfType == "cylinder") {
    CylinderGeometryPiece* cgp =
      dynamic_cast<CylinderGeometryPiece*>(bc.getSurface());
    out << "        "
        << "radius = " << cgp->radius() << " top = " << cgp->top()
        << " bottom = " << cgp->bottom() << std::endl;
  } else if (surfType == "sphere") {
    SphereGeometryPiece* sgp =
      dynamic_cast<SphereGeometryPiece*>(bc.getSurface());
    out << "        "
        << "radius = " << sgp->radius() << " origin = " << sgp->origin()
        << std::endl;
  }
  out << "    Time vs. Velocity = " << std::endl;
  LoadCurve<Vector>* lc = bc.getLoadCurve();
  int numPts            = lc->numberOfPointsOnLoadCurve();
  for (int ii = 0; ii < numPts; ++ii) {
    out << "        time = " << lc->getTime(ii)
        << " velocity = " << lc->getLoad(ii) << std::endl;
  }
  out << "End MPM Velocity BC # = " << bc.loadCurveID() << std::endl;
  return out;
}

} // end namespace Uintah
