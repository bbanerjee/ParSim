/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
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
PressureBC::PressureBC(ProblemSpecP& ps,
                       const GridP& grid,
                       const MPMFlags* flags)
{
  // First read the geometry information
  // d_surface is the geometry object containing the surface to be loaded.
  // The sign of the pressure load is +ve if applied in the direction
  // of the outward normal and -ve if applied in the direction of the
  // inward normal
  // **WARNING** Currently allows only for box, cylinder or sphere.
  if (flags->d_useCBDI) {
    ps->require("outward_normal", d_outwardNormal);
  }
  d_dxpp = Vector(1., 1., 1.); // Only needed for axisymmetric end, see below
  ProblemSpecP adult  = ps->findBlock("geom_object");
  ProblemSpecP child  = adult->findBlock();
  std::string go_type = child->getNodeName();
  // std::cerr << "PressureBC::go_type = " << go_type << endl;
  if (go_type == "box") {
    d_surface = scinew BoxGeometryPiece(child);
    // Box box = d_surface->getBoundingBox();
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
      "** ERROR ** No surface specified for pressure BC.", __FILE__, __LINE__);
  }

  d_volFracInsideDomain = 1.0;
  ps->get("volume_fraction_inside_domain", d_volFracInsideDomain);

  d_numMaterialPoints = 0; // this value is read in on a restart
  ps->get("numberOfParticlesOnLoadSurface", d_numMaterialPoints);

  // Read the scaling function for the load curve
  ps->getWithDefault(
    "load_curve_scaling_function", d_scaling_function_expr, "1.0");

  // Parse the expression
  d_time  = 0.0;
  d_pos_x = 0.0;
  d_pos_y = 0.0;
  d_pos_z = 0.0;
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
  d_loadCurve = scinew LoadCurve<double>(ps);

  //__________________________________
  //   Bulletproofing
  // user shouldn't specify a geometry object that is bigger than the domain
  Box boundingBox = d_surface->getBoundingBox();
  BBox compDomain;
  grid->getSpatialRange(compDomain);

  Point BB_min = boundingBox.lower();
  Point CD_min = compDomain.min();
  Point BB_max = boundingBox.upper();
  Point CD_max = compDomain.max();

  if ((BB_min.x() < CD_min.x()) || (BB_min.y() < CD_min.y()) ||
      (BB_min.z() < CD_min.z()) || (BB_max.x() > CD_max.x()) ||
      (BB_max.y() > CD_max.y()) || (BB_max.z() > CD_max.z())) {
    if (!d_axisymmetric_end && !d_axisymmetric_side) {
      proc0cout
        << "_________________________________________________________\n";
      proc0cout << "\n Input File WARNING: <PhysicalBC : MPM : Pressure> \n"
                << " The geometry Object [" << d_surface->getType()
                << "] exceeds the dimensions of the computational domain.\n"
                << " \n Please change the parameters so it doesn't. \n\n"
                << " There is a flaw in the surface area calculation for the "
                   "geometry object,\n"
                << " it does not take into account that the object exceeds the "
                   "domain\n";
      proc0cout
        << "_________________________________________________________\n";
    }
  }
}

// Destroy the pressure BCs
PressureBC::~PressureBC()
{
  delete d_surface;
  delete d_loadCurve;
}

void
PressureBC::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP press_ps = ps->appendChild("pressure");
  ProblemSpecP geom_ps  = press_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  press_ps->appendElement("volume_fraction_inside_domain",
                          d_volFracInsideDomain);
  press_ps->appendElement("numberOfParticlesOnLoadSurface",
                          d_numMaterialPoints);
  d_loadCurve->outputProblemSpec(press_ps);
  press_ps->appendElement("res", d_res);
  press_ps->appendElement("load_curve_scaling_function",
                          d_scaling_function_expr);
}

// Get the type of this object for BC application
std::string
PressureBC::getType() const
{
  return "Pressure";
}

// Locate and flag the material points to which this pressure BC is
// to be applied. Assumes that the "checkForSurface" function in
// ParticleCreator.cc has been used to identify this material point as being on
// the surface of the body. WARNING : For this logic to work, the surface object
// should be a box (zero volume), cylinder, sphere geometry piece that touches
// or contains the surface on which the pressure is to be applied.
bool
PressureBC::flagMaterialPoint(const Point& p, const Vector& dxpp)
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
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }

  return flag;
}

// Calculate the area of the surface on which the pressure BC
// is applied
double
PressureBC::getSurfaceArea() const
{
  double area = 0.0;
  if (d_surfaceType == "box") {

    BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
    area                 = gp->volume() / gp->smallestSide();
    // std::cout << "Surface area = " << area << std::endl;

  } else if (d_surfaceType == "cylinder") {
    CylinderGeometryPiece* gp = dynamic_cast<CylinderGeometryPiece*>(d_surface);
    if (!d_cylinder_end && !d_axisymmetric_end) { // Not a cylinder end
      area = gp->surfaceArea();
      if (d_axisymmetric_side) {
        area /= (2.0 * M_PI);
      }
    } else if (d_cylinder_end) {
      area = gp->surfaceAreaEndCaps() / 2.0;
    } else if (d_axisymmetric_end) {
      area = (gp->radius() * gp->radius()) / 2.0; // area of a 1 radian wedge
    }
  } else if (d_surfaceType == "sphere") {
    SphereGeometryPiece* gp = dynamic_cast<SphereGeometryPiece*>(d_surface);
    area                    = gp->surfaceArea();
  } else {
    throw ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }
  return area;
}

// Calculate the force per particle at a certain time
double
PressureBC::forcePerParticle(double time)
{
  if (d_numMaterialPoints < 1)
    return 0.0;

  // Get the area of the surface on which the pressure BC is applied
  double area = getSurfaceArea();
  area *= d_volFracInsideDomain;

  // Get the initial pressure that is applied ( t = 0.0 )
  double press = pressure(time);

  // Calculate the force per particle
  double force = (press * area) / static_cast<double>(d_numMaterialPoints);
  // std::cout << "Force = " << force << " pressure = " << press << " area = "
  // << area
  //           << " num. points = " << d_numMaterialPoints << std::endl;
  return force;
}

double
PressureBC::forcePerParticle(double time, const Point& pX)
{
  if (d_numMaterialPoints < 1)
    return 0.0;

  // Get the area of the surface on which the pressure BC is applied
  double area = getSurfaceArea();
  area *= d_volFracInsideDomain;

  // Get the initial pressure that is applied ( t = 0.0 )
  double press = pressure(time, pX);

  // Calculate the force per particle
  return (press * area) / static_cast<double>(d_numMaterialPoints);
}

double
PressureBC::pressure(double t)
{
  double load = d_loadCurve->getLoad(t);

  d_time  = t;
  d_pos_x = 0.0;
  d_pos_y = 0.0;
  d_pos_z = 0.0;

  double scale_factor = d_expression.value();
  // std::cout << "scale_factor = " << scale_factor << std::endl;

  return (load * scale_factor);
}

double
PressureBC::pressure(double t, const Point& pX)
{
  double load = d_loadCurve->getLoad(t);

  d_time  = t;
  d_pos_x = pX.x();
  d_pos_y = pX.y();
  d_pos_z = pX.z();

  double scale_factor = d_expression.value();
  // std::cout << "scale_factor = " << scale_factor << std::endl;

  return (load * scale_factor);
}

// Calculate the force vector to be applied to a particular
// material point location
Vector
PressureBC::getForceVector(const Point& px,
                           const Vector& pDisp,
                           double forcePerParticle,
                           const double time,
                           const Matrix3& defGrad)
{
  // Get reference position
  Point pX = px - pDisp;

  // Compute F^{-T} for Nanson's relation
  double JJ     = defGrad.Determinant();
  Matrix3 FinvT = (defGrad.Inverse()).Transpose();

  Vector force(0.0, 0.0, 0.0);
  if (d_surfaceType == "box") {

    BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
    Vector normalRef(0.0, 0.0, 0.0);
    normalRef[gp->thicknessDirection()] = 1.0;
    Vector scaledNormalCur              = (FinvT * normalRef) * JJ;
    force                               = scaledNormalCur * forcePerParticle;

    // std::cout << " Force vector = " << force << std::endl;
    //          << " Normal vector = " << scaledNormalCur << std::endl
    //          << " DefGrad = " << defGrad << std::endl
    //          << " J = " << JJ << " Area ratio = " << scaledNormalCur.length()
    //          << std::endl;

  } else if (d_surfaceType == "cylinder") {

    CylinderGeometryPiece* gp = dynamic_cast<CylinderGeometryPiece*>(d_surface);
    Vector normalRef          = gp->radialDirection(px);
    Vector scaledNormalCur    = (FinvT * normalRef) * JJ;
    force                     = scaledNormalCur * forcePerParticle;

    if (d_cylinder_end || d_axisymmetric_end) {

      normalRef =
        (gp->top() - gp->bottom()) / (gp->top() - gp->bottom()).length();
      scaledNormalCur = (FinvT * normalRef) * JJ;

      if (!d_axisymmetric_end) {

        force = scaledNormalCur * forcePerParticle;

      } else { // It IS on an axisymmetric end

        double pArea = px.x() * d_dxpp.x() * 1.0; /*(theta = 1 radian)*/
        double press = pressure(time, pX);
        double fpP   = pArea * press;
        force        = scaledNormalCur * fpP;
      }
    }

  } else if (d_surfaceType == "sphere") {

    SphereGeometryPiece* gp = dynamic_cast<SphereGeometryPiece*>(d_surface);
    Vector normalRef        = gp->radialDirection(px);
    Vector scaledNormalCur  = (FinvT * normalRef) * JJ;
    force                   = scaledNormalCur * forcePerParticle;

  } else {
    throw ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }
  return force;
}

// Calculate the force vector to be applied to a particular
// material point location
Vector
PressureBC::getForceVectorCBDI(const Point& px,
                               const Vector& pDisp,
                               const Matrix3& psize,
                               const Matrix3& pDeformationMeasure,
                               double forcePerParticle,
                               const double time,
                               Point& pExternalForceCorner1,
                               Point& pExternalForceCorner2,
                               Point& pExternalForceCorner3,
                               Point& pExternalForceCorner4,
                               const Vector& dxCell)
{
  // Get reference position
  Point pX = px - pDisp;

  // Compute force vector
  Vector force(0.0, 0.0, 0.0);
  Vector normal(0.0, 0.0, 0.0);
  if (d_surfaceType == "box") {
    BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
    normal[gp->thicknessDirection()] = 1.0;
    force                            = normal * forcePerParticle;
  } else if (d_surfaceType == "cylinder") {
    CylinderGeometryPiece* gp = dynamic_cast<CylinderGeometryPiece*>(d_surface);
    normal                    = gp->radialDirection(px);
    force                     = normal * forcePerParticle;
    // std::cout << "px = " << px << " force = " << force << " normal = " <<
    // normal
    //           << " fpp = " << forcePerParticle << "\n";
    if (d_cylinder_end || d_axisymmetric_end) {
      normal = (gp->top() - gp->bottom()) / (gp->top() - gp->bottom()).length();
      if (!d_axisymmetric_end) {
        force = normal * forcePerParticle;
      } else { // It IS on an axisymmetric end
        double pArea = px.x() * d_dxpp.x() * 1.0; /*(theta = 1 radian)*/
        double press = pressure(time, pX);
        double fpP   = pArea * press;
        force        = normal * fpP;
      }
    }
  } else if (d_surfaceType == "sphere") {
    SphereGeometryPiece* gp = dynamic_cast<SphereGeometryPiece*>(d_surface);
    normal                  = gp->radialDirection(px);
    force                   = normal * forcePerParticle;
  } else {
    throw ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }
  // 25% of total particle force goes to each corner
  force = force * 0.25;
  // modify the sign of force if outward normal is not correctly defined
  if (!d_outwardNormal) {
    force = force * (-1.0);
  }
  // determine four boundary-corners of the particle
  std::vector<Point> corners(4);
  Matrix3 pSize_new = pDeformationMeasure * psize;
  auto i1i2 =
    getParticleBoundaryCorners(px, pSize_new, normal, dxCell, corners);
  pExternalForceCorner1 = corners[0];
  pExternalForceCorner2 = corners[1];
  pExternalForceCorner3 = corners[2];
  pExternalForceCorner4 = corners[3];

  // Recalculate the force based on area changes (current vs. initial)
  int i1 = i1i2.first;
  int i2 = i1i2.second;
  if (i1 != i2) {
    Vector iniVec1(psize(0, i1), psize(1, i1), psize(2, i1));
    Vector iniVec2(psize(0, i2), psize(1, i2), psize(2, i2));
    Vector curVec1(pSize_new(0, i1), pSize_new(1, i1), pSize_new(2, i1));
    Vector curVec2(pSize_new(0, i2), pSize_new(1, i2), pSize_new(2, i2));
    Vector iniA    = Cross(iniVec1, iniVec2);
    Vector curA    = Cross(curVec1, curVec2);
    double iniArea = iniA.length();
    double curArea = curA.length();
    force          = force * (curArea / iniArea);
  }
  return force;
}

// Get the particle corners for CBDI
std::vector<Point>
PressureBC::getParticleCornersCBDI(const Point& pX,
                                   const Matrix3& pSize,
                                   const Matrix3& pDefGrad,
                                   const Vector& dxCell)
{
  Vector pNormal(0.0, 0.0, 0.0);
  if (d_surfaceType == "box") {
    BoxGeometryPiece* gp = dynamic_cast<BoxGeometryPiece*>(d_surface);
    pNormal[gp->thicknessDirection()] = 1.0;
  } else if (d_surfaceType == "cylinder") {
    CylinderGeometryPiece* gp = dynamic_cast<CylinderGeometryPiece*>(d_surface);
    pNormal                   = gp->radialDirection(pX);
  } else if (d_surfaceType == "sphere") {
    SphereGeometryPiece* gp = dynamic_cast<SphereGeometryPiece*>(d_surface);
    pNormal                 = gp->radialDirection(pX);
  } else {
    throw ParameterNotFound(
      "ERROR: Unknown surface specified for pressure BC", __FILE__, __LINE__);
  }
  // determine four boundary-corners of the particle
  std::vector<Point> corners(4);
  Matrix3 pSize_new = pDefGrad * pSize;
  getParticleBoundaryCorners(pX, pSize_new, pNormal, dxCell, corners);

  return corners;
}

// Determine four boundary-corners of the particle
std::pair<int, int>
PressureBC::getParticleBoundaryCorners(const Point& pX,
                                       const Matrix3& pSize_new,
                                       const Vector& pNormal,
                                       const Vector& dxCell,
                                       std::vector<Point>& corners) const
{
  int i1 = 0, i2 = 0;
  Point px1 = pX;
  for (int ii = 0; ii < 3; ++ii) {
    Vector dxDeformed(pSize_new(0, ii) * dxCell[0],
                      pSize_new(1, ii) * dxCell[1],
                      pSize_new(2, ii) * dxCell[2]);
    dxDeformed /= 2.0;

    double dxLength     = dxDeformed.length();
    double normalLength = pNormal.length();
    double cos_angle    = Dot(pNormal, dxDeformed) / (normalLength * dxLength);

    if (std::abs(cos_angle - 1.0) < 0.1) {
      px1 = pX + dxDeformed;
      i1  = (ii + 1) % 3;
      i2  = (ii + 2) % 3;
    } else if (std::abs(cos_angle + 1.0) < 0.1) {
      px1 = pX - dxDeformed;
      i1  = (ii + 1) % 3;
      i2  = (ii + 2) % 3;
    }
    // else {
    // std::cout << "cos_angle = " << cos_angle << " px1 = " << px1 << "\n";
    // std::cout << "dx = " << dxDeformed << " normal = " << pNormal << "\n";
    // }
  }
  // px1 is the position of the center of the boundary particle face that is on
  // the physical boundary.
  if (i1 == 0 && i2 == 0) {
    corners[0] = px1;
    corners[1] = px1;
    corners[2] = px1;
    corners[3] = px1;
  } else {
    corners[0] = Point(px1.x() - pSize_new(0, i1) * dxCell[0] / 2.0 -
                         pSize_new(0, i2) * dxCell[0] / 2.0,
                       px1.y() - pSize_new(1, i1) * dxCell[1] / 2.0 -
                         pSize_new(1, i2) * dxCell[1] / 2.0,
                       px1.z() - pSize_new(2, i1) * dxCell[2] / 2.0 -
                         pSize_new(2, i2) * dxCell[2] / 2.0);
    corners[1] = Point(px1.x() + pSize_new(0, i1) * dxCell[0] / 2.0 -
                         pSize_new(0, i2) * dxCell[0] / 2.0,
                       px1.y() + pSize_new(1, i1) * dxCell[1] / 2.0 -
                         pSize_new(1, i2) * dxCell[1] / 2.0,
                       px1.z() + pSize_new(2, i1) * dxCell[2] / 2.0 -
                         pSize_new(2, i2) * dxCell[2] / 2.0);
    corners[2] = Point(px1.x() - pSize_new(0, i1) * dxCell[0] / 2.0 +
                         pSize_new(0, i2) * dxCell[0] / 2.0,
                       px1.y() - pSize_new(1, i1) * dxCell[1] / 2.0 +
                         pSize_new(1, i2) * dxCell[1] / 2.0,
                       px1.z() - pSize_new(2, i1) * dxCell[2] / 2.0 +
                         pSize_new(2, i2) * dxCell[2] / 2.0);
    corners[3] = Point(px1.x() + pSize_new(0, i1) * dxCell[0] / 2.0 +
                         pSize_new(0, i2) * dxCell[0] / 2.0,
                       px1.y() + pSize_new(1, i1) * dxCell[1] / 2.0 +
                         pSize_new(1, i2) * dxCell[1] / 2.0,
                       px1.z() + pSize_new(2, i1) * dxCell[2] / 2.0 +
                         pSize_new(2, i2) * dxCell[2] / 2.0);
  }

  return std::make_pair(i1, i2);
}

// Update the load curve
void
PressureBC::updateLoadCurve(const std::vector<double>& time,
                            const std::vector<double>& pressure)
{
  d_loadCurve->setTimeLoad(time, pressure);
}

namespace Uintah {
// A method to print out the pressure bcs
ostream&
operator<<(std::ostream& out, const PressureBC& bc)
{
  out << "Begin MPM Pressure BC # = " << bc.loadCurveID() << endl;
  std::string surfType = bc.getSurfaceType();
  out << "    Surface of application = " << surfType << endl;
  if (surfType == "box") {
    Box box = (bc.getSurface())->getBoundingBox();
    out << "        " << box << endl;
  } else if (surfType == "cylinder") {
    CylinderGeometryPiece* cgp =
      dynamic_cast<CylinderGeometryPiece*>(bc.getSurface());
    out << "        "
        << "radius = " << cgp->radius() << " top = " << cgp->top()
        << " bottom = " << cgp->bottom() << endl;
  } else if (surfType == "sphere") {
    SphereGeometryPiece* sgp =
      dynamic_cast<SphereGeometryPiece*>(bc.getSurface());
    out << "        "
        << "radius = " << sgp->radius() << " origin = " << sgp->origin()
        << endl;
  }
  out << "    Time vs. Load = " << endl;
  LoadCurve<double>* lc = bc.getLoadCurve();
  int numPts            = lc->numberOfPointsOnLoadCurve();
  for (int ii = 0; ii < numPts; ++ii) {
    out << "        time = " << lc->getTime(ii)
        << " pressure = " << lc->getLoad(ii) << endl;
  }
  out << "End MPM Pressure BC # = " << bc.loadCurveID() << endl;
  return out;
}

} // end namespace Uintah
