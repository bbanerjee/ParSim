/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <Core/Grid/BoundaryConditions/BoundCondReader.h>

#include <Core/Grid/BoundaryConditions/AnnulusBCData.h>
#include <Core/Grid/BoundaryConditions/BCData.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BCGeomBase.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/BoundaryConditions/BoundCondBase.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/BoundaryConditions/CircleBCData.h>
#include <Core/Grid/BoundaryConditions/DifferenceBCData.h>
#include <Core/Grid/BoundaryConditions/EllipseBCData.h>
#include <Core/Grid/BoundaryConditions/RectangleBCData.h>
#include <Core/Grid/BoundaryConditions/RectangulusBCData.h>
#include <Core/Grid/BoundaryConditions/SideBCData.h>
#include <Core/Grid/BoundaryConditions/UnionBCData.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/StringUtil.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>

namespace {

// Usage: export SCI_DEBUG="BCR_dbg:+"
Uintah::Dout BCR_dbg{ "BCBCR_dbg",
                      "BoundaryCondReader",
                      "report info regarding the BC setup",
                      false };

} // namespace

namespace Uintah {

// given a set of lower or upper bounds for multiple boxes (points)
// this function checks if Point p (usually center of circle, ellipse, or
// annulus)
// is on a given face on any of the boxes.
auto
BoundCondReader::is_on_face(const int dir,
                            const Point p_in,
                            const std::vector<Point>& points) const -> bool
{
  for (const auto& pt : points) {
    if (p_in(dir) == pt(dir)) {
      return true;
    }
  }
  return false;
}

auto
BoundCondReader::isPtOnFace(const int dir,
                            const int plusMinusFaces,
                            const Point pt,
                            const std::vector<Point>& grid_LoPts,
                            const std::vector<Point>& grid_HiPts) const -> bool
{
  bool isOnFace = false;
  if (plusMinusFaces == -1) { // x-, y-, z- faces
    isOnFace = is_on_face(dir, pt, grid_LoPts);
  }

  if (plusMinusFaces == 1) { // x+, y+, z+ faces
    isOnFace = is_on_face(dir, pt, grid_HiPts);
  }
  return isOnFace;
}

void
BoundCondReader::whichPatchFace(const std::string fc,
                                Patch::FaceType& face_side,
                                int& plusMinusFaces,
                                int& p_dir) const
{
  if (fc == "x-") {
    plusMinusFaces = -1;
    p_dir          = 0;
    face_side      = Patch::xminus;
  }
  if (fc == "x+") {
    plusMinusFaces = 1;
    p_dir          = 0;
    face_side      = Patch::xplus;
  }
  if (fc == "y-") {
    plusMinusFaces = -1;
    p_dir          = 1;
    face_side      = Patch::yminus;
  }
  if (fc == "y+") {
    plusMinusFaces = 1;
    p_dir          = 1;
    face_side      = Patch::yplus;
  }
  if (fc == "z-") {
    plusMinusFaces = -1;
    p_dir          = 2;
    face_side      = Patch::zminus;
  }
  if (fc == "z+") {
    plusMinusFaces = 1;
    p_dir          = 2;
    face_side      = Patch::zplus;
  }
}

auto
BoundCondReader::createBoundaryConditionFace(ProblemSpecP& face_ps,
                                             const ProblemSpecP& grid_ps,
                                             Patch::FaceType& face_side)
  -> std::shared_ptr<BCGeomBase>
{

  // Determine the Level 0 grid high and low points, need by
  // the bullet proofing
  Point grid_LoPt(1e30, 1e30, 1e30);
  Point grid_HiPt(-1e30, -1e30, -1e30);

  std::vector<Point> grid_LoPts; // store the lower bounds of all boxes
  std::vector<Point> grid_HiPts; // store the upper bounds of all boxes

  for (ProblemSpecP level_ps = grid_ps->findBlock("Level"); level_ps != 0;
       level_ps              = level_ps->findNextBlock("Level")) {

    // find upper/lower corner
    for (ProblemSpecP box_ps = level_ps->findBlock("Box"); box_ps != 0;
         box_ps              = box_ps->findNextBlock("Box")) {
      Point lower;
      Point upper;
      box_ps->require("lower", lower);
      box_ps->require("upper", upper);
      grid_LoPts.push_back(lower);
      grid_HiPts.push_back(upper);
      grid_LoPt = Min(lower, grid_LoPt);
      grid_HiPt = Max(upper, grid_HiPt);
    }
  }

  std::map<std::string, std::string> values;
  face_ps->getAttributes(values);

  // Possible boundary condition types for a face:
  //    side (original -- entire side is one bc)
  //   Optional geometry objects on a side
  //    circle
  //    annulus
  //    rectangulus
  //    ellipse
  //    rectangle
  // This allows a user to specify variable boundary conditions on a given
  // side.  Will use the notion of a UnionBoundaryCondtion and Difference
  // BoundaryCondition.

  std::shared_ptr<BCGeomBase> bcGeom{ nullptr };

  if (values.find("side") != values.end()) {
    bcGeom = createSideBC(values, face_side);
  } else if (values.find("circle") != values.end()) {
    bcGeom = createCircleBC(values, grid_LoPts, grid_HiPts, face_side);
  } else if (values.find("annulus") != values.end()) {
    bcGeom = createAnnulusBC(values, grid_LoPts, grid_HiPts, face_side);
  } else if (values.find("rectangulus") != values.end()) {
    bcGeom = createRectangulusBC(values, grid_LoPts, grid_HiPts, face_side);
  } else if (values.find("ellipse") != values.end()) {
    bcGeom = createEllipseBC(values, grid_LoPts, grid_HiPts, face_side);
  } else if (values.find("rectangle") != values.end()) {
    bcGeom = createRectangleBC(values, grid_LoPts, grid_HiPts, face_side);
  } else {
    std::ostringstream warn;
    warn << "ERROR\n Boundary condition geometry not correctly specified "
            " Valid options (side, circle, rectangle, annulus";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // name the boundary condition object:
  std::string bcname;
  if (values.find("name") != values.end()) {
    std::string name = values["name"];
    DOUT(BCR_dbg, "Setting name to: " << name);
    bcGeom->setBCName(name);
  }

  // get the bctype - mainly used by wasatch:
  if (values.find("type") != values.end()) {
    string bndType = values["type"];
    DOUT(BCR_dbg, "Setting bc type to: " << bndType);
    bcGeom->setBndType(bndType);
  }

  if (face_ps->findBlock("ParticleBC")) {

    ProblemSpecP particleBCps = face_ps->findBlock("ParticleBC");
    ProblemSpecP pWallBC      = particleBCps->findBlock("Wall");
    ProblemSpecP pInletBC     = particleBCps->findBlock("Inlet");

    BCGeomBase::ParticleBndSpec pBndSpec;
    if (pWallBC) {
      pBndSpec.bndType = BCGeomBase::ParticleBndSpec::WALL;
      string wallType;
      pWallBC->getAttribute("walltype", wallType);

      if (wallType == "Elastic") {
        pBndSpec.wallType        = BCGeomBase::ParticleBndSpec::ELASTIC;
        pBndSpec.restitutionCoef = 1.0;
      } else if (wallType == "Inelastic") {
        pBndSpec.wallType        = BCGeomBase::ParticleBndSpec::INELASTIC;
        pBndSpec.restitutionCoef = 0.0;
      } else if (wallType == "PartiallyElastic") {
        pBndSpec.wallType = BCGeomBase::ParticleBndSpec::PARTIALLYELASTIC;
        pWallBC->get("Restitution", pBndSpec.restitutionCoef);
      }
    } else if (pInletBC) {
      pBndSpec.bndType = BCGeomBase::ParticleBndSpec::INLET;
      pInletBC->get("ParticlesPerSecond", pBndSpec.particlesPerSec);
    }
    bcGeom->setParticleBndSpec(pBndSpec);
  }

  return bcGeom;
}

auto
BoundCondReader::createSideBC(const std::map<std::string, std::string>& values,
                              Patch::FaceType& face_side) const
  -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("side");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);
  return std::make_shared<SideBCData>();
}

auto
BoundCondReader::createCircleBC(
  const std::map<std::string, std::string>& values,
  const std::vector<Point>& grid_LoPts,
  const std::vector<Point>& grid_HiPts,
  Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("circle");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);

  const std::string& origin = values.at("origin");
  const std::string& radius = values.at("radius");

  std::stringstream origin_stream(origin);
  std::stringstream radius_stream(radius);
  if (!radius_stream || !origin_stream) {
    std::cout << "WARNING: BoundCondReader.cc:  std::stringstream failed..."
              << std::endl;
  }

  double r, o[3];
  radius_stream >> r;
  origin_stream >> o[0] >> o[1] >> o[2];

  Point p(o[0], o[1], o[2]);

  //  bullet proofing-- origin must be on the same plane as the face
  bool isOnFace = false;
  if (plusMinusFaces == -1) { // x-, y-, z- faces
    isOnFace = is_on_face(p_dir, p, grid_LoPts);
  }

  if (plusMinusFaces == 1) { // x+, y+, z+ faces
    isOnFace = is_on_face(p_dir, p, grid_HiPts);
  }

  if (!isOnFace) {
    std::ostringstream warn;
    warn << "ERROR: Input file\n The Circle BC geometry is not correctly "
            "specified."
         << " The origin " << p << " must be on the same plane"
         << " as face (" << fc
         << "). Double check the origin and Level:box spec. \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (origin == "" || radius == "") {
    std::ostringstream warn;
    warn << "ERROR\n Circle BC geometry not correctly specified \n"
         << " you must specify origin [x,y,z] and radius [r] \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  return std::make_shared<CircleBCData>(p, r);
}

auto
BoundCondReader::createAnnulusBC(
  const std::map<std::string, std::string>& values,
  const std::vector<Point>& grid_LoPts,
  const std::vector<Point>& grid_HiPts,
  Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("annulus");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);

  const std::string& origin     = values.at("origin");
  const std::string& in_radius  = values.at("inner_radius");
  const std::string& out_radius = values.at("outer_radius");

  std::stringstream origin_stream(origin);
  std::stringstream in_radius_stream(in_radius);
  std::stringstream out_radius_stream(out_radius);

  double i_r, o_r, o[3];
  in_radius_stream >> i_r;
  out_radius_stream >> o_r;
  origin_stream >> o[0] >> o[1] >> o[2];

  Point p(o[0], o[1], o[2]);

  //  bullet proofing-- origin must be on the same plane as the face
  bool isOnFace = false;

  if (plusMinusFaces == -1) { // x-, y-, z- faces
    isOnFace = is_on_face(p_dir, p, grid_LoPts);
  }

  if (plusMinusFaces == 1) { // x+, y+, z+ faces
    isOnFace = is_on_face(p_dir, p, grid_HiPts);
  }

  if (!isOnFace) {
    std::ostringstream warn;
    warn << "ERROR: Input file\n The Annulus BC geometry is not correctly "
            "specified."
         << " The origin " << p << " must be on the same plane"
         << " as face (" << fc
         << "). Double check the origin and Level:box spec. \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (origin == "" || in_radius == "" || out_radius == "") {
    std::ostringstream warn;
    warn << "ERROR\n Annulus BC geometry not correctly specified \n"
         << " you must specify origin [x,y,z], inner_radius [r] outer_radius "
            "[r] \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  return std::make_shared<AnnulusBCData>(p, i_r, o_r);
}

auto
BoundCondReader::createEllipseBC(
  const std::map<std::string, std::string>& values,
  const std::vector<Point>& grid_LoPts,
  const std::vector<Point>& grid_HiPts,
  Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("ellipse");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);

  const std::string& str_origin       = values.at("origin");
  const std::string& str_minor_radius = values.at("minor_radius");
  const std::string& str_major_radius = values.at("major_radius");
  const std::string& str_angle        = values.at("angle");

  std::stringstream origin_stream(str_origin);
  std::stringstream minor_radius_stream(str_minor_radius);
  std::stringstream major_radius_stream(str_major_radius);
  std::stringstream angle_stream(str_angle);

  double minor_r, major_r, origin[3], angle;
  minor_radius_stream >> minor_r;
  major_radius_stream >> major_r;
  origin_stream >> origin[0] >> origin[1] >> origin[2];

  Point p(origin[0], origin[1], origin[2]);
  angle_stream >> angle;

  //  bullet proofing-- origin must be on the same plane as the face
  bool isOnFace = false;

  if (plusMinusFaces == -1) { // x-, y-, z- faces
    isOnFace = is_on_face(p_dir, p, grid_LoPts);
  }

  if (plusMinusFaces == 1) { // x+, y+, z+ faces
    isOnFace = is_on_face(p_dir, p, grid_HiPts);
  }

  if (!isOnFace) {
    std::ostringstream warn;
    warn << "ERROR: Input file\n The Ellipse BC geometry is not correctly "
            "specified."
         << " The origin " << p << " must be on the same plane"
         << " as face (" << fc
         << "). Double check the origin and Level:box spec. \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (major_r < minor_r) {
    std::ostringstream warn;
    warn << "ERROR\n Ellipse BC geometry not correctly specified \n"
         << " Major radius must be larger than minor radius \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (str_origin == "" || str_minor_radius == "" || str_major_radius == "") {
    std::ostringstream warn;
    warn << "ERROR\n Ellipse BC geometry not correctly specified \n"
         << " you must specify origin [x,y,z], inner_radius [r] outer_radius "
            "[r] \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  return std::make_shared<EllipseBCData>(p, minor_r, major_r, fc, angle);
}

auto
BoundCondReader::createRectangleBC(
  const std::map<std::string, std::string>& values,
  const std::vector<Point>& grid_LoPts,
  const std::vector<Point>& grid_HiPts,
  Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("rectangle");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);

  const std::string& low = values.at("lower");
  const std::string& up  = values.at("upper");

  std::stringstream low_stream(low), up_stream(up);
  double lower[3], upper[3];
  low_stream >> lower[0] >> lower[1] >> lower[2];
  up_stream >> upper[0] >> upper[1] >> upper[2];

  Point l(lower[0], lower[1], lower[2]), u(upper[0], upper[1], upper[2]);

  //  bullet proofing-- rectangle must be on the same plane as the face
  bool isOnFace = false;

  if (plusMinusFaces == -1) { // x-, y-, z- faces
    isOnFace =
      is_on_face(p_dir, l, grid_LoPts) && is_on_face(p_dir, u, grid_LoPts);
  }

  if (plusMinusFaces == 1) { // x+, y+, z+ faces
    isOnFace =
      is_on_face(p_dir, l, grid_HiPts) && is_on_face(p_dir, u, grid_HiPts);
  }

  if (!isOnFace) {
    std::ostringstream warn;
    warn << "ERROR: Input file\n The rectangle BC geometry is not correctly "
            "specified."
         << " The low " << l << " high " << u
         << " points must be on the same plane"
         << " as face (" << fc
         << "). Double check against and Level:box spec. \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (low == "" || up == "") {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry not correctly specified \n"
         << " you must specify lower [x,y,z] and upper[x,y,z] \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if ((l.x() > u.x() || l.y() > u.y() || l.z() > u.z()) ||
      (l.x() == u.x() && l.y() == u.y() && l.z() == u.z())) {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry not correctly specified \n"
         << " lower pt " << l << " upper pt " << u;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  return std::make_shared<RectangleBCData>(l, u);
}

auto
BoundCondReader::createRectangulusBC(
  const std::map<std::string, std::string>& values,
  const std::vector<Point>& grid_LoPts,
  const std::vector<Point>& grid_HiPts,
  Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>
{
  const auto& fc = values.at("rectangulus");
  DOUT(BCR_dbg, "Face = " << fc);

  int plusMinusFaces = 0;
  int p_dir          = 0;
  whichPatchFace(fc, face_side, plusMinusFaces, p_dir);

  const std::string& in_low_str  = values.at("inner_lower");
  const std::string& in_up_str   = values.at("inner_upper");
  const std::string& out_low_str = values.at("outer_lower");
  const std::string& out_up_str  = values.at("outer_upper");

  if (in_low_str == "" || in_up_str == "") {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry on face " << fc
         << " was not correctly specified \n"
         << " you must specify inner_lower \"x,y,z\" and "
            "inner_upper\"x,y,z\" \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (out_low_str == "" || out_up_str == "") {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry on face " << fc
         << " was not correctly specified \n"
         << " you must specify outer_lower \"x,y,z\" and outer_upper "
            "\"x,y,z\" \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  std::vector<char> delimiters = { ' ', ',' }; // used to parse input strings
  Point in_low                 = string_to_Point(in_low_str, delimiters);
  Point in_up                  = string_to_Point(in_up_str, delimiters);
  Point out_low                = string_to_Point(out_low_str, delimiters);
  Point out_up                 = string_to_Point(out_up_str, delimiters);

  //  bullet proofing-- both rectangles must be on the same plane as the
  //  face
  bool isOnFace_inLow =
    isPtOnFace(p_dir, plusMinusFaces, in_low, grid_LoPts, grid_HiPts);
  bool isOnFace_inUp =
    isPtOnFace(p_dir, plusMinusFaces, in_up, grid_LoPts, grid_HiPts);
  bool isOnFace_outLow =
    isPtOnFace(p_dir, plusMinusFaces, out_low, grid_LoPts, grid_HiPts);
  bool isOnFace_outUp =
    isPtOnFace(p_dir, plusMinusFaces, out_up, grid_LoPts, grid_HiPts);

  if (!isOnFace_inLow || !isOnFace_inUp || !isOnFace_outLow ||
      !isOnFace_outUp) {
    std::ostringstream warn;
    warn << "ERROR: Input file\n The rectangle BC geometry on face " << fc
         << " was not correctly specified."
         << " The low " << in_low << " high " << in_up
         << " points must be on the same plane"
         << " The low " << out_low << " high " << out_up
         << " points must be on the same plane"
         << " as face (" << fc
         << "). Double check against and Level:box spec. \n\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if ((in_low.x() > in_up.x() || in_low.y() > in_up.y() ||
       in_low.z() > in_up.z()) ||
      (in_low.x() == in_up.x() && in_low.y() == in_up.y() &&
       in_low.z() == in_up.z())) {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry on face " << fc
         << " was not correctly specified \n"
         << " inner_lower " << in_low << " inner_upper " << in_up;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if ((out_low.x() > out_up.x() || out_low.y() > out_up.y() ||
       out_low.z() > out_up.z()) ||
      (out_low.x() == out_up.x() && out_low.y() == out_up.y() &&
       out_low.z() == out_up.z())) {
    std::ostringstream warn;
    warn << "ERROR\n Rectangle BC geometry on face " << fc
         << " was not correctly specified \n"
         << " outer_lower " << out_low << " outer_upper " << out_up;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  return std::make_shared<RectangulusBCData>(in_low, in_up, out_low, out_up);
}

auto
BoundCondReader::createInteriorBndBoundaryConditionFace(
  ProblemSpecP& face_ps,
  const ProblemSpecP& grid_ps,
  Patch::FaceType& face_side,
  const LevelP level) -> std::shared_ptr<BCGeomBase>
{

  // Determine the Level 0 grid high and low points, need by
  // the bullet proofing
  Point grid_LoPt(1e30, 1e30, 1e30);
  Point grid_HiPt(-1e30, -1e30, -1e30);

  std::vector<Point> grid_LoPts; // store the lower bounds of all boxes
  std::vector<Point> grid_HiPts; // store the upper bounds of all boxes

  for (ProblemSpecP level_ps = grid_ps->findBlock("Level"); level_ps != nullptr;
       level_ps              = level_ps->findNextBlock("Level")) {

    // Find upper/lower corner:
    for (ProblemSpecP box_ps = level_ps->findBlock("Box"); box_ps != nullptr;
         box_ps              = box_ps->findNextBlock("Box")) {
      Point lower;
      Point upper;
      box_ps->require("lower", lower);
      box_ps->require("upper", upper);
      grid_LoPts.push_back(lower);
      grid_HiPts.push_back(upper);
      grid_LoPt = Min(lower, grid_LoPt);
      grid_HiPt = Max(upper, grid_HiPt);
    }
  }

  std::map<std::string, std::string> values;
  face_ps->getAttributes(values);

  // Have three possible types for the boundary condition face:
  // a. side (original -- entire side is one bc)
  // b. cirle (part of the side consists of a circle)
  // c. rectangle (part of the side consists of a rectangle)
  // This allows us to specify variable boundary conditions on a given
  // side.  Will use the notion of a UnionBoundaryCondtion and Difference
  // BoundaryCondition.

  std::string fc;
  int plusMinusFaces, p_dir;
  std::shared_ptr<BCGeomBase> bcGeom;

  if (values.find("side") != values.end()) {
    std::ostringstream warn;
    warn << "ERROR: You cannot specify an internal side boundary condition.";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  } else if (values.find("circle") != values.end()) {
    fc = values["circle"];
    whichPatchFace(fc, face_side, plusMinusFaces, p_dir);
    std::string origin = values["origin"];
    std::string radius = values["radius"];
    std::stringstream origin_stream(origin);
    std::stringstream radius_stream(radius);
    double r, o[3];
    radius_stream >> r;
    origin_stream >> o[0] >> o[1] >> o[2];
    Point p0(o[0], o[1], o[2]);
    Point p = Uintah::BCReaderUtils::moveToClosestNode(
      level, p_dir, plusMinusFaces, p0);
    if (!radius_stream || !origin_stream) {
      std::cout << "WARNING: BoundCondReader.cc:  std::stringstream failed..."
                << std::endl;
    }
    bcGeom = std::make_shared<CircleBCData>(p, r);
  } else if (values.find("annulus") != values.end()) {
    fc = values["annulus"];
    whichPatchFace(fc, face_side, plusMinusFaces, p_dir);
    std::string origin     = values["origin"];
    std::string in_radius  = values["inner_radius"];
    std::string out_radius = values["outer_radius"];
    std::stringstream origin_stream(origin);
    std::stringstream in_radius_stream(in_radius);
    std::stringstream out_radius_stream(out_radius);
    double i_r, o_r, o[3];
    in_radius_stream >> i_r;
    out_radius_stream >> o_r;
    origin_stream >> o[0] >> o[1] >> o[2];
    Point p0(o[0], o[1], o[2]);
    if (origin == "" || in_radius == "" || out_radius == "") {
      std::ostringstream warn;
      warn << "ERROR\n Annulus BC geometry not correctly specified \n"
           << " you must specify origin [x,y,z], inner_radius [r] outer_radius "
              "[r] \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
    Point p = Uintah::BCReaderUtils::moveToClosestNode(
      level, p_dir, plusMinusFaces, p0);
    bcGeom = std::make_shared<AnnulusBCData>(p, i_r, o_r);
  } else if (values.find("ellipse") != values.end()) {
    fc = values["ellipse"];
    whichPatchFace(fc, face_side, plusMinusFaces, p_dir);
    std::string str_origin       = values["origin"];
    std::string str_minor_radius = values["minor_radius"];
    std::string str_major_radius = values["major_radius"];
    std::string str_angle        = values["angle"];
    std::stringstream origin_stream(str_origin);
    std::stringstream minor_radius_stream(str_minor_radius);
    std::stringstream major_radius_stream(str_major_radius);
    std::stringstream angle_stream(str_angle);
    double minor_r, major_r, origin[3], angle;
    minor_radius_stream >> minor_r;
    major_radius_stream >> major_r;
    origin_stream >> origin[0] >> origin[1] >> origin[2];
    Point p0(origin[0], origin[1], origin[2]);
    Point p = Uintah::BCReaderUtils::moveToClosestNode(
      level, p_dir, plusMinusFaces, p0);
    angle_stream >> angle;

    if (major_r < minor_r) {
      std::ostringstream warn;
      warn << "ERROR\n Ellipse BC geometry not correctly specified \n"
           << " Major radius must be larger than minor radius \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    if (str_origin == "" || str_minor_radius == "" || str_major_radius == "") {
      std::ostringstream warn;
      warn << "ERROR\n Ellipse BC geometry not correctly specified \n"
           << " you must specify origin [x,y,z], inner_radius [r] outer_radius "
              "[r] \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    bcGeom = std::make_shared<EllipseBCData>(p, minor_r, major_r, fc, angle);
  }

  else if (values.find("rectangle") != values.end()) {
    fc = values["rectangle"];
    whichPatchFace(fc, face_side, plusMinusFaces, p_dir);
    std::string low = values["lower"];
    std::string up  = values["upper"];
    std::stringstream low_stream(low), up_stream(up);
    double lower[3], upper[3];
    low_stream >> lower[0] >> lower[1] >> lower[2];
    up_stream >> upper[0] >> upper[1] >> upper[2];
    Point l0(lower[0], lower[1], lower[2]), u0(upper[0], upper[1], upper[2]);
    Point l = Uintah::BCReaderUtils::moveToClosestNode(
      level, p_dir, plusMinusFaces, l0);
    Point u = Uintah::BCReaderUtils::moveToClosestNode(
      level, p_dir, plusMinusFaces, u0);

    if (low == "" || up == "") {
      std::ostringstream warn;
      warn << "ERROR\n Rectangle BC geometry not correctly specified \n"
           << " you must specify lower [x,y,z] and upper[x,y,z] \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
    if ((l.x() > u.x() || l.y() > u.y() || l.z() > u.z()) ||
        (l.x() == u.x() && l.y() == u.y() && l.z() == u.z())) {
      std::ostringstream warn;
      warn << "ERROR\n Rectangle BC geometry not correctly specified \n"
           << " lower pt " << l << " upper pt " << u;
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    bcGeom = std::make_shared<RectangleBCData>(l, u);
  }

  else {
    std::ostringstream warn;
    warn << "ERROR\n Boundary condition geometry not correctly specified "
            " Valid options (side, circle, rectangle, annulus";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // name the boundary condition object:
  std::string bcname;
  if (values.find("name") != values.end()) {
    std::string name = values["name"];
    DOUT(BCR_dbg, "Setting name to: " << name);
    bcGeom->setBCName(name);
  }

  // get the bctype - mainly used by wasatch:
  if (values.find("type") != values.end()) {
    std::string bndType = values["type"];
    DOUT(BCR_dbg, "Setting bc type to: " << bndType);
    bcGeom->setBndType(bndType);
  }

  if (face_ps->findBlock("ParticleBC")) {
    ProblemSpecP particleBCps = face_ps->findBlock("ParticleBC");
    ProblemSpecP pWallBC      = particleBCps->findBlock("Wall");
    ProblemSpecP pInletBC     = particleBCps->findBlock("Inlet");
    BCGeomBase::ParticleBndSpec pBndSpec;
    if (pWallBC) {
      pBndSpec.bndType = BCGeomBase::ParticleBndSpec::WALL;
      std::string wallType;
      pWallBC->getAttribute("walltype", wallType);
      if (wallType == "Elastic") {
        pBndSpec.wallType        = BCGeomBase::ParticleBndSpec::ELASTIC;
        pBndSpec.restitutionCoef = 1.0;
      } else if (wallType == "Inelastic") {
        pBndSpec.wallType        = BCGeomBase::ParticleBndSpec::INELASTIC;
        pBndSpec.restitutionCoef = 0.0;
      } else if (wallType == "PartiallyElastic") {
        pBndSpec.wallType = BCGeomBase::ParticleBndSpec::PARTIALLYELASTIC;
        pWallBC->get("Restitution", pBndSpec.restitutionCoef);
      }
    } else if (pInletBC) {
      pBndSpec.bndType = BCGeomBase::ParticleBndSpec::INLET;
      pInletBC->get("ParticlesPerSecond", pBndSpec.particlesPerSec);
    }
    bcGeom->setParticleBndSpec(pBndSpec);
  }

  DOUT(BCR_dbg, "Face = " << fc);
  return bcGeom;
}

void
BoundCondReader::read(ProblemSpecP& bc_ps,
                      const ProblemSpecP& grid_ps,
                      const Uintah::LevelP level)
{
  readDomainBCs(bc_ps, grid_ps);
  readInteriorBndBCs(bc_ps, grid_ps, level);
}

void
BoundCondReader::readDomainBCs(ProblemSpecP& bc_ps, const ProblemSpecP& grid_ps)
{
  // This function first looks for the geometric specification for the
  // boundary condition which includes the tags side, circle and rectangle.
  // The function createBoundaryConditionFace parses the tag and creates the
  // appropriate class.  Once this class is created, then the actual boundary
  // conditions are parsed from the input file (Pressure, Density, etc.).
  // Boundary conditions can be specified for various materials within a
  // face.  This complicates things, so we have to check for this and then
  // separate them out.  Once things are separated out, we then must take
  // all the boundary conditions for a given face and material id and combine
  // them so that any circle or rectangles that are specified can be combined
  // appropriately with the side case.  Multiple circle/rectangles are added
  // together and stored in a Union class.  This union class is then subtracted
  // off from the side class resulting in a difference class.  The difference
  // class represents the region of the side minus any circles/rectangle.

  for (ProblemSpecP face_ps = bc_ps->findBlock("Face"); face_ps != 0;
       face_ps              = face_ps->findNextBlock("Face")) {

    Patch::FaceType face_side;
    auto bcGeom = createBoundaryConditionFace(face_ps, grid_ps, face_side);

    std::string face_label = "none";
    face_ps->getAttribute("name", face_label);
    DOUT(BCR_dbg, "DomainBCs:: Face Label = " << face_label);

    DOUT(BCR_dbg,
         "DomainBCs:: Face = " << face_side << " Geometry type = "
                               << typeid(*bcGeom).name() << " " << bcGeom);

    std::multimap<int, BoundCondBaseSP> bctype_data;

    for (ProblemSpecP child = face_ps->findBlock("BCType"); child != 0;
         child              = child->findNextBlock("BCType")) {
      int mat_id;
      BoundCondBaseSP bc = BoundCondFactory::create(child, mat_id, face_label);
      DOUT(BCR_dbg,
           "DomainBCs:: Inserting into mat_id = "
             << mat_id << " bc = " << bc->getBCVariable()
             << " bctype = " << bc->getBCType() << " " << bc);

      bctype_data.insert(std::pair<int, BoundCondBaseSP>(mat_id, bc->clone()));
    }

    // Print out all of the bcs just created
    for (auto& [mat_id, bc] : bctype_data) {
      DOUT(BCR_dbg,
           "DomainBCs:: Getting out mat_id = "
             << mat_id << " bc = " << bc->getBCVariable()
             << " bctype = " << bc->getBCType());
    }

    // Search through the newly created boundary conditions and create
    // new std::shared_ptr<BCGeomBase> clones if there are multi materials
    // specified in the give <Face>.  This is usually a problem when Pressure is
    // specified for material id = 0, and other bcs such as velocity,
    // temperature, etc. for material_id != 0.

    std::map<int, std::shared_ptr<BCGeomBase>> bcgeom_data;

    // Search through the bctype_data and make sure that there are
    // enough bcGeom clones for each material.
    for (const auto& [mat_id, bc] : bctype_data) {
      if (bcgeom_data.find(mat_id) == bcgeom_data.end()) {
        bcgeom_data[mat_id] = bcGeom->clone();
      }

      DOUT(BCR_dbg,
           "DomainBCs:: Storing in  = "
             << typeid(bcgeom_data[mat_id]).name() << " " << bcgeom_data[mat_id]
             << " " << typeid(bc.get()).name() << " " << bc);

      bcgeom_data[mat_id]->addBC(bc);
    }

    if (d_BCReaderData[face_side] == nullptr) {
      d_BCReaderData[face_side] = std::make_shared<BCDataArray>();
    }

    for (auto& [mat_id, bc_geom] : bcgeom_data) {
      DOUT(BCR_dbg,
           "DomainBCs:: Adding BC data ... " << mat_id << ":" << bc_geom);
      d_BCReaderData[face_side]->addBCData(mat_id, bc_geom->clone());
    }

    DOUT(BCR_dbg, "DomainBCs:: Printing out bcDataArray ... ");
    d_BCReaderData[face_side]->print();

    // Delete stuff in bctype_data
    bctype_data.clear();
    bcgeom_data.clear();

  } // loop over faces

  // Find the mat_id = "all" (-1) information and store it in each
  // materials boundary condition section.
  DOUT(BCR_dbg, "DomainBCs:: Add 'all' boundary condition information");
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {

    // If no domain BCs have been read in for this face then skip
    if (d_BCReaderData[face] == nullptr) {
      continue;
    }

    for (auto& [mat_id, bcgeom_vec] : d_BCReaderData[face]->d_BCDataArray) {
      if (mat_id == -1) {
        for (auto& bc_geom : bcgeom_vec) {
          d_BCReaderData[face]->addBCData(mat_id, bc_geom->clone());
        }
      }
    }

    DOUT(BCR_dbg, std::endl << "Combining BCGeometryTypes for face " << face);
    for (auto& [mat_id, bcgeom_vec] : d_BCReaderData[face]->d_BCDataArray) {
      DOUT(BCR_dbg, "mat_id = " << mat_id);
      d_BCReaderData[face]->combineBCGeometryTypes_NEW(mat_id);
    }

    DOUT(BCR_dbg,
         std::endl
           << "DomainBCs:: Printing out bcDataArray for face " << face
           << " after adding 'all' ... ");
    d_BCReaderData[face]->print();
  } // face loop

  // Need to take the individual boundary conditions and combine them into
  // a single different (side and the union of any holes (circles or
  // rectangles.  This only happens if there are more than 1 bc_data per
  // face.

  DOUT(BCR_dbg,
       std::endl
         << "DomainBCs:: Before combineBCS() . . ." << std::endl);
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {

    // If no domain BCs have been read in for this face then skip
    if (d_BCReaderData[face] == nullptr) {
      continue;
    }

    DOUT(BCR_dbg, "DomainBCs:: Before Face . . ." << face);
    d_BCReaderData[face]->print();
  }

  bulletProofing();

  combineBCS();

  DOUT(BCR_dbg,
       std::endl
         << "DomainBCs:: After combineBCS() . . ." << std::endl);
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {

    // If no domain BCs have been read in for this face then skip
    if (d_BCReaderData[face] == nullptr) {
      continue;
    }

    DOUT(BCR_dbg, "DomainBCs:: After Face . . .  " << face);
    d_BCReaderData[face]->print();
  }
}

void
BoundCondReader::readInteriorBndBCs(ProblemSpecP& bc_ps,
                                    const ProblemSpecP& grid_ps,
                                    const Uintah::LevelP level)
{
  // This function first looks for the geometric specification for the
  // boundary condition which includes the tags side, circle and rectangle.
  // The function createBoundaryConditionFace parses the tag and creates the
  // appropriate class.  Once this class is created, then the actual boundary
  // conditions are parsed from the input file (Pressure, Density, etc.).
  // Boundary conditions can be specified for various materials within a
  // face.  This complicates things, so we have to check for this and then
  // separate them out.  Once things are separated out, we then must take
  // all the boundary conditions for a given face and material id and combine
  // them so that any circle or rectangles that are specified can be combined
  // appropriately with the side case.  Multiple circle/rectangles are added
  // together and stored in a Union class.  This union class is then subtracted
  // off from the side class resulting in a difference class.  The difference
  // class represents the region of the side minus any circles/rectangle.

  std::string defaultMat      = "";
  ProblemSpecP defaultMatSpec = bc_ps->findBlock("DefaultMaterial");
  if (defaultMatSpec) {
    bc_ps->get("DefaultMaterial", defaultMat);
  }

  for (ProblemSpecP face_ps = bc_ps->findBlock("InteriorFace");
       face_ps != nullptr;
       face_ps = face_ps->findNextBlock("InteriorFace")) {

    Patch::FaceType face_side;
    std::shared_ptr<BCGeomBase> bcGeom = createInteriorBndBoundaryConditionFace(
      face_ps, grid_ps, face_side, level);

    std::string face_label = "none";
    face_ps->getAttribute("name", face_label);
    DOUT(BCR_dbg, "InteriorBCs:: Face Label = " << face_label);

    DOUT(BCR_dbg,
         "InteriorBCs:: Face = " << face_side << " Geometry type = "
                                 << typeid(*bcGeom).name() << " " << bcGeom);

    std::multimap<int, BoundCondBaseSP> bctype_data;

    for (ProblemSpecP child = face_ps->findBlock("BCType"); child != nullptr;
         child              = child->findNextBlock("BCType")) {
      int mat_id;

      std::map<std::string, std::string> bc_attr;
      child->getAttributes(bc_attr);
      bool foundMatlID = (bc_attr.find("id") != bc_attr.end());
      if (!foundMatlID) {
        if (defaultMat == "") {
          SCI_THROW(ProblemSetupException(
            "ERROR: No material id was specified in the BCType tag and I could "
            "not find a DefaulMaterial to use! Please revise your input file.",
            __FILE__,
            __LINE__));
        } else {
          mat_id = (defaultMat == "all") ? -1 : atoi(defaultMat.c_str());
        }
      } else {
        std::string id = bc_attr["id"];
        mat_id         = (id == "all") ? -1 : atoi(id.c_str());
      }

      BoundCondBaseSP bc = BoundCondFactory::create(child, mat_id, face_label);
      DOUT(BCR_dbg,
           "InteriorBCs:: Inserting into mat_id = "
             << mat_id << " bc = " << bc->getBCVariable()
             << " bctype = " << bc->getBCType() << " " << bc);

      bctype_data.insert(std::pair<int, BoundCondBaseSP>(mat_id, bc->clone()));
    }

    // Print out all of the bcs just created
    for (auto& [mat_id, bc] : bctype_data) {
      DOUT(BCR_dbg,
           "InteriorBCs:: Getting out mat_id = "
             << mat_id << " bc = " << bc->getBCVariable()
             << " bctype = " << bc->getBCType());
    }

    // Search through the newly created boundary conditions and create
    // new std::shared_ptr<BCGeomBase> clones if there are multi materials
    // specified in the give <Face>.  This is usually a problem when Pressure is
    // specified for material id = 0, and other bcs such as velocity,
    // temperature, etc. for material_id != 0.

    std::map<int, std::shared_ptr<BCGeomBase>> bcgeom_data;

    // Search through the bctype_data and make sure that there are
    // enough bcGeom clones for each material.
    for (auto& [mat_id, bc] : bctype_data) {
      if (bcgeom_data.find(mat_id) == bcgeom_data.end()) {
        bcgeom_data[mat_id] = bcGeom->clone();
      }

      DOUT(BCR_dbg,
           "InteriorBCs:: Storing in  = " << typeid(bcgeom_data[mat_id]).name()
                                          << " " << bcgeom_data[mat_id] << " "
                                          << typeid(*(bc)).name() << " " << bc);

      bcgeom_data[mat_id]->addBC(bc);
    }

    //____________________________________________________________________
    // CAUTION! tsaad: If NO BCs have been specified on this boundary, then NO
    // iterators for that boundary will be added. The next if-statement
    // circumvents that problem for lack of a better design. This is done to
    // reduce input-file clutter. For example, for a constant density flow
    // problem a stationary-wall boundary is well defined and there's no reason
    // for the user to input any BCs there. To be able to set BCs through the
    // code, we still need access to the iterator for that boundary.
    if (bctype_data.size() == 0) {
      bcgeom_data[-1] = bcGeom->clone();
    }
    //-------------------------------------------------------------------

    if (d_interiorBndBCReaderData[face_side] == nullptr) {
      d_interiorBndBCReaderData[face_side] = std::make_shared<BCDataArray>();
    }

    for (auto& [mat_id, bcgeom_sp] : bcgeom_data) {
      d_interiorBndBCReaderData[face_side]->addBCData(mat_id,
                                                      bcgeom_sp->clone());
    }

    DOUT(
      BCR_dbg,
      "InteriorBCs:: Printing out bcDataArray for face_side = " << face_side);
    d_interiorBndBCReaderData[face_side]->print();

    // Delete stuff in bctype_data
    bctype_data.clear();
    bcgeom_data.clear();

  } // loop over faces

  // Find the mat_id = "all" (-1) information and store it in each
  // materials boundary condition section.
  DOUT(BCR_dbg, "InteriorBCs:: Add 'all' boundary condition information");
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {

    // If no interior BCs have been read in for this face then skip
    if (d_interiorBndBCReaderData[face] == nullptr) {
      continue;
    }

    for (auto& [mat_id, bcgeom_vec] :
         d_interiorBndBCReaderData[face]->d_BCDataArray) {
      if (mat_id == -1) {
        for (auto& bc_geom : bcgeom_vec) {
          d_interiorBndBCReaderData[face]->addBCData(mat_id, bc_geom->clone());
        }
      }
    }

    DOUT(BCR_dbg, "InteriorBCs:: Combining BCGeometryTypes for face " << face);
    for (auto& [id, bc_geom] : d_interiorBndBCReaderData[face]->d_BCDataArray) {
      DOUT(BCR_dbg, "mat_id = " << id);
      // d_interiorBndBCReaderData[face].combineBCGeometryTypes_NEW(bc_geom_itr->first);
    }

    DOUT(BCR_dbg,
         "InteriorBCs:: Printing out bcDataArray for face "
           << face << " after adding 'all' . . ");
    d_interiorBndBCReaderData[face]->print();
  } // face loop

  // Need to take the individual boundary conditions and combine them into
  // a single different (side and the union of any holes (circles or
  // rectangles.  This only happens if there are more than 1 bc_data per
  // face.

  DOUT(BCR_dbg, "InteriorBCs:: Before combineBCS() . . .");
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {
    DOUT(BCR_dbg, "InteriorBCs:: Before Face . . ." << face);

    // If no interior BCs have been read in for this face then skip
    if (d_interiorBndBCReaderData[face] == nullptr) {
      continue;
    }

    d_interiorBndBCReaderData[face]->print();
  }

  // bulletProofing();

  // combineBCS();

  DOUT(BCR_dbg,
       std::endl
         << "InteriorBCs:: After combineBCS() . . ." << std::endl);
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {
    DOUT(BCR_dbg, "InteriorBCs:: After Face . . .  " << face);

    // If no interior BCs have been read in for this face then skip
    if (d_interiorBndBCReaderData[face] == nullptr) {
      continue;
    }

    d_interiorBndBCReaderData[face]->print();
  }
}

auto
BoundCondReader::getBCDataArray(Patch::FaceType& face) const
  -> const BCDataArray
{
  auto m = this->d_BCReaderData;
  return *m[face];
}

void
BoundCondReader::combineBCS()
{
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {
    DOUT(BCR_dbg, std::endl << "Working on Face = " << face);
    DOUT(BCR_dbg, std::endl << "Original inputs");

    auto rearranged = std::make_shared<BCDataArray>();
    auto original   = d_BCReaderData[face];
    if (original == nullptr) {
      continue;
    }

    original->print();
    DOUT(BCR_dbg, std::endl);

    for (auto& [mat_id, bcgeom_vec] : original->d_BCDataArray) {

      DOUT(BCR_dbg, "Mat ID = " << mat_id);

      // Don't do anything if the only BCGeomBase element is a SideBC
      if ((bcgeom_vec.size() == 1) &&
          (std::count_if(bcgeom_vec.begin(),
                         bcgeom_vec.end(),
                         cmp_type<SideBCData>()) == 1)) {

        rearranged->addBCData(mat_id, bcgeom_vec[0]->clone());
      }

      // If there is more than one BCGeomBase element find the SideBC
      if (bcgeom_vec.size() > 1) {

        int num_other = std::count_if(
          bcgeom_vec.begin(), bcgeom_vec.end(), not_type<SideBCData>());

        DOUT(BCR_dbg, "num_other = " << num_other << std::endl);

        if (num_other == 1) {

          auto side_index = std::find_if(
            bcgeom_vec.begin(), bcgeom_vec.end(), cmp_type<SideBCData>());
          auto other_index = std::find_if(
            bcgeom_vec.begin(), bcgeom_vec.end(), not_type<SideBCData>());

          auto side_bc  = (*side_index)->clone();
          auto other_bc = (*other_index)->clone();

          auto diff_bc = std::make_shared<DifferenceBCData>(side_bc, other_bc);

          diff_bc->setBCName(
            side_bc->getBCName()); // make sure the new piece has the right name
          diff_bc->setBndType(
            side_bc->getBndType()); // make sure the new piece has the correct
                                    // boundary type
          diff_bc->setParticleBndSpec(side_bc->getParticleBndSpec());

          rearranged->addBCData(mat_id, diff_bc->clone());
          rearranged->addBCData(mat_id, other_bc->clone());

        } else {

          auto union_bc = std::make_shared<UnionBCData>();
          // Need to clone the RectangleBC that are being inserted
          // into the UnionBC
          for (auto bc_geom : bcgeom_vec) {
            if (!cmp_type<SideBCData>()(bc_geom)) {
              union_bc->child.push_back(bc_geom->clone());
            }
          }

          auto side_index = std::find_if(
            bcgeom_vec.begin(), bcgeom_vec.end(), cmp_type<SideBCData>());

          auto side_bc = (*side_index)->clone();

          auto union_bc_clone = union_bc->clone();
          auto union_bc_ptr   = static_cast<UnionBCData*>(union_bc_clone.get());

          auto diff_bc =
            std::make_shared<DifferenceBCData>(side_bc, union_bc_clone);

          diff_bc->setBCName(
            side_bc->getBCName()); // make sure the new piece has the right name
          diff_bc->setBndType(
            side_bc->getBndType()); // make sure the new piece has the correct
                                    // boundary type
          diff_bc->setParticleBndSpec(side_bc->getParticleBndSpec());

          rearranged->addBCData(mat_id, diff_bc->clone());
          for (auto& bc : union_bc_ptr->child) {
            rearranged->addBCData(mat_id, bc->clone());
          }
        }
      }

      std::for_each(
        bcgeom_vec.begin(), bcgeom_vec.end(), delete_object<BCGeomBase>());
      bcgeom_vec.clear();
    }
    DOUT(BCR_dbg, std::endl << "Printing out rearranged list");
    rearranged->print();

    d_BCReaderData[face] = rearranged;

    DOUT(BCR_dbg,
         std::endl
           << "Printing out rearranged from d_BCReaderData list");

    d_BCReaderData[face]->print();
  }

  DOUT(BCR_dbg, std::endl << "Printing out in combineBCS()");
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {
    DOUT(BCR_dbg, "After Face . . .  " << face);
    const auto& data = d_BCReaderData[face];
    if (data == nullptr) {
      continue;
    }
    data->print();
  }
}

auto
BoundCondReader::compareBCData([[maybe_unused]] std::shared_ptr<BCGeomBase> b1,
                               [[maybe_unused]] std::shared_ptr<BCGeomBase> b2)
  -> bool
{
  return false;
}

void
BoundCondReader::bulletProofing()
{
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {

    auto& original = d_BCReaderData[face];
    if (original == nullptr) {
      continue;
    }

    for (auto& [id, bcgeom_vec] : original->d_BCDataArray) {

      // There must be 1 and only 1 side BC specified
      int nSides = std::count_if(
        bcgeom_vec.begin(), bcgeom_vec.end(), cmp_type<SideBCData>());

      if (nSides != 1) {
        std::ostringstream warn;
        warn << "ERROR: <BoundaryConditions> <" << Patch::getFaceName(face)
             << ">: SideBCData # = " << nSides << "\n"
             << "There must be at least 1 and only 1 side boundary condition "
                "specified "
                "\n\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
    } // BCDataArray
  }   // patch faces
}

} // end namespace Uintah

namespace Uintah::BCReaderUtils {

void
print(std::shared_ptr<BCGeomBase> p)
{
  DOUT(BCR_dbg, "type = " << typeid(*p).name());
}

auto
moveToClosestNode(const LevelP level,
                  const int facedir,
                  const int plusMinusFaces,
                  const Point& p0) -> Point
{
  //
  // now find the closest node
  Vector halfdx = level->dCell() / 2.0;
  Point newPos  = p0;

  // find the closest cell
  IntVector cellIdx = level->getCellIndex(
    p0); // find the closest cell to the center of this circle
  Point closestCell = level->getCellPosition(cellIdx);

  // move the appropriate coordinate to the closest cell
  newPos(facedir) = closestCell(facedir);

  Point leftNode  = newPos;
  Point rightNode = newPos;
  leftNode(facedir) -= halfdx[facedir];
  rightNode(facedir) += halfdx[facedir];

  Vector diffLeft  = p0 - leftNode;
  Vector diffRight = rightNode - p0;

  if (diffRight.length() >
      diffLeft.length()) { // object is closer to the left node
    newPos(facedir) -=
      halfdx[facedir]; // move the circle to the closest layer of nodes
  } else if (diffRight.length() < diffLeft.length()) {
    newPos(facedir) +=
      halfdx[facedir]; // move the circle to the closest layer of nodes
  } else {
    newPos(facedir) +=
      plusMinusFaces *
      halfdx[facedir]; // move the circle to the closest layer of nodes
  }
  return newPos;
}

} // namespace Uintah::BCReaderUtils
