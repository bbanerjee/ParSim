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

/********************************************************************************
    Crack.cc
    PART ONE: CONSTRUCTOR, DECONSTRUCTOR, READ IN AND DISCRETIZE CRACKS

    Created by Yajun Guo in 2002-2005.
********************************************************************************/

#include "Crack.h"
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Crack/CrackGeometry.h>
#include <CCA/Components/MPM/Crack/CrackGeometryFactory.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <format>
#include <fstream>
#include <iostream>
#include <vector>

using namespace Uintah;

using std::string;
using std::vector;

#define MAX_BASIS 27

Crack::Crack(const ProblemSpecP& ps,
             MaterialManagerP& mat_manager,
             Output* data_archiver,
             MPMLabel* mpm_labels,
             MPMFlags* mpm_flags)
{
  MPI_Comm_dup(MPI_COMM_WORLD, &d_mpi_crack_comm);

  // Task 1: Initialization of fracture analysis

  d_mat_manager  = mat_manager;
  d_dataArchiver = data_archiver;
  lb             = mpm_labels;
  d_flag         = mpm_flags;
  d_n8or27       = d_flag->d_8or27;

  if (d_n8or27 == 8) {
    d_NGP = 1;
    d_NGN = 1;
  } else if (d_n8or27 == MAX_BASIS) {
    d_NGP = 2;
    d_NGN = 2;
  }

  d_calFractParameters = false;
  d_doCrackPropagation = false;
  d_useVolumeIntegral  = false;
  d_smoothCrackFront   = false;
  d_saveCrackGeometry  = true;

  d_rdadx     = 1.;  // Ratio of crack incremental to cell-size
  d_rJ        = -1.; // Radius of J-integral contour
  d_NJ        = 2;   // J-integral contour size
  d_CODOption = 0;   // Calculate COD at a fixed location by default

  d_computeJKInterval = 0.; // Intervals of calculating J & K
  d_growCrackInterval = 0.; // Interval of doing crack propagation

  d_GLP = Point(-9e99, -9e99, -9e99); // Highest global grid
  d_GHP = Point(9e99, 9e99, 9e99);    // Lowest global grid

  // Initialize boundary type
  for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
       face                 = Patch::nextFace(face)) {
    d_GridBCType[face] = "None";
  }

  // Task 2: Read in MPM parameters related to fracture analysis

  ProblemSpecP mpm_soln_ps = ps->findBlock("MPM");
  if (mpm_soln_ps) {
    mpm_soln_ps->get("calculate_fracture_parameters", d_calFractParameters);
    mpm_soln_ps->get("do_crack_propagation", d_doCrackPropagation);
    mpm_soln_ps->get("use_volume_integral", d_useVolumeIntegral);
    mpm_soln_ps->get("smooth_crack_front", d_smoothCrackFront);
    mpm_soln_ps->get("J_radius", d_rJ);
    mpm_soln_ps->get("dadx", d_rdadx);
    mpm_soln_ps->get("CODOption", d_CODOption);
  }

  // Get .uda directory
  ProblemSpecP uda_ps = ps->findBlock("DataArchiver");
  uda_ps->get("filebase", d_udaDir);
  uda_ps->get("save_crack_geometry", d_saveCrackGeometry);
  uda_ps->get("computeJKInterval", d_computeJKInterval);
  uda_ps->get("growCrackInterval", d_growCrackInterval);

  // Read in extent of the global grid
  ProblemSpecP grid_level_ps =
    ps->findBlock("Grid")->findBlock("Level")->findBlock("Box");
  grid_level_ps->get("lower", d_GLP);
  grid_level_ps->get("upper", d_GHP);

  // Read in boundary-condition types
  ProblemSpecP grid_bc_ps =
    ps->findBlock("Grid")->findBlock("BoundaryConditions");
  for (ProblemSpecP face_ps = grid_bc_ps->findBlock("Face"); face_ps != 0;
       face_ps              = face_ps->findNextBlock("Face")) {
    std::map<std::string, std::string> values;
    face_ps->getAttributes(values);
    ProblemSpecP bcType_ps = face_ps->findBlock("BCType");
    std::map<std::string, std::string> bc_attr;
    bcType_ps->getAttributes(bc_attr);
    if (values["side"] == "x-") {
      d_GridBCType[Patch::xminus] = bc_attr["var"];
    } else if (values["side"] == "x+") {
      d_GridBCType[Patch::xplus] = bc_attr["var"];
    } else if (values["side"] == "y-") {
      d_GridBCType[Patch::yminus] = bc_attr["var"];
    } else if (values["side"] == "y+") {
      d_GridBCType[Patch::yplus] = bc_attr["var"];
    } else if (values["side"] == "z-") {
      d_GridBCType[Patch::zminus] = bc_attr["var"];
    } else if (values["side"] == "z+") {
      d_GridBCType[Patch::zplus] = bc_attr["var"];
    }
  }

  // Task 3: Allocate memory for crack geometrical data

  int numMPMMatls = 0;
  ProblemSpecP mpm_ps =
    ps->findBlockWithOutAttribute("MaterialProperties")->findBlock("MPM");
  for (ProblemSpecP mat_ps = mpm_ps->findBlock("material"); mat_ps != 0;
       mat_ps              = mat_ps->findNextBlock("material")) {
    numMPMMatls++;
  }
  // Physical properties of cracks
  d_stressState.resize(numMPMMatls);
  d_crackType.resize(numMPMMatls);
  d_cmu.resize(numMPMMatls);

  // General quad crack segments
  d_quads.resize(numMPMMatls);
  d_quadN12.resize(numMPMMatls);
  d_quadN23.resize(numMPMMatls);
  d_quadCrackSidesAtFront.resize(numMPMMatls);
  d_quadRepetition.resize(numMPMMatls);
  d_quadOffset.resize(numMPMMatls);

  // Curved quad crack segments
  d_cquads.resize(numMPMMatls);
  d_cquadNStraightSides.resize(numMPMMatls);
  d_cquadPtsSide2.resize(numMPMMatls);
  d_cquadPtsSide4.resize(numMPMMatls);
  d_cquadCrackSidesAtFront.resize(numMPMMatls);
  d_cquadRepetition.resize(numMPMMatls);
  d_cquadOffset.resize(numMPMMatls);

  // Rriangular crack segments
  d_triangles.resize(numMPMMatls);
  d_triNCells.resize(numMPMMatls);
  d_triCrackSidesAtFront.resize(numMPMMatls);
  d_triRepetition.resize(numMPMMatls);
  d_triOffset.resize(numMPMMatls);

  // Arc crack segments
  d_arcs.resize(numMPMMatls);
  d_arcNCells.resize(numMPMMatls);
  d_arcCrkFrtSegID.resize(numMPMMatls);

  // Elliptic or partial alliptic crack segments
  d_ellipses.resize(numMPMMatls);
  d_ellipseNCells.resize(numMPMMatls);
  d_ellipseCrkFrtSegID.resize(numMPMMatls);
  d_pellipses.resize(numMPMMatls);
  d_pellipseNCells.resize(numMPMMatls);
  d_pellipseCrkFrtSegID.resize(numMPMMatls);
  d_pellipseExtent.resize(numMPMMatls);

  // Crack extent
  d_cmin.resize(numMPMMatls);
  d_cmax.resize(numMPMMatls);

  // Task 4:  Read in cracks

  int m = 0;
  for (ProblemSpecP mat_ps = mpm_ps->findBlock("material"); mat_ps != 0;
       mat_ps              = mat_ps->findNextBlock("material")) {
    ProblemSpecP crk_ps = mat_ps->findBlock("crack");
    if (crk_ps == 0) {
      d_crackType[m] = "NO_CRACK";
    }
    if (crk_ps != 0) {
      // Crack surface contact type, either "friction", "stick" or "null"
      crk_ps->require("type", d_crackType[m]);
      if (d_crackType[m] != "friction" && d_crackType[m] != "stick" &&
          d_crackType[m] != "null") {
        std::cout << "Error: unknown crack type: " << d_crackType[m] << "\n";
        exit(1);
      }

      // Friction coefficient needed for friction contact, zero by default
      d_cmu[m] = 0.0;
      if (d_crackType[m] == "friction") {
        crk_ps->get("mu", d_cmu[m]);
      }

      // Stress state at crack front
      d_stressState[m] = "planeStress";
      crk_ps->get("stress_state", d_stressState[m]);

// Read in crack segments. Presently seven kinds of basic shapes are available.
// More complicated crack plane can be input by combining the basic shapes.
#if 0
       CrackGeometry* cg = CrackGeometryFactory::create(crk_ps);
       d_crackGeometry.push_back(cg);
#endif

      ProblemSpecP geom_ps = crk_ps->findBlock("crack_segments");
      ReadQuadCracks(m, geom_ps);
      ReadCurvedQuadCracks(m, geom_ps);
      ReadTriangularCracks(m, geom_ps);
      ReadArcCracks(m, geom_ps);
      ReadEllipticCracks(m, geom_ps);
      ReadPartialEllipticCracks(m, geom_ps);
    }
    m++;
  }

  OutputInitialCrackPlane(numMPMMatls);
}

void
Crack::outputProblemSpec(ProblemSpecP& mpm_ps, ProblemSpecP& uda_ps)
{
  if (mpm_ps) {
    mpm_ps->appendElement("calculate_fracture_parameters",
                          d_calFractParameters);
    mpm_ps->appendElement("do_crack_propagation", d_doCrackPropagation);
    mpm_ps->appendElement("use_volume_integral", d_useVolumeIntegral);
    mpm_ps->appendElement("smooth_crack_front", d_smoothCrackFront);
    mpm_ps->appendElement("J_radius", d_rJ);
    mpm_ps->appendElement("dadx", d_rdadx);
    mpm_ps->appendElement("CODOption", d_CODOption);
  }

  if (uda_ps) {
    uda_ps->appendElement("save_crack_geometry", d_saveCrackGeometry);
    uda_ps->appendElement("computeJKInterval", d_computeJKInterval);
    uda_ps->appendElement("growCrackInterval", d_growCrackInterval);
  }

  int mat = 0;
  for (ProblemSpecP mat_ps = mpm_ps->findBlock("material"); mat_ps != nullptr;
       mat_ps              = mat_ps->findNextBlock("material")) {
    ProblemSpecP crk_ps = mat_ps->appendChild("crack");
    if (crk_ps) {
      crk_ps->appendElement("type", d_crackType[mat]);
      if (d_crackType[mat] == "friction") {
        crk_ps->appendElement("mu", d_cmu[mat]);
      }
    }

    crk_ps->appendElement("stress_state", d_stressState[mat]);
    ProblemSpecP geom_ps = crk_ps->appendChild("crack_segments");

    // *TODO* Output crack geometry inputs
    // outputQuadCrackProblemSpec(m, geom_ps);
    // outputReadCurvedQuadCrackProblemSpec(m, geom_ps);
    // outputReadTriangularCrackProblemSpec(m, geom_ps);
    // outputReadArcCrackProblemSpec(m, geom_ps);
    // outputReadEllipticCrackProblemSpec(m, geom_ps);
    // outputReadPartialEllipticCrackProblemSpec(m, geom_ps);

    mat++;
  }
}

void
Crack::ReadCurvedQuadCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP cquad_ps = geom_ps->findBlock("curved_quad"); cquad_ps != 0;
       cquad_ps              = cquad_ps->findNextBlock("curved_quad")) {

    // Four vertices of the curved quad
    Point p;
    std::vector<Point> vertices;
    cquad_ps->require("p1", p);
    vertices.push_back(p);
    cquad_ps->require("p2", p);
    vertices.push_back(p);
    cquad_ps->require("p3", p);
    vertices.push_back(p);
    cquad_ps->require("p4", p);
    vertices.push_back(p);
    d_cquads[m].push_back(vertices);
    vertices.clear();

    // Mesh resolution on two opposite straight sides
    int n = 1;
    cquad_ps->get("resolution_straight_sides", n);
    d_cquadNStraightSides[m].push_back(n);

    // Characteristic points on two opposite cuvered sides
    std::vector<Point> ptsSide2, ptsSide4;
    ProblemSpecP side2_ps = cquad_ps->findBlock("points_curved_side2");
    for (ProblemSpecP pt_ps = side2_ps->findBlock("point"); pt_ps != 0;
         pt_ps              = pt_ps->findNextBlock("point")) {
      pt_ps->get("val", p);
      ptsSide2.push_back(p);
    }
    d_cquadPtsSide2[m].push_back(ptsSide2);

    ProblemSpecP side4_ps = cquad_ps->findBlock("points_curved_side4");
    for (ProblemSpecP pt_ps = side4_ps->findBlock("point"); pt_ps != 0;
         pt_ps              = pt_ps->findNextBlock("point")) {
      pt_ps->get("val", p);
      ptsSide4.push_back(p);
    }
    d_cquadPtsSide4[m].push_back(ptsSide4);

    if (ptsSide4.size() != ptsSide2.size()) {
      std::cout << "Error: The points on curved side 2 and side 4 "
                << "should appear in pairs."
                << "\n";
    }
    ptsSide2.clear();
    ptsSide4.clear();

    // Crack front
    std::vector<short> crackSidesAtFront;
    string cfsides("");
    cquad_ps->get("crack_front_sides", cfsides);
    if (cfsides.length() == 4) {
      for (string::const_iterator iter = cfsides.begin(); iter != cfsides.end();
           iter++) {
        short atFront = NO;
        if (*iter == 'Y' || *iter == 'y') {
          atFront = YES;
        }
        crackSidesAtFront.push_back(atFront);
      }
    } else if (cfsides.length() == 0) {
      for (int i = 0; i < 4; i++) {
        crackSidesAtFront.push_back(NO);
      }
    }
    d_cquadCrackSidesAtFront[m].push_back(crackSidesAtFront);
    crackSidesAtFront.clear();

    // Repetition information
    n = 1;
    cquad_ps->get("repetition", n);
    d_cquadRepetition[m].push_back(n);

    Vector offset = Vector(0., 0., 0.);
    if (n > 1) {
      cquad_ps->require("offset", offset);
    }
    d_cquadOffset[m].push_back(offset);
  }
}

void
Crack::ReadQuadCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP quad_ps = geom_ps->findBlock("quad"); quad_ps != 0;
       quad_ps              = quad_ps->findNextBlock("quad")) {

    // Four vertices (p1-p4) of the quad
    Point p1, p2, p3, p4, p5, p6, p7, p8;
    std::vector<Point> vertices;
    quad_ps->require("p1", p1);
    vertices.push_back(p1);
    quad_ps->require("p2", p2);
    vertices.push_back(p2);
    quad_ps->require("p3", p3);
    vertices.push_back(p3);
    quad_ps->require("p4", p4);
    vertices.push_back(p4);

    // Four middle points (p5-p8) of the quad
    if (!quad_ps->get("p5", p5)) {
      p5 = p1 + 0.5 * (p2 - p1);
    }
    vertices.push_back(p5);
    if (!quad_ps->get("p6", p6)) {
      p6 = p2 + 0.5 * (p3 - p2);
    }
    vertices.push_back(p6);
    if (!quad_ps->get("p7", p7)) {
      p7 = p3 + 0.5 * (p4 - p3);
    }
    vertices.push_back(p7);
    if (!quad_ps->get("p8", p8)) {
      p8 = p1 + 0.5 * (p4 - p1);
    }
    vertices.push_back(p8);

    d_quads[m].push_back(vertices);
    vertices.clear();

    // Mesh resolutions
    int n12 = 1, n23 = 1;
    quad_ps->get("resolution_p1_p2", n12);
    d_quadN12[m].push_back(n12);
    quad_ps->get("resolution_p2_p3", n23);
    d_quadN23[m].push_back(n23);

    // Crack front
    std::vector<short> crackSidesAtFront;
    string cfsides("");
    quad_ps->get("crack_front_sides", cfsides);
    if (cfsides.length() == 4) {
      for (auto side : cfsides) {
        short atFront = NO;
        if (std::tolower(side) == 'y') {
          atFront = YES;
        }
        crackSidesAtFront.push_back(atFront);
      }
    } else if (cfsides.length() == 0) {
      for (int i = 0; i < 4; i++) {
        crackSidesAtFront.push_back(NO);
      }
    }
    d_quadCrackSidesAtFront[m].push_back(crackSidesAtFront);
    crackSidesAtFront.clear();

    // Repetition information
    int n = 1;
    quad_ps->get("repetition", n);
    d_quadRepetition[m].push_back(n);

    Vector offset = Vector(0., 0., 0.);
    if (n > 1) {
      quad_ps->require("offset", offset);
    }
    d_quadOffset[m].push_back(offset);
  }
}

void
Crack::ReadTriangularCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP tri_ps = geom_ps->findBlock("triangle"); tri_ps != 0;
       tri_ps              = tri_ps->findNextBlock("triangle")) {

    // Three vertices (p1-p3) of the triangle
    Point p1, p2, p3, p4, p5, p6;
    std::vector<Point> vertices;
    tri_ps->require("p1", p1);
    vertices.push_back(p1);
    tri_ps->require("p2", p2);
    vertices.push_back(p2);
    tri_ps->require("p3", p3);
    vertices.push_back(p3);

    // Three middle points (p4-p6) of the triangle
    if (!tri_ps->get("p4", p4)) {
      p4 = p1 + 0.5 * (p2 - p1);
    }
    vertices.push_back(p4);
    if (!tri_ps->get("p5", p5)) {
      p5 = p2 + 0.5 * (p3 - p2);
    }
    vertices.push_back(p5);
    if (!tri_ps->get("p6", p6)) {
      p6 = p3 + 0.5 * (p1 - p3);
    }
    vertices.push_back(p6);

    d_triangles[m].push_back(vertices);
    vertices.clear();

    // Mesh resolution
    int n = 1;
    tri_ps->get("resolution", n);
    d_triNCells[m].push_back(n);

    // Crack front
    string cfsides("");
    std::vector<short> crackSidesAtFront;
    tri_ps->get("crack_front_sides", cfsides);
    if (cfsides.length() == 3) {
      for (string::const_iterator iter = cfsides.begin(); iter != cfsides.end();
           iter++) {
        short atFront = NO;
        if (*iter == 'Y' || *iter == 'n') {
          atFront = YES;
        }
        crackSidesAtFront.push_back(atFront);
      }
    } else if (cfsides.length() == 0) {
      for (int i = 0; i < 3; i++) {
        crackSidesAtFront.push_back(NO);
      }
    }
    d_triCrackSidesAtFront[m].push_back(crackSidesAtFront);
    crackSidesAtFront.clear();

    // Repetition information
    n = 1;
    tri_ps->get("repetition", n);
    d_triRepetition[m].push_back(n);

    Vector offset = Vector(0., 0., 0.);
    if (n > 1) {
      tri_ps->require("offset", offset);
    }
    d_triOffset[m].push_back(offset);
  }
}

void
Crack::ReadArcCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP arc_ps = geom_ps->findBlock("arc"); arc_ps != 0;
       arc_ps              = arc_ps->findNextBlock("arc")) {

    // Three points on the arc
    Point p;
    std::vector<Point> thisArcPts;
    arc_ps->require("start_point", p);
    thisArcPts.push_back(p);
    arc_ps->require("middle_point", p);
    thisArcPts.push_back(p);
    arc_ps->require("end_point", p);
    thisArcPts.push_back(p);
    d_arcs[m].push_back(thisArcPts);
    thisArcPts.clear();

    // Resolution on circumference
    int n = 1;
    arc_ps->require("resolution_circumference", n);
    d_arcNCells[m].push_back(n);

    // Crack front segment ID, -1 by default which means all segments are crack
    // front
    int cfsID;
    if (!arc_ps->get("crack_front_segment_ID", cfsID)) {
      cfsID = -1;
    }
    d_arcCrkFrtSegID[m].push_back(cfsID);
  }
}

void
Crack::ReadEllipticCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP ellipse_ps = geom_ps->findBlock("ellipse"); ellipse_ps != 0;
       ellipse_ps              = ellipse_ps->findNextBlock("ellipse")) {

    // Three points on the arc
    Point p;
    std::vector<Point> thisEllipsePts;
    ellipse_ps->require("point1_axis1", p);
    thisEllipsePts.push_back(p);
    ellipse_ps->require("point_axis2", p);
    thisEllipsePts.push_back(p);
    ellipse_ps->require("point2_axis1", p);
    thisEllipsePts.push_back(p);
    d_ellipses[m].push_back(thisEllipsePts);
    thisEllipsePts.clear();

    // Resolution on circumference
    int n = 1;
    ellipse_ps->require("resolution_circumference", n);
    d_ellipseNCells[m].push_back(n);

    // Crack front segment ID, -1 by default which means all segments are crack
    // front
    int cfsID;
    if (!ellipse_ps->get("crack_front_segment_ID", cfsID)) {
      cfsID = -1;
    }
    d_ellipseCrkFrtSegID[m].push_back(cfsID);
  }
}

void
Crack::ReadPartialEllipticCracks(const int& m, const ProblemSpecP& geom_ps)
{
  for (ProblemSpecP pellipse_ps = geom_ps->findBlock("partial_ellipse");
       pellipse_ps != 0;
       pellipse_ps = pellipse_ps->findNextBlock("partial_ellipse")) {

    // Center,two points on major and minor axes
    Point p;
    std::vector<Point> thispEllipsePts;
    pellipse_ps->require("center", p);
    thispEllipsePts.push_back(p);
    pellipse_ps->require("point_axis1", p);
    thispEllipsePts.push_back(p);
    pellipse_ps->require("point_axis2", p);
    thispEllipsePts.push_back(p);
    d_pellipses[m].push_back(thispEllipsePts);
    thispEllipsePts.clear();

    // Extent in degree of the partial ellipse rotating from axis1
    double extent = 360.;
    pellipse_ps->get("extent", extent);
    d_pellipseExtent[m].push_back(extent);

    // Resolution on circumference
    int n = 1;
    pellipse_ps->require("resolution_circumference", n);
    d_pellipseNCells[m].push_back(n);

    // Crack front segment ID, -1 by default which means all segments are crack
    // front
    int cfsID;
    if (!pellipse_ps->get("crack_front_segment_ID", cfsID)) {
      cfsID = -1;
    }
    d_pellipseCrkFrtSegID[m].push_back(cfsID);
  }
}

void
Crack::OutputInitialCrackPlane(const int& numMatls)
{
  int pid;
  MPI_Comm_rank(d_mpi_crack_comm, &pid);
  if (pid == 0) { // output from the first rank
    for (int m = 0; m < numMatls; m++) {
      if (d_crackType[m] == "NO_CRACK") {
        std::cout << "\nMaterial " << m << ": no crack exists"
                  << "\n";
      } else {
        std::cout << "\nMaterial " << m << ":\n"
                  << "  Crack contact type: " << d_crackType[m] << "\n";
        if (d_crackType[m] == "friction") {
          std::cout << "    Frictional coefficient: " << d_cmu[m] << "\n";
        }

        std::cout << "  Crack geometry:"
                  << "\n";
        // general quad cracks
        for (int i = 0; i < (int)d_quads[m].size(); i++) {
          std::cout << "  * Quad " << i + 1 << ": meshed by ["
                    << d_quadN12[m][i] << ", " << d_quadN23[m][i] << ", "
                    << d_quadN12[m][i] << ", " << d_quadN23[m][i] << "]"
                    << "\n";
          for (int j = 0; j < 8; j++) {
            std::cout << "    p" << j + 1 << ": " << d_quads[m][i][j] << "\n";
          }
          for (int j = 0; j < 4; j++) {
            if (d_quadCrackSidesAtFront[m][i][j]) {
              int j2 = (j + 2 < 5 ? j + 2 : 1);
              std::cout << "    Side " << j + 1 << " (p" << j + 1 << "-"
                        << "p" << j2 << ") is a crack front."
                        << "\n";
            }
          }
          // repetition information
          if (d_quadRepetition[m][i] > 1) {
            std::cout << "    The quad is repeated " << d_quadRepetition[m][i]
                      << " times with the offset " << d_quadOffset[m][i] << "."
                      << "\n";
          }
        }

        // curved quad cracks
        for (int i = 0; i < (int)d_cquads[m].size(); i++) {
          std::cout << "  * Curved quad " << i + 1 << ":"
                    << "\n";
          std::cout << "    Four vertices:"
                    << "\n";
          // four vertices
          for (int j = 0; j < 4; j++) {
            std::cout << "      p" << j + 1 << ": " << d_cquads[m][i][j]
                      << "\n";
          }
          // resolution on straight sides 1 & 3
          std::cout
            << "    Resolution on straight sides (sides p1-p2 and p3-p4):"
            << d_cquadNStraightSides[m][i] << "\n";
          // points on curved egde 2
          std::cout << "    Points on curved side 2 (p2-p3): "
                    << "\n";
          for (int j = 0; j < (int)d_cquadPtsSide2[m][i].size(); j++) {
            std::cout << "      p" << j + 1 << ": " << d_cquadPtsSide2[m][i][j]
                      << "\n";
          }
          // points on curved side 3
          std::cout << "    Points on curved side 4 (p1-p4): "
                    << "\n";
          for (int j = 0; j < (int)d_cquadPtsSide4[m][i].size(); j++) {
            std::cout << "      p" << j + 1 << ": " << d_cquadPtsSide4[m][i][j]
                      << "\n";
          }
          // crack-front sides
          for (int j = 0; j < 4; j++) {
            if (d_cquadCrackSidesAtFront[m][i][j]) {
              int j2 = (j + 2 < 5 ? j + 2 : 1);
              std::cout << "    Side " << j + 1 << " (p" << j + 1 << "-"
                        << "p" << j2 << ") is a crack front."
                        << "\n";
            }
          }
          // repetition information
          if (d_cquadRepetition[m][i] > 1) {
            std::cout << "    The quad is repeated " << d_cquadRepetition[m][i]
                      << " times with the offset " << d_cquadOffset[m][i] << "."
                      << "\n";
          }
        }

        // Triangular cracks
        for (int i = 0; i < (int)d_triangles[m].size(); i++) {
          std::cout << "  * Triangle " << i + 1 << ": meshed by ["
                    << d_triNCells[m][i] << ", " << d_triNCells[m][i] << ", "
                    << d_triNCells[m][i] << "]"
                    << "\n";
          for (int j = 0; j < 6; j++) {
            std::cout << "    p" << j + 1 << ": " << d_triangles[m][i][j]
                      << "\n";
          }
          for (int j = 0; j < 3; j++) {
            if (d_triCrackSidesAtFront[m][i][j]) {
              int j2 = (j + 2 < 4 ? j + 2 : 1);
              std::cout << "    side " << j + 1 << " (p" << j + 1 << "-"
                        << "p" << j2 << ") is a crack front."
                        << "\n";
            }
          }
          // repetition information
          if (d_triRepetition[m][i] > 1) {
            std::cout << "    The triangle is repeated "
                      << d_triRepetition[m][i] << " times with the offset "
                      << d_triOffset[m][i] << "."
                      << "\n";
          }
        }

        // Arc cracks
        for (int i = 0; i < (int)d_arcs[m].size(); i++) {
          std::cout << "  * Arc " << i + 1 << ": meshed by "
                    << d_arcNCells[m][i] << " cells on the circumference."
                    << "\n";
          if (d_arcCrkFrtSegID[m][i] == -1) {
            std::cout << "   crack front: on the arc"
                      << "\n";
          } else {
            std::cout << "   crack front segment ID: " << d_arcCrkFrtSegID[m][i]
                      << "\n";
          }
          std::cout << "\n    start, middle and end points of the arc:"
                    << "\n";
          for (int j = 0; j < 3; j++) {
            std::cout << "    p" << j + 1 << ": " << d_arcs[m][i][j] << "\n";
          }
        }

        // Elliptic cracks
        for (int i = 0; i < (int)d_ellipses[m].size(); i++) {
          std::cout << "  * Ellipse " << i + 1 << ": meshed by "
                    << d_ellipseNCells[m][i] << " cells on the circumference."
                    << "\n";
          if (d_ellipseCrkFrtSegID[m][i] == -1) {
            std::cout << "    crack front: on the ellipse circumference"
                      << "\n";
          } else {
            std::cout << "    crack front segment ID: "
                      << d_ellipseCrkFrtSegID[m][i] << "\n";
          }
          std::cout << "    end point on axis2: " << d_ellipses[m][i][1]
                    << "\n";
          std::cout << "    another end point on axis1: " << d_ellipses[m][i][2]
                    << "\n";
        }

        // Partial elliptic cracks
        for (int i = 0; i < (int)d_pellipses[m].size(); i++) {
          std::cout << "  * Partial ellipse " << i + 1 << " ("
                    << d_pellipseExtent[m][i] << " degree): meshed by "
                    << d_pellipseNCells[m][i] << " cells on the circumference."
                    << "\n";
          if (d_pellipseCrkFrtSegID[m][i] == -1) {
            std::cout << "    crack front: on the ellipse circumference"
                      << "\n";
          } else {
            std::cout << "    crack front segment ID: "
                      << d_pellipseCrkFrtSegID[m][i] << "\n";
          }
          std::cout << "    center: " << d_pellipses[m][i][0] << "\n";
          std::cout << "    end point on axis1: " << d_pellipses[m][i][1]
                    << "\n";
          std::cout << "    end point on axis2: " << d_pellipses[m][i][2]
                    << "\n";
        }
      }
    } // End of loop over materials

    // Ratio of crack propagation incremental to cell-size
    if (d_doCrackPropagation) {
      std::cout << "  Ratio of crack increment to cell size (dadx) = "
                << d_rdadx << "."
                << "\n"
                << "\n";
    }
  }
}

void
Crack::addComputesAndRequiresCrackDiscretization(
  Task* /*t*/,
  const PatchSet* /*patches*/,
  const MaterialSet* /*matls*/) const
{
  // Do nothing currently
}

void
Crack::CrackDiscretization(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* /*matls*/,
                           DataWarehouse* /*old_dw*/,
                           DataWarehouse* /*new_dw*/)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int pid, rankSize;
    MPI_Comm_rank(d_mpi_crack_comm, &pid);
    MPI_Comm_size(d_mpi_crack_comm, &rankSize);

    // Set radius (rJ) of J-integral contour or number of cells
    Vector dx     = patch->dCell();
    double dx_min = Min(dx.x(), dx.y(), dx.z());

    if (d_rJ < 0.) { // Input NJ, and calculate rJ
      d_rJ = d_NJ * dx_min;
    } else { // Input rJ, and calculate NJ
      d_NJ = (int)(d_rJ / dx_min);
    }

    // Allocate memories for crack mesh
    int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    d_css.resize(numMPMMatls);
    d_csa.resize(numMPMMatls);
    d_cx.resize(numMPMMatls);
    d_ce.resize(numMPMMatls);
    d_cfSegNodes.resize(numMPMMatls);
    d_cfSegTime.resize(numMPMMatls);
    d_cfSegDis.resize(numMPMMatls);
    d_cfSegVel.resize(numMPMMatls);
    d_cfSegPreIdx.resize(numMPMMatls);
    d_cfSegMinIdx.resize(numMPMMatls);
    d_cfSegMaxIdx.resize(numMPMMatls);
    d_cfSegPtsT.resize(numMPMMatls);
    d_cfSegV1.resize(numMPMMatls);
    d_cfSegV2.resize(numMPMMatls);
    d_cfSegV3.resize(numMPMMatls);
    d_cfSegJ.resize(numMPMMatls);
    d_cfSegK.resize(numMPMMatls);
    d_cnset.resize(numMPMMatls);
    d_cfnset.resize(numMPMMatls);
    d_cfsset.resize(numMPMMatls);

    for (int m = 0; m < numMPMMatls; m++) {
      d_cnset[m].resize(rankSize);
      d_cfnset[m].resize(rankSize);
      d_cfsset[m].resize(rankSize);

      // Initialize crack extent
      d_cmin[m] = Point(1.e200, 1.e200, 1.e200);
      d_cmax[m] = Point(-1.e200, -1.e200, -1.e200);

      if (d_crackType[m] != "NO_CRACK") {
        // Discretize crack segments
        int nnode0 = 0;
        DiscretizeQuadCracks(m, nnode0);
        DiscretizeCurvedQuadCracks(m, nnode0);
        DiscretizeTriangularCracks(m, nnode0);
        DiscretizeArcCracks(m, nnode0);
        DiscretizeEllipticCracks(m, nnode0);
        DiscretizePartialEllipticCracks(m, nnode0);

        // Combine crack segments
        CombineCrackSegments(m);

        // Determine crack extent
        for (int i = 0; i < (int)d_cx[m].size(); i++) {
          d_cmin[m] = Min(d_cmin[m], d_cx[m][i]);
          d_cmax[m] = Max(d_cmax[m], d_cx[m][i]);
        }

        // Controlling parameters for fracture parameter calculation and crack
        // propagation
        if (d_calFractParameters || d_doCrackPropagation) {
          // Initialize parameters of crack-front nodes
          int num = (int)d_cfSegNodes[m].size();
          d_cfSegVel[m].resize(num);
          d_cfSegTime[m].resize(num);
          d_cfSegDis[m].resize(num);
          for (int i = 0; i < num; i++) {
            d_cfSegVel[m][i]  = 0.0;
            d_cfSegTime[m][i] = 0.0;
            d_cfSegDis[m][i]  = 0.0;
          }

          // Get average length of crack-front segments
          d_css[m]    = 0.;
          int ncfSegs = num / 2;
          for (int i = 0; i < ncfSegs; i++) {
            int n1 = d_cfSegNodes[m][2 * i];
            int n2 = d_cfSegNodes[m][2 * i + 1];
            d_css[m] += (d_cx[m][n1] - d_cx[m][n2]).length();
          }
          if (ncfSegs > 0) {
            d_css[m] /= ncfSegs;
          }

          // Determine connectivity of crack-front nodes
          FindCrackFrontNodeIndexes(m);

          // Get average angle (in degree) of crack-front segments
          d_csa[m]  = 0.;
          int count = 0;
          for (int i = 0; i < num; i++) {
            int preIdx = d_cfSegPreIdx[m][i];
            /* if (preIdx > 0) */
            // 1. Ensure a previous occurrence was found (preIdx != -1).
            // 2. Ensure (i - 2) is a valid non-negative index (i >= 2).
            // 3. Ensure (i + 1) is a valid index within bounds (i < num - 1).
            if (preIdx != -1 && i >= 2 && i < num - 1)
            {
              //std::cout << std::format("m = {}, i = {}, num = {}", m, i, num)
              //          << std::format("preIdx = {}, cfSegNodes = {}",
              //                         preIdx,
              //                         d_cfSegNodes[m][i])
              //          << std::endl;
              Point p   = d_cx[m][d_cfSegNodes[m][i]];
              Point p1  = d_cx[m][d_cfSegNodes[m][i - 2]];
              Point p2  = d_cx[m][d_cfSegNodes[m][i + 1]];
              Vector v1 = TwoPtsDirCos(p1, p);
              Vector v2 = TwoPtsDirCos(p, p2);
              d_csa[m] += std::abs(acos(Dot(v1, v2))) * 180 / 3.141592654;
              count++;
            }
          }
          if (count != 0) {
            d_csa[m] /= count;
          } else {
            d_csa[m] = 180;
          }

          // Calculate normals of crack plane at crack-front nodes
          if (d_smoothCrackFront) {
            if (!SmoothCrackFrontAndCalculateNormals(m)) {
              CalculateCrackFrontNormals(m);
            }
          } else {
            CalculateCrackFrontNormals(m);
          }
        }
#if 0
        OutputInitialCrackMesh(m);
#endif
      }
    } // End of loop over matls
  }
}

void
Crack::DiscretizeCurvedQuadCracks(const int& m, int& nnode0)
{
  int k, k1, i, j, ni, nj, n1, n2, n3, num;
  int nstart1, nstart2, nstart3;
  Point p1, p2, p3, p4, pt;

  for (k = 0; k < (int)d_cquads[m].size(); k++) {
    Vector offset = d_cquadOffset[m][k];
    for (k1 = 0; k1 < (int)d_cquadRepetition[m][k]; k1++) {
      // Four vertices of the curved quad
      p1 = d_cquads[m][k][0] + k1 * offset;
      p2 = d_cquads[m][k][1] + k1 * offset;
      p3 = d_cquads[m][k][2] + k1 * offset;
      p4 = d_cquads[m][k][3] + k1 * offset;

      // Mesh resolutions on curved sides (ni) & straight sides (nj)
      ni = d_cquadNStraightSides[m][k];
      nj = (int)(d_cquadPtsSide2[m][k].size()) + 1;

      // total number of nodes of the quad
      num = (ni + 1) * (nj + 1) + ni * nj;

      // Flag if node i is on edge j, initialized by NO
      short** nodeOnEdge = scinew short* [num];
      for (i = 0; i < num; i++) {
        nodeOnEdge[i] = scinew short[4];
      }
      for (i = 0; i < num; i++) {
        for (j = 0; j < 4; j++) {
          nodeOnEdge[i][j] = NO;
        }
      }

      // Nodes on curved sides 2 (p2-p3) & 4 (p1-p4) - "j" direction
      Point* p_s2  = new Point[2 * nj + 1];
      Point* p_s4  = new Point[2 * nj + 1];
      p_s2[0]      = p2;
      p_s2[2 * nj] = p3;
      p_s4[0]      = p1;
      p_s4[2 * nj] = p4;
      for (int l = 2; l < 2 * nj; l += 2) {
        p_s2[l] = d_cquadPtsSide2[m][k][l / 2 - 1] + k1 * offset;
        p_s4[l] = d_cquadPtsSide4[m][k][l / 2 - 1] + k1 * offset;
      }
      for (int l = 1; l < 2 * nj; l += 2) {
        p_s2[l] = p_s2[l - 1] + (p_s2[l + 1] - p_s2[l - 1]) / 2.;
        p_s4[l] = p_s4[l - 1] + (p_s4[l + 1] - p_s4[l - 1]) / 2.;
      }

      // Generate crack nodes
      int count = -1;
      for (j = 0; j <= nj; j++) {
        for (i = 0; i <= ni; i++) {
          // Detect edge nodes
          count++;
          if (j == 0) {
            nodeOnEdge[count][0] = YES;
          }
          if (i == ni) {
            nodeOnEdge[count][1] = YES;
          }
          if (j == nj) {
            nodeOnEdge[count][2] = YES;
          }
          if (i == 0) {
            nodeOnEdge[count][3] = YES;
          }
          pt = p_s4[2 * j] + (p_s2[2 * j] - p_s4[2 * j]) * (float)i / ni;
          d_cx[m].push_back(pt);
        }
        if (j != nj) {
          for (i = 0; i < ni; i++) {
            count++;
            int jj = 2 * j + 1;
            pt =
              p_s4[jj] + (p_s2[jj] - p_s4[jj]) * (float)(2 * i + 1) / (2 * ni);
            d_cx[m].push_back(pt);
          }
        }
      }
      delete[] p_s2;
      delete[] p_s4;

      // Generate crack elements
      for (j = 0; j < nj; j++) {
        nstart1 = nnode0 + (2 * ni + 1) * j;
        nstart2 = nstart1 + (ni + 1);
        nstart3 = nstart2 + ni;
        for (i = 0; i < ni; i++) {
          // the 1st element
          n1 = nstart2 + i;
          n2 = nstart1 + i;
          n3 = nstart1 + (i + 1);
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 2nd element
          n1 = nstart2 + i;
          n2 = nstart3 + i;
          n3 = nstart1 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 3rd element
          n1 = nstart2 + i;
          n2 = nstart1 + (i + 1);
          n3 = nstart3 + (i + 1);
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 4th element
          n1 = nstart2 + i;
          n2 = nstart3 + (i + 1);
          n3 = nstart3 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
        } // End of loop over j
      } // End of loop over i

      // Collect crack-front nodes
      for (int j = 0; j < 4; j++) { // Loop over sides of the quad
        if (d_cquadCrackSidesAtFront[m][k][j]) {
          for (i = 0; i < (int)d_ce[m].size(); i++) {
            // three element nodes
            n1 = d_ce[m][i].x();
            n2 = d_ce[m][i].y();
            n3 = d_ce[m][i].z();
            if (n1 < nnode0 || n2 < nnode0 || n3 < nnode0) {
              continue;
            }
            for (int s = 0; s < 3; s++) { // Loop over sides of the element
              int sn = n1, en = n2;
              if (s == 1) {
                sn = n2;
                en = n3;
              }
              if (s == 2) {
                sn = n3;
                en = n1;
              }
              if (nodeOnEdge[sn - nnode0][j] && nodeOnEdge[en - nnode0][j]) {
                d_cfSegNodes[m].push_back(sn);
                d_cfSegNodes[m].push_back(en);
              }
            }
          } // End of loop over i
        }
      } // End of loop over j
      nnode0 += num;
      delete[] nodeOnEdge;
    }
  } // End of loop over k
}

void
Crack::DiscretizeQuadCracks(const int& m, int& nnode0)
{
  int k, l, i, j, ni, nj, n1, n2, n3, num;
  int nstart1, nstart2, nstart3;
  double ksi, eta;
  Point pt;

  for (k = 0; k < (int)d_quads[m].size(); k++) {
    for (l = 0; l < (int)d_quadRepetition[m][k]; l++) {
      // Mesh resolutions of the quad
      ni = d_quadN12[m][k];
      nj = d_quadN23[m][k];

      // total number of nodes of the quad
      num = (ni + 1) * (nj + 1) + ni * nj;

      // Flag if node i is on edge j, initialized by NO
      short** nodeOnEdge = scinew short* [num];
      for (i = 0; i < num; i++) {
        nodeOnEdge[i] = scinew short[4];
      }
      for (i = 0; i < num; i++) {
        for (j = 0; j < 4; j++) {
          nodeOnEdge[i][j] = NO;
        }
      }

      // Generate crack nodes
      int count = -1;
      for (j = 0; j <= nj; j++) {
        for (i = 0; i <= ni; i++) {
          // Detect edge nodes
          count++;
          if (j == 0) {
            nodeOnEdge[count][0] = YES;
          }
          if (i == ni) {
            nodeOnEdge[count][1] = YES;
          }
          if (j == nj) {
            nodeOnEdge[count][2] = YES;
          }
          if (i == 0) {
            nodeOnEdge[count][3] = YES;
          }
          // Intrinsic coordinates
          ksi = -1.0 + (float)(2 * i) / ni;
          eta = -1.0 + (float)(2 * j) / nj;
          // Global coordinates by interpolation with shape function
          GetGlobalCoordinatesQuad(m, k, l, ksi, eta, pt);
          d_cx[m].push_back(pt);
        }
        if (j != nj) {
          for (i = 0; i < ni; i++) {
            count++;
            // intrinsic coordinates
            ksi = -1.0 + (float)(2 * i + 1) / ni;
            eta = -1.0 + (float)(2 * j + 1) / nj;
            // Global coordinates
            GetGlobalCoordinatesQuad(m, k, l, ksi, eta, pt);
            d_cx[m].push_back(pt);
          }
        }
      }

      // Generate crack elements
      for (j = 0; j < nj; j++) {
        nstart1 = nnode0 + (2 * ni + 1) * j;
        nstart2 = nstart1 + (ni + 1);
        nstart3 = nstart2 + ni;
        for (i = 0; i < ni; i++) {
          // the 1st element
          n1 = nstart2 + i;
          n2 = nstart1 + i;
          n3 = nstart1 + (i + 1);
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 2nd element
          n1 = nstart2 + i;
          n2 = nstart3 + i;
          n3 = nstart1 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 3rd element
          n1 = nstart2 + i;
          n2 = nstart1 + (i + 1);
          n3 = nstart3 + (i + 1);
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // the 4th element
          n1 = nstart2 + i;
          n2 = nstart3 + (i + 1);
          n3 = nstart3 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
        }
      }

      // Collect crack-front nodes
      for (int j = 0; j < 4; j++) { // Loop over sides of the quad
        if (d_quadCrackSidesAtFront[m][k][j]) {
          for (i = 0; i < (int)d_ce[m].size(); i++) {
            // three element nodes
            n1 = d_ce[m][i].x();
            n2 = d_ce[m][i].y();
            n3 = d_ce[m][i].z();
            if (n1 < nnode0 || n2 < nnode0 || n3 < nnode0) {
              continue;
            }
            for (int s = 0; s < 3; s++) { // Loop over sides of the element
              int sn = n1, en = n2;
              if (s == 1) {
                sn = n2;
                en = n3;
              }
              if (s == 2) {
                sn = n3;
                en = n1;
              }
              if (nodeOnEdge[sn - nnode0][j] && nodeOnEdge[en - nnode0][j]) {
                d_cfSegNodes[m].push_back(sn);
                d_cfSegNodes[m].push_back(en);
              }
            }
          } // End of loop over i
        }
      } // End of loop over j
      nnode0 += num;
      for (int i = 0; i < num; i++) {
        delete nodeOnEdge[i];
      }
      delete[] nodeOnEdge;
    }
  } // End of loop over quads
}

void
Crack::GetGlobalCoordinatesQuad(const int& m,
                                const int& k,
                                const int& l,
                                const double& x,
                                const double& y,
                                Point& pt)
{
  // (x,y): intrinsic coordinates of point "pt".

  // Shape functions of the serendipity eight-noded quadrilateral element
  double sf[8];
  sf[0] = (1. - x) * (1. - y) * (-1. - x - y) / 4.;
  sf[1] = (1. + x) * (1. - y) * (-1. + x - y) / 4.;
  sf[2] = (1. + x) * (1. + y) * (-1. + x + y) / 4.;
  sf[3] = (1. - x) * (1. + y) * (-1. - x + y) / 4.;
  sf[4] = (1. - x * x) * (1. - y) / 2.;
  sf[5] = (1. + x) * (1. - y * y) / 2.;
  sf[6] = (1. - x * x) * (1. + y) / 2.;
  sf[7] = (1. - x) * (1. - y * y) / 2.;

  // Global coordinates of (x,y)
  double px = 0., py = 0., pz = 0.;
  for (int j = 0; j < 8; j++) {
    px += sf[j] * (d_quads[m][k][j].x() + l * d_quadOffset[m][k].x());
    py += sf[j] * (d_quads[m][k][j].y() + l * d_quadOffset[m][k].y());
    pz += sf[j] * (d_quads[m][k][j].z() + l * d_quadOffset[m][k].z());
  }
  pt = Point(px, py, pz);
}

void
Crack::DiscretizeTriangularCracks(const int& m, int& nnode0)
{
  int k, l, i, j;
  int neq, num, nstart1, nstart2, n1 = 0, n2 = 0, n3 = 0;
  Point pt;

  for (k = 0; k < (int)d_triangles[m].size(); k++) {
    for (l = 0; l < (int)d_triRepetition[m][k]; l++) {
      // Mesh resolution of the triangle
      neq = d_triNCells[m][k];

      // total number of nodes of the triangle
      num = (neq + 1) * (neq + 2) / 2;

      // Flag if node 'i' is on edge 'j', initialized by NO
      short** nodeOnEdge = scinew short* [num];
      for (i = 0; i < num; i++) {
        nodeOnEdge[i] = scinew short[3];
      }
      for (i = 0; i < num; i++) {
        for (j = 0; j < 3; j++) {
          nodeOnEdge[i][j] = NO;
        }
      }

      // Generate crack nodes
      int count = -1;
      for (j = 0; j <= neq; j++) {
        for (i = 0; i <= neq - j; i++) {
          // Detect edge nodes
          count++;
          if (j == 0) {
            nodeOnEdge[count][0] = YES;
          }
          if (i + j == neq) {
            nodeOnEdge[count][1] = YES;
          }
          if (i == 0) {
            nodeOnEdge[count][2] = YES;
          }
          // Intrinsic coordinates
          double ksi = (float)i / neq;
          double eta = (float)j / neq;
          // Global coordinates by interpolation with shape function
          GetGlobalCoordinatesTriangle(m, k, l, ksi, eta, pt);
          d_cx[m].push_back(pt);
        }
      }

      // Generate crack elements
      nstart2 = nnode0;
      for (j = 0; j < neq - 1; j++) {
        nstart2 += (neq + 1 - j);
        nstart1 = nstart2 - (neq + 1 - j);
        for (i = 0; i < neq - (j + 1); i++) {
          // left element
          n1 = nstart1 + i;
          n2 = n1 + 1;
          n3 = nstart2 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
          // right element
          n1 = nstart1 + (i + 1);
          n2 = nstart2 + (i + 1);
          n3 = nstart2 + i;
          d_ce[m].push_back(IntVector(n1, n2, n3));
        }
        d_ce[m].push_back(IntVector(n1, n1 + 1, n2));
      }
      d_ce[m].push_back(IntVector(nstart2, nstart2 + 1, nstart2 + 2));

      // Collect crack-front nodes
      for (int j = 0; j < 3; j++) { // Loop over sides of the triangle
        if (d_triCrackSidesAtFront[m][k][j]) {
          for (i = 0; i < (int)d_ce[m].size(); i++) {
            // three nodes of the element
            n1 = d_ce[m][i].x();
            n2 = d_ce[m][i].y();
            n3 = d_ce[m][i].z();
            if (n1 < nnode0 || n2 < nnode0 || n3 < nnode0) {
              continue;
            }
            for (int s = 0; s < 3; s++) { // Loop over sides of the element
              int sn = n1, en = n2;
              if (s == 1) {
                sn = n2;
                en = n3;
              }
              if (s == 2) {
                sn = n3;
                en = n1;
              }
              if (nodeOnEdge[sn - nnode0][j] && nodeOnEdge[en - nnode0][j]) {
                d_cfSegNodes[m].push_back(sn);
                d_cfSegNodes[m].push_back(en);
              }
            }
          } // End of loop over i
        }
      } // End of loop over j
      nnode0 += num;
      delete[] nodeOnEdge;
    } // End of loop over l
  }
}

void
Crack::GetGlobalCoordinatesTriangle(const int& m,
                                    const int& k,
                                    const int& l,
                                    const double& r,
                                    const double& s,
                                    Point& pt)
{
  // (r,s): intrinsic coordinates of point "pt".

  // Shape functions of the serendipity six-noded triangular element
  double sf[6];
  sf[5] = 4. * s * (1. - r - s);
  sf[4] = 4. * r * s;
  sf[3] = 4. * r * (1. - r - s);
  sf[2] = s - 0.5 * (sf[4] + sf[5]);
  sf[1] = r - 0.5 * (sf[3] + sf[4]);
  sf[0] = (1. - r - s) - 0.5 * (sf[3] + sf[5]);

  // Global coordinates of (r,s)
  double px = 0., py = 0., pz = 0.;
  for (int j = 0; j < 6; j++) {
    px += sf[j] * (d_triangles[m][k][j].x() + l * d_triOffset[m][k].x());
    py += sf[j] * (d_triangles[m][k][j].y() + l * d_triOffset[m][k].y());
    pz += sf[j] * (d_triangles[m][k][j].z() + l * d_triOffset[m][k].z());
  }
  pt = Point(px, py, pz);
}

void
Crack::DiscretizeArcCracks(const int& m, int& nnode0)
{
  for (int k = 0; k < (int)d_arcs[m].size(); k++) {
    // Three points of the arc
    Point p1 = d_arcs[m][k][0];
    Point p2 = d_arcs[m][k][1];
    Point p3 = d_arcs[m][k][2];
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    x1 = p1.x();
    y1 = p1.y();
    z1 = p1.z();
    x2 = p2.x();
    y2 = p2.y();
    z2 = p2.z();
    x3 = p3.x();
    y3 = p3.y();
    z3 = p3.z();

    // Find center of the arc
    double a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
    a1 = 2 * (x2 - x1);
    b1 = 2 * (y2 - y1);
    c1 = 2 * (z2 - z1);
    d1 = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2 + z1 * z1 - z2 * z2;
    a2 = 2 * (x3 - x1);
    b2 = 2 * (y3 - y1);
    c2 = 2 * (z3 - z1);
    d2 = x1 * x1 - x3 * x3 + y1 * y1 - y3 * y3 + z1 * z1 - z3 * z3;
    FindPlaneEquation(p1, p2, p3, a3, b3, c3, d3);

    double delt, deltx, delty, deltz;
    delt  = Matrix3(a1, b1, c1, a2, b2, c2, a3, b3, c3).Determinant();
    deltx = Matrix3(-d1, b1, c1, -d2, b2, c2, -d3, b3, c3).Determinant();
    delty = Matrix3(a1, -d1, c1, a2, -d2, c2, a3, -d3, c3).Determinant();
    deltz = Matrix3(a1, b1, -d1, a2, b2, -d2, a3, b3, -d3).Determinant();
    double x0, y0, z0;
    x0            = deltx / delt;
    y0            = delty / delt;
    z0            = deltz / delt;
    Point origin  = Point(x0, y0, z0);
    double radius = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) +
                         (z1 - z0) * (z1 - z0));

    // Define local coordinates
    Vector v1, v2, v3;
    double temp = sqrt(a3 * a3 + b3 * b3 + c3 * c3);
    v3          = Vector(a3 / temp, b3 / temp, c3 / temp);
    v1          = TwoPtsDirCos(origin, p1);
    Vector v31  = Cross(v3, v1);
    v2          = v31 / v31.length();
    double lx, mx, nx, ly, my, ny;
    lx = v1.x();
    mx = v1.y();
    nx = v1.z();
    ly = v2.x();
    my = v2.y();
    ny = v2.z();

    // Angle of the arc
    double angleOfArc;
    double PI = 3.141592654;
    double x3prime, y3prime;
    x3prime         = lx * (x3 - x0) + mx * (y3 - y0) + nx * (z3 - z0);
    y3prime         = ly * (x3 - x0) + my * (y3 - y0) + ny * (z3 - z0);
    double cosTheta = x3prime / radius;
    double sinTheta = y3prime / radius;
    double thetaQ   = std::abs(asin(y3prime / radius));
    if (sinTheta >= 0.) {
      if (cosTheta >= 0) {
        angleOfArc = thetaQ;
      } else {
        angleOfArc = PI - thetaQ;
      }
    } else {
      if (cosTheta <= 0.) {
        angleOfArc = PI + thetaQ;
      } else {
        angleOfArc = 2 * PI - thetaQ;
      }
    }

    // Generate crack nodes
    d_cx[m].push_back(origin);
    for (int j = 0; j <= d_arcNCells[m][k]; j++) {
      double thetai  = angleOfArc * j / d_arcNCells[m][k];
      double xiprime = radius * cos(thetai);
      double yiprime = radius * sin(thetai);
      double xi      = lx * xiprime + ly * yiprime + x0;
      double yi      = mx * xiprime + my * yiprime + y0;
      double zi      = nx * xiprime + ny * yiprime + z0;
      d_cx[m].push_back(Point(xi, yi, zi));
    }

    // Generate crack elements
    for (int j = 1; j <= d_arcNCells[m][k]; j++) {
      int n1 = nnode0;
      int n2 = nnode0 + j;
      int n3 = nnode0 + (j + 1);
      d_ce[m].push_back(IntVector(n1, n2, n3));
      // Crack front nodes
      if (d_arcCrkFrtSegID[m][k] == -1 || d_arcCrkFrtSegID[m][k] == j) {
        d_cfSegNodes[m].push_back(n2);
        d_cfSegNodes[m].push_back(n3);
      }
    }
    nnode0 += d_arcNCells[m][k] + 2;
  }
}

void
Crack::DiscretizeEllipticCracks(const int& m, int& nnode0)
{
  for (int k = 0; k < (int)d_ellipses[m].size(); k++) {
    // Three points of the ellipse
    Point p1 = d_ellipses[m][k][0];
    Point p2 = d_ellipses[m][k][1];
    Point p3 = d_ellipses[m][k][2];

    // Center and half axial lengths of the ellipse
    double x0, y0, z0, a, b;
    Point origin = p3 + (p1 - p3) * 0.5;
    x0           = origin.x();
    y0           = origin.y();
    z0           = origin.z();
    a            = (p1 - origin).length();
    b            = (p2 - origin).length();

    // Local coordinates
    Vector v1, v2, v3;
    v1         = TwoPtsDirCos(origin, p1);
    v2         = TwoPtsDirCos(origin, p2);
    Vector v12 = Cross(v1, v2);
    v3         = v12 / v12.length();
    double lx, mx, nx, ly, my, ny;
    lx = v1.x();
    mx = v1.y();
    nx = v1.z();
    ly = v2.x();
    my = v2.y();
    ny = v2.z();

    // Generate crack nodes
    d_cx[m].push_back(origin);
    for (int j = 0; j < d_ellipseNCells[m][k]; j++) {
      double PI      = 3.141592654;
      double thetai  = j * (2 * PI) / d_ellipseNCells[m][k];
      double xiprime = a * cos(thetai);
      double yiprime = b * sin(thetai);
      double xi      = lx * xiprime + ly * yiprime + x0;
      double yi      = mx * xiprime + my * yiprime + y0;
      double zi      = nx * xiprime + ny * yiprime + z0;
      d_cx[m].push_back(Point(xi, yi, zi));
    }

    // Generate crack elements
    for (int j = 1; j <= d_ellipseNCells[m][k]; j++) {
      int j1 = (j == d_ellipseNCells[m][k] ? 1 : j + 1);
      int n1 = nnode0;
      int n2 = nnode0 + j;
      int n3 = nnode0 + j1;
      d_ce[m].push_back(IntVector(n1, n2, n3));
      // Collect crack-front nodes
      if (d_ellipseCrkFrtSegID[m][k] == -1 || d_ellipseCrkFrtSegID[m][k] == j) {
        d_cfSegNodes[m].push_back(n2);
        d_cfSegNodes[m].push_back(n3);
      }
    }
    nnode0 += d_ellipseNCells[m][k] + 1;
  }
}

void
Crack::DiscretizePartialEllipticCracks(const int& m, int& nnode0)
{
  for (int k = 0; k < (int)d_pellipses[m].size(); k++) {
    double extent = d_pellipseExtent[m][k] / 360.;

    // Center, end points on major and minor axes
    Point origin  = d_pellipses[m][k][0];
    Point major_p = d_pellipses[m][k][1];
    Point minor_p = d_pellipses[m][k][2];
    double x0, y0, z0, a, b;
    x0 = origin.x();
    y0 = origin.y();
    z0 = origin.z();
    a  = (major_p - origin).length();
    b  = (minor_p - origin).length();

    // Local coordinates
    Vector v1, v2, v3;
    v1         = TwoPtsDirCos(origin, major_p);
    v2         = TwoPtsDirCos(origin, minor_p);
    Vector v12 = Cross(v1, v2);
    v3         = v12 / v12.length();
    double lx, mx, nx, ly, my, ny;
    lx = v1.x();
    mx = v1.y();
    nx = v1.z();
    ly = v2.x();
    my = v2.y();
    ny = v2.z();

    // Generate crack nodes
    d_cx[m].push_back(origin);
    for (int j = 0; j <= d_pellipseNCells[m][k]; j++) {
      double PI      = 3.141592654;
      double thetai  = j * (2 * PI * extent) / d_pellipseNCells[m][k];
      double xiprime = a * cos(thetai);
      double yiprime = b * sin(thetai);
      double xi      = lx * xiprime + ly * yiprime + x0;
      double yi      = mx * xiprime + my * yiprime + y0;
      double zi      = nx * xiprime + ny * yiprime + z0;
      d_cx[m].push_back(Point(xi, yi, zi));
    }

    // Generate crack elements
    for (int j = 1; j <= d_pellipseNCells[m][k]; j++) {
      int n1 = nnode0;
      int n2 = nnode0 + j;
      int n3 = nnode0 + j + 1;
      d_ce[m].push_back(IntVector(n1, n2, n3));
      // Collect crack-front nodes
      if (d_pellipseCrkFrtSegID[m][k] == -1 ||
          d_pellipseCrkFrtSegID[m][k] == j) {
        d_cfSegNodes[m].push_back(n2);
        d_cfSegNodes[m].push_back(n3);
      }
    }
    nnode0 += d_pellipseNCells[m][k] + 2;
  }
}

void
Crack::CombineCrackSegments(const int& m)
{
  int num = (int)d_cx[m].size();

  for (int i = 0; i < num; i++) {
    for (int j = 0; j < i; j++) {
      if (TwoPointsCoincide(d_cx[m][i], d_cx[m][j])) {
        ResetCrackNodes(m, i, j);
        ResetCrackElements(m, i, j);
        ResetCrackFrontNodes(m, i, j);
        // After dropping node i, d_cx[m].size() decreses by 1
        num--;
        // The components (>i) in d_cx[m] move back a record.
        // It's tricky, but necessary for multiple coincidences.
        i--;
      }
    }
  }

  // Reorder crack-front nodes
  ReorderCrackFrontNodes(m);

  // Check for errors
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < num; j++) {
      if (i != j && TwoPointsCoincide(d_cx[m][i], d_cx[m][j])) {
        std::cout << "Error: duplicate crack nodes are found: d_cx[" << m
                  << "][" << i << "] = "
                  << "d_cx[" << m << "][" << j << "]. Program terminated."
                  << "\n";
        exit(1);
      }
    }
  }
  for (int i = 0; i < (int)d_ce[m].size(); i++) {
    int n1 = d_ce[m][i].x(), n2 = d_ce[m][i].y(), n3 = d_ce[m][i].z();
    if (n1 == n2 || n1 == n3 || n2 == n3) {
      std::cout << "Error: crack element d_ce[m][" << i << "] = " << d_ce[m][i]
                << " has two same nodes. Program terminted."
                << "\n";
      exit(1);
    }
  }
}

short
Crack::TwoPointsCoincide(const Point& p1, const Point& p2)
{
  double t = 1.e-6;
  return TwoDoublesEqual(p1.x(), p2.x(), t) &&
         TwoDoublesEqual(p1.y(), p2.y(), t) &&
         TwoDoublesEqual(p1.z(), p2.z(), t);
}

short
Crack::TwoDoublesEqual(const double& db1,
                       const double& db2,
                       const double& tolerance)
{
  double ab1 = std::abs(db1);
  double ab2 = std::abs(db2);
  double change;

  if (db1 == db2) {
    return (YES);
  } else if (ab1 > ab2) {
    change = std::abs(db1 - db2) / ab1;
  } else {
    change = std::abs(db1 - db2) / ab2;
  }

  // Equal if different by less than 100 ppm
  if (change < tolerance) {
    return (YES);
  } else {
    if (ab1 < tolerance / 100. && ab2 < tolerance / 100.) {
      return (YES);
    } else {
      return (NO);
    }
  }
}

void
Crack::ResetCrackNodes(const int& m, const int& n1, const int& /*n2*/)
{
  // If d_cx[m][n1]=d_cx[m][n2], drop node n1

  int num    = (int)d_cx[m].size();
  Point* tmp = scinew Point[num];
  for (int i = 0; i < num; i++) {
    tmp[i] = d_cx[m][i];
  }

  d_cx[m].clear();
  d_cx[m].resize(num - 1);
  for (int i = 0; i < num - 1; i++) {
    if (i < n1) {
      d_cx[m][i] = tmp[i];
    } else {
      d_cx[m][i] = tmp[i + 1];
    }
  }
  delete[] tmp;
}

void
Crack::ResetCrackElements(const int& m, const int& n1, const int& n2)
{
  int num        = (int)d_ce[m].size();
  IntVector* tmp = scinew IntVector[num];
  for (int i = 0; i < num; i++) {
    tmp[i] = d_ce[m][i];
  }

  for (int i = 0; i < num; i++) {
    int n = -1;
    // first node of the element
    n = tmp[i].x();
    if (n < n1) {
      d_ce[m][i][0] = n;
    } else if (n == n1) {
      d_ce[m][i][0] = n2;
    } else {
      d_ce[m][i][0] = n - 1;
    }
    // second node of the element
    n = tmp[i].y();
    if (n < n1) {
      d_ce[m][i][1] = n;
    } else if (n == n1) {
      d_ce[m][i][1] = n2;
    } else {
      d_ce[m][i][1] = n - 1;
    }
    // third node of the element
    n = tmp[i].z();
    if (n < n1) {
      d_ce[m][i][2] = n;
    } else if (n == n1) {
      d_ce[m][i][2] = n2;
    } else {
      d_ce[m][i][2] = n - 1;
    }
  }
  delete[] tmp;
}

void
Crack::ResetCrackFrontNodes(const int& m, const int& n1, const int& n2)
{
  int num  = (int)d_cfSegNodes[m].size();
  int* tmp = scinew int[num];
  for (int i = 0; i < num; i++) {
    tmp[i] = d_cfSegNodes[m][i];
  }

  for (int i = 0; i < num; i++) {
    int n = tmp[i];
    if (n < n1) {
      d_cfSegNodes[m][i] = n;
    } else if (n == n1) {
      d_cfSegNodes[m][i] = n2;
    } else {
      d_cfSegNodes[m][i] = n - 1;
    }
  }
  delete[] tmp;
}

void
Crack::ReorderCrackFrontNodes(const int& m)
{
  int k1 = -1, k2 = -1, segs[2];
  std::vector<int> tmp;

  int num = (int)d_cfSegNodes[m].size();
  for (int i = 0; i < num / 2; i++) {
    // two nodes of the crack-front element
    k1 = d_cfSegNodes[m][2 * i];
    k2 = d_cfSegNodes[m][2 * i + 1];
    // element(s) connected by the first node k1
    FindSegsFromNode(m, k1, segs);
    if (segs[R] < 0) { // a right edge element
      tmp.push_back(k1);
      tmp.push_back(k2);
      // the rest elements of the sub crack front
      FindSegsFromNode(m, k2, segs);
      while (segs[L] >= 0) { // Node k2 is connected by segs[L]
        k1 = d_cfSegNodes[m][2 * segs[L]];
        k2 = d_cfSegNodes[m][2 * segs[L] + 1];
        tmp.push_back(k1);
        tmp.push_back(k2);
        FindSegsFromNode(m, k2, segs);
      }
    } // End of if(segs[R]<0)
  }

  if ((int)tmp.size() == 0) { // for enclosed cracks
    // start from the first element
    k1 = d_cfSegNodes[m][0];
    k2 = d_cfSegNodes[m][1];
    tmp.push_back(k1);
    tmp.push_back(k2);
    // the rest elements of the sub crack front
    FindSegsFromNode(m, k2, segs);
    while (segs[L] >= 0 && (int)tmp.size() < num - 1) {
      k1 = d_cfSegNodes[m][2 * segs[L]];
      k2 = d_cfSegNodes[m][2 * segs[L] + 1];
      tmp.push_back(k1);
      tmp.push_back(k2);
      FindSegsFromNode(m, k2, segs);
    }
  }

  // Save the reordered crack-front nodes
  for (int i = 0; i < num; i++) {
    d_cfSegNodes[m][i] = tmp[i];
  }
}

/**
 * @brief Ensures the outer vectors (d_cfSegNodes, d_cfSegPreIdx, etc.) are
 * sized appropriately to contain the given segment index `m_idx`.
 * This prevents out-of-bounds access when accessing d_cfSegNodes[m_idx], etc.
 * @param m_idx The segment index to ensure capacity for.
 */
void
Crack::ensure_segment_vectors_sized(int m_idx)
{
  if (m_idx < 0) {
    // Negative index is invalid; handle as an error or assert
    return;
  }
  // Resize all relevant member vectors if m_idx is beyond their current
  // capacity
  if (static_cast<size_t>(m_idx) >= d_cfSegNodes.size()) {
    d_cfSegNodes.resize(m_idx + 1);
  }
  if (static_cast<size_t>(m_idx) >= d_cfSegPreIdx.size()) {
    d_cfSegPreIdx.resize(m_idx + 1);
  }
  if (static_cast<size_t>(m_idx) >= d_cfSegMinIdx.size()) {
    d_cfSegMinIdx.resize(m_idx + 1);
  }
  if (static_cast<size_t>(m_idx) >= d_cfSegMaxIdx.size()) {
    d_cfSegMaxIdx.resize(m_idx + 1);
  }
}

/**
 * @brief Finds crack-front node indexes for a specific segment `m_idx`.
 * This includes:
 * 1. `d_cfSegPreIdx`: The previous index of a crack-front node
 * with the same value (if any, otherwise -1).
 * 2. `d_cfSegMinIdx` and `d_cfSegMaxIdx`: The minimum and maximum
 * node indexes of the crack-front segment to which the node belongs.
 * @param m_idx The index of the crack-front segment to process.
 */
void
Crack::FindCrackFrontNodeIndexes(const int& m_idx)
{
  // Ensure the outer vectors are sufficiently sized for m_idx
  ensure_segment_vectors_sized(m_idx);
  if (m_idx < 0) {
    return;
  }

  // Get the number of nodes in the current segment
  int num = static_cast<int>(d_cfSegNodes[m_idx].size());

  // --- Part 1: Calculate d_cfSegPreIdx (Previous Node Index) ---
  // For each node, find the index of its last previous occurrence in the
  // segment. Optimized from O(N^2) to average O(N) using a hash map.
  // Initialize d_cfSegPreIdx[m_idx] with -1 and resize to `num`.
  d_cfSegPreIdx[m_idx].assign(num, -1);
  // Map to store the most recent index where each node value was encountered.
  std::unordered_map<int, int> last_seen_index;

  for (int i = 0; i < num; ++i) {
    int current_node_value = d_cfSegNodes[m_idx][i];
    // Check if this node value has been seen before (C++20 `contains`)
    if (last_seen_index.contains(current_node_value)) {
      // If yes, set d_cfSegPreIdx[m_idx][i] to the last seen index.
      d_cfSegPreIdx[m_idx][i] = last_seen_index[current_node_value];
    }
    // Always update the last seen index for the current node value to its
    // current position `i`.
    last_seen_index[current_node_value] = i;
  }

  // --- Part 2: Calculate d_cfSegMinIdx and d_cfSegMaxIdx (Segment Min/Max
  // Indexes) --- These identify the boundaries of logical crack-front
  // segments within d_cfSegNodes[m_idx]. Initialize both vectors with -1 and
  // resize to `num`.
  d_cfSegMinIdx[m_idx].assign(num, -1);
  d_cfSegMaxIdx[m_idx].assign(num, -1);

  // If the segment is empty, there's nothing more to do.
  if (num == 0) {
    return;
  }

  // `current_segment_min_idx` tracks the starting index of the segment
  // currently being defined.
  int current_segment_min_idx = 0;

  // Iterate through each node in the segment.
  // We determine segments block by block. Once a segment's min/max are found,
  // we fill values for all nodes within that segment.
  for (int i = 0; i < num; ++i) {
    // Only search for a new segment's boundaries if `i` is the actual start
    // of a new segment. This prevents redundant calculations for nodes
    // already assigned to a segment.
    if (i == current_segment_min_idx) {
      // `current_segment_max_idx` will store the end index of the current
      // segment. Initialize it to the end of the array, assuming it's a
      // single segment for now.
      int current_segment_max_idx = num - 1;

      // Determine `j_start_for_boundary_check` as per original logic:
      // If `i` is odd, `j_start` is `i`. If `i` is even, `j_start` is `i +
      // 1`. This means `j_start` will always be an odd index, or `num` if `i`
      // is `num-1` (and even).
      int j_start_for_boundary_check = (i % 2 != 0) ? i : (i + 1);

      // Iterate from `j_start_for_boundary_check` with a step of 2 to find
      // segment boundaries. This specific iteration pattern suggests that
      // segment breaks are determined by comparing nodes at `odd_index` and
      // `odd_index + 1`.
      for (int j = j_start_for_boundary_check; j < num; j += 2) {
        // Check if `j` is the last element. If so, the segment ends here.
        if (j == num - 1) {
          current_segment_max_idx = j;
          break; // End of array reached, segment ends here.
        }
        // Check if the node at `j` is different from the node at `j + 1`.
        // If they are different, it signifies a segment boundary.
        if (d_cfSegNodes[m_idx][j] != d_cfSegNodes[m_idx][j + 1]) {
          current_segment_max_idx = j;
          break; // Segment break detected.
        }
      }

      // After the inner loop, `current_segment_max_idx` holds the true end
      // index of the segment that started at `current_segment_min_idx`.

      // Now, assign these determined `minIdx` and `maxIdx` values to all
      // nodes within this identified segment.
      for (int k = current_segment_min_idx; k <= current_segment_max_idx; ++k) {
        d_cfSegMinIdx[m_idx][k] = current_segment_min_idx;
        d_cfSegMaxIdx[m_idx][k] = current_segment_max_idx;
      }

      // Advance `current_segment_min_idx` to the start of the next potential
      // segment.
      current_segment_min_idx = current_segment_max_idx + 1;
    }
    // The outer loop `i` will eventually catch up to
    // `current_segment_min_idx` to start processing the next block.
  }
}

/*
// Determine how the crack-font nodes are connected
void
Crack::FindCrackFrontNodeIndexes(const int& m)
{
// The previous node index of a crack-front node (d_cfSegPreIdx)
// for which node[i]=node[preIdx] (preIdx<i)
d_cfSegPreIdx[m].clear();
int num = (int)d_cfSegNodes[m].size();
d_cfSegPreIdx[m].resize(num);

for (int i = 0; i < num; i++) {
  int preIdx   = -1;
  int thisNode = d_cfSegNodes[m][i];
  for (int j = i - 1; j >= 0; j--) {
    int preNode = d_cfSegNodes[m][j];
    if (thisNode == preNode) {
      preIdx = j;
      break;
    }
  }
  d_cfSegPreIdx[m][i] = preIdx;
}

// The minimum and maximum node indexes of the crack-front
// on which the node resides: d_cfSegMaxIdx and d_cfSegMinIdx
d_cfSegMaxIdx[m].clear();
d_cfSegMinIdx[m].clear();
d_cfSegMaxIdx[m].resize(num);
d_cfSegMinIdx[m].resize(num);

int maxIdx = -1, minIdx = 0;
for (int i = 0; i < num; i++) {
  if (!(i >= minIdx && i <= maxIdx)) {
    for (int j = ((i % 2) != 0 ? i : i + 1); j < num; j += 2) {
      if (j == num - 1 ||
          (j < num - 1 && d_cfSegNodes[m][j] != d_cfSegNodes[m][j + 1])) {
        maxIdx = j;
        break;
      }
    }
  }
  d_cfSegMinIdx[m][i] = minIdx;
  d_cfSegMaxIdx[m][i] = maxIdx;
  if (i == maxIdx) {
    minIdx = maxIdx + 1;
  }
}
}
*/

// Calculate direction cosines of line p1->p2
Vector
Crack::TwoPtsDirCos(const Point& p1, const Point& p2)
{
  Vector v = Vector(0., 0., 0.);

  if (p1 != p2) {
    double l12 = (p1 - p2).length();
    double dx, dy, dz;
    dx = p2.x() - p1.x();
    dy = p2.y() - p1.y();
    dz = p2.z() - p1.z();
    v  = Vector(dx / l12, dy / l12, dz / l12);
  }

  return v;
}

// Smoothe crack front by cubic-spline fit, and then calculate the normals
// Usually, it is not used.
short
Crack::SmoothCrackFrontAndCalculateNormals(const int& mm)
{
  int i = -1, l = -1, k = -1;
  int cfNodeSize = (int)d_cfSegNodes[mm].size();

  // Task 1: Calculate tangential normals at crack-front nodes
  //         by cubic spline fitting, and/or smooth crack front

  short flag = 1;     // Smooth successfully
  double ep  = 1.e-6; // Tolerance

  d_cfSegV3[mm].clear();
  d_cfSegV3[mm].resize(cfNodeSize);

  // Minimum and maximum index of each sub-crack
  int minIdx = -1, maxIdx = -1;
  int minNode = -1, maxNode = -1;
  int numSegs = -1, numPts = -1;
  std::vector<Point> pts;  // Crack-front point subset of the sub-crack
  std::vector<Vector> V3;  // Crack-front point tangential vector
  std::vector<double> dis; // Arc length from the starting point
  std::vector<int> idx;

  for (k = 0; k < cfNodeSize; k++) {
    // Step a: Collect crack points for current sub-crack
    if (k > maxIdx) { // The next sub-crack
      maxIdx = d_cfSegMaxIdx[mm][k];
      minIdx = d_cfSegMinIdx[mm][k];

      // numbers of segments and points of this sub-crack
      minNode = d_cfSegNodes[mm][minIdx];
      maxNode = d_cfSegNodes[mm][maxIdx];
      numSegs = (maxIdx - minIdx + 1) / 2;
      numPts  = numSegs + 1;

      // Allocate memories for the sub-crack
      pts.resize(numPts);
      V3.resize(numPts);
      dis.resize(numPts);
      idx.resize(maxIdx + 1);
    }

    if (k >= minIdx && k <= maxIdx) { // For the sub-crack
      short preIdx = d_cfSegPreIdx[mm][k];
      int ki       = (k - minIdx + 1) / 2;
      if (preIdx < 0 || preIdx == minIdx) {
        pts[ki] = d_cx[mm][d_cfSegNodes[mm][k]];
        // Arc length
        if (k == minIdx) {
          dis[ki] = 0.;
        } else {
          dis[ki] = dis[ki - 1] + (pts[ki] - pts[ki - 1]).length();
        }
      }
      idx[k] = ki;
      if (k < maxIdx) {
        continue; // Collect next points
      }
    }

    // Step b: Define how to smooth the sub-crack
    int n = numPts; // number of points (>=2)
    // int m=(int)(numSegs/2)+2; // number of intervals (>=2)
    int m  = 2; // just two segments
    int n1 = 7 * m - 3;

    // Arries starting from 1
    double* S = new double[n + 1]; // arc-length to the first point
    double* X = new double[n + 1]; // x indexed from 1
    double* Y = new double[n + 1]; // y indexed from 1
    double* Z = new double[n + 1]; // z indexed from 1
    for (i = 1; i <= n; i++) {
      S[i] = dis[i - 1];
      X[i] = pts[i - 1].x();
      Y[i] = pts[i - 1].y();
      Z[i] = pts[i - 1].z();
    }

    int* g     = new int[n + 1];    // segID
    int* j     = new int[m + 1];    // number of points
    double* s  = new double[m + 1]; // positions of intervals
    double* ex = new double[n1 + 1];
    double* ey = new double[n1 + 1];
    double* ez = new double[n1 + 1];

    // Positins of the intervals
    s[1] = S[1] - (S[2] - S[1]) / 50.;
    for (l = 2; l <= m; l++) {
      s[l] = s[1] + (S[n] - s[1]) / m * (l - 1);
    }

    // Number of crack-front nodes of each seg & the segs to which
    // the points belongs
    for (l = 1; l <= m; l++) { // Loop over segs
      j[l] = 0;                // Number of points in the seg
      for (i = 1; i <= n; i++) {
        if ((l < m && S[i] > s[l] && S[i] <= s[l + 1]) ||
            (l == m && S[i] > s[l] && S[i] <= S[n])) {
          j[l]++;   // Number of points in seg l
          g[i] = l; // Seg ID of point i
        }
      }
    }

    // Step c: Smooth the sub-crack points
    if (CubicSpline(n, m, n1, S, X, s, j, ex, ep) &&
        CubicSpline(n, m, n1, S, Y, s, j, ey, ep) &&
        CubicSpline(n, m, n1, S, Z, s, j, ez, ep)) { // Smooth successfully
      for (i = 1; i <= n; i++) {
        l        = g[i];
        double t = 0., dtdS = 0.;
        if (l < m) {
          t    = 2 * (S[i] - s[l]) / (s[l + 1] - s[l]) - 1.;
          dtdS = 2. / (s[l + 1] - s[l]);
        }
        if (l == m) {
          t    = 2 * (S[i] - s[l]) / (S[n] - s[l]) - 1.;
          dtdS = 2. / (S[n] - s[l]);
        }

        double Xv0, Xv1, Xv2, Xv3, Yv0, Yv1, Yv2, Yv3, Zv0, Zv1, Zv2, Zv3;
        Xv0 = ex[7 * l - 6];
        Xv1 = ex[7 * l - 5];
        Xv2 = ex[7 * l - 4];
        Xv3 = ex[7 * l - 3];
        Yv0 = ey[7 * l - 6];
        Yv1 = ey[7 * l - 5];
        Yv2 = ey[7 * l - 4];
        Yv3 = ey[7 * l - 3];
        Zv0 = ez[7 * l - 6];
        Zv1 = ez[7 * l - 5];
        Zv2 = ez[7 * l - 4];
        Zv3 = ez[7 * l - 3];

        double t0 = 1.;
        double t1 = t;
        double t2 = 2 * t * t - 1.;
        double t3 = 4 * t * t * t - 3 * t;
        // double t0p = 0.;
        double t1p = dtdS;
        double t2p = 4 * t * dtdS;
        double t3p = (12. * t * t - 3.) * dtdS;

        V3[i - 1].x(Xv1 * t1p + Xv2 * t2p + Xv3 * t3p);
        V3[i - 1].y(Yv1 * t1p + Yv2 * t2p + Yv3 * t3p);
        V3[i - 1].z(Zv1 * t1p + Zv2 * t2p + Zv3 * t3p);
        pts[i - 1].x(Xv0 * t0 + Xv1 * t1 + Xv2 * t2 + Xv3 * t3);
        pts[i - 1].y(Yv0 * t0 + Yv1 * t1 + Yv2 * t2 + Yv3 * t3);
        pts[i - 1].z(Zv0 * t0 + Zv1 * t1 + Zv2 * t2 + Zv3 * t3);
      }
    } else { // Not smooth successfully, use the raw data
      flag = 0;
      for (i = 0; i < n; i++) {
        Point pt1 = (i == 0 ? pts[i] : pts[i - 1]);
        Point pt2 = (i == n - 1 ? pts[i] : pts[i + 1]);
        V3[i]     = TwoPtsDirCos(pt1, pt2);
      }
    }

    delete[] g;
    delete[] j;
    delete[] s;
    delete[] ex;
    delete[] ey;
    delete[] ez;
    delete[] S;
    delete[] X;
    delete[] Y;
    delete[] Z;

    // Step d: Smooth crack-front points and store tangential vectors
    for (i = minIdx; i <= maxIdx; i++) { // Loop over all nodes on the sub-crack
      int ki = idx[i];
      // Smooth crack-front points
      int ni       = d_cfSegNodes[mm][i];
      d_cx[mm][ni] = pts[ki];

      // Store tangential vectors
      if (minNode == maxNode && (i == minIdx || i == maxIdx)) {
        // for the first and last points (They coincide) of enclosed cracks
        int k1           = idx[minIdx];
        int k2           = idx[maxIdx];
        Vector averageV3 = (V3[k1] + V3[k2]) / 2.;
        d_cfSegV3[mm][i] = -averageV3 / averageV3.length();
      } else {
        d_cfSegV3[mm][i] = -V3[ki] / V3[ki].length();
      }
    }
    pts.clear();
    idx.clear();
    dis.clear();
    V3.clear();
  } // End of loop over k

  // Task 2: Calculate normals of crack plane at crack-front nodes
  d_cfSegV2[mm].clear();
  d_cfSegV2[mm].resize(cfNodeSize);
  for (k = 0; k < cfNodeSize; k++) {
    int node   = d_cfSegNodes[mm][k];
    int preIdx = d_cfSegPreIdx[mm][k];

    if (preIdx < 0) { // Not operated
      Vector v2T       = Vector(0., 0., 0.);
      double totalArea = 0.;
      for (i = 0; i < (int)d_ce[mm].size(); i++) { // Loop over crack elems
        // Three nodes of the elems
        int n1 = d_ce[mm][i].x();
        int n2 = d_ce[mm][i].y();
        int n3 = d_ce[mm][i].z();
        if (node == n1 || node == n2 || node == n3) {
          // Three points of the triangle
          Point p1 = d_cx[mm][n1];
          Point p2 = d_cx[mm][n2];
          Point p3 = d_cx[mm][n3];
          // Lengths of sides of the triangle
          double a = (p1 - p2).length();
          double b = (p1 - p3).length();
          double c = (p2 - p3).length();
          // Half of perimeter of the triangle
          double s = (a + b + c) / 2.;
          // Area of the triangle
          double thisArea = sqrt(s * (s - a) * (s - b) * (s - c));
          // Normal of the triangle
          Vector thisNorm = TriangleNormal(p1, p2, p3);
          // Area-weighted normal vector
          v2T += thisNorm * thisArea;
          // Total area of crack plane related to the node
          totalArea += thisArea;
        }
      } // End of loop over crack elems
      v2T /= totalArea;
      d_cfSegV2[mm][k] = v2T / v2T.length();
    } else { // Calculated
      d_cfSegV2[mm][k] = d_cfSegV2[mm][preIdx];
    }
  } // End of loop over crack-front nodes

  // Task 3: Calculate bi-normals of crack plane at crack-front nodes
  //         and adjust crack-plane normals to make sure the three axes
  //         are perpendicular to each other.
  d_cfSegV1[mm].clear();
  d_cfSegV1[mm].resize(cfNodeSize);
  for (k = 0; k < cfNodeSize; k++) {
    Vector V1        = Cross(d_cfSegV2[mm][k], d_cfSegV3[mm][k]);
    d_cfSegV1[mm][k] = V1 / V1.length();
    Vector V2        = Cross(d_cfSegV3[mm][k], d_cfSegV1[mm][k]);
    d_cfSegV2[mm][k] = V2 / V2.length();
  }

  return flag;
}

short
Crack::CubicSpline(const int& n,
                   const int& m,
                   const int& n1,
                   double x[],
                   double y[],
                   double z[],
                   int j[],
                   double e[],
                   const double& ep)
{
  short flag = 1;
  int i, k, n3, l, j1, nk, lk, llk, jj, lly, nnj, mmi, nn, ii, my, jm, ni, nij;
  double h1, h2, xlk, xlk1, a1, a2, a3, a4, t;

  double** f = new double*[n1 + 1];
  for (i = 0; i < n1 + 1; i++) {
    f[i] = new double[14];
  }

  for (i = 1; i <= n1; i++) {
    e[i] = 0.;
    for (k = 1; k <= 13; k++) {
      f[i][k] = 0.;
    }
  }

  n3 = 0;
  for (l = 1; l <= m; l++) {
    if (l < m) {
      h1 = 1. / (z[l + 1] - z[l]);
    } else {
      h1 = 1. / (x[n] - z[m]);
    }

    j1 = j[l];
    for (k = 1; k <= j1; k++) {
      nk   = n3 + k;
      xlk  = 2. * (x[nk] - z[l]) * h1 - 1.;
      xlk1 = xlk * xlk;
      a1   = 1.;
      a2   = xlk;
      a3   = 2. * xlk1 - 1.;
      a4   = (4. * xlk1 - 3.) * xlk;
      e[7 * l - 6] += a1 * y[nk];
      e[7 * l - 5] += a2 * y[nk];
      e[7 * l - 4] += a3 * y[nk];
      e[7 * l - 3] += a4 * y[nk];
      f[7 * l - 6][7] += a1 * a1;
      f[7 * l - 5][7] += a2 * a2;
      f[7 * l - 4][7] += a3 * a3;
      f[7 * l - 3][7] += a4 * a4;
      f[7 * l - 6][8] += a1 * a2;
      f[7 * l - 5][8] += a2 * a3;
      f[7 * l - 4][8] += a3 * a4;
      f[7 * l - 6][9] += a1 * a3;
      f[7 * l - 5][9] += a2 * a4;
      f[7 * l - 6][10] += a1 * a4;
    }

    f[7 * l - 5][6] = f[7 * l - 6][8];
    f[7 * l - 4][5] = f[7 * l - 6][9];
    f[7 * l - 3][4] = f[7 * l - 6][10];

    f[7 * l - 4][6] = f[7 * l - 5][8];
    f[7 * l - 3][5] = f[7 * l - 5][9];

    f[7 * l - 3][6] = f[7 * l - 4][8];

    f[7 * l - 6][4]  = -0.5;
    f[7 * l - 4][2]  = -0.5;
    f[7 * l - 5][3]  = 0.5;
    f[7 * l - 3][1]  = 0.5;
    f[7 * l - 6][11] = 0.5;
    f[7 * l - 5][10] = 0.5;
    f[7 * l - 4][9]  = 0.5;
    f[7 * l - 3][8]  = 0.5;
    f[7 * l - 5][4]  = -h1;
    f[7 * l - 5][11] = h1;
    f[7 * l - 4][3]  = 4. * h1;
    f[7 * l - 4][10] = 4. * h1;
    f[7 * l - 4][11] = 8. * h1 * h1;
    f[7 * l - 4][4]  = -8. * h1 * h1;
    f[7 * l - 3][2]  = -9. * h1;
    f[7 * l - 3][9]  = 9. * h1;
    f[7 * l - 3][3]  = 48. * h1 * h1;
    f[7 * l - 3][10] = 48. * h1 * h1;

    if (l <= m - 1) {
      if (l < m - 1) {
        h2 = 1. / (z[l + 2] - z[l + 1]);
      } else {
        h2 = 1. / (x[n] - z[m]);
      }

      f[7 * l - 2][3]  = 1.;
      f[7 * l - 2][4]  = 1.;
      f[7 * l - 2][5]  = 1.;
      f[7 * l - 2][6]  = 1.;
      f[7 * l - 2][11] = 1.;
      f[7 * l - 2][13] = 1.;
      f[7 * l - 2][10] = -1.;
      f[7 * l - 2][12] = -1.;
      f[7 * l - 1][3]  = 2. * h1;
      f[7 * l - 1][4]  = 8. * h1;
      f[7 * l - 1][5]  = 18. * h1;
      f[7 * l - 1][10] = -2. * h2;
      f[7 * l - 1][11] = 8. * h2;
      f[7 * l - 1][12] = -18. * h2;
      f[7 * l][3]      = 16. * h1 * h1;
      f[7 * l][4]      = 96. * h1 * h1;
      f[7 * l][10]     = -16. * h2 * h2;
      f[7 * l][11]     = 96. * h2 * h2;
    }
    n3 += j[l];
  }

  lk  = 7;
  llk = lk - 1;
  for (jj = 1; jj <= llk; jj++) {
    lly = lk - jj;
    nnj = n1 + 1 - jj;
    for (i = 1; i <= lly; i++) {
      for (k = 2; k <= 13; k++) {
        f[jj][k - 1] = f[jj][k];
      }
      f[jj][13]   = 0.;
      mmi         = 14 - i;
      f[nnj][mmi] = 0.;
    }
  }

  nn = n1 - 1;
  for (i = 1; i <= nn; i++) {
    k  = i;
    ii = i + 1;
    for (my = ii; my <= lk; my++) {
      if (std::abs(f[my][1]) <= std::abs(f[k][1])) {
        continue;
      }
      k = my;
    }

    if (k != i) {
      t    = e[i];
      e[i] = e[k];
      e[k] = t;
      for (jj = 1; jj <= 13; jj++) {
        t        = f[i][jj];
        f[i][jj] = f[k][jj];
        f[k][jj] = t;
      }
    }

    if (ep >= std::abs(f[i][1])) {
      flag = 0;
      return flag; // unsuccessful
    } else {
      e[i] /= f[i][1];
      for (jj = 2; jj <= 13; jj++) {
        f[i][jj] /= f[i][1];
      }

      ii = i + 1;
      for (my = ii; my <= lk; my++) {
        t = f[my][1];
        e[my] -= t * e[i];
        for (jj = 2; jj <= 13; jj++) {
          f[my][jj - 1] = f[my][jj] - t * f[i][jj];
        }
        f[my][13] = 0.;
      }

      if (lk == n1) {
        continue;
      }
      lk++;
    }
  }

  e[n1] /= f[n1][1];
  jm = 2;
  nn = n1 - 1;
  for (i = 1; i <= nn; i++) {
    ni = n1 - i;
    for (jj = 2; jj <= jm; jj++) {
      nij = ni - 1 + jj;
      e[ni] -= f[ni][jj] * e[nij];
    }
    if (jm == 13) {
      continue;
    }
    jm++;
  }

  return flag;
}

// Calculate crack-front normals by weighted average method,
// which is used by defualt
void
Crack::CalculateCrackFrontNormals(const int& mm)
{
  // Task 1: Calculate crack-front tangential normals
  int num = d_cfSegNodes[mm].size();
  d_cfSegV3[mm].clear();
  d_cfSegV3[mm].resize(num);
  for (int k = 0; k < num; k++) {
    int node = d_cfSegNodes[mm][k];
    Point pt = d_cx[mm][node];

    int preIdx = d_cfSegPreIdx[mm][k];
    if (preIdx < 0) { // a duplicate node, not operated
      int minIdx  = d_cfSegMinIdx[mm][k];
      int maxIdx  = d_cfSegMaxIdx[mm][k];
      int minNode = d_cfSegNodes[mm][minIdx];
      int maxNode = d_cfSegNodes[mm][maxIdx];

      // The node to the right of pt
      int node1 = -1;
      if (minNode == maxNode && (k == minIdx || k == maxIdx)) {
        // for the ends of enclosed crack
        node1 = d_cfSegNodes[mm][maxIdx - 2];
      } else { // for the nodes of non-enclosd cracks or
               // non-end-nodes of enclosed cracks
        int k1 = -1;
        if ((maxIdx - minIdx + 1) / 2 > 2) { // three or more segments
          k1 = (k - 2) < minIdx + 1 ? minIdx + 1 : k - 2;
        } else { // one or two segments
          k1 = (k - 2) < minIdx ? minIdx : k - 2;
        }
        node1 = d_cfSegNodes[mm][k1];
      }
      Point pt1 = d_cx[mm][node1];

      // The node to the left of pt
      int node2 = -1;
      if (minNode == maxNode && (k == minIdx || k == maxIdx)) {
        // for the ends of enclosed crack
        node2 = d_cfSegNodes[mm][minIdx + 2];
      } else { // for the nodes of non-enclosd cracks or
               // non-end-nodes of enclosed cracks
        int k2 = -1;
        if ((maxIdx - minIdx + 1) > 2) { // Three or more segments
          k2 = (k + 2) > maxIdx - 1 ? maxIdx - 1 : k + 2;
        } else { // one or two segments
          k2 = (k + 2) > maxIdx ? maxIdx : k + 2;
        }
        node2 = d_cfSegNodes[mm][k2];
      }
      Point pt2 = d_cx[mm][node2];

      // Weighted tangential vector between pt1->pt->pt2
      double l1  = (pt1 - pt).length();
      double l2  = (pt - pt2).length();
      Vector v1  = (l1 == 0. ? Vector(0., 0., 0.) : TwoPtsDirCos(pt1, pt));
      Vector v2  = (l2 == 0. ? Vector(0., 0., 0.) : TwoPtsDirCos(pt, pt2));
      Vector v3T = (l1 * v1 + l2 * v2) / (l1 + l2);
      d_cfSegV3[mm][k] = -v3T / v3T.length();
    } else { // calculated
      d_cfSegV3[mm][k] = d_cfSegV3[mm][preIdx];
    }
  } // End of loop over k

  // Reset tangential vectors for edge nodes outside material
  // to the values of the nodes next to them.
  // This way will eliminate the effect of edge nodes.
  for (int k = 0; k < num; k++) {
    int node = d_cfSegNodes[mm][k];
    int segs[2];
    FindSegsFromNode(mm, node, segs);
    if (segs[R] < 0) {
      d_cfSegV3[mm][k] = d_cfSegV3[mm][k + 1];
    }
    if (segs[L] < 0) {
      d_cfSegV3[mm][k] = d_cfSegV3[mm][k - 1];
    }
  }

  // Task 2: Calculate normals of crack plane at crack-front nodes
  d_cfSegV2[mm].clear();
  d_cfSegV2[mm].resize(num);
  for (int k = 0; k < num; k++) {
    int node = d_cfSegNodes[mm][k];

    // Detect if it is an edge node
    int segs[2];
    FindSegsFromNode(mm, node, segs);
    short edgeNode = NO;
    if (segs[L] < 0 || segs[R] < 0) {
      edgeNode = YES;
    }

    // End information of the sub crack
    int minNode = d_cfSegNodes[mm][d_cfSegMinIdx[mm][k]];
    int maxNode = d_cfSegNodes[mm][d_cfSegMaxIdx[mm][k]];

    int preIdx = d_cfSegPreIdx[mm][k];
    if (preIdx < 0) { // a duplicate node, not operated
      Vector v2T       = Vector(0., 0., 0.);
      double totalArea = 0.;
      for (int i = 0; i < (int)d_ce[mm].size(); i++) {
        // Three nodes of the elems
        int n1 = d_ce[mm][i].x();
        int n2 = d_ce[mm][i].y();
        int n3 = d_ce[mm][i].z();

        // Detect if the elem is connected to the node
        short elemRelatedToNode = NO;
        if (node == n1 || node == n2 || node == n3) {
          elemRelatedToNode = YES;
        }

        // Detect if the elem is an inner elem
        short innerElem = YES;
        if (minNode != maxNode &&
            (n1 == minNode || n2 == minNode || n3 == minNode || n1 == maxNode ||
             n2 == maxNode || n3 == maxNode)) {
          innerElem = NO;
        }

        // The elem will be used if it is connected to the node AND
        // if the node is an edge node or the elem is an interior element.
        if (elemRelatedToNode && (innerElem || edgeNode)) {
          // Three points of the triangle
          Point p1 = d_cx[mm][n1];
          Point p2 = d_cx[mm][n2];
          Point p3 = d_cx[mm][n3];
          // Lengths of sides of the triangle
          double a = (p1 - p2).length();
          double b = (p1 - p3).length();
          double c = (p2 - p3).length();
          // Half of perimeter of the triangle
          double s = (a + b + c) / 2.;
          // Area of the triangle
          double thisArea = sqrt(s * (s - a) * (s - b) * (s - c));
          // Normal of the triangle
          Vector thisNorm = TriangleNormal(p1, p2, p3);
          // Area-weighted normal vector
          v2T += thisNorm * thisArea;
          // Total area of crack plane related to the node
          totalArea += thisArea;
        }
      } // End of loop over crack elems

      if (totalArea != 0.) {
        v2T /= totalArea;
      } else {
        std::cout << "Error: divided by zero in calculating outer normal"
                  << " at crack front node, cfSegNodes[" << mm << "][" << k
                  << "] = " << d_cfSegNodes[mm][k] << "\n";
        exit(1);
      }

      d_cfSegV2[mm][k] = v2T / v2T.length();
    } else { // Calculated
      d_cfSegV2[mm][k] = d_cfSegV2[mm][preIdx];
    }
  } // End of loop over crack-front nodes

  // Reset normals for edge nodes outside material
  // to the values of the nodes next to them.
  // This way will eliminate the effect of edge nodes.
  for (int k = 0; k < num; k++) {
    int node = d_cfSegNodes[mm][k];
    int segs[2];
    FindSegsFromNode(mm, node, segs);
    if (segs[R] < 0) {
      d_cfSegV2[mm][k] = d_cfSegV2[mm][k + 1];
    }
    if (segs[L] < 0) {
      d_cfSegV2[mm][k] = d_cfSegV2[mm][k - 1];
    }
  }

  // Task 3: Calculate bi-normals of crack plane at crack-front nodes
  // and adjust tangential normals to make sure the three axes
  // are perpendicular to each other. Keep V2 unchanged.
  d_cfSegV1[mm].clear();
  d_cfSegV1[mm].resize(num);
  for (int k = 0; k < num; k++) {
    Vector V1        = Cross(d_cfSegV2[mm][k], d_cfSegV3[mm][k]);
    d_cfSegV1[mm][k] = V1 / V1.length();
    Vector V3        = Cross(d_cfSegV1[mm][k], d_cfSegV2[mm][k]);
    d_cfSegV3[mm][k] = V3 / V3.length();
  }
}

// Find the segment numbers which are connected by the same node
void
Crack::FindSegsFromNode(const int& m, const int& node, int segs[])
{
  // segs[R] -- the segment on the right of the node
  // segs[L] -- the segment on the left of the node
  segs[R] = segs[L] = -1;

  int ncfSegs = (int)d_cfSegNodes[m].size() / 2;
  for (int j = 0; j < ncfSegs; j++) {
    int node0 = d_cfSegNodes[m][2 * j];
    int node1 = d_cfSegNodes[m][2 * j + 1];
    if (node == node1) { // the right seg
      segs[R] = j;
    }
    if (node == node0) { // the left seg
      segs[L] = j;
    }
  } // End of loop over j

  // See if reasonable
  if (segs[R] < 0 && segs[L] < 0) {
    std::cout << "Error: failure to find the crack-front segments for node "
              << node << ". Program terminated."
              << "\n";
    exit(1);
  }
}

// Calculate outer normal of a triangle
Vector
Crack::TriangleNormal(const Point& p1, const Point& p2, const Point& p3)
{
  double x21, x31, y21, y31, z21, z31;
  double a, b, c;
  Vector norm;

  x21 = p2.x() - p1.x();
  x31 = p3.x() - p1.x();
  y21 = p2.y() - p1.y();
  y31 = p3.y() - p1.y();
  z21 = p2.z() - p1.z();
  z31 = p3.z() - p1.z();

  a = y21 * z31 - z21 * y31;
  b = x31 * z21 - z31 * x21;
  c = x21 * y31 - y21 * x31;

  if (Vector(a, b, c).length() > 1.e-16) {
    norm = Vector(a, b, c) / Vector(a, b, c).length();
  } else {
    norm = Vector(a, b, c);
  }

  return norm;
}

void
Crack::OutputInitialCrackMesh(const int& m)
{
  int pid;
  MPI_Comm_rank(d_mpi_crack_comm, &pid);
  if (pid == 0) { // Output from the first rank
    std::cout << "\n---Initial Crack mesh---"
              << "\n";
    std::cout << "MatID: " << m << "\n";
    std::cout << "  Number of crack elements: " << (int)d_ce[m].size()
              << "\n  Number of crack nodes: " << (int)d_cx[m].size()
              << "\n  Number of crack-front elements: "
              << (int)d_cfSegNodes[m].size() / 2 << "\n";

    std::cout << "  Crack elements (" << (int)d_ce[m].size()
              << " elements in total):"
              << "\n";
    for (int i = 0; i < (int)d_ce[m].size(); i++) {
      std::cout << "     Elem " << i << ": " << d_ce[m][i] << "\n";
    }

    std::cout << "  Crack nodes (" << (int)d_cx[m].size() << " nodes in total):"
              << "\n";
    for (int i = 0; i < (int)d_cx[m].size(); i++) {
      std::cout << "     Node " << i << ": " << d_cx[m][i] << "\n";
    }

    std::cout << "  Crack-front elements and normals ("
              << (int)d_cfSegNodes[m].size() / 2 << " elements in total)"
              << "\n";
    std::cout << "     V1: bi-normal; V2: outer normal; V3: tangential normal."
              << "\n";
    for (int i = 0; i < (int)d_cfSegNodes[m].size(); i++) {
      std::cout << "     Seg " << i / 2 << ": node " << d_cfSegNodes[m][i]
                << d_cx[m][d_cfSegNodes[m][i]] << ", V1: " << d_cfSegV1[m][i]
                << ", V2: " << d_cfSegV2[m][i] << ", V3: " << d_cfSegV3[m][i]
                << "\n";
      if (i % 2 != 0) {
        std::cout << "\n";
      }
    }

    std::cout << "  Average length of crack-front segments, css[m]=" << d_css[m]
              << "\n";
    std::cout << "  Average angle of crack-front segments, csa[m]=" << d_csa[m]
              << " degree."
              << "\n";
    std::cout << "  Crack extent: " << d_cmin[m] << "-->" << d_cmax[m] << "\n"
              << "\n";
  }
}
