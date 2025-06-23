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
    Crack.h
    Created by Yajun Guo in 2002-2005.
********************************************************************************/
#ifndef UINTAH_HOMEBREW_CRACK_H
#define UINTAH_HOMEBREW_CRACK_H

#include <CCA/Components/MPM/Crack/CrackGeometry.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Math/Matrix3.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iomanip>
#include <iostream>

using std::string;
using std::vector;
namespace Uintah {

using CrackElements = std::vector<Uintah::IntVector>;
using CrackPolygon = std::vector<Uintah::Point>;
using CrackOffsets = std::vector<Uintah::Vector>;
using CrackEdgeIndices = std::vector<short>;
using CrackCellIndices = std::vector<int>;
using CrackNodeIndices = std::vector<int>;
using CrackSegmentIndices = std::vector<int>;
using CrackScalarArray = std::vector<double>;

using CrackCoordsVec = std::vector<CrackPolygon>;
using CrackEdgeIndicesVec = std::vector<CrackEdgeIndices>;
using CrackNodeIndicesVec = std::vector<CrackNodeIndices>;

class DataWarehouse;
class MPMLabel;
class MPMFlags;
class ProcessorGroup;
class Patch;
class VarLabel;
class Task;

class Crack
{
public:
  // Constructor
  Crack(const ProblemSpecP& ps,
        MaterialManagerP& d_sS,
        Output* dataArchiver,
        MPMLabel* lb,
        MPMFlags* MFlag);

  // Destructor
  ~Crack() = default;

  void
  outputProblemSpec(ProblemSpecP& mpm_ps, ProblemSpecP& uda_ps);

  // Public methods in Crack.cc
  void
  addComputesAndRequiresCrackDiscretization(Task* task,
                                            const PatchSet* patches,
                                            const MaterialSet* matls) const;
  void
  CrackDiscretization(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  // Public methods in ParticleNodePairVelocityField.cc
  void
  addComputesAndRequiresParticleVelocityField(Task* task,
                                              const PatchSet* patches,
                                              const MaterialSet* matls) const;
  void
  ParticleVelocityField(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  // Public methods in CrackSurfaceContact.cc
  void
  addComputesAndRequiresAdjustCrackContactInterpolated(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  AdjustCrackContactInterpolated(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw);
  void
  addComputesAndRequiresAdjustCrackContactIntegrated(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  AdjustCrackContactIntegrated(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw);

  // Public methods in FractureParametersCalculation.cc
  void
  addComputesAndRequiresGetNodalSolutions(Task* task,
                                          const PatchSet* patches,
                                          const MaterialSet* matls) const;
  void
  GetNodalSolutions(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);
  void
  addComputesAndRequiresCalculateFractureParameters(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  CalculateFractureParameters(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  // Public methods in CrackPropagation.cc
  void
  addComputesAndRequiresPropagateCrackFrontPoints(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  PropagateCrackFrontPoints(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* matls,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);
  void
  addComputesAndRequiresConstructNewCrackFrontElems(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  ConstructNewCrackFrontElems(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  // Public methods in MoveCracks.cc
  void
  addComputesAndRequiresCrackPointSubset(Task* task,
                                         const PatchSet* patches,
                                         const MaterialSet* matls) const;
  void
  CrackPointSubset(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);
  void
  addComputesAndRequiresMoveCracks(Task* task,
                                   const PatchSet* patches,
                                   const MaterialSet* matls) const;
  void
  MoveCracks(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse* old_dw,
             DataWarehouse* new_dw);

  // Public methods in UpdateCrackFront.cc
  void
  addComputesAndRequiresCrackFrontNodeSubset(Task* task,
                                             const PatchSet* patches,
                                             const MaterialSet* matls) const;
  void
  CrackFrontNodeSubset(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);
  void
  addComputesAndRequiresRecollectCrackFrontSegments(
    Task* task,
    const PatchSet* patches,
    const MaterialSet* matls) const;
  void
  RecollectCrackFrontSegments(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

private:

  // PRIVATE METHODS
  // Private methods in Crack.cc
  void
  ReadQuadCracks(const int&, const ProblemSpecP&);
  void
  ReadCurvedQuadCracks(const int&, const ProblemSpecP&);
  void
  ReadTriangularCracks(const int&, const ProblemSpecP&);
  void
  ReadArcCracks(const int&, const ProblemSpecP&);
  void
  ReadEllipticCracks(const int&, const ProblemSpecP&);
  void
  ReadPartialEllipticCracks(const int&, const ProblemSpecP&);
  void
  OutputInitialCrackPlane(const int&);
  void
  DiscretizeQuadCracks(const int&, int&);
  void
  DiscretizeCurvedQuadCracks(const int&, int&);
  void
  GetGlobalCoordinatesQuad(const int&,
                           const int&,
                           const int&,
                           const double&,
                           const double&,
                           Point&);
  void
  GetGlobalCoordinatesTriangle(const int&,
                               const int&,
                               const int&,
                               const double&,
                               const double&,
                               Point&);
  void
  DiscretizeTriangularCracks(const int&, int&);
  void
  DiscretizeArcCracks(const int&, int&);
  void
  DiscretizeEllipticCracks(const int&, int&);
  void
  DiscretizePartialEllipticCracks(const int&, int&);
  void
  CombineCrackSegments(const int&);
  short
  TwoPointsCoincide(const Point&, const Point&);
  short
  TwoDoublesEqual(const double&, const double&, const double&);
  void
  ResetCrackNodes(const int&, const int&, const int&);
  void
  ResetCrackElements(const int&, const int&, const int&);
  void
  ResetCrackFrontNodes(const int&, const int&, const int&);
  void
  ReorderCrackFrontNodes(const int&);
  void
  FindCrackFrontNodeIndexes(const int&);
  short
  SmoothCrackFrontAndCalculateNormals(const int& m);
  short
  CubicSpline(const int& n,
              const int& m,
              const int& n1,
              double[],
              double[],
              double[],
              int[],
              double[],
              const double&);
  void
  CalculateCrackFrontNormals(const int& m);
  void
  FindSegsFromNode(const int&, const int&, int[]);
  void
  OutputInitialCrackMesh(const int&);
  Vector
  TwoPtsDirCos(const Point&, const Point&);
  Vector
  TriangleNormal(const Point&, const Point&, const Point&);

  // Private methods in ParticleNodePairVelocityField.cc
  IntVector
  CellOffset(const Point&, const Point&, Vector);
  short
  ParticleNodeCrackPLaneRelation(const Point&,
                                 const Point&,
                                 const Point&,
                                 const Point&,
                                 const Point&);
  double
  Volume(const Point&, const Point&, const Point&, const Point&);

  // Private methods in FractureParametersCalculation.cc
  void
  DetectIfDoingFractureAnalysisAtThisTimestep(double);
  void
  FindJIntegralPath(const Point&,
                    const Vector&,
                    const Vector&,
                    const Vector&,
                    double[]);
  bool
  FindIntersectionJPathAndCrackPlane(const int&,
                                     const double& r,
                                     const double[],
                                     Point&);
  void
  FindPlaneEquation(const Point&,
                    const Point&,
                    const Point&,
                    double&,
                    double&,
                    double&,
                    double&);
  short
  PointInTriangle(const Point&, const Point&, const Point&, const Point&);
  void
  GetPositionToComputeCOD(const int&, const Point&, const Matrix3&, double&);
  void
  OutputCrackFrontResults(const int&, double time, double delT);

  // Private methods in CrackPropagation.cc
  void
  TrimLineSegmentWithBox(const Point&, Point&, const Point&, const Point&);
  void
  PruneCrackFrontAfterPropagation(const int& m, const double& ca);

  // Private methods in MoveCracks.cc
  short
  PhysicalGlobalGridContainsPoint(const double&, const Point&);
  void
  ApplySymmetricBCsToCrackPoints(const Vector&, const Point&, Point&);

  // Private methods in UpdateCrackFront.cc
  void
  OutputCrackGeometry(const int&, const int&);

protected:
  MPMLabel* lb;

private:

  // PRIVATE DATA MEMBERS
  MPI_Comm d_mpi_crack_comm;
  MaterialManagerP d_mat_manager;
  MPMFlags* d_flag;
  int d_n8or27;
  int d_NGP;
  int d_NGN;
  enum
  {
    NO = 0,
    YES
  }; // No (NO=0) or Yes (YES=1)
  enum
  {
    R = 0,
    L
  };                                  // Right (R=0) or left (L=1)
  Output* d_dataArchiver;               // Data archiving information
  string d_udaDir;                      // Base file directory
  string d_GridBCType[Patch::numFaces]; // BC types of global grid
  Point d_GLP, d_GHP;            // Lowest and highest pt of real global grid
  int d_NJ;                    // rJ = NJ*min_dCell
  double d_rJ;                 // NJ = rJ/min_dCell
  double d_rdadx;              // Ratio of crack incremental to cell-size
  double d_computeJKInterval;  // Interval of calculating fracture parameters
                             // zero by default (every time step)
  double d_growCrackInterval;  // Interval of crack propagation
                             // zero by default (every time step)
  bool d_useVolumeIntegral;    // Use volume integral in J, "no" by default
  bool d_saveCrackGeometry;    // Save crack geometry, "yes" by default
  bool d_smoothCrackFront;     // Smoothe crack-front, "no" by default
  bool d_calFractParameters; // Calculate J or K, "no" by default
  bool d_doCrackPropagation; // Do crack propagation, "no" by default
  bool d_calFractParametersStep;   // Calculate J or K at this step
  bool d_doCrackPropagationStep;   // Do crack propagation at this step
  int d_CODOption;             // CODOption=0 (by default):
                             //   calculate COD at a fixed location;
                 // CODOption=1: calculate COD at the farthest position
                 //   on the crack element at crack-front;
                 // CODOption=2: calculate COD at the intersection
                 //   between J-integral contour and crack plane;

  // Physical parameters of cracks
  std::vector<std::string> d_stressState; // Crack front stress state
  std::vector<std::string> d_crackType;   // Crack contact type
  std::vector<double> d_cmu;              // Crack surface friction coefficient

  std::vector<CrackGeometry*> d_crackGeometry;

  // Geometrical parameters of crack segments
  std::vector<CrackCoordsVec> d_quads;
  std::vector<CrackCellIndices> d_quadN12, d_quadN23;
  std::vector<CrackEdgeIndicesVec> d_quadCrackSidesAtFront;
  std::vector<CrackSegmentIndices> d_quadRepetition;
  std::vector<CrackOffsets> d_quadOffset;
  std::vector<CrackCoordsVec> d_cquads;
  std::vector<CrackCellIndices> d_cquadNStraightSides;
  std::vector<CrackCoordsVec> d_cquadPtsSide2;
  std::vector<CrackCoordsVec> d_cquadPtsSide4;
  std::vector<CrackEdgeIndicesVec> d_cquadCrackSidesAtFront;
  std::vector<CrackSegmentIndices> d_cquadRepetition;
  std::vector<CrackOffsets> d_cquadOffset;
  std::vector<CrackCoordsVec> d_triangles;
  std::vector<CrackCellIndices> d_triNCells;
  std::vector<CrackEdgeIndicesVec> d_triCrackSidesAtFront;
  std::vector<CrackSegmentIndices> d_triRepetition;
  std::vector<CrackOffsets> d_triOffset;
  std::vector<CrackCoordsVec> d_arcs;
  std::vector<CrackCellIndices> d_arcNCells;
  std::vector<CrackSegmentIndices> d_arcCrkFrtSegID;
  std::vector<CrackCoordsVec> d_ellipses;
  std::vector<CrackCellIndices> d_ellipseNCells;
  std::vector<CrackSegmentIndices> d_ellipseCrkFrtSegID;
  std::vector<CrackCoordsVec> d_pellipses;
  std::vector<CrackCellIndices> d_pellipseNCells;
  std::vector<CrackSegmentIndices> d_pellipseCrkFrtSegID;
  std::vector<CrackScalarArray> d_pellipseExtent;
  CrackPolygon d_cmin, d_cmax;

  // Crack data after mesh
  CrackScalarArray d_css;             // Average length of crack-front segments
  CrackScalarArray d_csa;             // Average angle of crack-front segments
  std::vector<CrackPolygon> d_cx;       // Coordinates of crack nodes
  std::vector<CrackElements> d_ce;   // Crack elements
  std::vector<CrackNodeIndices> d_cfSegNodes; // Crack-front nodes
  std::vector<CrackScalarArray> d_cfSegVel;  // Velocity of crack-front nodes
  std::vector<CrackScalarArray> d_cfSegTime; // Time instant of crack propagation
  std::vector<CrackScalarArray> d_cfSegDis;  // Crack incremental
  std::vector<CrackNodeIndices> d_cfSegPreIdx;  // node[i]=node[preIdx]
  std::vector<CrackNodeIndices> d_cfSegMinIdx;  // Minimum node-index of the sub-crack
  std::vector<CrackNodeIndices> d_cfSegMaxIdx;  // Maximum node-index of the sub-crack
  std::vector<CrackPolygon> d_cfSegPtsT;  // Crack-front points after propagation
  std::vector<CrackOffsets> d_cfSegV1;   // Bi-normals at crack-front nodes
  std::vector<CrackOffsets> d_cfSegV2;   // Outer normals at crack-front nodes
  std::vector<CrackOffsets> d_cfSegV3;   // Tangential normals at crack-front nodes
  std::vector<CrackOffsets> d_cfSegJ; // J-integral at crack-front nodes
  std::vector<CrackOffsets> d_cfSegK; // SIF at crack-front nodes
  std::vector<CrackNodeIndicesVec> d_cnset;  // Crack-node subset in each patch
  std::vector<CrackNodeIndicesVec> d_cfnset; // Crack-front-node index subset
  std::vector<CrackNodeIndicesVec> d_cfsset; // Crack-front-seg subset in each patch
};

} // namespace Uintah

#endif /* __CRACK_H__*/
