/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

/********************************************************************************
    Crack.h
    Created by Yajun Guo in 2002-2005.
********************************************************************************/
#ifndef UINTAH_HOMEBREW_CRACK_H
#define UINTAH_HOMEBREW_CRACK_H

#include <iomanip>
#include <iostream>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Math/Matrix3.h>
#include <CCA/Ports/Output.h>
#include <CCA/Components/MPM/Crack/CrackGeometry.h>

namespace Uintah {
   using namespace Uintah;
   using std::vector;
   using std::string;

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
    Crack(const ProblemSpecP& ps, MaterialManagerP& d_sS,
          Output* dataArchiver,
          MPMLabel* lb,MPMFlags* MFlag);

    // Destructor
     ~Crack();

    // Public methods in ReadAndDiscretizeCracks.cc
    void addComputesAndRequiresCrackDiscretization(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;  
    void CrackDiscretization(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

    // Public methods in ParticleNodePairVelocityField.cc
    void addComputesAndRequiresParticleVelocityField(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void ParticleVelocityField(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

    // Public methods in CrackSurfaceContact.cc
    void addComputesAndRequiresAdjustCrackContactInterpolated(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void AdjustCrackContactInterpolated(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    void addComputesAndRequiresAdjustCrackContactIntegrated(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void AdjustCrackContactIntegrated(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

    // Public methods in FractureParametersCalculation.cc
    void addComputesAndRequiresGetNodalSolutions(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void GetNodalSolutions(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    void addComputesAndRequiresCalculateFractureParameters(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void CalculateFractureParameters(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

    // Public methods in CrackPropagation.cc
    void addComputesAndRequiresPropagateCrackFrontPoints(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void PropagateCrackFrontPoints(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    void addComputesAndRequiresConstructNewCrackFrontElems(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void ConstructNewCrackFrontElems(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

    // Public methods in MoveCracks.cc
    void addComputesAndRequiresCrackPointSubset(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void CrackPointSubset(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    void addComputesAndRequiresMoveCracks(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void MoveCracks(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    
    // Public methods in UpdateCrackFront.cc
    void addComputesAndRequiresCrackFrontNodeSubset(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void CrackFrontNodeSubset(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);
    void addComputesAndRequiresRecollectCrackFrontSegments(Task* task,
                                const PatchSet* patches,
                                const MaterialSet* matls) const;
    void RecollectCrackFrontSegments(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

 private:

    // PRIVATE DATA MEMBERS
    MPI_Comm mpi_crack_comm;
    MaterialManagerP 
 d_mat_manager;
    MPMFlags* flag;
    int n8or27;
    int NGP;
    int NGN;
    enum {NO=0,YES};               // No (NO=0) or Yes (YES=1)
    enum {R=0, L};                 // Right (R=0) or left (L=1)
    Output* dataArchiver;          // Data archiving information
    string udaDir;                 // Base file directory
    string 
    GridBCType[Patch::numFaces];   // BC types of global grid
    Point GLP, GHP;                // Lowest and highest pt of real global grid 
    int    NJ;                     // rJ = NJ*min_dCell
    double rJ;                     // NJ = rJ/min_dCell
    double rdadx;                  // Ratio of crack incremental to cell-size
    double computeJKInterval;      // Interval of calculating fracture parameters
                                   // zero by default (every time step) 
    double growCrackInterval;      // Interval of crack propagation
                                   // zero by default (every time step)     
    bool useVolumeIntegral;        // Use volume integral in J, "no" by default
    bool saveCrackGeometry;        // Save crack geometry, "yes" by default
    bool smoothCrackFront;         // Smoothe crack-front, "no" by default
    bool d_calFractParameters;     // Calculate J or K, "no" by default
    bool d_doCrackPropagation;     // Do crack propagation, "no" by default
    bool calFractParameters;       // Calculate J or K at this step
    bool doCrackPropagation;       // Do crack propagation at this step
    int  CODOption;                // CODOption=0 (by default): 
                                   //   calculate COD at a fixed location; 
                                   // CODOption=1: calculate COD at the farthest position 
                                   //   on the crack element at crack-front;
                                   // CODOption=2: calculate COD at the intersection 
                                   //   between J-integral contour and crack plane; 

    // Physical parameters of cracks
    std::vector<std::string> stressState;    // Crack front stress state 
    std::vector<std::string> crackType;      // Crack contact type
    std::vector<double> cmu;            // Crack surface friction coefficient

    std::vector<CrackGeometry*> d_crackGeometry;

    // Geometrical parameters of crack segments
    std::vector<vector<vector<Point> > > quads;
    std::vector<vector<int> >            quadN12,quadN23;       
    std::vector<vector<vector<short> > > quadCrackSidesAtFront;
    std::vector<vector<int> >            quadRepetition;
    std::vector<vector<Vector> >         quadOffset;    
    std::vector<vector<vector<Point> > > cquads;
    std::vector<vector<int> >            cquadNStraightSides;
    std::vector<vector<vector<Point> > > cquadPtsSide2;
    std::vector<vector<vector<Point> > > cquadPtsSide4;
    std::vector<vector<vector<short> > > cquadCrackSidesAtFront;
    std::vector<vector<int> >            cquadRepetition;
    std::vector<vector<Vector> >         cquadOffset;
    std::vector<vector<vector<Point> > > triangles;
    std::vector<vector<int> >            triNCells;
    std::vector<vector<vector<short> > > triCrackSidesAtFront;
    std::vector<vector<int> >            triRepetition;
    std::vector<vector<Vector> >         triOffset;    
    std::vector<vector<vector<Point> > > arcs;
    std::vector<vector<int> >            arcNCells;
    std::vector<vector<int> >            arcCrkFrtSegID;
    std::vector<vector<vector<Point> > > ellipses;
    std::vector<vector<int> >            ellipseNCells;
    std::vector<vector<int> >            ellipseCrkFrtSegID;
    std::vector<vector<vector<Point> > > pellipses;
    std::vector<vector<int> >            pellipseNCells;
    std::vector<vector<int> >            pellipseCrkFrtSegID;
    std::vector<vector<double> >         pellipseExtent;
    std::vector<Point>                   cmin,cmax;  

    // Crack data after mesh  
    std::vector<double>               css;  // Average length of crack-front segments
    std::vector<double>               csa;  // Average angle of crack-front segments  
    std::vector<vector<Point> >        cx;  // Coordinates of crack nodes
    std::vector<vector<IntVector> >    ce;  // Crack elements
    std::vector<vector<int> >  cfSegNodes;  // Crack-front nodes 
    std::vector<vector<double> > cfSegVel;  // Velocity of crack-front nodes
    std::vector<vector<double> >cfSegTime;  // Time instant of crack propagation
    std::vector<vector<double> >cfSegDis;   // Crack incremental 
    std::vector<vector<int> > cfSegPreIdx;  // node[i]=node[preIdx]
    std::vector<vector<int> > cfSegMinIdx;  // Minimum node-index of the sub-crack
    std::vector<vector<int> > cfSegMaxIdx;  // Maximum node-index of the sub-crack
    std::vector<vector<Point> > cfSegPtsT;  // Crack-front points after propagation
    std::vector<vector<Vector> >  cfSegV1;  // Bi-normals at crack-front nodes
    std::vector<vector<Vector> >  cfSegV2;  // Outer normals at crack-front nodes
    std::vector<vector<Vector> >  cfSegV3;  // Tangential normals at crack-front nodes
    std::vector<vector<Vector> >   cfSegJ;  // J-integral at crack-front nodes
    std::vector<vector<Vector> >   cfSegK;  // SIF at crack-front nodes
    std::vector<vector<vector<int> > > cnset;  // Crack-node subset in each patch
    std::vector<vector<vector<int> > > cfnset; // Crack-front-node index subset
    std::vector<vector<vector<int> > > cfsset; // Crack-front-seg subset in each patch


    // PRIVATE METHODS
           // Private methods in ReadAndDiscretizeCracks.cc
    void   ReadQuadCracks(const int&,const ProblemSpecP&);
    void   ReadCurvedQuadCracks(const int&,const ProblemSpecP&);
    void   ReadTriangularCracks(const int&,const ProblemSpecP&);
    void   ReadArcCracks(const int&,const ProblemSpecP&);
    void   ReadEllipticCracks(const int&,const ProblemSpecP&);
    void   ReadPartialEllipticCracks(const int&,const ProblemSpecP&);
    void   OutputInitialCrackPlane(const int&);
    void   DiscretizeQuadCracks(const int&,int&);
    void   DiscretizeCurvedQuadCracks(const int&,int&);
    void   GetGlobalCoordinatesQuad(const int&, const int&,const int&, 
                                    const double&, const double&, Point&);
    void   GetGlobalCoordinatesTriangle(const int&, const int&,const int&,
                                        const double&, const double&, Point&);
    void   DiscretizeTriangularCracks(const int&,int&);
    void   DiscretizeArcCracks(const int&,int&);
    void   DiscretizeEllipticCracks(const int&,int&);
    void   DiscretizePartialEllipticCracks(const int&,int&);
    void   CombineCrackSegments(const int&);
    short  TwoPointsCoincide(const Point&,const Point&);
    short  TwoDoublesEqual(const double&, const double&, const double&);
    void   ResetCrackNodes(const int&, const int&, const int&);
    void   ResetCrackElements(const int&, const int&, const int&);
    void   ResetCrackFrontNodes(const int&, const int&, const int&);
    void   ReorderCrackFrontNodes(const int&);
    void   FindCrackFrontNodeIndexes(const int&);   
    short  SmoothCrackFrontAndCalculateNormals(const int& m);    
    short  CubicSpline(const int& n, const int& m, const int& n1,
                       double [], double [], double [],
                       int [], double [], const double&);    
    void   CalculateCrackFrontNormals(const int& m);
    void   FindSegsFromNode(const int&,const int&, int []);
    void   OutputInitialCrackMesh(const int&);    
    Vector TwoPtsDirCos(const Point&,const Point&);
    Vector TriangleNormal(const Point&,const Point&,const Point&);
    
           // Private methods in ParticleNodePairVelocityField.cc
    IntVector CellOffset(const Point&, const Point&, Vector);    
    short  ParticleNodeCrackPLaneRelation(const Point&,const Point&,
                            const Point&,const Point&,const Point&);        
    double Volume(const Point&,const Point&,const Point&,const Point&);
    
           // Private methods in FractureParametersCalculation.cc
    void   DetectIfDoingFractureAnalysisAtThisTimestep(double);
    void   FindJIntegralPath(const Point&,const Vector&,const Vector&,
                             const Vector&,double []);    
    bool   FindIntersectionJPathAndCrackPlane(const int&,const double& r,
                                              const double [],Point&);    
    void   FindPlaneEquation(const Point&,const Point&,const Point&,
                             double&,double&,double&,double&);
    short  PointInTriangle(const Point&,const Point&,const Point&,const Point&); 
    void   GetPositionToComputeCOD(const int&,const Point&,const Matrix3&,double&);
    void   OutputCrackFrontResults(const int&, double time, double delT);    

           // Private methods in CrackPropagation.cc
    void   TrimLineSegmentWithBox(const Point&,Point&,const Point&,const Point&);
    void   PruneCrackFrontAfterPropagation(const int& m,const double& ca);
    
           // Private methods in MoveCracks.cc
    short  PhysicalGlobalGridContainsPoint(const double&,const Point&);
    void   ApplySymmetricBCsToCrackPoints(const Vector&,const Point&,Point&); 

           // Private methods in UpdateCrackFront.cc
    void   OutputCrackGeometry(const int&,const int&);

 protected:
     MPMLabel* lb;
};

}

#endif  /* __CRACK_H__*/
