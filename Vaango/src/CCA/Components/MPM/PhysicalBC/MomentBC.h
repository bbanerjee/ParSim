/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef UINTAH_MPM_MOMENTBC_H
#define UINTAH_MPM_MOMENTBC_H

#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/MPM/PhysicalBC/LoadCurve.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/Grid.h>
#include <Core/Math/Matrix3.h>
#include <iosfwd>

namespace Uintah {

  using namespace Uintah;

  class GeometryPiece;
  class ParticleCreator;

/**************************************

CLASS
	MomentBC

	Moment Boundary Conditions for MPM

GENERAL INFORMATION

	MomentBC.h

	Danny Cheng
	Sensing & Automation
	Callaghan Innovation

KEYWORDS
	MomentBC

DESCRIPTION
   Stores the moment load curves and boundary information for
   moment boundary conditions that can be applied to surfaces
   with simple geometry -  planes, cylinders and spheres.

WARNING

****************************************/

  class MomentBC : public MPMPhysicalBC  {

  public:

    // Construct a MomentBC object that contains
    // the area over which moment is to be applied
    // and the value of that moment (in the form
    // of a load curve)
    MomentBC(ProblemSpecP& ps, const GridP& grid, const MPMFlags* flags);
    ~MomentBC();
    virtual std::string getType() const;

    virtual void outputProblemSpec(ProblemSpecP& ps);

    // Locate and flag the material points to which this moment BC is
    // to be applied.
    bool flagMaterialPoint(const Point& p, const Vector& dxpp);

    // Get the load curve number for this moment BC
    inline int loadCurveID() const {return d_loadCurve->getID();}

    // Get the surface
    inline GeometryPiece* getSurface() const {return d_surface;}

    // Get the surface type
    inline std::string getSurfaceType() const {return d_surfaceType;}

    // Set the number of material points on the surface
    inline void numMaterialPoints(long num) {d_numMaterialPoints = num;}

    // Get the number of material points on the surface
    inline long numMaterialPoints() const {return d_numMaterialPoints;}

    // Get the load curve
    inline LoadCurve<double>* getLoadCurve() const {return d_loadCurve;}

    // Get the applied moment at time t
    inline double moment(double t) const {return d_loadCurve->getLoad(t);}

    // Get the force per particle at time t
    double forcePerParticle(double time) const;

    // Get the force vector to be applied at a point
    Vector getForceVector(const Point& px, double forcePerParticle,
			  const double time,
                          const Matrix3& defGrad) const;

    // Get the force vector to be applied at 4 corners of the point
    Vector getForceVectorCBDI(const Point& px, const Matrix3& pSize,
			      const Matrix3& pDeformationMeasure,
			      double forcePerParticle, const double time,
			      Point& pExternalForceCorner1,
			      Point& pExternalForceCorner2,
			      Point& pExternalForceCorner3,
			      Point& pExternalForceCorner4,
			      const Vector& dxCell) const;

  private:

    // Prevent empty constructor
    MomentBC();

    // Prevent copying
    MomentBC(const MomentBC&);
    MomentBC& operator=(const MomentBC&);

    // Private Data
    // Surface information
    GeometryPiece* d_surface;
    std::string d_surfaceType;
    long d_numMaterialPoints;
    bool d_cylinder_end;
    bool d_axisymmetric_end;
    bool d_axisymmetric_side;
    bool d_outwardNormal;

    // Load curve information (Moment and time)
    LoadCurve<double>* d_loadCurve;

    // Normal plane information
    Vector d_norm_norm;	// Normal plane vector.
    double d_norm_d;	// Normal plane d.
    double d_norm_L1L2; // L1^3 + L2^3, where L1 and L2 are the distance from the normal plane to the edges of the boundary.

  public:
    Vector d_dxpp;
    IntVector d_res;

    friend std::ostream& operator<<(std::ostream& out, const Uintah::MomentBC& bc);
  };
} // End namespace Uintah

#endif
