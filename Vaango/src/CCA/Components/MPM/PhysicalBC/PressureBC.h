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

#ifndef UINTAH_MPM_PRESSUREBC_H
#define UINTAH_MPM_PRESSUREBC_H

#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/MPM/PhysicalBC/LoadCurve.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/Grid.h>
#include <Core/Math/Matrix3.h>

#include <CCA/Components/MPM/PhysicalBC/exprtk/exprtk.hpp>

#include <iosfwd>

namespace Uintah {

using namespace Uintah;

class GeometryPiece;
class ParticleCreator;
   
/**************************************

CLASS
   PressureBC
   
   Pressure Boundary Conditions for MPM
 
GENERAL INFORMATION

   PressureBC.h

   Biswajit Banerjee
   Department of Mechanical Engineering, University of Utah
   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

KEYWORDS
   PressureBC

DESCRIPTION
   Stores the pressure load curves and boundary imformation for
   pressure boundary conditions that can be applied to surfaces
   with simple geometry -  planes, cylinders and spheres.

WARNING
  
****************************************/

   class PressureBC : public MPMPhysicalBC  {

   public:

      // Construct a PressureBC object that contains
      // the area over which pressure is to be applied
      // and the value of that pressure (in the form
      // of a load curve)
      PressureBC(ProblemSpecP& ps, const GridP& grid, const MPMFlags* flags);
      ~PressureBC();
      virtual std::string getType() const override;

      virtual void outputProblemSpec(ProblemSpecP& ps) override;

      // Locate and flag the material points to which this pressure BC is
      // to be applied. 
      bool flagMaterialPoint(const Point& p, const Vector& dxpp);
      
      // Get the load curve number for this pressure BC
      inline int loadCurveID() const {return d_loadCurve->getID();}

      // Get the surface 
      inline GeometryPiece* getSurface() const {return d_surface;}

      // Get the surface type
      inline std::string getSurfaceType() const {return d_surfaceType;}

      // Set the number of material points on the surface
      inline void numMaterialPoints(long num) {d_numMaterialPoints = num;}

      // Get the number of material points on the surface
      inline long numMaterialPoints() const {return d_numMaterialPoints;}

      // Get the area of the surface
      double getSurfaceArea() const;

      // Get the load curve 
      inline LoadCurve<double>* getLoadCurve() const {return d_loadCurve;}

      // Update the load curve
      void updateLoadCurve(const std::vector<double>& time,
                           const std::vector<double>& pressure);

      // Get the applied pressure at time t
      double pressure(double t);
      double pressure(double t, const Point& pX);

      // Get the force per particle at time t
      double forcePerParticle(double time);
      double forcePerParticle(double time, const Point& pX);

      // Get the force vector to be applied at a point 
      Vector getForceVector(const Point& px, 
                            const Vector& pDisp,
                            double forcePerParticle,
                            const double time,
                            const Matrix3& defGrad);

      // Get the force vector to be applied at 4 corners of the point 
      Vector getForceVectorCBDI(const Point& px, 
                                const Vector& pDisp,
                                const Matrix3& pSize,
                                const Matrix3& pDeformationMeasure,
                                double forcePerParticle, const double time,
                                Point& pExternalForceCorner1,
                                Point& pExternalForceCorner2,
                                Point& pExternalForceCorner3,
                                Point& pExternalForceCorner4,
                                const Vector& dxCell);

      std::vector<Point> getParticleCornersCBDI(const Point& pX, 
                                                const Matrix3& pSize,
                                                const Matrix3& pDefGrad,
                                                const Vector& dxCell);
   private:

      // Prevent empty constructor
      PressureBC();

      // Prevent copying
      PressureBC(const PressureBC&);
      PressureBC& operator=(const PressureBC&);
      
      // Private Data
      // Surface information
      GeometryPiece* d_surface;
      std::string d_surfaceType;
      long d_numMaterialPoints;
      bool d_cylinder_end;
      bool d_axisymmetric_end;
      bool d_axisymmetric_side;
      bool d_outwardNormal;
      double d_volFracInsideDomain;

      // Load curve information (Pressure and time)
      LoadCurve<double>* d_loadCurve;

      // Scaling function for load curve
      std::string d_scaling_function_expr;

      // Typedefs for expression parser
      typedef exprtk::symbol_table<double> symbol_table_t;
      typedef exprtk::expression<double>   expression_t;
      typedef exprtk::parser<double>       parser_t;

      // Storage for parsed expression
      symbol_table_t d_symbol_table;
      expression_t   d_expression;
      parser_t       d_parser;
      double         d_time;
      double         d_pos_x;
      double         d_pos_y;
      double         d_pos_z;

      std::pair<int, int> 
      getParticleBoundaryCorners(const Point& pX,
                                 const Matrix3& pSize_new,
                                 const Vector& pNormal,
                                 const Vector& dxCell,
                                 std::vector<Point>& corners) const;

    public:
      Vector d_dxpp;
      IntVector d_res;



      friend std::ostream& operator<<(std::ostream& out, const Uintah::PressureBC& bc);
   };
} // End namespace Uintah

#endif
