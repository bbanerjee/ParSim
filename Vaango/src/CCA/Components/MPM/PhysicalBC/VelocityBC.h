/*
 * The MIT License
 *
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

#ifndef VAANGO_MPM_PHYSICALBC_VELOCITYBC_H
#define VAANGO_MPM_PHYSICALBC_VELOCITYBC_H

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
   VelocityBC
   
   Velocity Boundary Conditions for MPM particles
 
GENERAL INFORMATION

   VelocityBC.h

   Biswajit Banerjee

KEYWORDS
   VelocityBC

DESCRIPTION
   Stores the velocity load curves and boundary imformation for
   velocity boundary conditions that can be applied to surfaces
   with simple geometry -  planes, cylinders and spheres.

WARNING
  
****************************************/

   class VelocityBC : public MPMPhysicalBC  {

   public:

      // Construct a VelocityBC object that contains
      // the area over which velocity is to be applied
      // and the value of that velocity (in the form
      // of a load curve)
      VelocityBC(ProblemSpecP& ps, const GridP& grid, const MPMFlags* flags);
      ~VelocityBC();
      virtual std::string getType() const;

      virtual void outputProblemSpec(ProblemSpecP& ps);

      // Locate and flag the material points to which this velocity BC is
      // to be applied. 
      bool flagMaterialPoint(const Point& p, const Vector& dxpp);
      
      // Get the load curve number for this velocity BC
      inline int loadCurveID() const {return d_loadCurve->getID();}

      // Get the surface 
      inline GeometryPiece* getSurface() const {return d_surface;}

      // Get the surface type
      inline std::string getSurfaceType() const {return d_surfaceType;}

      // Get the load curve 
      inline LoadCurve<Vector>* getLoadCurve() const {return d_loadCurve;}

      // Update the load curve
      void updateLoadCurve(const std::vector<double>& time,
                           const std::vector<Vector>& velocity);

      // Get the applied velocity at time t
      Vector velocity(double t, const Point& pX);

      // Get the force vector to be applied at a point 
      Vector getVelocityVector(const Point& px, 
                               const Vector& pDisp,
                               const double time);

   private:

      // Prevent empty constructor
      VelocityBC();

      // Prevent copying
      VelocityBC(const VelocityBC&);
      VelocityBC& operator=(const VelocityBC&);
      
      // Private Data
      // Surface information
      GeometryPiece* d_surface;
      std::string d_surfaceType;
      bool d_cylinder_end;
      bool d_axisymmetric_end;
      bool d_axisymmetric_side;

      // Load curve information (Velocity and time)
      LoadCurve<Vector>* d_loadCurve;

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

    public:
      Vector d_dxpp;
      IntVector d_res;

      friend std::ostream& operator<<(std::ostream& out, const Uintah::VelocityBC& bc);
   };
} // End namespace Uintah

#endif
