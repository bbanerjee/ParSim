/*
 * The MIT License
 *
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

#include <GeometryPiece/SphereGeometryPiece.h>
//#include <CellInterpolationVector.h>
#include <Exception.h>
#include <GeometryMath/Box3D.h>
#include <GeometryMath/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

using namespace BrMPM;

 SphereGeometryPiece::SphereGeometryPiece()
    :
      d_ghost(0.0, 0.0, 0.0),
      d_num_particles_per_cell(2, 2, 2),
      d_radius(1.e-4),
      d_center(0.0, 0.0, 0.0),
      d_outerlower(0.0, 0.0, 0.0), d_outerupper(0.0, 0.0, 0.0),
     // d_num_grids(0.0, 0.0, 0.0),
      d_cellsize(0.0, 0.0, 0.0)
 {
 }

 SphereGeometryPiece::~SphereGeometryPiece()
 {
 }


  void
  SphereGeometryPiece::initialise(MPMPatch& Patch)
   {
	  d_num_particles_per_cell = Patch.particlesperelement();
      d_ghost = Patch.nGhost();
      //Vector3D Grid;
      //Grid = Patch.d();
     // d_num_grids(Grid);
     // Vector3D Size;
     // Size = Patch.dX();
      d_cellsize = Patch.dX();
   }
    

 SphereGeometryPiece::SphereGeometryPiece(Uintah::ProblemSpecP& ps,
                                          Point3DParticleData& setPointsOfSphere)
 {
  initialise(d_patch);
  d_name = "sphere";
  Uintah::Vector center(0.0, 0.0, 0.0);
  double radius = 0.0;
  ps->require("center", center);
  ps->require("radius", radius);
  d_center[0] = center[0];
  d_center[1] = center[1];
  d_center[2] = center[2];
  d_radius = radius;
  if (radius<=0) {
     std::ostringstream out;
     out << "**ERROR** The radius of sphere should be a positive number. Radius = " << d_radius
     << std::endl;
     throw Exception(out.str(), __FILE__, __LINE__);
     }
  
  //create particles
  createParticles(setPointsOfSphere);
  }

 void
 SphereGeometryPiece::outerBox(Box3D& Box)
 
 {

  d_outerlower.x(d_center.x() - d_radius);
  d_outerlower.y(d_center.y() - d_radius);
  d_outerlower.z(d_center.z() - d_radius);

  d_outerupper.x(d_center.x() + d_radius);
  d_outerupper.y(d_center.y() + d_radius);
  d_outerupper.z(d_center.z() + d_radius);

  Box = Box3D(d_outerlower, d_outerupper);

 }

 bool 
 SphereGeometryPiece::inside (const Point3D& pt) const

 {
  Vector3D relativePosition(d_center ,pt);
  double distance = relativePosition.lengthSq();
 /* if (distance <= d_radius)
     {return true;}
    else
     {return false;} */
  return (distance <= d_radius);
 }

 void
 SphereGeometryPiece::createParticles(Point3DParticleData& setPoints)

 {
  const int tolerance = 1.e-14;
  //int num = d_num_particles_per_cell;
  Vector3D numberParticles;
  numberParticles.x(ceil((2*d_radius/(d_num_particles_per_cell.x()*d_cellsize.x())) - tolerance));
  numberParticles.y(ceil((2*d_radius/(d_num_particles_per_cell.y()*d_cellsize.y())) - tolerance));
  numberParticles.z(ceil((2*d_radius/(d_num_particles_per_cell.z()*d_cellsize.z())) - tolerance));

  Vector3D sizeElement;
  sizeElement.x(2*d_radius/numberParticles.x());
  sizeElement.y(2*d_radius/numberParticles.y());
  sizeElement.z(2*d_radius/numberParticles.z());

  const double location = 0.5;
  Point3D point;
  Box3D outBox;
 
  outerBox(outBox);

  for (int kk = 0; kk < numberParticles.z(); kk++) {
       for (int jj = 0; jj < numberParticles.y(); jj++) {
            for (int ii = 0; ii < numberParticles.x(); ii++) {
                 point.x(outBox.lower().x() + sizeElement.x()*(ii + location));
                 point.y(outBox.lower().x() + sizeElement.y()*(jj + location));
                 point.z(outBox.lower().z() + sizeElement.z()*(kk + location));
                 if (d_patch.insidePatch(point) && inside(point)) {
                     setPoints.emplace_back(point);
                 }
            }
        }
   }
 }
                     
                 
                 
  
  

 
 
     


  

  

  
