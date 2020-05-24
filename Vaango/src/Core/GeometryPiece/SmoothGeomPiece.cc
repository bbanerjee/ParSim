/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <Core/GeometryPiece/SmoothGeomPiece.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace Uintah;
using std::vector;
using std::string;


const string SmoothGeomPiece::TYPE_NAME = "smooth_geom";

SmoothGeomPiece::SmoothGeomPiece()
{
  d_dx = 1.0;
}

SmoothGeomPiece::~SmoothGeomPiece()
{
}

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle locations */
//////////////////////////////////////////////////////////////////////
vector<Point>* 
SmoothGeomPiece::getPoints()
{
  return &d_points;
}

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle volumes */
//////////////////////////////////////////////////////////////////////
vector<double>* 
SmoothGeomPiece::getVolume()
{
  return &d_volume;
}

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle temperatures */
//////////////////////////////////////////////////////////////////////
vector<double>* 
SmoothGeomPiece::getTemperature()
{
  return &d_temperature;
}
//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle color*/
//////////////////////////////////////////////////////////////////////
vector<double>* 
SmoothGeomPiece::getColors()
{
  return &d_color;
}

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle forces */
//////////////////////////////////////////////////////////////////////
vector<Vector>* 
SmoothGeomPiece::getForces()
{
  return &d_forces;
}

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle fiber directions */
//////////////////////////////////////////////////////////////////////
vector<Vector>* 
SmoothGeomPiece::getFiberDirs()
{
  return &d_fiberdirs;
}

//////////////////////////////////////////////////////////////////////
/* Returns the vectors containing the CPDI or CPTI R-vectors        */
//////////////////////////////////////////////////////////////////////
vector<Vector>* 
SmoothGeomPiece::getRvec1()
{
  return &d_rvec1;
}
vector<Vector>* 
SmoothGeomPiece::getRvec2()
{
  return &d_rvec2;
}
vector<Vector>* 
SmoothGeomPiece::getRvec3()
{
  return &d_rvec3;
}

////////////////////////////////////////////////////////////////////// // gcd adds
/* Returns the vector containing the set of particle velocity components */
//////////////////////////////////////////////////////////////////////
vector<Vector>* 
SmoothGeomPiece::getVelocity()
{
  return &d_velocity;
}                                                  // end gcd add

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set ofparticle size tensors    */
//////////////////////////////////////////////////////////////////////
vector<Matrix3>* 
SmoothGeomPiece::getSize()
{
  return &d_size;
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle locations */
//////////////////////////////////////////////////////////////////////
void 
SmoothGeomPiece::deletePoints()
{
  d_points.clear();
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle volumes */
//////////////////////////////////////////////////////////////////////
void 
SmoothGeomPiece::deleteVolume()
{
  d_volume.clear();
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle sizes          */
//////////////////////////////////////////////////////////////////////
void 
SmoothGeomPiece::deleteSizes()
{
  d_size.clear();
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle temperatures */
//////////////////////////////////////////////////////////////////////
void 
SmoothGeomPiece::deleteTemperature()
{
  d_temperature.clear();
}

void 
SmoothGeomPiece::writePoints(const std::string& f_name, const std::string& var)
{
  if (var == "vol") {
    std::ofstream file(f_name.c_str());
    file.setf(std::ios::scientific, std::ios::floatfield);
    file.precision(8);
    file << "x_coord " << " y_coord " << " z_coord " << var << "\n";
    for (unsigned int jj = 0; jj < d_points.size(); ++jj) {
      file << d_points[jj].x() << " " << d_points[jj].y() << " " << d_points[jj].z() 
           << " " << d_volume[jj] << "\n";
    }
    file.close();
    std::cout << "Wrote output file " << f_name << std::endl;
  } else {
    std::cout << "** No output file " << f_name << " written." << std::endl;
  }
}

int 
SmoothGeomPiece::returnPointCount() const
{
  return d_points.size();
}

void 
SmoothGeomPiece::setParticleSpacing(double dx)
{
  d_dx = dx;
}

void 
SmoothGeomPiece::setCellSize(Vector DX)
{
  d_DX = DX;
}
