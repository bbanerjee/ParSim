/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __CORE_GEOMETRYPIECE_SPECIAL_GEOMETRYPIECE_H__
#define __CORE_GEOMETRYPIECE_SPECIAL_GEOMETRYPIECE_H__

#include <Core/GeometryPiece/GeometryPiece.h>

#include <Core/Geometry/Point.h>
#include <Core/Grid/Box.h>
#include <Core/Math/Matrix3.h>

#include <cmath>
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!

\class SpecialGeomPiece

\brief Abstract base class for smooth geometry pieces

\warning Does not allow for correct application of symmetry
boundary conditions.  Use symmetry at your own risk.
The end caps are exactly the same diameter as the outer
diameter of the cylinder and are welded perfectly to the
cylinder.

\author  Biswajit Banerjee \n
C-SAFE and Department of Mechanical Engineering \n
University of Utah \n
*/
/////////////////////////////////////////////////////////////////////////////

class SpecialGeomPiece : public GeometryPiece
{
public:
  using ParticleScalarData = std::vector<double>;
  using ParticleVectorData = std::vector<Vector>;
  using ParticleTensorData = std::vector<Matrix3>;

  //////////////////////////////////////////////////////////////////////
  /*! Constructor */
  //////////////////////////////////////////////////////////////////////
  SpecialGeomPiece() = default;

  //////////////////////////////////////////////////////////////////////
  /*! Destructor */
  //////////////////////////////////////////////////////////////////////
  virtual ~SpecialGeomPiece() = default;

  /// Make a clone
  virtual GeometryPieceP
  clone() const = 0;

  static const std::string TYPE_NAME;
  virtual std::string
  getType() const
  {
    return TYPE_NAME;
  }

  //////////////////////////////////////////////////////////////////////
  /*! Determines whether a point is inside the cylinder. */
  //////////////////////////////////////////////////////////////////////
  virtual bool
  inside(const Point& p) const = 0;

  //////////////////////////////////////////////////////////////////////
  /*! Returns the bounding box surrounding the box. */
  //////////////////////////////////////////////////////////////////////
  virtual Box
  getBoundingBox() const = 0;

  //////////////////////////////////////////////////////////////////////
  /*! Creates points and returns count of points */
  //////////////////////////////////////////////////////////////////////
  virtual unsigned int
  createPoints() = 0;

  //////////////////////////////////////////////////////////////////////
  /*! Returns the vector containing the set of particle locations */
  //////////////////////////////////////////////////////////////////////
  std::vector<Point>*
  getPoints();

  //////////////////////////////////////////////////////////////////////
  /*! Returns the vector containing the set of particle data using name */
  //////////////////////////////////////////////////////////////////////
  const ParticleScalarData*
  getScalar(const std::string& scalar_name) const;

  const ParticleVectorData*
  getVector(const std::string& vector_name) const;

  const ParticleTensorData*
  getTensor(const std::string& tensor_name) const;

  //////////////////////////////////////////////////////////////////////
  /*! Deletes the vector containing the set of particle locations */
  //////////////////////////////////////////////////////////////////////
  void
  deletePoints();

  //////////////////////////////////////////////////////////////////////
  /* Deletes the vector containing the set of particle scalars */
  //////////////////////////////////////////////////////////////////////
  void
  deleteScalar(const std::string& name);

  //////////////////////////////////////////////////////////////////////
  /* Deletes the vector containing the set of particle vectors          */
  //////////////////////////////////////////////////////////////////////
  void
  deleteVector(const std::string& name);

  //////////////////////////////////////////////////////////////////////
  /* Deletes the vector containing the set of particle tensors          */
  //////////////////////////////////////////////////////////////////////
  void
  deleteTensor(const std::string& name);

  //////////////////////////////////////////////////////////////////////
  /*! Returns the number of particles */
  //////////////////////////////////////////////////////////////////////
  int
  returnPointCount() const;

  //////////////////////////////////////////////////////////////////////
  /*! Set the particle spacing */
  //////////////////////////////////////////////////////////////////////
  void
  setParticleSpacing(double dx);

  //////////////////////////////////////////////////////////////////////
  /*! Set the grid cell size                                          */
  //////////////////////////////////////////////////////////////////////
  void
  setCellSize(Vector DX);

protected:
  //////////////////////////////////////////////////////////////////////
  /*! Writes the particle locations to a file that can be read by
      the FileGeometryPiece */
  //////////////////////////////////////////////////////////////////////
  void
  writePoints(const std::string& f_name, const std::string& var);

protected:
  std::vector<Point> d_points;
  std::map<std::string, ParticleScalarData> d_scalars;
  std::map<std::string, ParticleVectorData> d_vectors;
  std::map<std::string, ParticleTensorData> d_tensors;

  double d_dx{ 1.0 };
  Vector d_DX{ 0.0, 0.0, 0.0 };
};
} // End namespace Uintah

#endif // __CORE_GEOMETRYPIECE_SPECIAL_GEOMETRYPIECE_H__
