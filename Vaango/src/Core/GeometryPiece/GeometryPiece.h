/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __GEOMETRY_PIECE_H__
#define __GEOMETRY_PIECE_H__

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/RefCounted.h>

#include <memory>
#include <string>

namespace Uintah {

class Box;
class GeometryPiece;

template<class T>
class Handle;
using GeometryPieceHP = Handle<GeometryPiece>;
using GeometryPieceUP = std::unique_ptr<GeometryPiece>;
using GeometryPieceP  = std::shared_ptr<GeometryPiece>;

/**************************************

CLASS
   GeometryPiece

   Short description...

GENERAL INFORMATION

   GeometryPiece.h

   John A. Schmidt
   Department of Mechanical Engineering
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

KEYWORDS
   GeometryPiece

DESCRIPTION
   Long description...

WARNING

****************************************/

class GeometryPiece : public RefCounted
{

public:
  GeometryPiece()          = default;
  virtual ~GeometryPiece() = default;

  /// Clone a geometry piece
  virtual GeometryPieceP
  clone() const = 0;

  void
  outputProblemSpec(ProblemSpecP& ps) const;

  virtual Box
  getBoundingBox() const = 0;

  virtual bool
  inside(const Point& p) const = 0;

  virtual bool
  inside(const Point& p, bool default_val) const 
  {
    return inside(p);
  }

  std::string
  getName() const
  {
    return d_name;
  }

  // Returns the type as a string (eg: sphere, box, etc).
  // The string must/will match the string checks found in
  // GeometryPieceFactory.cc::create()
  virtual std::string
  getType() const = 0;

  void
  setName(const std::string& name)
  {
    d_nameSet = true;
    d_name    = name;
  }

  // Call at the beginning of outputing (ProblemSpec) so that this
  // object will output the full spec the first time, and only a
  // reference subsequently.
  void
  resetOutput() const
  {
    d_firstOutput = true;
  }

protected:
  virtual void
  outputHelper(ProblemSpecP& ps) const = 0;

  bool d_nameSet{ false }; // defaults to false
  std::string d_name;

  // Used for outputing the problem spec... on the 1st output, the
  // entire object is output, on the 2nd, only a reference is output.
  // Must be 'mutable' as most of these objects are (mostly) 'const'
  mutable bool d_firstOutput{ true };
};

} // End namespace Uintah

#endif // __GEOMETRY_PIECE_H__
