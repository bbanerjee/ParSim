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

#ifndef __CORE_GEOMETRY_GEOMETRY_PIECE_FACTORY_H__
#define __CORE_GEOMETRY_GEOMETRY_PIECE_FACTORY_H__

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/Patch.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <map>
#include <string>
#include <vector>

namespace Uintah {

class GeometryPieceFactory
{
public:
  // this function has a switch for all known go_types
  // and calls the proper class' readParameters()
  // addMaterial() calls this
  static void
  create(const ProblemSpecP& ps,
         std::vector<GeometryPieceP>& objs,
         const std::string& input_ups_dir = "");

  // Clears out the saved geometry piece information...  In theory, this should
  // only be called by (and necessary for) the Switcher component (and only if
  // a component that is being switched to happens to be a 'copy' of a previous
  // component).
  static void
  resetFactory();

  // Runs through all the GeometryPiece that have been created and
  // sets their flag for first time output.  This should be done at
  // the beginning of any output of a problemspec.
  static void
  resetGeometryPiecesOutput();

  // does the geometry piece already exist in the map
  static int
  geometryPieceExists(const ProblemSpecP& ps, const bool isTopLevel = true);

  /*
   *  \brief Returns a map of all the named geometry pieces. This maps the
   geometry name to the GeometryPiece.
   */
  const static std::map<std::string, GeometryPieceP>&
  getNamedGeometryPieces();

  /*
   *  \brief Finds the inside points for ALL geometry pieces. This can be called
   for example in a pre-processing step.
   */
  static void
  findInsidePoints(const Uintah::Patch* const patch);

  /*
   *  \brief Finds and returns a reference to a vector that contains the points
   inside the geomPieceName on the specified patch. The points are of type
   Uintah::Point. If no points are found, then an empty vector is store. If
   points were already found by previous call, then this returns a reference to
   the points that were already found (which avoids additional computation).
   */
  const static std::vector<Point>&
  getInsidePoints(const std::string geomPieceName,
                  const Uintah::Patch* const patch);

  /*
   *  \brief Finds and returns a reference to a vector that contains the points
  inside ALL the named geometries that live on patch.
   */
  const static std::vector<Point>&
  getInsidePoints(const Uintah::Patch* const patch);

private:
  // This variable records all named GeometryPieces, so that if they
  // are referenced a 2nd time, they don't have to be rebuilt.
  //
  // Assuming multiple GeometryPieceFactory's will not exist and if
  // they do, they won't be executing at the same time (in different
  // threads)... If this is not the case, then this variable should
  // be locked...
  inline static std::map<std::string, GeometryPieceP> s_namedPieces;
  inline static std::vector<GeometryPieceP> s_unnamedPieces;

  // geompiece name -> (patchID, inside points vector)
  inline static std::map<std::string, std::map<int, std::vector<Point>>>
    s_insidePointsMap;

  // patchID -> insidePoints. returns ALL points inside geometries for a given
  // patch
  inline static std::map<int, std::vector<Point>> s_allInsidePointsMap;

  /*
   *  \brief A private helper function to check whether we already looked for
   the inside points for this patch and geometry or not.
   */
  static bool
  foundInsidePoints(const std::string geomName, const int patchID);
};

} // End namespace Uintah

#endif /* __CORE_GEOMETRY_GEOMETRY_PIECE_FACTORY_H__ */
