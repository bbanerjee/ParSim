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

#include <GeometryPiece/GeometryPieceFactory.h>
#include <GeometryPiece/GeometryPiece.h>
#include <GeometryPiece/BoxGeometryPiece.h>
#include <GeometryPiece/GeometryReader.h>
#include <GeometryPiece/PlaneGeometryReader.h>
#include <Core/Exception.h>

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

GeometryPiece* 
GeometryPieceFactory::create(Uintah::ProblemSpecP& ps,
                             NodePArray& nodes,
                             ElementPArray& elements, Vector3D& gridSize)
{
  // Get the geometry 
  Uintah::ProblemSpecP geom_ps = ps->findBlock("Geometry");
  if (!geom_ps) {
    std::ostringstream out;
    out << "**ERROR** No geometry information found for body";
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
  
  // Get the attribute of the geometry
  std::string geom_type;
  if (!geom_ps->getAttribute("type", geom_type)) {
    std::ostringstream out;
    out << "**ERROR** Geometry does not have type information" << geom_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (geom_type == "box") {
    return (new BoxGeometryPiece(geom_ps, nodes, elements, gridSize));
  } else if (geom_type == "file") {
    // This is to read 3d files generated with Abaqus
    return (new GeometryReader(geom_ps, nodes, elements));
  } else if (geom_type == "plane_file") {
    // This is to read 2d files used in EMUNE
    return (new PlaneGeometryReader(geom_ps, nodes, elements));
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown geometry type" << geom_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
