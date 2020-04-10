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

#ifndef __MATITI_BOX_GEOMETRY_PIECE_H__
#define __MATITI_BOX_GEOMETRY_PIECE_H__

#include <GeometryPiece/GeometryPiece.h>
#include <Types/Types.h>
#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>
#include <map>

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>

namespace Matiti 
{

  class BoxGeometryPiece : public GeometryPiece
  {
  public:

    BoxGeometryPiece(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elements, Vector3D& gridSize);

    BoxGeometryPiece(const Point3D& lower, 
                     const Point3D& upper, 
                     const Uintah::IntVector& numElem,
                     NodePArray& nodes, 
                     ElementPArray& elements, 
                     Vector3D& gridSize);

    virtual ~BoxGeometryPiece();

    Box3D boundingBox() const;

    bool inside (const Point3D& pt) const;

    std::string name() const;


  protected:

    void createNodes(NodePArray& nodes, Vector3D& gridSize);

    void createElements(ElementPArray& elem);

    void findNodalAdjacentElements(ElementPArray& elements);

  private:

    Box3D d_box;
    IntArray3 d_num_elements;
    double d_dx;
    double d_dy;
    double d_dz;

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

  }; // end class

} // end namespace

#endif
