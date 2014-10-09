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

#ifndef MATITI_PLANE_GEOMETRY_READER_H
#define MATITI_PLANE_GEOMETRY_READER_H

#include <GeometryPiece/GeometryPiece.h>
#include <Types/Types.h>
#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>
#include <Geometry/Point3D.h>
#include <Geometry/Box3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <map>

namespace Matiti 
{
  class PlaneGeometryReader : public GeometryPiece
  {
  public: 

    PlaneGeometryReader(Uintah::ProblemSpecP& ps,
                        NodePArray& nodes,
                        ElementPArray& elems);

    virtual ~PlaneGeometryReader();

    Box3D boundingBox() const;

    bool inside (const Point3D& pt) const;

    std::string name() const;

  protected:

    void readGeometryInputFiles(Uintah::ProblemSpecP& ps,
                                NodePArray& nodes,
                                ElementPArray& elems);

    void readMeshNodesAndElements(const std::string& nodeFileName,
                                  const std::string& elemFileName,
                                  NodePArray& nodes,
                                  ElementPArray& elements);

    int readNodeFile(const std::string& inputLine,
                     NodePArray& nodes);

    void readElementFile(const std::string& inputLine,
                         ElementPArray& elements,
                         int numNodes);

    void findNodalAdjacentElements(ElementPArray& elements);

  private:

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

    double d_xmax, d_ymax, d_zmax, d_xmin, d_ymin, d_zmin;

  }; // end class

} // end namespace
#endif
