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

#ifndef MATITI_GEOMETRY_READER_H
#define MATITI_GEOMETRY_READER_H

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
  class GeometryReader : public GeometryPiece
  {
  public: 

    GeometryReader(Uintah::ProblemSpecP& ps,
                   NodePArray& nodes,
                   ElementPArray& elems);

    virtual ~GeometryReader();

    Box3D boundingBox() const;

    bool inside (const Point3D& pt) const;

    std::string name() const;

  protected:

    void readGeometryInputFiles(Uintah::ProblemSpecP& ps,
                                NodePArray& nodes,
                                ElementPArray& elems);
    void readSurfaceMeshNodes(const std::string& fileName);

    void readVolumeMeshNodesAndElements(const std::string& fileName,
                                        NodePArray& nodes,
                                        ElementPArray& elements);

    void readVolumeMeshNode(const std::string& inputLine,
                            NodePArray& nodes);

    void readVolumeMeshElement(const std::string& inputLine,
                               ElementPArray& elements);

    void findNodalAdjacentElements(ElementPArray& elements);

    void findSurfaceNodes(NodePArray& nodes);

  private:

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

    typedef std::multimap<long64, int> BucketIDNodeIDMap;
    BucketIDNodeIDMap d_bucket_to_node_map;

    std::vector<Point3D> d_surf_pts; 
    double d_xmax, d_ymax, d_zmax, d_xmin, d_ymin, d_zmin;
    int d_num_buckets_x;

  }; // end class

} // end namespace
#endif
