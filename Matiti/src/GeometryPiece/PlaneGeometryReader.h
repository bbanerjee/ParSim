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

    void readNodeFile(const std::string& inputLine,
                      NodePArray& nodes);

    void readElementFile(const std::string& inputLine,
                         ElementPArray& elements);

    void findNodalAdjacentElements(ElementPArray& elements);

  private:

    typedef std::map<int, NodeP> NodeIDMap;
    NodeIDMap d_id_ptr_map;

    double d_xmax, d_ymax, d_zmax, d_xmin, d_ymin, d_zmin;

  }; // end class

} // end namespace
#endif
