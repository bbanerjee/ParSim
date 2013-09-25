#ifndef MATITI_GEOMETRY_READER_H
#define MATITI_GEOMETRY_READER_H

#include <GeometryPiece/GeometryPiece.h>
#include <Types.h>
#include <NodeP.h>
#include <NodePArray.h>
#include <ElementPArray.h>
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
