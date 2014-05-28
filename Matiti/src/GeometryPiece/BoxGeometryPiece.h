#ifndef __MATITI_BOX_GEOMETRY_PIECE_H__
#define __MATITI_BOX_GEOMETRY_PIECE_H__

#include <GeometryPiece/GeometryPiece.h>
#include <Types/Types.h>
#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>
#include <map>

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti 
{

  class BoxGeometryPiece : public GeometryPiece
  {
  public:

    BoxGeometryPiece(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elements, Vector3D& gridSize);
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
