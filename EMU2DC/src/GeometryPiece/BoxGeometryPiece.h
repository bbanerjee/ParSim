#ifndef __EMU2DC_BOX_GEOMETRY_PIECE_H__
#define __EMU2DC_BOX_GEOMETRY_PIECE_H__

#include <GeometryPiece/GeometryPiece.h>
#include <NodePArray.h>
#include <ElementPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Emu2DC 
{

  class BoxGeometryPiece : public GeometryPiece
  {
  public:

    BoxGeometryPiece(Uintah::ProblemSpecP& ps, NodePArray& nodes, ElementPArray& elements);
    virtual ~BoxGeometryPiece();

    Box3D boundingBox() const;

    bool inside (const Point3D& pt) const;

    std::string name() const;

  protected:

    void createNodes(NodePArray& nodes, ElementPArray& elements);

  private:

    Box3D d_box;
    int d_num_edge_nodes;

  }; // end class

} // end namespace

#endif
