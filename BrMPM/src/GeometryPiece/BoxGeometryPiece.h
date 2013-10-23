#ifndef __MATITI_BOX_GEOMETRY_PIECE_H__
#define __MATITI_BOX_GEOMETRY_PIECE_H__

#include <GeometryPiece/GeometryPiece.h>
#include <Types.h>
#include <MPMDatawarehouseP.h>
#include <map>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace BrMPM 
{

  class BoxGeometryPiece : public GeometryPiece
  {
  public:

    BoxGeometryPiece(Uintah::ProblemSpecP& ps, MPMDatawarehouseP& dw);
    virtual ~BoxGeometryPiece();

    Box3D boundingBox() const;

    bool inside (const Point3D& pt) const;

    std::string name() const;

  protected:

    void createParticles(MPMDatawarehouseP& dw);

  private:

    Box3D d_box;

  }; // end class

} // end namespace

#endif
