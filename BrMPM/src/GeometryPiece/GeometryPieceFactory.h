#ifndef __MATITI_GEOMETRY_PIECE_FACTORY_H__
#define __MATITI_GEOMETRY_PIECE_FACTORY_H__

#include <MPMDatawarehouseP.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace BrMPM {

  class GeometryPiece;

  class GeometryPieceFactory
  {
  public:

    static GeometryPiece* create(Uintah::ProblemSpecP& ps,
                                 MPMDatawarehouseP& dw);
  }; // end class

}  // end namespace

#endif
