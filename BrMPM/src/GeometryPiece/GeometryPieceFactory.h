#ifndef __MATITI_GEOMETRY_PIECE_FACTORY_H__
#define __MATITI_GEOMETRY_PIECE_FACTORY_H__

#include <NodePArray.h>
#include <ElementPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace BrMPM {

  class GeometryPiece;

  class GeometryPieceFactory
  {
  public:

    static GeometryPiece* create(Uintah::ProblemSpecP& ps,
                                 NodePArray& nodes,
                                 ElementPArray& elements);
  }; // end class

}  // end namespace

#endif
