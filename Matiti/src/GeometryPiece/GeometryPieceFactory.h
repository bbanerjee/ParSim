#ifndef __MATITI_GEOMETRY_PIECE_FACTORY_H__
#define __MATITI_GEOMETRY_PIECE_FACTORY_H__

#include <Containers/NodePArray.h>
#include <Containers/ElementPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {

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
