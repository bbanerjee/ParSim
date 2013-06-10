#ifndef __EMU2DC_GEOMETRY_PIECE_FACTORY_H__
#define __EMU2DC_GEOMETRY_PIECE_FACTORY_H__

#include <NodePArray.h>
#include <ElementPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Emu2DC {

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
