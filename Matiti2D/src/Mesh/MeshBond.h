#ifndef MATITI_MESHBOND_H
#define MATITI_MESHBOND_H

#include <vector>

namespace Matiti {

  class MeshNode;

  class MeshBond
  {
    public:
      MeshBond();
      ~MeshBond();

    private:
      MeshNode* d_startNode;
      MeshNode* d_endNode;
      BondMaterial* d_bondMat;

      // prevent blank creation and copying
      MeshBond(const MeshBond& bond);

  };
} // End namespace Matiti

#endif // MATITI_MESHBOND_H
