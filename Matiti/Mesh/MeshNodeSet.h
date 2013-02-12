#ifndef MESHNODE_SET_H
#define MESHNODE_SET_H

#include <Common/RefCounted.h>

namespace Matiti {
  typedef int meshNodeIndex;
  typedef int meshNodeId;

  class MeshNodeVariableBase;

  class MeshNodeSet : public RefCounted {
  public:
    MeshNodeSet(int num_meshNodes, int matlIndex);
    MeshNodeSet();
    ~MeshNodeSet();
    
    bool operator==(const MeshNodeSet& ps) const {
      return d_numMeshNodes == ps.d_numMeshNodes && d_matlIndex == ps.d_matlIndex ;
    }
      
    void addMeshNode(meshNodeIndex idx) {
      if(d_numMeshNodes >= d_allocatedSize)
        expand(1);
      d_meshNodes[d_numMeshNodes++] = idx;
    }
    meshNodeIndex addMeshNodes(meshNodeIndex count);

    void resize(meshNodeIndex idx);

    typedef meshNodeIndex* iterator;
      
    iterator begin() {
      return d_meshNodes;
    }
      
    iterator end() {
      return d_meshNodes+d_numMeshNodes;
    }
      
    meshNodeIndex* getPointer()
    {
      return d_meshNodes;
    }
      
    meshNodeIndex numMeshNodes() {
      return d_numMeshNodes;
    }
      
    void set(meshNodeIndex idx, meshNodeIndex value) {
      d_meshNodes[idx] = value;
    }

    int getMatlIndex() const {
      return d_matlIndex;
    }

    void expand(meshNodeIndex minSizeIncrement);

    // sort the set by meshNode IDs
    void sort(MeshNodeVariableBase* meshNodeIDs);

    friend ostream& operator<<(ostream& out, Matiti::MeshNodeSet& nodeset);

   private:

    meshNodeIndex* d_meshNodes;
    meshNodeIndex d_numMeshNodes;
    meshNodeIndex d_allocatedSize;
    int d_numExpansions;

    int d_matlIndex;

    void fillset();

    void init();
    MeshNodeSet(const MeshNodeSet& copy);
    MeshNodeSet& operator=(const MeshNodeSet&);
  };
} // End namespace Matiti

#endif
