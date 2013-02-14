#ifndef MATITI_MESHNODESUBSET_H
#define MATITI_MESHNODESUBSET_H

#include <Common/RefCounted.h>
#include <Core/Geometry/IntVector.h>

#include <vector>
#include <iostream>


using std::ostream;
using SCIRun::IntVector;

namespace Matiti {
  typedef int meshNodeIndex;
  typedef int meshNodeId;

  class MeshNodeVariableBase;

  class MeshNodeSubset : public RefCounted {
  public:
    MeshNodeSubset(int num_meshNodes, int matlIndex);
    MeshNodeSubset(int num_meshNodes, int matlIndex, 
                   const std::vector<MeshNodeSubset*>& subsets);
    MeshNodeSubset();
    ~MeshNodeSubset();
    
    bool operator==(const MeshNodeSubset& ps) const {
      return d_numMeshNodes == ps.d_numMeshNodes && 
        d_matlIndex == ps.d_matlIndex;
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

    const std::vector<MeshNodeSubset*>& getNeighborSubsets() const {
      return neighbor_subsets;
    }
    
    friend ostream& operator<<(ostream& out, Matiti::MeshNodeSubset& pset);

   private:
    
    meshNodeIndex* d_meshNodes;
    meshNodeIndex d_numMeshNodes;
    meshNodeIndex d_allocatedSize;
    int d_numExpansions;

    int d_matlIndex;

    std::vector<MeshNodeSubset*> neighbor_subsets;

    void fillset();

    void init();
    MeshNodeSubset(const MeshNodeSubset& copy);
    MeshNodeSubset& operator=(const MeshNodeSubset&);
  };
} // End namespace Matiti

#endif
