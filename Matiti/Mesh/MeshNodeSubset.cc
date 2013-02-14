#include <Mesh/MeshNodeSubset.h>
#include <Mesh/MeshNodeVariable.h>
#include <Core/Disclosure/TypeUtils.h>

#include <algorithm>
#include <iostream>

using namespace Matiti;
using namespace SCIRun;


MeshNodeSubset::~MeshNodeSubset()
{
  for(int i=0;i<(int)neighbor_subsets.size();i++)
    if(neighbor_subsets[i]->removeReference())
      delete neighbor_subsets[i];
  if(d_meshNodes)
    delete[] d_meshNodes;
}

MeshNodeSubset::MeshNodeSubset() : d_numMeshNodes(0)
{
  init();
}

MeshNodeSubset::MeshNodeSubset(int num_meshNodes, int matlIndex)
    : d_numMeshNodes(num_meshNodes), d_matlIndex(matlIndex)
{
  init();
  fillset();
}

MeshNodeSubset::MeshNodeSubset(int num_meshNodes, int matlIndex, 
                               const vector<MeshNodeSubset*>& neighbor_subsets)
  : d_numMeshNodes(num_meshNodes), d_matlIndex(matlIndex), 
    neighbor_subsets(neighbor_subsets)
{
  init();
  for(int i=0;i<(int)neighbor_subsets.size();i++)
    neighbor_subsets[i]->addReference();
  fillset();
}

void
MeshNodeSubset::fillset()
{
  if (d_numMeshNodes > 0) {
    d_meshNodes = new meshNodeIndex[d_numMeshNodes];
    for(int i=0;i<d_numMeshNodes;i++)
      d_meshNodes[i]=i;
    d_allocatedSize = d_numMeshNodes;
  }
}


class compareIDFunctor
{
public:
  compareIDFunctor(MeshNodeVariable<long64>* meshNodeIDs)
    : meshNodeIDs_(meshNodeIDs) {}
  
  bool operator()(meshNodeIndex x, meshNodeIndex y)
  {
    return (*meshNodeIDs_)[x] < (*meshNodeIDs_)[y];
  }

private:
  MeshNodeVariable<long64>* meshNodeIDs_;
};

void
MeshNodeSubset::sort(MeshNodeVariableBase* meshNodeIDs)
{
  MeshNodeVariable<long64>* pIDs =
    dynamic_cast<MeshNodeVariable<long64>*>(meshNodeIDs);
  compareIDFunctor comp(pIDs);
  std::sort(d_meshNodes, d_meshNodes+d_numMeshNodes, comp);
}

void
MeshNodeSubset::init()
{
  d_meshNodes = 0;
  d_allocatedSize = 0;
  d_numExpansions = 0;
}

void
MeshNodeSubset::resize(meshNodeIndex newSize)
{
  // Check for spurious resizes
  d_allocatedSize = d_numMeshNodes = newSize;
  d_meshNodes = new meshNodeIndex[newSize];
}

void
MeshNodeSubset::expand(meshNodeIndex amount)
{
  meshNodeIndex minAmount = d_numMeshNodes>>2;
  if(minAmount < 10)
    minAmount = 10;
  if(amount < minAmount)
    amount = minAmount;
  d_allocatedSize += amount;
  meshNodeIndex* newmeshNodes = new meshNodeIndex[d_allocatedSize];
  if(d_meshNodes){
    for(meshNodeIndex i = 0; i < d_numMeshNodes; i++)
      newmeshNodes[i] = d_meshNodes[i];
    delete[] d_meshNodes;
  }
  d_meshNodes = newmeshNodes;
}

meshNodeIndex MeshNodeSubset::addMeshNodes(meshNodeIndex count)
{
  if(d_numMeshNodes + count > d_allocatedSize)
    expand(count);


  meshNodeIndex oldsize = d_numMeshNodes;
  d_numMeshNodes += count;

  for(meshNodeIndex idx = oldsize; idx < d_numMeshNodes; idx++)
    d_meshNodes[idx] = idx;
  return oldsize;  // The beginning of the new index range
}

namespace Matiti {
ostream& operator<<(ostream& out, MeshNodeSubset& pset)
{
    out << "pset, matl "
        << pset.getMatlIndex() << 
        << pset.numMeshNodes() << " meshNodes, " ;
    return out;
}
} // end namespace Matiti
