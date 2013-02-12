#include <Mesh/MeshNodeSet.h>
#include <Mesh/MeshNodeVariable.h>

#include <algorithm>
#include <iostream>

using namespace Matiti;


MeshNodeSet::~MeshNodeSet()
{
  if(d_meshNodes)
    delete[] d_meshNodes;
}

MeshNodeSet::MeshNodeSet() : d_numMeshNodes(0)
{
  init();
}

MeshNodeSet::MeshNodeSet(int num_meshNodes, int matlIndex);
    : d_numMeshNodes(num_meshNodes), d_matlIndex(matlIndex)
{
  init();
  fillset();
}

void
MeshNodeSet::fillset()
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
MeshNodeSet::sort(MeshNodeVariableBase* meshNodeIDs)
{
  MeshNodeVariable<long64>* pIDs =
    dynamic_cast<MeshNodeVariable<long64>*>(meshNodeIDs);
  if (pIDs == 0) std::cerr << "meshNodeID variable must be MeshNodeVariable<long64>" << endl;
  compareIDFunctor comp(pIDs);
  std::sort(d_meshNodes, d_meshNodes+d_numMeshNodes, comp);
}

void
MeshNodeSet::init()
{
  d_meshNodes = 0;
  d_allocatedSize = 0;
  d_numExpansions = 0;
}

void
MeshNodeSet::resize(meshNodeIndex newSize)
{
  // Check for spurious resizes
  if(d_meshNodes) std::cerr << "MeshNodeSets should not be resized after creation" << endl;
  d_allocatedSize = d_numMeshNodes = newSize;
  d_meshNodes = new meshNodeIndex[newSize];
}

void
MeshNodeSet::expand(meshNodeIndex amount)
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

meshNodeIndex MeshNodeSet::addMeshNodes(meshNodeIndex count)
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
  ostream& operator<<(ostream& out, MeshNodeSet& nodeset)
  {
    out << "nodeset: " <<  ", matl "
        << nodeset.getMatlIndex() << << nodeset.numMeshNodes() << " meshNodes, " ;
    return out;
  }
} // end namespace Matiti
