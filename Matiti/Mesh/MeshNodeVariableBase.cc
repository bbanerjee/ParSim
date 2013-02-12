#include <Mesh/MeshNodeVariableBase.h>

using namespace Matiti;

MeshNodeVariableBase::~MeshNodeVariableBase()
{       
   if(d_nodeset && d_nodeset->removeReference())
      delete d_nodeset;
}

MeshNodeVariableBase::MeshNodeVariableBase(MeshNodeSet* nodeset)
   : d_nodeset(nodeset)
{
   if(d_nodeset)
      d_nodeset->addReference();
}

MeshNodeVariableBase::MeshNodeVariableBase(const MeshNodeVariableBase& copy)
   : d_nodeset(copy.d_nodeset)
{
   if(d_nodeset)
      d_nodeset->addReference();
}   

MeshNodeVariableBase& MeshNodeVariableBase::operator=(const MeshNodeVariableBase& copy)
{
   if(this != &copy){
      if(d_nodeset && d_nodeset->removeReference())
         delete d_nodeset;
      d_nodeset = copy.d_nodeset;
      if(d_nodeset)
         d_nodeset->addReference();
   }
   return *this;
}

void MeshNodeVariableBase::setMeshNodeSet(MeshNodeSet* subset)
{
  if(d_nodeset && d_nodeset->removeReference())
    delete d_nodeset;
  d_nodeset = subset;
  if(d_nodeset)
    d_nodeset->addReference();
}
