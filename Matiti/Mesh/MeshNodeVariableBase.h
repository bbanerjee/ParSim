#ifndef MATITI_MESHNODE_VARIABLE_BASE_H
#define MATITI_MESHNODE_VARIABLE_BASE_H

#include <Mesh/MeshNodeSet.h>
#include <Mesh/Variable.h>
#include <Mesh/constVariable.h>

#include <vector>


namespace Matiti {

  class MeshNodeSet;

  typedef constVariableBase<MeshNodeVariableBase> constMeshNodeVariableBase;

  class MeshNodeVariableBase : public Variable {

    public:
      
      virtual ~MeshNodeVariableBase();

      virtual MeshNodeVariableBase* clone() = 0;     
      virtual MeshNodeVariableBase* cloneSet(MeshNodeSet*) = 0;

      // Make a new default object of the base class.
      virtual MeshNodeVariableBase* cloneType() const = 0;
      virtual constMeshNodeVariableBase* cloneConstType() const = 0;

      virtual void copyData(const MeshNodeVariableBase* src) = 0;
      
      virtual void allocate(MeshNodeSet*) = 0;
      virtual void allocate(int totalParticles) = 0;
      virtual void gather(MeshNodeSet* dest,
                          const std::vector<MeshNodeSet*> &subsets,
                          const std::vector<MeshNodeVariableBase*> &srcs,
                          particleIndex extra = 0) = 0;
      virtual int size() = 0;

      MeshNodeSet* getMeshNodeSet() const {
        //ASSERT(!isForeign());
        return d_nodeset;
      }

      virtual void* getBasePointer() const = 0;
      virtual RefCounted* getRefCounted() = 0;
      virtual void getSizeInfo(std::string& elems, unsigned long& totsize,
                               void*& ptr) const = 0;
   protected:

      MeshNodeVariableBase(const MeshNodeVariableBase&);
      MeshNodeVariableBase(MeshNodeSet* nodeset);
      MeshNodeVariableBase& operator=(const MeshNodeVariableBase&);
      
      MeshNodeSet*  d_nodeset;

   private:
   };

} // End namespace Matiti

#endif
