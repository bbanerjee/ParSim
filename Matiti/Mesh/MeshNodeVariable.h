#ifndef MATITI_MESHNODE_VARIABLE_H
#define MATITI_MESHNODE_VARIABLE_H


#include <Mesh/MeshNodeVariableBase.h>
#include <Mesh/MeshNodeData.h>
#include <Mesh/MeshNodeSet.h>

#include <iostream>


namespace Matiti {

  template<class T>
  class MeshNodeVariable : public MeshNodeVariableBase {

    friend class constVariable<MeshNodeVariableBase, MeshNodeVariable<T>, T, meshNodeIndex>;

    public:

      MeshNodeVariable();
      virtual ~MeshNodeVariable();
      MeshNodeVariable(MeshNodeSet* nodeset);
      MeshNodeVariable(MeshNodeData<T>*, MeshNodeSet* nodeset);
      
      void resync() {
        d_nodedata->resize(getMeshNodeSet()->numNodes());
      }
      
      virtual MeshNodeVariableBase* clone();
      virtual const MeshNodeVariableBase* clone() const;
      virtual MeshNodeVariableBase* cloneSet(MeshNodeSet*);
      virtual const MeshNodeVariableBase* cloneSet(MeshNodeSet*) const;

      virtual MeshNodeVariableBase* cloneType() const { 
        return new MeshNodeVariable<T>(); 
      }

      virtual constMeshNodeVariableBase* cloneConstType() const { 
        return new constVariable<MeshNodeVariableBase, MeshNodeVariable<T>, T, meshNodeIndex>();
      }
  
      void copyData(const MeshNodeVariable<T>& src);

      virtual void copyData(const MeshNodeVariableBase* src) { 
        copyData(castFromBase(src)); 
      }
  
      inline T& operator[](meshNodeIndex idx) {
        //ASSERTRANGE(idx, 0, (meshNodeIndex)d_nodedata->size);
        return d_nodedata->data[idx];
      }
      
      inline const T& operator[](meshNodeIndex idx) const {
        // ASSERTRANGE(idx, 0, (meshNodeIndex)d_nodedata->size);
        return d_nodedata->data[idx];
      }

      virtual void copyPointer(MeshNodeVariable<T>&);
      virtual void copyPointer(Variable&);
      virtual void allocate(MeshNodeSet*);
      virtual void allocate(int totalMeshNodes);

      virtual int size() { return d_nodedata->size; }

      // specialized for T=Point
      virtual void gather(MeshNodeSet* dest,
                          const std::vector<MeshNodeSet*> &subsets,
                          const std::vector<MeshNodeVariableBase*> &srcs,
                          meshNodeIndex extra = 0);
  
      virtual void emitNormal(ostream& out); 
  
      virtual void readNormal(istream& in); 
  
      virtual void* getBasePointer() const;

      virtual RefCounted* getRefCounted() {
        return d_nodedata;
      }

      virtual void getSizeInfo(std::string& elems, unsigned long& totsize,
                           void*& ptr) const {
        std::ostringstream str;
        str << getMeshNodeSet()->numNodes();
        elems=str.str();
        totsize = getMeshNodeSet()->numNodes()*sizeof(T);
        ptr = getBasePointer();
      }

    protected:

      MeshNodeVariable(const MeshNodeVariable<T>&);
      MeshNodeVariable<T>& operator=(const MeshNodeVariable<T>&);

    private:

      MeshNodeData<T>* d_nodedata;
      Vector offset_; // only used when T is Point

      static const MeshNodeVariable<T>& castFromBase(const MeshNodeVariableBase* srcptr);
      static Variable* maker();
  };

  template<class T>
  Variable*
  MeshNodeVariable<T>::maker()
  {
    return new MeshNodeVariable<T>();
  }
   
  template<class T>
  MeshNodeVariable<T>::MeshNodeVariable()
    : MeshNodeVariableBase(0), d_nodedata(0)
  {
  }
   
  template<class T>
  MeshNodeVariable<T>::~MeshNodeVariable()
  {
    if(d_nodedata && d_nodedata->removeReference())
      delete d_nodedata;
  }
   
  template<class T>
  MeshNodeVariable<T>::MeshNodeVariable(MeshNodeSet* nodeset)
    : MeshNodeVariableBase(nodeset)
  {
    d_nodedata=new MeshNodeData<T>(nodeset->numNodes());
    d_nodedata->addReference();
  }
   
  template<class T>
  void MeshNodeVariable<T>::allocate(int totalMeshNodes)
  {
    //ASSERT(isForeign());
    //ASSERT(d_nodeset == 0);

    // this is a nodeset-less storage as it could have several.  Should be used for
    // foreign data only.  To iterate over mesh nodes in this nodeset, use gather
    d_nodedata=new MeshNodeData<T>(totalMeshNodes);
    d_nodedata->addReference();
  }

  template<class T>
  void MeshNodeVariable<T>::allocate(MeshNodeSet* nodeset)
  {
    //TAU_PROFILE_TIMER(t1, "Release old MeshNodeVariable<T>::allocate()", "", TAU_USER3);
    //TAU_PROFILE_TIMER(t2, "Allocate Data MeshNodeVariable<T>::allocate()", "", TAU_USER3);

    //TAU_PROFILE_START(t1);
    if(d_nodedata && d_nodedata->removeReference())
      delete d_nodedata;
    if(d_nodeset && d_nodeset->removeReference())
      delete d_nodeset;
    //TAU_PROFILE_STOP(t1);

    d_nodeset=nodeset;
    d_nodeset->addReference();

    //TAU_PROFILE_START(t2);
    d_nodedata=new MeshNodeData<T>(nodeset->numNodes());
    //TAU_PROFILE_STOP(t2);

    d_nodedata->addReference();
  }
   
  template<class T>
  MeshNodeVariableBase*
  MeshNodeVariable<T>::clone()
  { return new MeshNodeVariable<T>(*this); }

  template<class T>
  const MeshNodeVariableBase*
  MeshNodeVariable<T>::clone() const
  { return new MeshNodeVariable<T>(*this); }
   
  template<class T>
  MeshNodeVariableBase*
  MeshNodeVariable<T>::cloneSet(MeshNodeSet* nodeset)
  { return new MeshNodeVariable<T>(d_nodedata, nodeset); }

  template<class T>
  const MeshNodeVariableBase*
  MeshNodeVariable<T>::cloneSet(MeshNodeSet* nodeset) const
  { return new MeshNodeVariable<T>(d_nodedata, nodeset); }

  template<class T>
  const MeshNodeVariable<T>& MeshNodeVariable<T>::castFromBase(const MeshNodeVariableBase* srcptr)
  {
    const MeshNodeVariable<T>* c = dynamic_cast<const MeshNodeVariable<T>* >(srcptr);
    if(!c)
      std::cerr << "Type mismatch in MeshNode variable" << std::endl;
    return *c;
  }

  template<class T>
    void MeshNodeVariable<T>::copyData(const MeshNodeVariable<T>& src)
  {
    //ASSERT(*d_nodeset == *src.d_nodeset);
    *d_nodedata = *src.d_nodedata;
  }


  template<class T>
  MeshNodeVariable<T>::MeshNodeVariable(MeshNodeData<T>* nodedata,
                                        MeshNodeSet* nodeset)
    : MeshNodeVariableBase(nodeset), d_nodedata(nodedata)
  {
    if(d_nodedata)
      d_nodedata->addReference();
  }
   
  template<class T>
  MeshNodeVariable<T>::MeshNodeVariable(const MeshNodeVariable<T>& copy)
    : MeshNodeVariableBase(copy), d_nodedata(copy.d_nodedata)
  {
    if(d_nodedata)
      d_nodedata->addReference();
  }
   
  template<class T>
  void
  MeshNodeVariable<T>::copyPointer(MeshNodeVariable<T>& copy)
  {
    if(this != &copy){
      MeshNodeVariableBase::operator=(copy);
      if(d_nodedata && d_nodedata->removeReference())
        delete d_nodedata;
      d_nodedata = copy.d_nodedata;
      if(d_nodedata)
        d_nodedata->addReference();
    }
  }
   
  template<class T>
  void
  MeshNodeVariable<T>::copyPointer(Variable& copy)
  {
    MeshNodeVariable<T>* c = dynamic_cast<MeshNodeVariable<T>* >(&copy);
    if(!c) std::cerr << "Type mismatch in mesh node variable" << std::endl;
    copyPointer(*c);
  }
  
  // specialization for T=Point
  template<class T>
  void
  MeshNodeVariable<T>::gather(MeshNodeSet* nodeset,
                              const std::vector<MeshNodeSet*> &subsets,
                              const std::vector<MeshNodeVariableBase*> &srcs,
                              meshNodeIndex extra)
  {
    if(d_nodedata && d_nodedata->removeReference())
      delete d_nodedata;
    if(d_nodeset && d_nodeset->removeReference())
      delete d_nodeset;
    d_nodeset = nodeset;
    nodeset->addReference();
    d_nodedata=new MeshNodeData<T>(nodeset->numNodes());
    d_nodedata->addReference();
    //ASSERTEQ(subsets.size(), srcs.size());
    MeshNodeSet::iterator dstiter = nodeset->begin();
    for(int i=0;i<(int)subsets.size();i++){
      MeshNodeVariable<T>* srcptr = dynamic_cast<MeshNodeVariable<T>*>(srcs[i]);
      if(!srcptr) std::cerr << "Type mismatch in MeshNodeVariable::gather" << std::endl;
      MeshNodeVariable<T>& src = *srcptr;
      MeshNodeSet* subset = subsets[i];
      for(MeshNodeSet::iterator srciter = subset->begin();
          srciter != subset->end(); srciter++){
        (*this)[*dstiter] = src[*srciter];
        dstiter++;
      }
    }
    //ASSERT(dstiter+extra == nodeset->end());
    extra = extra;   // This is to shut up the REMARKS from the MIPS compiler
  }
  
  template<class T>
  void*
  MeshNodeVariable<T>::getBasePointer() const
  {
    return &d_nodedata->data[0];
  }
  
  // Specialized in MeshNodeVariable_special.cc
  template<>
   void
   MeshNodeVariable<double>::emitNormal(ostream& out); 

  template<class T>
  void
  MeshNodeVariable<T>::emitNormal(ostream& out)
  {
    // This could be optimized...
    MeshNodeSet::iterator iter = d_nodeset->begin();
    while(iter != d_nodeset->end()){
      meshNodeIndex start = *iter;
      iter++;
      meshNodeIndex end = start+1;
      while(iter != d_nodeset->end() && *iter == end) {
        end++;
        iter++;
      }
      ssize_t size = (ssize_t)(sizeof(T)*(end-start));
      out.write((char*)&(*this)[start], size);
    }
  }

  template<class T>
  void
  MeshNodeVariable<T>::readNormal(istream& in)
  {
    // This could be optimized...
    MeshNodeSet::iterator iter = d_nodeset->begin();
    while(iter != d_nodeset->end()){
      meshNodeIndex start = *iter;
      iter++;
      meshNodeIndex end = start+1;
      while(iter != d_nodeset->end() && *iter == end) {
        end++;
        iter++;
      }
      ssize_t size = (ssize_t)(sizeof(T)*(end-start));
      in.read((char*)&(*this)[start], size);
    }
  }

  template <class T>
  class constMeshNodeVariable : public constVariable<MeshNodeVariableBase, MeshNodeVariable<T>, T, meshNodeIndex>
  {
  public:
    constMeshNodeVariable()
      : constVariable<MeshNodeVariableBase, MeshNodeVariable<T>, T, meshNodeIndex>() {}
    
    constMeshNodeVariable(const MeshNodeVariable<T>& copy)
      : constVariable<MeshNodeVariableBase, MeshNodeVariable<T>, T, meshNodeIndex>(copy) {}

    MeshNodeSet* getMeshNodeSet() const {
      return this->rep_.getMeshNodeSet();
    }
  };

} // End namespace Matiti

#endif
