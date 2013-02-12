#ifndef MATITI_MESHNODE_DATA_H
#define MATITI_MESHNODE_DATA_H

#include <Common/RefCounted.h>

namespace Matiti {

template<class T>
   class MeshNodeVariable;

   template<class T> class MeshNodeData : public RefCounted {
   public:
      MeshNodeData();
      MeshNodeData(meshNodeIndex size);
      virtual ~MeshNodeData();

      void resize(int newSize) {
        T* newdata = new T[newSize];
        if(data){
          int smaller = ((newSize < size ) ? newSize:size);
          for(int i = 0; i < smaller; i++)
            newdata[i] = data[i];
          delete[] data;
        }
        data = newdata;
        size = newSize;
      }

   private:
      MeshNodeData(const MeshNodeData<T>&);
      MeshNodeData<T>& operator=(const MeshNodeData<T>&);
      friend class MeshNodeVariable<T>;
      
      T* data;
      meshNodeIndex size;
   };
   
   template<class T>
      MeshNodeData<T>::MeshNodeData()
      {
        data=0;
      }
   
   template<class T>
     MeshNodeData<T>::MeshNodeData(meshNodeIndex size)
     : size(size)
      {
        data = new T[size];
      }
      
   template<class T>
      MeshNodeData<T>::~MeshNodeData()
      {
        if(data)
          delete[] data;
      }

   template<class T>
     MeshNodeData<T>& MeshNodeData<T>::operator=(const MeshNodeData<T>& copy)
     {
       for(meshNodeIndex i=0;i<size;i++)
         data[i] = copy.data[i];
       return *this;
     }
} // End namespace Matiti
   
#endif
