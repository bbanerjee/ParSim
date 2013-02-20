#ifndef MATITI_MESHELEMENTITERATOR_H
#define MATITI_MESHELEMENTITERATOR_H

#include <Core/Grid/Variables/BaseIterator.h>
#include <Core/Geometry/IntVector.h>
#include <iterator>

namespace Matiti {

 using SCIRun::IntVector;
 using std::ostream;


 class MeshElementIterator : public BaseIterator 
 {
   public:
     inline ~MeshElementIterator() {}

     inline void operator++(int) {
       this->operator++();
     }

     inline MeshElementIterator& operator++() {
       if(++d_cur.modifiable_x() >= d_e.x()){
         d_cur.modifiable_x() = d_s.x();
         if(++d_cur.modifiable_y() >= d_e.y()){
           d_cur.modifiable_y() = d_s.y();
           ++d_cur.modifiable_z();
           if(d_cur.modifiable_z() >= d_e.z())
             d_done=true;
         }
       }
       return *this;
     }

     inline MeshElementIterator operator+=(int step) {
       MeshElementIterator old(*this);

       for (int i = 0; i < step; i++) {
         if(++d_cur.modifiable_x() >= d_e.x()){
           d_cur.modifiable_x() = d_s.x();
           if(++d_cur.modifiable_y() >= d_e.y()){
             d_cur.modifiable_y() = d_s.y();
             ++d_cur.modifiable_z();
             if(d_cur.modifiable_z() >= d_e.z())
               d_done=true;
           }
         }
         if (done())
           break;
       }
       return old;
     }

     inline bool done() const {
       return d_done;
     }

     IntVector operator*() const {
       //ASSERT(!d_done);
       return d_cur;
     }
     inline MeshElementIterator(const IntVector& s, const IntVector& e)
       : d_s(s), d_e(e){
       reset();
     }
     inline IntVector begin() const {
       return d_s;
     }
     inline IntVector end() const {
       return d_e;
     }
     /**
     * Return the number of cells in the iterator
     */
     inline unsigned int size() const
     {
       IntVector size=d_e-d_s;
       if(size.x()<=0 || size.y()<=0 || size.z()<=0)
         return 0;
       else
         return size.x()*size.y()*size.z();
     };
     inline MeshElementIterator(const MeshElementIterator& copy)
       : d_s(copy.d_s), d_e(copy.d_e), d_cur(copy.d_cur), d_done(copy.d_done) {
       }

     inline MeshElementIterator& operator=( const MeshElementIterator& copy ) {
       d_s    = copy.d_s;
       d_e    = copy.d_e;
       d_cur  = copy.d_cur;
       d_done = copy.d_done;
       return *this;
     }
     bool operator==(const MeshElementIterator& o) const
     {
       return begin()==o.begin() && end()==o.end() && d_cur==o.d_cur;
     }

     bool operator!=(const MeshElementIterator& o) const
     {
       return begin()!=o.begin() || end()!=o.end() || d_cur!=o.d_cur;
     }
     friend std::ostream& operator<<(std::ostream& out, const Matiti::MeshElementIterator& b);

     friend class MeshNodeIterator;

     inline void reset()
     {
       d_cur=d_s;
       d_done=d_s.x() >= d_e.x() || d_s.y() >= d_e.y() || d_s.z() >= d_e.z();
     }

     ostream& limits(ostream& out) const
     {
       out << begin() << " " << end() - IntVector(1,1,1);
       return out;
     }
   private:
     MeshElementIterator();

     MeshElementIterator* clone() const
     {
       return new MeshElementIterator(*this);
     }

     ostream& put(ostream& out) const
     {
       out << *this;
       return out;
     }

     IntVector d_s,d_e;
     IntVector d_cur;
     bool d_done;

}; // end class MeshElementIterator

} // End namespace Matiti
  
#endif
