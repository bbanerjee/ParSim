#ifndef MATITI_MESHNODEITERATOR_H
#define MATITI_MESHNODEITERATOR_H

#include <Core/Grid/Variables/BaseIterator.h>
#include <Core/Geometry/IntVector.h>

namespace Matiti {

  using SCIRun::IntVector;

  class MeshNodeIterator  :
    public BaseIterator {
      public:
        inline ~MeshNodeIterator() {}

        //////////
        // Insert Documentation Here:
        inline MeshNodeIterator& operator++() {

          if(++d_ix >= d_e.x()){
            d_ix = d_s.x();
            if(++d_iy >= d_e.y()){
              d_iy = d_s.y();
              ++d_iz;
            }
          }
          return *this;
        }

        //////////
        // Insert Documentation Here:
        inline void operator++(int) {
          this->operator++();
        }

        //////////
        // Insert Documentation Here:
        inline MeshNodeIterator operator+=(int step) {
          MeshNodeIterator old(*this);

          for (int i = 0; i < step; i++) {
            if(++d_ix >= d_e.x()){
              d_ix = d_s.x();
              if(++d_iy >= d_e.y()){
                d_iy = d_s.y();
                ++d_iz;
              }
            }
            if (done())
              break;
          }
          return old;
        }

        inline bool done() const {
          return d_ix >= d_e.x() || d_iy >= d_e.y() || d_iz >= d_e.z();
        }

        IntVector operator*() const {
          return IntVector(d_ix, d_iy, d_iz);
        }
        IntVector index() const {
          return IntVector(d_ix, d_iy, d_iz);
        }
        inline MeshNodeIterator(const IntVector& s, const IntVector& e)
          : d_s(s), d_e(e) {
            d_ix = s.x();
            d_iy = s.y();
            d_iz = s.z();
          }
        inline IntVector current() const {
          return IntVector(d_ix, d_iy, d_iz);
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
        inline MeshNodeIterator(const MeshNodeIterator& copy)
          : d_s(copy.d_s), d_e(copy.d_e),
          d_ix(copy.d_ix), d_iy(copy.d_iy), d_iz(copy.d_iz) {
          }

        friend class GridIterator;
        
        inline MeshNodeIterator& operator=( const MeshNodeIterator& copy ) {
          d_s   = copy.d_s;
          d_e   = copy.d_e;
          d_ix  = copy.d_ix;
          d_iy  = copy.d_iy;
          d_iz  = copy.d_iz;
          return *this;
        }

        bool operator==(const MeshNodeIterator& o) const
        {
          return begin()==o.begin() && end()==o.end() && index()==o.index();
        }

        bool operator!=(const MeshNodeIterator& o) const
        {
          return begin()!=o.begin() || end()!=o.end() || index()!=o.index();
        }

        inline void reset()
        {
          d_ix=d_s.x();
          d_iy=d_s.y();
          d_iz=d_s.z();
        }

        std::ostream& limits(std::ostream& out) const
        {
          out << begin() << " " << end() - IntVector(1,1,1);
          return out;
        }

      private:
        MeshNodeIterator();

        MeshNodeIterator* clone() const
        {
          return new MeshNodeIterator(*this);
        }
        
        std::ostream& put(std::ostream& out) const
        {
          out << *this;
          return out;
        }

        IntVector d_s,d_e;
        int d_ix, d_iy, d_iz;
        
        friend std::ostream& operator<<(std::ostream& out, const Matiti::MeshNodeIterator& b);
    }; // end class MeshNodeIterator

} // End namespace Matiti

#endif
