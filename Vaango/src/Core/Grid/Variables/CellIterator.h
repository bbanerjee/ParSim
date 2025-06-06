/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */


#ifndef UINTAH_HOMEBREW_CellIterator_H
#define UINTAH_HOMEBREW_CellIterator_H

#include <Core/Geometry/IntVector.h>
#include <Core/Util/Assert.h>
#include <iterator>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Variables/BaseIterator.h>

namespace Uintah {

using Uintah::IntVector;
 using std::ostream;

/**************************************

CLASS
   CellIterator
   
   Short Description...

GENERAL INFORMATION

   CellIterator.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  

KEYWORDS
   CellIterator

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

 class CellIterator : 
 public BaseIterator {
   public:
     inline ~CellIterator() {}

     //////////
     // Insert Documentation Here:
     inline void operator++(int) override {
       this->operator++();
     }

     inline CellIterator& operator++() override {
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

     //////////
     // Insert Documentation Here:
     inline CellIterator operator+=(int step) {
       CellIterator old(*this);

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

     //////////
     // Insert Documentation Here:
     inline bool done() const override {
       return d_done;
     }

     IntVector operator*() const override {
       ASSERT(!d_done);
       return d_cur;
     }
     inline CellIterator(const IntVector& s, const IntVector& e)
       : d_s(s), d_e(e){
       reset();
     }
     inline IntVector begin() const override {
       return d_s;
     }
     inline IntVector end() const override {
       return d_e;
     }
     /**
     * Return the number of cells in the iterator
     */
     inline unsigned int size() const override
     {
       IntVector size=d_e-d_s;
       if(size.x()<=0 || size.y()<=0 || size.z()<=0)
         return 0;
       else
         return size.x()*size.y()*size.z();
     };
     inline CellIterator(const CellIterator& copy)
       : d_s(copy.d_s), d_e(copy.d_e), d_cur(copy.d_cur), d_done(copy.d_done) {
       }

     inline CellIterator& operator=( const CellIterator& copy ) {
       d_s    = copy.d_s;
       d_e    = copy.d_e;
       d_cur  = copy.d_cur;
       d_done = copy.d_done;
       return *this;
     }
     bool operator==(const CellIterator& o) const
     {
       return begin()==o.begin() && end()==o.end() && d_cur==o.d_cur;
     }

     bool operator!=(const CellIterator& o) const
     {
       return begin()!=o.begin() || end()!=o.end() || d_cur!=o.d_cur;
     }
     friend std::ostream& operator<<(std::ostream& out, const Uintah::CellIterator& b);

     friend class GridIterator;

     inline void reset() override
     {
       d_cur=d_s;
       d_done=d_s.x() >= d_e.x() || d_s.y() >= d_e.y() || d_s.z() >= d_e.z();
     }

     ostream& limits(std::ostream& out) const override
     {
       out << begin() << " " << end() - IntVector(1,1,1);
       return out;
     }
   private:
     CellIterator();

     CellIterator* clone() const override
     {
       return scinew CellIterator(*this);
     }

     ostream& put(std::ostream& out) const override
     {
       out << *this;
       return out;
     }
     //////////
     // Insert Documentation Here:
     IntVector d_s,d_e;
     IntVector d_cur;
     bool d_done;

}; // end class CellIterator

} // End namespace Uintah
  
#endif
