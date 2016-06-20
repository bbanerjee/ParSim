/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __VAANGO_NEIGHBOR_LIST_H__
#define __VAANGO_NEIGHBOR_LIST_H__

#include <Core/Util/Assert.h>
#include <Core/Datatypes/TypeName.h>
#include <Core/Disclosure/TypeUtils.h>  // Contains long64 and ParticleID

#include <cmath>
#include <iosfwd>
#include <vector>
#include <string>

namespace Uintah {
  class TypeDescription;
  class Piostream;
}


namespace Uintah {

  class NeighborList 
  {
  public: 

    friend std::ostream& operator<<(std::ostream& out, 
                                    const Uintah::NeighborList& family);

  private:

    Uintah::ParticleID d_family[216];  // 6 x 6 x 6 
    
  public: 

    /**
     * Constructors/Destructors
     **/
    inline NeighborList();
    inline NeighborList(std::vector<long64>& family);
    inline ~NeighborList();

    /**
     * Access methods
     */
    inline Uintah::ParticleID operator[](int ii) const;
    inline Uintah::ParticleID& operator[](int ii);
    static const std::string& get_h_file_path();

  };  // end class

  inline NeighborList::NeighborList()
  {
    for (int ii = 0; ii < 216; ii++) {
      d_family[ii] = (long64) 0;
    }
  }
  
  inline NeighborList::NeighborList(std::vector<long64>& family)
  {
    for (int ii = 0; ii < 216; ii++) {
      d_family[ii] = (long64) 0;
    }
  
    int ii = 0;
    for (auto iter = family.begin(); iter != family.end(); ++iter) {
      d_family[ii] = *iter;
      ii++;
    }
  }
  
  inline NeighborList::~NeighborList()
  {
  }
  
  inline Uintah::ParticleID 
  NeighborList::operator[](int ii) const
  {
    return d_family[ii];
  }
  
  inline Uintah::ParticleID& 
  NeighborList::operator[](int ii) 
  {
    return d_family[ii];
  }

} // end namespace

// Added for compatibility with core types
#include <Core/Datatypes/TypeName.h>
#include <string>
namespace Uintah {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborList& broken);
  template<>  const std::string find_type_name(Uintah::NeighborList*);
  const TypeDescription* get_type_description(Uintah::NeighborList*);
  void Pio( Piostream&, Uintah::NeighborList& );
} // namespace Uintah


#endif
