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

#ifndef __VAANGO_NEIGHBOR_CONNECTIVITY_H__
#define __VAANGO_NEIGHBOR_CONNECTIVITY_H__

#include <Core/Util/Assert.h>
#include <Core/Datatypes/TypeName.h>
#include <Core/Disclosure/TypeUtils.h>  // Contains long64 and ParticleID

#include <cmath>
#include <iosfwd>
#include <vector>
#include <string>
#include <iostream>

namespace Uintah {

  class NeighborConnectivity 
  {
  public: 

    friend std::ostream& operator<<(std::ostream& out, 
                                    const Uintah::NeighborConnectivity& conn);

  private:

    bool d_connected[216];  // 6 x 6 x 6 
    
  public: 

    /**
     * Constructors/Destructors
     **/
    inline NeighborConnectivity();
    inline ~NeighborConnectivity();

    /**
     * Access methods
     */
    inline bool operator[](int ii) const;
    inline bool& operator[](int ii);
    static const std::string& get_h_file_path();

  };  // end class

  
  inline NeighborConnectivity::NeighborConnectivity()
  {
    for (int ii = 0; ii < 216; ii++) {
      d_connected[ii] = false; 
    }
  }

  inline NeighborConnectivity::~NeighborConnectivity()
  {
  }

  inline bool NeighborConnectivity::operator[](int ii) const
  {
    return d_connected[ii];
  }

  inline bool& NeighborConnectivity::operator[](int ii) 
  {
    return d_connected[ii];
  }
} // end namespace

// Added for compatibility with core types
#include <Core/Datatypes/TypeName.h>
#include <string>
namespace Uintah {

  void swapbytes(Uintah::NeighborConnectivity& broken);
  template<>  const std::string find_type_name(Uintah::NeighborConnectivity*);
  
} // namespace Uintah


#endif
