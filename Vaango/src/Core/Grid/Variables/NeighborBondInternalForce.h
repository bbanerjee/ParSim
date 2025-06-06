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

#ifndef __VAANGO_NEIGHBOR_BOND_INTERNAL_FORCE_H__
#define __VAANGO_NEIGHBOR_BOND_INTERNAL_FORCE_H__

#include <Core/Util/Assert.h>
#include <Core/Datatypes/TypeName.h>
#include <Core/Disclosure/TypeUtils.h>  // Contains double and ParticleID

#include <cmath>
#include <iosfwd>
#include <vector>
#include <string>

namespace Uintah {
  class TypeDescription;
  class Piostream;
}


namespace Uintah {

  using Uintah::Vector;

  class NeighborBondInternalForce 
  {
  public: 

    friend std::ostream& operator<<(std::ostream& out, 
                                    const Uintah::NeighborBondInternalForce& force);

  private:

    Vector d_bondInternalForce[216];  // 6 x 6 x 6 
    
  public: 

    /**
     * Constructors/Destructors
     **/
    inline NeighborBondInternalForce();
    inline NeighborBondInternalForce(std::vector<Vector>& bondInternalForce);
    inline ~NeighborBondInternalForce();

    /**
     * Access methods
     */
    inline Vector operator[](int ii) const;
    inline Vector& operator[](int ii);
    static const std::string& get_h_file_path();

  };  // end class

  
  inline NeighborBondInternalForce::NeighborBondInternalForce()
  {
    for (int ii = 0; ii < 216; ii++) {
      d_bondInternalForce[ii] = Vector(0.0, 0.0, 0.0);
    }
  }

  inline NeighborBondInternalForce::NeighborBondInternalForce(std::vector<Vector>& bondInternalForce)
  {
    int ii = 0;
    for (auto iter = bondInternalForce.begin(); iter != bondInternalForce.end(); ++iter) {
      d_bondInternalForce[ii++] = *iter;
    }
  }

  inline NeighborBondInternalForce::~NeighborBondInternalForce()
  {
  }

  inline Vector NeighborBondInternalForce::operator[](int ii) const
  {
    return d_bondInternalForce[ii];
  }

  inline Vector& NeighborBondInternalForce::operator[](int ii) 
  {
    return d_bondInternalForce[ii];
  }
} // end namespace

// Added for compatibility with core types
#include <Core/Datatypes/TypeName.h>
#include <string>
namespace Uintah {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborBondInternalForce& force);
  template<>  const std::string find_type_name(Uintah::NeighborBondInternalForce*);
  const FETypeDescription* get_fetype_description(Uintah::NeighborBondInternalForce*);
  void Pio( Piostream&, Uintah::NeighborBondInternalForce& );
} // namespace Uintah


#endif
