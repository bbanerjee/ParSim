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

#ifndef __VAANGO_NEIGHBOR_BOND_ENERGY_H__
#define __VAANGO_NEIGHBOR_BOND_ENERGY_H__

#include <Core/Util/Assert.h>
#include <Core/Datatypes/TypeName.h>
#include <Core/Disclosure/TypeUtils.h>  // Contains double and ParticleID

#include <cmath>
#include <iosfwd>
#include <vector>
#include <string>

namespace SCIRun {
  class TypeDescription;
  class Piostream;
}


namespace Uintah {


  class NeighborBondEnergy 
  {
  public: 

    friend std::ostream& operator<<(std::ostream& out, 
                                    const Uintah::NeighborBondEnergy& energy);

  private:

    double d_bondEnergy[216];  // 6 x 6 x 6 
    
  public: 

    /**
     * Constructors/Destructors
     **/
    inline NeighborBondEnergy();
    inline NeighborBondEnergy(std::vector<double>& bondEnergy);
    inline ~NeighborBondEnergy();

    /**
     * Access methods
     */
    inline double operator[](int ii) const;
    inline double& operator[](int ii);
    static const std::string& get_h_file_path();

  };  // end class

  
  inline NeighborBondEnergy::NeighborBondEnergy()
  {
    for (int ii = 0; ii < 216; ii++) {
      d_bondEnergy[ii] = -1.0;
    }
  }

  inline NeighborBondEnergy::NeighborBondEnergy(std::vector<double>& bondEnergy)
  {
    int ii = 0;
    for (auto iter = bondEnergy.begin(); iter != bondEnergy.end(); ++iter) {
      d_bondEnergy[ii++] = *iter;
    }
  }

  inline NeighborBondEnergy::~NeighborBondEnergy()
  {
  }

  inline double NeighborBondEnergy::operator[](int ii) const
  {
    return d_bondEnergy[ii];
  }

  inline double& NeighborBondEnergy::operator[](int ii) 
  {
    return d_bondEnergy[ii];
  }
} // end namespace

// Added for compatibility with core types
#include <Core/Datatypes/TypeName.h>
#include <string>
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborBondEnergy& broken);
  template<>  const std::string find_type_name(Uintah::NeighborBondEnergy*);
  const TypeDescription* get_type_description(Uintah::NeighborBondEnergy*);
  void Pio( Piostream&, Uintah::NeighborBondEnergy& );
} // namespace SCIRun


#endif
