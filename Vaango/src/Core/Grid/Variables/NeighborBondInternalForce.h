#ifndef __VAANGO_NEIGHBOR_BOND_INTERNAL_FORCE_H__
#define __VAANGO_NEIGHBOR_BOND_INTERNAL_FORCE_H__

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

  using SCIRun::Vector;

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
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborBondInternalForce& force);
  template<>  const std::string find_type_name(Uintah::NeighborBondInternalForce*);
  const TypeDescription* get_type_description(Uintah::NeighborBondInternalForce*);
  void Pio( Piostream&, Uintah::NeighborBondInternalForce& );
} // namespace SCIRun


#endif
