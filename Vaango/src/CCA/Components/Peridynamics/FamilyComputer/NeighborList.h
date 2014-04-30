#ifndef __VAANGO_NEIGHBOR_LIST_H__
#define __VAANGO_NEIGHBOR_LIST_H__

#include <Core/Util/Assert.h>
#include <Core/Datatypes/TypeName.h>
#include <Core/Disclosure/TypeUtils.h>  // Contains long64 and ParticleID

#include <cmath>
#include <iosfwd>
#include <vector>
#include <string>

namespace SCIRun {
  class TypeDescription;
  class Piostream;
}


namespace Vaango {

  class NeighborList 
  {
  private:

    Uintah::ParticleID d_family[216];  // 6 x 6 x 6 
    
  public: 

    /**
     * Constructors/Destructors
     **/
    inline NeighborList();
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
      d_family[ii] = false; 
    }
  }

  inline NeighborList::~NeighborList()
  {
  }

  inline Uintah::ParticleID NeighborList::operator[](int ii) const
  {
    return d_family[ii];
  }

  inline Uintah::ParticleID& NeighborList::operator[](int ii) const
  {
    return d_family[ii];
  }
} // end namespace

// Added for compatibility with core types
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Vaango::NeighborList& broken);
  template<>  const std::string find_type_name(Vaango::NeighborList*);
  const TypeDescription* get_type_description(Vaango::NeighborList*);
  void Pio( Piostream&, Vaango::NeighborList& );
} // namespace SCIRun


#endif
