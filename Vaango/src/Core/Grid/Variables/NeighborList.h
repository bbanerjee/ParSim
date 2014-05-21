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


namespace Uintah {

  class NeighborList 
  {
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
      d_family[ii] = 0;
    }
  }

  inline NeighborList::NeighborList(std::vector<long64>& family)
  {
    int ii = 0;
    for (auto iter = family.begin(); iter != family.end(); ++iter) {
      d_family[ii++] = *iter;
    }
  }

  inline NeighborList::~NeighborList()
  {
  }

  inline Uintah::ParticleID NeighborList::operator[](int ii) const
  {
    return d_family[ii];
  }

  inline Uintah::ParticleID& NeighborList::operator[](int ii) 
  {
    return d_family[ii];
  }
} // end namespace

// Added for compatibility with core types
#include <Core/Datatypes/TypeName.h>
#include <string>
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborList& broken);
  template<>  const std::string find_type_name(Uintah::NeighborList*);
  const TypeDescription* get_type_description(Uintah::NeighborList*);
  void Pio( Piostream&, Uintah::NeighborList& );
} // namespace SCIRun


#endif
