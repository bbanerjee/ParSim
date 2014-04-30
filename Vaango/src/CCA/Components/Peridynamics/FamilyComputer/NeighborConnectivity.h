#ifndef __VAANGO_NEIGHBOR_CONNECTIVITY_H__
#define __VAANGO_NEIGHBOR_CONNECTIVITY_H__

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

  class NeighborConnectivity 
  {
  private:

    bool d_broken[216];  // 6 x 6 x 6 
    
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
      d_broken[ii] = false; 
    }
  }

  inline NeighborConnectivity::~NeighborConnectivity()
  {
  }

  inline bool NeighborConnectivity::operator[](int ii) const
  {
    return d_broken[ii];
  }

  inline bool& NeighborConnectivity::operator[](int ii) const
  {
    return d_broken[ii];
  }
} // end namespace

// Added for compatibility with core types
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Vaango::NeighborConnectivity& broken);
  template<>  const std::string find_type_name(Vaango::NeighborConnectivity*);
  const TypeDescription* get_type_description(Vaango::NeighborConnectivity*);
  void Pio( Piostream&, Vaango::NeighborConnectivity& );
} // namespace SCIRun


#endif
