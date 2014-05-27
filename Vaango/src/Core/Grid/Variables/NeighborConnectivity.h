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
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Uintah::NeighborConnectivity& broken);
  template<>  const std::string find_type_name(Uintah::NeighborConnectivity*);
  const TypeDescription* get_type_description(Uintah::NeighborConnectivity*);
  void Pio( Piostream&, Uintah::NeighborConnectivity& );
} // namespace SCIRun


#endif
