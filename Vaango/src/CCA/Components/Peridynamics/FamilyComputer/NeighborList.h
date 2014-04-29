#ifndef __VAANGO_NEIGHBORLIST_H__
#define __VAANGO_NEIGHBORLIST_H__

#include <CCA/Components/Peridynamics/FamilyComputer/Bond.h>
#include <Core/Datatypes/TypeName.h>

#include <vector>
#include <string>

namespace Vaango {

  class NeighborList 
  {
  public: 

    NeighborList();
    virtual ~NeighborList();

    /**
     * Set methods
     */
    void addNeighbor(const Bond* bond) {
      d_bonds.push_back(bond);
    }

    void removeNeighbor(const Bond* bond) {
      std::vector<Bond*>::iterator iter = d_bonds.begin();
      for (; iter != d_bonds.end(); ++iter) {
        if (*iter == bond) {
          (*iter)->isBroken(true);
          return; 
        }
      }
    }

    /**
     * Get methods
     */
    std::vector<Bond*>::iterator begin() {
      return d_bonds.begin();
    }
    
    std::vector<Bond*>::iterator end() {
      return d_bonds.end();
    }

    std::vector<Bond*>::const_iterator begin() const {
      return d_bonds.begin();
    }
    
    std::vector<Bond*>::const_iterator end() const {
      return d_bonds.end();
    }

  private:

    std::vector<Bond*> d_bonds;

  };  // end class

} // end namespace

// Added for compatibility with core types
namespace SCIRun {

  class TypeDescription;
  class Piostream;

  void swapbytes(Vaango::NeighborList& family);
  template<>  const std::string find_type_name(Vaango::NeighborList*);
  const TypeDescription* get_type_description(Vaango::NeighborList*);
  void Pio( Piostream&, Vaango::NeighborList& );
} // namespace SCIRun


#endif
