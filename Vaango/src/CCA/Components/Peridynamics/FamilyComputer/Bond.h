#ifndef VAANGO_BOND_H
#define VAANGO_BOND_H

#include <Core/Disclosure/TypeUtils.h>  // Contains long64 and ParticleID

#include <iostream>

namespace Vaango {

  class Bond 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Vaango::Bond& bond);

  public: 

    Bond();
    Bond(const Uintah::ParticleID& start, const Uintah::ParticleID& end);
    virtual ~Bond();

    /**
     * Set methods
     */
    void first(const Uintah::ParticleID& particle) { d_start = particle; }
    void second(const Uintah::ParticleID& particle) { d_end = particle; }
    void isBroken(bool broken) { d_broken = broken;}

    /**
     * Get methods
     */
    const Uintah::ParticleID& first() const { return d_start; }
    const Uintah::ParticleID& second() const { return d_end; }
    bool isBroken() const { return d_broken; }
    
    /**
     * Check if two bonds are identical
     *   True if the start and end points are the same 
     */
    bool operator==(const Bond& bond) const;

  private:

    Uintah::ParticleID d_start;
    Uintah::ParticleID d_end;
    bool d_broken;

  };  // end class

} // end namespace

#endif
