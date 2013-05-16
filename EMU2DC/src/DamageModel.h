#ifndef EMU2DC_DAMAGE_MODEL_H
#define EMU2DC_DAMAGE_MODEL_H

#include <Types.h>
#include <NodeP.h>
#include <DamageModelUP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace Emu2DC {

  class DamageModel 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::DamageModel& dam);

  public:
  
    DamageModel();
    DamageModel(const DamageModel& dam);
    virtual ~DamageModel();

    void clone(const DamageModelUP& dam);

    void initialize(const Uintah::ProblemSpecP& ps);

    /**
     *  Compute the damage factor from the damage stretch coefficients
     *  given the damage index at a node.
     *  **WARNING** Don't ubderstand this.
     */
    double computeDamageFactor(const double& damage_index) const;

    const Array3& damageViscosity() const {return d_damage_viscosity;}
    const Array3& damageStretch() const {return d_damage_stretch;}

  protected:

    Array3 d_damage_viscosity;
    Array3 d_damage_stretch;
    double d_damage_index_max; // Not used yet

  }; // end class

} // end namespace

#endif
