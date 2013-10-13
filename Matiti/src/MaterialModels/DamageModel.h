#ifndef MATITI_DAMAGE_MODEL_H
#define MATITI_DAMAGE_MODEL_H

#include <Pointers/NodeP.h>
#include <Pointers/DamageModelUP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Geometry/Vector3D.h>
#include <iostream>

namespace Matiti {

  class DamageModel 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::DamageModel& dam);

  public:
  
    DamageModel();
    DamageModel(const DamageModel& dam);
    virtual ~DamageModel();

    void clone(const DamageModelUP& dam);
    void clone(const DamageModelUP& dam,
               double randomNumber,
               double coeffOfVar);
    void cloneAverage(const DamageModelUP& dam1, const DamageModelUP& dam2);

    void initialize(const Uintah::ProblemSpecP& ps);

    /**
     *  Compute the damage factor from the damage stretch coefficients
     *  given the damage index at a node.
     *  **WARNING** Don't ubderstand this.
     */
    double computeDamageFactor(const double& damage_index) const;

    const Vector3D& damageViscosity() const {return d_damage_viscosity;}
    const Vector3D& damageStretch() const {return d_damage_stretch;}

  protected:

    Vector3D d_damage_viscosity;
    Vector3D d_damage_stretch;
    double d_damage_index_max; // Not used yet

  }; // end class

} // end namespace

#endif
