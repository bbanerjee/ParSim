#ifndef __MATITI_DAMAGE_MODEL_BASE_H__
#define __MATITI_DAMAGE_MODEL_BASE_H__

#include <Pointers/DamageModelSP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {
  
  class DamageModelBase {
  
  public:

    DamageModelBase();
    ~DamageModelBase();

    /**
     *  Initialize the damage model
     */
    virtual void initialize(Uintah::ProblemSpecP& ps) = 0;

    /* Make copies without having to read in stuff */
    virtual DamageModelSP clone() = 0;
    virtual void setVariation(double randomNumber, double coeffOfVar) = 0;
    virtual void setAverage(const DamageModelBase* dam1, const DamageModelBase* dam2) = 0;

    //virtual void clone(const DamageModelBase* dam);
    //virtual void clone(const DamageModelBase* dam, double randomNumber, double coeffOfVar);
    //virtual void cloneAverage(const DamageModelBase* dam1, const DamageModelBase* dam2);

    /**
     *  Compute the damage factor  given the damage index at a node.
     */
    virtual double computeDamageFactor(const double& damage_index) const;


  private:

    // prevent copying
    DamageModelBase(const DamageModelBase& bc);
    DamageModelBase& operator=(const DamageModelBase& bc);

  }; // end class

} // end namespace
#endif

