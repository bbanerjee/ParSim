#ifndef MATITI_SIMULATION_STATE_H
#define MATITI_SIMULATION_STATE_H

#include <Types.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <iostream>

namespace Matiti {

  class SimulationState {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::SimulationState& state);

  public:
 
    SimulationState();
    ~SimulationState();

    void initialize(const Uintah::ProblemSpecP& ps);

    enum class ModulusType {
      Constant=0,
      Conical=1
    };
    inline bool isDynamic() const {return d_is_dynamic;}
    inline ModulusType modulusType() const {return d_modulus_type;}
    inline double horizonFactor() const {return d_horizon_factor;}

  private:


    // Peridynamics flags
    bool d_is_dynamic;          // static or dyanamic simulation
    Array3 d_damping_coeff;     // Artificial viscosity coeffs
    ModulusType d_modulus_type; // constant = 0 or conical = 1
    double d_horizon_factor;    // Factor with which to multiply the nodal volume to get
                                // horizon size 
  }; // end class
   
} // end namespace

#endif
