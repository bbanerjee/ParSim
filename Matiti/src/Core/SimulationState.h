/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MATITI_SIMULATION_STATE_H
#define MATITI_SIMULATION_STATE_H

#include <Types/Types.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <iostream>

namespace Matiti {

  class SimulationState {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::SimulationState& state);

    enum class ModulusType {
      Constant=0,
      Conical=1
    };

  public:
 
    SimulationState();
    SimulationState(bool isDynamic, Array3 dampingCoeff, 
                    ModulusType modulusType, double horizonFactor);
    ~SimulationState();
    void clone(const SimulationState& state);

    void initialize(const Uintah::ProblemSpecP& ps);

    inline bool isDynamic() const {return d_is_dynamic;}
    inline ModulusType modulusType() const {return d_modulus_type;}
    inline double horizonFactor() const {return d_horizon_factor;}
    inline Array3 dampingCoeff() const {return d_damping_coeff;}

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
