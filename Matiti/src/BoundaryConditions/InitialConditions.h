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

#ifndef MATITI_INITIALCONDITIONS_H
#define MATITI_INITIALCONDITIONS_H

#include <Containers/NodePArray.h>
#include <Containers/CrackSPArray.h>
#include <Geometry/Vector3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {
  
  class InitialConditions {
  
  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::InitialConditions& ic);

  public:

    InitialConditions();
    InitialConditions(const Vector3D& initialVel,
                      const Vector3D& gravity);
    ~InitialConditions();

    void readInitialConditions(Uintah::ProblemSpecP& ps);

    void applyInitialVelocity(NodePArray& nodes);
    void applyBodyForce(NodePArray& nodes);
    void removeBondsIntersectedByCracks(NodePArray& nodes);

    /**
     * Get methods
     */
    const Vector3D& initialVelocity() const {return d_initial_velocity;}
    const Vector3D& bodyForce() const {return d_body_force;}
    const CrackSPArray& cracks() const {return d_cracks;}

    /**
     * Set methods
     */
    void initialVelocity(const Vector3D& vel) {d_initial_velocity = vel;}
    void bodyForce(const Vector3D& gravity) {d_body_force = gravity;}

  private:

    Vector3D d_initial_velocity; // Initial velocity
    Vector3D d_body_force;       // Gravity (essentially)
    CrackSPArray d_cracks;       // Initial crack geometries

    // prevent copying
    InitialConditions(const InitialConditions& bc);
    InitialConditions& operator=(const InitialConditions& bc);

  }; // end class

} // end namespace
#endif

