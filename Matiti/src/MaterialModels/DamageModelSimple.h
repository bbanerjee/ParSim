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

#ifndef MATITI_DAMAGE_MODEL_SIMPLE_H
#define MATITI_DAMAGE_MODEL_SIMPLE_H

#include <MaterialModels/DamageModelBase.h>
#include <Pointers/NodeP.h>
#include <Pointers/DamageModelSimpleSP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Geometry/Vector3D.h>
#include <iostream>

namespace Matiti {

  class DamageModelSimple : public DamageModelBase
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::DamageModelSimple& dam);

  public:
  
    DamageModelSimple();
    DamageModelSimple(const DamageModelSimple& dam);
    virtual ~DamageModelSimple();

    DamageModelSP clone();
    void setVariation(double randomNumber, double coeffOfVar);
    //void setAverage(const DamageModelSimple* dam1, const DamageModelSimple* dam2);
    void setAverage(const DamageModelBase* dam1, const DamageModelBase* dam2);

    void clone(const DamageModelSimple* dam);
    void clone(const DamageModelSimple* dam,
               double randomNumber,
               double coeffOfVar);
    void cloneAverage(const DamageModelSimple* dam1, 
                      const DamageModelSimple* dam2);

    void initialize(Uintah::ProblemSpecP& ps);

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
