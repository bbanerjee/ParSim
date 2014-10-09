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

#ifndef __MATITI_DAMAGE_MODEL_BASE_H__
#define __MATITI_DAMAGE_MODEL_BASE_H__

#include <Pointers/DamageModelSP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Geometry/Vector3D.h>

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
    virtual const Vector3D& damageStretch() const = 0;

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

