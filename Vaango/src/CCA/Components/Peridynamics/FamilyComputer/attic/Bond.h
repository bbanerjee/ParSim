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
