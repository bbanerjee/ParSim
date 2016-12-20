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

#ifndef VAANGO_PARTICLE_FORCE_BC_H
#define VAANGO_PARTICLE_FORCE_BC_H

#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCBase.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadCurve.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   ParticleForceBC
    \brief   The input to this type of BC is a force vector and a 
             box that encapsulates the surface particles to which the BC should
             be applied.
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class ParticleForceBC : public ParticleLoadBCBase  {

  public:

    ParticleForceBC(Uintah::ProblemSpecP& ps);
    virtual ~ParticleForceBC();

    virtual std::string getType() const;
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    // Get the load curve
    inline ParticleLoadCurve<Uintah::Vector>* getLoadCurve() const {return d_loadCurve;}

    // Get the load curve number for this pressure BC
    inline int loadCurveID() const {return d_loadCurve->getID();}

    // Get the applied normal force at time t
    inline Uintah::Vector getLoad(double t) const {return d_loadCurve->getLoad(t);}

  private:

    // Load curve information (Force density and time)
    ParticleLoadCurve<Uintah::Vector>* d_loadCurve;

    // Prevent empty constructor
    ParticleForceBC();

    // No copying allowed
    ParticleForceBC(const ParticleForceBC&);
    ParticleForceBC& operator=(const ParticleForceBC&);
      
  };
} // End namespace Vaango

#endif
