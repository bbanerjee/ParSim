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

#ifndef VAANGO_PARTICLE_PRESSURE_BC_H
#define VAANGO_PARTICLE_PRESSURE_BC_H

#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCBase.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadCurve.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Math/Matrix3.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>

#include <iosfwd>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   ParticlePressureBC
    \brief   Stores the pressure load curves and boundary imformation for
             pressure boundary conditions that can be applied to surfaces
             with simple geometry -  planes, cylinders and spheres.
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class GeometryPiece;
  class ParticleCreator;
   
  class ParticlePressureBC : public ParticleLoadBCBase  {

  public:

    // Construct a ParticlePressureBC object that contains
    // the area over which pressure is to be applied
    // and the value of that pressure (in the form
    // of a load curve)
    ParticlePressureBC(Uintah::ProblemSpecP& ps);
    virtual ~ParticlePressureBC();

    virtual std::string getType() const;

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    // Get the load curve 
    inline ParticleLoadCurve<double>* getLoadCurve() const {return d_loadCurve;}

    // Get the load curve number for this pressure BC
    inline int loadCurveID() const {return d_loadCurve->getID();}

    // Get the applied pressure at time t
    inline double pressure(double t) const {return d_loadCurve->getLoad(t);}

    // Get the force per particle at time t
    double forcePerParticle(double time) const;

    // Get the force vector to be applied at a point 
    Uintah::Vector getForceVector(const Uintah::Point& px, double forcePerParticle,
                                  const double time) const;

  private:

    // Load curve information (Pressure and time)
    ParticleLoadCurve<double>* d_loadCurve;

    // Prevent empty constructor
    ParticlePressureBC();

    // Prevent copying
    ParticlePressureBC(const ParticlePressureBC&);
    ParticlePressureBC& operator=(const ParticlePressureBC&);
      
  public:

    friend std::ostream& operator<<(std::ostream& out, const Vaango::ParticlePressureBC& bc);
  };
} // End namespace Vaango

#endif
