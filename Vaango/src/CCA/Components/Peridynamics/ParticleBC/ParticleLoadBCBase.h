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

#ifndef VAANGO_PARTICLE_LOAD_BC_BASE_H
#define VAANGO_PARTICLE_LOAD_BC_BASE_H

#include <string>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>

namespace Uintah {
  class GeometryPiece;
}

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   ParticleLoadBCBase
    \brief   Abstract base class for various types of load boundary conditions
             that are applied to particle sets
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class ParticleLoadBCBase  {
  public:
    ParticleLoadBCBase() {};
    ParticleLoadBCBase(Uintah::ProblemSpecP& ps);
    virtual ~ParticleLoadBCBase();

    virtual std::string getType() const = 0;
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;
    virtual int loadCurveID() const = 0;

    // Locate and flag the material points to which this pressure BC is
    // to be applied. 
    bool flagSurfaceParticle(const Uintah::Point& p, const Uintah::Vector& dxpp);
      
    // Get the surface 
    inline Uintah::GeometryPiece* getSurface() const {return d_surface;}

    // Get the surface type
    inline std::string getSurfaceType() const {return d_surfaceType;}

    // Get the area of the surface
    double getSurfaceArea() const;

    // Set the number of material points on the surface
    inline void numParticlesOnLoadSurface(long num) {d_numParticlesOnLoadSurface = num;}

    // Get the number of material points on the surface
    inline long numParticlesOnLoadSurface() const {return d_numParticlesOnLoadSurface;}

  protected:

    // Private Data for this and derived classes
    // Surface information
    Uintah::GeometryPiece* d_surface;
    std::string d_surfaceType;
    long d_numParticlesOnLoadSurface;

  private:

    // Do not allow copying
    ParticleLoadBCBase(const ParticleLoadBCBase&);
    ParticleLoadBCBase& operator=(const ParticleLoadBCBase&);
  };

} // End namespace Uintah
   
#endif
