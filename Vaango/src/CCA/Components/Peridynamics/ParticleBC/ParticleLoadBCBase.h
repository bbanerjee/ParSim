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
        
    // Locate and flag the material points to which this pressure BC is
    // to be applied. 
    bool flagSurfaceParticle(const SCIRun::Point& p, const SCIRun::Vector& dxpp);
      
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
