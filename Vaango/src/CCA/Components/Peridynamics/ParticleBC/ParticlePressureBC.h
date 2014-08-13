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
    SCIRun::Vector getForceVector(const SCIRun::Point& px, double forcePerParticle,
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
