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
    inline ParticleLoadCurve<SCIRun::Vector>* getLoadCurve() const {return d_loadCurve;}

    // Get the load curve number for this pressure BC
    inline int loadCurveID() const {return d_loadCurve->getID();}

    // Get the applied normal force at time t
    inline SCIRun::Vector getLoad(double t) const {return d_loadCurve->getLoad(t);}

  private:

    // Load curve information (Force density and time)
    ParticleLoadCurve<SCIRun::Vector>* d_loadCurve;

    // Prevent empty constructor
    ParticleForceBC();

    // No copying allowed
    ParticleForceBC(const ParticleForceBC&);
    ParticleForceBC& operator=(const ParticleForceBC&);
      
  };
} // End namespace Vaango

#endif
