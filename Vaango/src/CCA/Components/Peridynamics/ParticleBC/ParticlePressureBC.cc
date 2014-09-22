#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Geometry/BBox.h>
#include <Core/Math/Matrix3.h>
#include <iostream>

using namespace Vaango;

// Store the geometry object and the load curve
ParticlePressureBC::ParticlePressureBC(Uintah::ProblemSpecP& ps)
   : ParticleLoadBCBase(ps)
{
  // Read and save the load curve information
  // The sign of the pressure load is +ve if applied in the direction
  // of the outward normal and -ve if applied in the direction of the
  // inward normal
  d_loadCurve = scinew ParticleLoadCurve<double>(ps);
}

// Destroy the pressure BCs
ParticlePressureBC::~ParticlePressureBC()
{
  delete d_loadCurve;
}

void ParticlePressureBC::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP press_ps = ps->appendChild("pressure");
  Uintah::ProblemSpecP geom_ps = press_ps->appendChild("geom_object");
  d_surface->outputProblemSpec(geom_ps);
  press_ps->appendElement("numberOfParticlesOnLoadSurface",d_numParticlesOnLoadSurface);
  d_loadCurve->outputProblemSpec(press_ps);
}

// Get the type of this object for BC application
std::string 
ParticlePressureBC::getType() const
{
  return "Pressure";
}

// Calculate the force per particle at a certain time
double 
ParticlePressureBC::forcePerParticle(double time) const
{
  if (d_numParticlesOnLoadSurface < 1) return 0.0;

  // Get the area of the surface on which the pressure BC is applied
  double area = getSurfaceArea();

  // Get the initial pressure that is applied ( t = 0.0 )
  double press = pressure(time);

  // Calculate the force per particle
  return (press*area)/static_cast<double>(d_numParticlesOnLoadSurface);
}

// Calculate the force vector to be applied to a particular
// particle location
SCIRun::Vector
ParticlePressureBC::getForceVector(const SCIRun::Point& px, double forcePerParticle,
                           const double time) const
{
  SCIRun::Vector force(0.0,0.0,0.0);
  if (d_surfaceType == "box") {
    Uintah::BoxGeometryPiece* gp = dynamic_cast<Uintah::BoxGeometryPiece*>(d_surface);
    SCIRun::Vector normal(0.0, 0.0, 0.0);
    normal[gp->thicknessDirection()] = 1.0;
    force = normal*forcePerParticle;
  } else if (d_surfaceType == "cylinder") {
    Uintah::CylinderGeometryPiece* gp = dynamic_cast<Uintah::CylinderGeometryPiece*>(d_surface);
    SCIRun::Vector normal = gp->radialDirection(px);
    force = normal*forcePerParticle;
  } else if (d_surfaceType == "sphere") {
    Uintah::SphereGeometryPiece* gp = dynamic_cast<Uintah::SphereGeometryPiece*>(d_surface);
    SCIRun::Vector normal = gp->radialDirection(px);
    force = normal*forcePerParticle;
  } else {
    throw Uintah::ParameterNotFound("ERROR: Unknown surface specified for pressure BC",
                            __FILE__, __LINE__);
  }
  return force;
}

namespace Vaango {

// A method to print out the pressure bcs
std::ostream& operator<<(std::ostream& out, const ParticlePressureBC& bc) 
{
   out << "Begin Particle Pressure BC # = " << bc.loadCurveID() << std::endl;
   std::string surfType = bc.getSurfaceType();
   out << "    Surface of application = " << surfType << std::endl;
   if (surfType == "box") {
      Uintah::Box box = (bc.getSurface())->getBoundingBox();
      out << "        " << box << std::endl;
   } else if (surfType == "cylinder") {
      Uintah::CylinderGeometryPiece* cgp = 
         dynamic_cast<Uintah::CylinderGeometryPiece*>(bc.getSurface());
      out << "        " << "radius = " << cgp->radius() 
                        << " top = " << cgp->top() 
                        << " bottom = " << cgp->bottom() << std::endl;
   } else if (surfType == "sphere") {
      Uintah::SphereGeometryPiece* sgp = 
         dynamic_cast<Uintah::SphereGeometryPiece*>(bc.getSurface());
      out << "        " << "radius = " << sgp->radius() 
                        << " origin = " << sgp->origin() << std::endl;
   }
   out << "    Time vs. Load = " << std::endl;
   ParticleLoadCurve<double>* lc = bc.getLoadCurve();
   int numPts = lc->numberOfPointsOnLoadCurve();
   for (int ii = 0; ii < numPts; ++ii) {
     out << "        time = " << lc->getTime(ii) 
         << " pressure = " << lc->getLoad(ii) << std::endl;
   }
   out << "End Particle Pressure BC # = " << bc.loadCurveID() << std::endl;
   return out;
}

} // end namespace Vaango
