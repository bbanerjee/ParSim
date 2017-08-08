// Function Definitions
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>
#include <SmoothParticleHydro/SPHParticle.h>

using namespace sph;

SPHParticle::SPHParticle()
  : d_id(0)
  , d_mass(0)
  , d_density(0)
  , d_volume(0)
  , d_pressure(0)
  , d_mu(0)
  , d_densityDot(0)
  , d_initialPos(0)
  , d_currPos(0)
  , d_velocity(0)
  , d_acceleration(0)
  , d_velocityCorrection(0)
  , d_localCoords(0)
  , d_type(SPHParticleType::NONE)
{
  d_demParticle = nullptr;
}

SPHParticle::SPHParticle(ParticleID id, REAL mass, REAL rho, REAL x, REAL y,
                         REAL z, const dem::Vec& local, SPHParticleType type)
{
  // Check whether the correct constructor has been called
  if (type != SPHParticleType::FREE && type != SPHParticleType::GHOST &&
      type != SPHParticleType::BOUNDARY && type != SPHParticleType::NONE) {
    std::cout << "Type of current SPH particle is " << static_cast<int>(type)
              << "\n"
              << "Type should be one of: 1: Free particle, "
              << " 2: Ghost particle, 3: Boundary particle, 4: None."
              << " Error in creating SPH free/ghost/boundary particle!"
              << std::endl;
    exit(-1);
  }

  d_id = id;
  d_mass = mass;
  d_density = rho;
  d_initialPos = dem::Vec(x, y, z);

  initialize();

  if (type == SPHParticleType::GHOST) {
    d_localCoords = local; // Ghost SPH point
  } else {
    d_localCoords = 0; // free or boundary SPH point
  }

  d_type = type;
  d_demParticle = nullptr;

} // end SPHParticle()

void
SPHParticle::initialize()
{
  d_volume = calculateVolume();
  d_pressure = calculatePressure();
  d_mu = calculateViscosity();
  d_densityDot = 0;

  d_currPos = d_initialPos;
  d_velocity = 0;
  d_acceleration = 0;
  d_velocityCorrection = 0;
} // end initial()

REAL
SPHParticle::calculatePressure()
{
  // Material constants
  REAL initialPressure = util::getParam<REAL>("P0");
  REAL initialDensity = util::getParam<REAL>("SPHInitialDensity");
  REAL gamma = util::getParam<REAL>("gamma");
  return initialPressure * (std::pow(d_density / initialDensity, gamma) - 1);
}

REAL
SPHParticle::calculateViscosity()
{
  // Material constants
  REAL kinematicViscosity = util::getParam<REAL>("nu");
  return d_density * kinematicViscosity;
}

// return the trial position, used for the boundary handling model
// based on momentum conservation
dem::Vec
SPHParticle::getTrialPosition() const
{
  REAL timeStep = util::getParam<REAL>("timeStep");
  return d_currPos + d_velocity * timeStep;
}

// here the forward Euler time integration is used
void
SPHParticle::update()
{
  REAL delT = util::getParam<REAL>("timeStep");
  REAL dampingCoeff = util::getParam<REAL>("sphDamping");
  updateDensity(delT);
  updateVelocity(delT);
  d_velocity -= d_velocity * dampingCoeff;
  updatePosition(delT);
} // end update()

// here the forward Euler time integration is used
void
SPHParticle::updateDensity(const REAL& delT)
{
  d_density += d_densityDot * delT;
}

void
SPHParticle::updatePosition(const REAL& delT)
{
  d_currPos += (d_velocity + d_velocityCorrection) * delT;
}

// update position and d_density based on equation (4.1)
void
SPHParticle::updatePositionDensityLeapFrog(const REAL& delT)
{
  updatePosition(delT);
  updateDensity(delT);
}

// update d_velocity based on equation (4.2)
void
SPHParticle::updateVelocity(const REAL& delT)
{
  d_velocity += d_acceleration * delT;
}

// update d_velocity based on equation (4.3)
void
SPHParticle::initialVelocityLeapFrog(const REAL& delT)
{
  updateVelocity(delT * 0.5);
}

template <>
bool
SPHParticle::isInsideDEMParticle<2>(const REAL& kernelSize,
                                    const dem::DEMParticleP& dem_particle,
                                    dem::Vec& localCoord,
                                    bool& insideGhostLayer)
{
  dem::Vec sph_pos = d_initialPos;
  sph_pos.setY(0);
  dem::Vec dem_pos = dem_particle->currentPosition();
  dem_pos.setY(0);
  return dem_particle->containsPoint(sph_pos, dem_pos, kernelSize, localCoord,
                                     insideGhostLayer);
}

template <>
bool
SPHParticle::isInsideDEMParticle<3>(const REAL& kernelSize,
                                    const dem::DEMParticleP& dem_particle,
                                    dem::Vec& localCoord,
                                    bool& insideGhostLayer)
{
  dem::Vec sph_pos = d_initialPos;
  dem::Vec dem_pos = dem_particle->currentPosition();
  return dem_particle->containsPoint(sph_pos, dem_pos, kernelSize, localCoord,
                                     insideGhostLayer);
}

template <>
bool
SPHParticle::isOutsideDomain<2>(const REAL& bufferLength,
                                const dem::Vec& minCorner,
                                const dem::Vec& maxCorner)
{
  REAL xminBoundary = minCorner.x() - bufferLength;
  REAL xmaxBoundary = maxCorner.x() + bufferLength;
  REAL zminBoundary = minCorner.z() - bufferLength;
  REAL zmaxBoundary = maxCorner.z() + bufferLength;
  dem::Box domain(xminBoundary, 0, zminBoundary, xmaxBoundary, 0, zmaxBoundary);

  if (domain.outside(d_initialPos))
    return true;
  return false;
}

template <>
bool
SPHParticle::isOutsideDomain<3>(const REAL& bufferLength,
                                const dem::Vec& minCorner,
                                const dem::Vec& maxCorner)
{
  REAL xminBoundary = minCorner.x() - bufferLength;
  REAL xmaxBoundary = maxCorner.x() + bufferLength;
  REAL yminBoundary = minCorner.y() - bufferLength;
  REAL ymaxBoundary = maxCorner.y() + bufferLength;
  REAL zminBoundary = minCorner.z() - bufferLength;
  REAL zmaxBoundary = maxCorner.z() + bufferLength;
  dem::Box domain(xminBoundary, yminBoundary, zminBoundary, xmaxBoundary,
                  ymaxBoundary, zmaxBoundary);

  if (domain.outside(d_initialPos))
    return true;
  return false;
}

template <>
bool
SPHParticle::isOutsideDomainWithoutZBottom<2>(const REAL& bufferLength,
                                              const dem::Vec& minCorner,
                                              const dem::Vec& maxCorner)
{
  REAL xminBoundary = minCorner.x() - bufferLength;
  REAL xmaxBoundary = maxCorner.x() + bufferLength;
  REAL zminBoundary = -1.0e20;
  REAL zmaxBoundary = maxCorner.z() + bufferLength;
  dem::Box domain(xminBoundary, 0, zminBoundary, xmaxBoundary, 0, zmaxBoundary);

  if (domain.outside(d_initialPos))
    return true;
  return false;
}

template <>
bool
SPHParticle::isOutsideDomainWithoutZBottom<3>(const REAL& bufferLength,
                                const dem::Vec& minCorner,
                                const dem::Vec& maxCorner)
{
  REAL xminBoundary = minCorner.x() - bufferLength;
  REAL xmaxBoundary = maxCorner.x() + bufferLength;
  REAL yminBoundary = minCorner.y() - bufferLength;
  REAL ymaxBoundary = maxCorner.y() + bufferLength;
  REAL zminBoundary = -1.0e20;
  REAL zmaxBoundary = maxCorner.z() + bufferLength;
  dem::Box domain(xminBoundary, yminBoundary, zminBoundary, xmaxBoundary,
                  ymaxBoundary, zmaxBoundary);

  if (domain.outside(d_initialPos))
    return true;
  return false;
}
