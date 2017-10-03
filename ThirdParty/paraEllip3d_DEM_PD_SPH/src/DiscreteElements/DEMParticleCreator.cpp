#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Util/Utility.h>
#include <Core/Math/Vec.h>

using Vec = dem::Vec;

using namespace dem;

/**
 * Create a set of DEM particles:
 *    layerFlag = 0 => A single DEM particle
 *    layerFlag = 1 => A single horizontal layer of DEM particles
 *    layerFlag = 2 => Multiple horizontal layers of DEM particles
 */
template <std::size_t layerFlag>
DEMParticlePArray
DEMParticleCreator::generateDEMParticles(const dem::Box& allContainer, 
  dem::Gradation& gradation)
{
  // Get particle parameters
  DEMParticleCreator::ParticleParameters params = getParticleParameters(gradation);

  // Create the particles
  DEMParticlePArray particles = createParticles<layerFlag>(allContainer,
    gradation, params);

  return particles;

} // generateDEMParticle


/**
 * This structure is used to pass parameters used in particle creation.
 * The spacing and offset is different for spheres.  It's just based 
 * on the maximum particle radius for artibitrary ellipsoids.
 */
DEMParticleCreator::ParticleParameters
DEMParticleCreator::getParticleParameters(const dem::Gradation& gradation)
{
  DEMParticleCreator::ParticleParameters params;

  // Get material parameters
  params.youngModulus = util::getParam<REAL>("young");
  params.poissonRatio = util::getParam<REAL>("poisson");

  // Get max particle size
  params.maxDiameter = gradation.getPtclMaxRadius()*2.0;

  // Default edge and offset
  params.edge = params.maxDiameter;
  params.offset = 0;

  // Updated edge and offset  for spherical particles
  if (gradation.getSize().size() == 1 && gradation.getPtclRatioBA() == 1.0 &&
      gradation.getPtclRatioCA() == 1.0) {
    params.edge = params.maxDiameter * 2.0;
    params.offset = params.maxDiameter * 0.25;
  }

  return params;
}

/** 
 * Create a single particle at the center of the container box
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<0>(const dem::Box& allContainer, 
  dem::Gradation& gradation, 
  const DEMParticleCreator::ParticleParameters& params)
{
  DEMParticlePArray particles;
  std::size_t particleNum = 0;

  DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 0,
    allContainer.getCenter(), gradation, params.youngModulus,
    params.poissonRatio);

  particles.push_back(particle);

  return particles;
}

/** 
 * Create a single horizontal layer of particles at the z-location
 * of the center of the container
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<1>(const dem::Box& allContainer, 
  dem::Gradation& gradation, 
  const DEMParticleCreator::ParticleParameters& params)
{
  auto z_cen = allContainer.getCenter().z();

  auto x_min = allContainer.getMinCorner().x() + params.edge;
  auto x_max = allContainer.getMinCorner().x() - params.edge;
  auto xCoords = util::linspaceApprox<REAL>(x_min, x_max, params.maxDiameter);

  auto y_min = allContainer.getMinCorner().y() + params.edge;
  auto y_max = allContainer.getMinCorner().y() - params.edge;
  auto yCoords = util::linspaceApprox<REAL>(y_min, y_max, params.maxDiameter);

  DEMParticlePArray particles;
  std::size_t particleNum = 0;

  for (auto xCoord : xCoords) {
    for (auto yCoord : yCoords) {
      DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 0,
        Vec(xCoord, yCoord, z_cen), gradation, params.youngModulus,
        params.poissonRatio);
      particles.push_back(particle);
    }
  }

  return particles;
}

/** 
 * Create a multiple horizontal layers of particles from a starting 
 * z-position to an ending z-position inside the container.
 * Particles in each layer are offset by a small amount in a zig-zag
 * manner.
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<2>(const dem::Box& allContainer, 
  dem::Gradation& gradation, 
  const DEMParticleCreator::ParticleParameters& params)
{
  // Get the z-location parameters
  REAL z_min = 0;
  try {
    z_min = util::getParam<REAL>("floatMinZ");
  } catch (const std::out_of_range& err) {
    z_min = allContainer.getMinCorner().z();
  }
  REAL z_max = 0;
  try {
    z_max = util::getParam<REAL>("floatMaxZ");
  } catch (const std::out_of_range& err) {
    z_max = allContainer.getMaxCorner().z();
  }
  z_min += params.maxDiameter;
  z_max -= params.maxDiameter;
  auto zCoords = util::linspaceApprox<REAL>(z_min, z_max, params.maxDiameter);

  auto x_min = allContainer.getMinCorner().x() + params.edge;
  auto x_max = allContainer.getMinCorner().x() - params.edge;
  auto xCoords = util::linspaceApprox<REAL>(x_min, x_max, params.maxDiameter);

  auto y_min = allContainer.getMinCorner().y() + params.edge;
  auto y_max = allContainer.getMinCorner().y() - params.edge;
  auto yCoords = util::linspaceApprox<REAL>(y_min, y_max, params.maxDiameter);

  DEMParticlePArray particles;
  std::size_t particleNum = 0;
  auto offset = params.offset;

  for (auto zCoord : zCoords) {
    for (auto xCoord : xCoords) {
      xCoord += offset;
      for (auto yCoord : yCoords) {
        yCoord += offset;
        DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 0,
          Vec(xCoord, yCoord, zCoord), gradation, params.youngModulus,
          params.poissonRatio);
        particles.push_back(particle);
      }
    }
    offset *= -1;
  }

  return particles;
}

namespace dem {
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<0>(
  const dem::Box& allContainer, dem::Gradation& gradation);
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<1>(
  const dem::Box& allContainer, dem::Gradation& gradation);
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<2>(
  const dem::Box& allContainer, dem::Gradation& gradation);
}