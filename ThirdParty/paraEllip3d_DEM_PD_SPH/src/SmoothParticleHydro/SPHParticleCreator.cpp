#include <SmoothParticleHydro/SPHParticleCreator.h>
#include <Core/Util/Utility.h>
#include <Core/Math/Vec.h>

using Vec = dem::Vec;

namespace sph {

// here the free particles, ghost particles and boundary particles
// will all stored in the same vector
//  if (getMPIRank() != 0) return;  // make only primary cpu generate particles
template <int dim>
SPHParticlePArray
SPHParticleCreator::generateSPHParticleDam(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles)
{
  // Determine length parameter and SPH point spacing
  auto sphLength = util::getParam<REAL>("waterLength");
  auto numSPHPoint = util::getParam<REAL>("nSPHPoint");
  auto spaceInterval = sphLength / (numSPHPoint - 1);
  dem::InputParameter::get().addParameter("spaceInterval", spaceInterval);
  auto small_value = 0.01 * spaceInterval;
  auto smoothLength = 1.5 * spaceInterval;
  auto kernelSize = 3 * smoothLength;
  //std::cout << "Space interval = " << spaceInterval << "\n";

  // sph parameters
  auto gravAccel = util::getParam<REAL>("gravAccel");
  auto gravScale = util::getParam<REAL>("gravScale");
  auto numLayers = util::getParam<int>("numLayers");
  auto gamma = util::getParam<REAL>("gamma");
  auto P0 = util::getParam<REAL>("P0");
  auto sphInitialDensity = util::getParam<REAL>("SPHInitialDensity");

  // Compute particle mass
  auto sphMass = computeMass<dim>(sphInitialDensity, sphLength, numSPHPoint);

  // get the dimensions of the sph domain
  Vec vmin = spatialDomain.minCorner();
  Vec vmax = spatialDomain.maxCorner();

  // Create an linearly spaced arrays of x/y/zcoords
  std::vector<REAL> xCoords, yCoords, zCoords;
  createCoords<dim>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                    zCoords);

  /*
  std::cout << "x:(" << xCoords.size() << ")";
  for (auto x: xCoords) { std::cout << x << " "; }
  std::cout << std::endl;
  std::cout << "y:(" << yCoords.size() << ")";
  for (auto y: yCoords) { std::cout << y << " "; }
  std::cout << std::endl;
  std::cout << "z:(" << zCoords.size() << ")";
  for (auto z: zCoords) { std::cout << z << " "; }
  std::cout << std::endl;
  */

  // Create the initial allSPHParticleArray
  SPHParticlePArray allSPHParticles =
    createParticleArray<dim>(sphMass, sphInitialDensity, xCoords, yCoords,
                           zCoords);
  // std::cout << "Created particle array: " << d_allSPHParticleVec.size() <<
  // "\n";

  // Identify type of particle and change density based on position
  for (auto& sph_particle : allSPHParticles) {
    dem::Vec local_pos = 0;
    dem::Vec sph_pos = sph_particle->getInitPosition();
    bool isInside = false;
    bool insideGhostLayer = false;

    for (auto& dem_particle : allDEMParticles) {

      if (sph_particle->isInsideDEMParticle<dim>(kernelSize, dem_particle,
                                                 local_pos, insideGhostLayer)) {
        isInside = true;
        if (insideGhostLayer) {
          sph_particle->setType(SPHParticleType::GHOST);
          sph_particle->setLocalPosition(local_pos);
          sph_particle->setDEMParticle(dem_particle.get());
        } else {
          // std::cout << sph_particle->getId() << " Inside but not in ghost
          // layer\n";
          sph_particle->setType(SPHParticleType::NONE);
          sph_particle->setLocalPosition(local_pos);
          sph_particle->setDEMParticle(dem_particle.get());
        }
        break;
      }
    } // end dem particles

    // free/boundary particles
    if (!isInside) {
      REAL density = sphInitialDensity;
      REAL bufferLength = spaceInterval - small_value;
      if (sph_particle->isOutsideDomain<dim>(bufferLength, vmin, vmax)) {
        if (sph_pos.z() < sphLength + spaceInterval) {
          REAL rhoGH = density * gravAccel * (sphLength - sph_pos.z());
          density *= pow(1 + rhoGH * gravScale / P0, 1.0 / gamma);
        }
        sph_particle->setType(SPHParticleType::BOUNDARY);
        sph_particle->setDensity(density);
        sph_particle->setLocalPosition(local_pos);
      } else {
        if (sph_pos.x() <= sphLength && sph_pos.z() <= sphLength) {
          REAL rhoGH = density * gravAccel * (sphLength - sph_pos.z());
          density *= pow(1 + rhoGH * gravScale / P0, 1.0 / gamma);
        }
        sph_particle->setType(SPHParticleType::FREE);
        sph_particle->setDensity(density);
        sph_particle->setLocalPosition(local_pos);
      } // End domain if
    }   // End if not isGhost
  }     // end d_allSPHParticleVec

  // Remove the SPH particle that are inside the DEM particles
  // but not in the ghost layer of each particle
  removeRedundantSPHParticles(allSPHParticles);

  // Compute volume, pressure, viscosity of each particle
  for (auto& sph_particle : allSPHParticles) {
    sph_particle->initialize();
  }

  return allSPHParticles;

} // generateSPHParticle2D

// For drainage problem with open bottomed domain
// make only primary cpu generate particles
// if (getMPIRank() != 0) return; 
template <int dim>
SPHParticlePArray
SPHParticleCreator::generateSPHParticleNoBottom(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles)
{
  // Determine length parameter and SPH point spacing
  REAL xMin = spatialDomain.minCorner().x();
  REAL xMax = spatialDomain.maxCorner().x();
  auto sphLength = xMin - xMax;
  auto spaceInterval = util::getParam<REAL>("spaceInterval");
  auto numSPHPoint = sphLength/spaceInterval + 1;

  auto small_value = 0.01 * spaceInterval;
  //auto smoothLength = 0.0;
  auto kernelSize = 0.0;

  // SPH parameters
  auto gravAccel = util::getParam<REAL>("gravAccel");
  auto gravScale = util::getParam<REAL>("gravScale");
  auto numLayers = util::getParam<int>("numLayers");
  auto gamma = util::getParam<REAL>("gamma");
  auto P0 = util::getParam<REAL>("P0");
  auto sphInitialDensity = util::getParam<REAL>("SPHInitialDensity");

  auto sphMass = computeMass<dim>(sphInitialDensity, sphLength, numSPHPoint);

  // get the dimensions of the sph domain
  Vec vmin = spatialDomain.minCorner();
  Vec vmax = spatialDomain.maxCorner();

  // Create an linearly spaced arrays of x/y/zcoords
  std::vector<REAL> xCoords, yCoords, zCoords;
  createCoords<dim>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                  zCoords);

  /*
  std::cout << "x:(" << xCoords.size() << ")";
  for (auto x: xCoords) { std::cout << x << " "; }
  std::cout << std::endl;
  std::cout << "y:(" << yCoords.size() << ")";
  for (auto y: yCoords) { std::cout << y << " "; }
  std::cout << std::endl;
  std::cout << "z:(" << zCoords.size() << ")";
  for (auto z: zCoords) { std::cout << z << " "; }
  std::cout << std::endl;
  */

  // Create the initial allSPHParticleVec
  SPHParticlePArray allSPHParticles =
    createParticleArray<dim>(sphMass, sphInitialDensity, xCoords, yCoords,
                         zCoords);

  //std::cout << "Created particle array: " << d_allSPHParticleVec.size() << "\n";

  // create and store SPHParticle objects into d_sphParticleVec
  for (auto& sph_particle : allSPHParticles) {
    dem::Vec local_pos = 0;
    dem::Vec sph_pos = sph_particle->getInitPosition();
    bool isInside = false;
    bool insideGhostLayer = false;

    for (auto& dem_particle : allDEMParticles) {

      if (sph_particle->isInsideDEMParticle<dim>(kernelSize, dem_particle,
                                               local_pos, insideGhostLayer)) {
        isInside = true;
        sph_particle->setType(SPHParticleType::GHOST);
        sph_particle->setLocalPosition(local_pos);
        sph_particle->setDEMParticle(dem_particle.get());
        break;
      }
    } // end dem particles

    // free/boundary particles
    if (!isInside) {
      REAL density = sphInitialDensity;
      REAL bufferLength = spaceInterval - small_value;
      if (sph_particle->isOutsideDomainWithoutZBottom<dim>(bufferLength, vmin, vmax)) {
        if (sph_pos.z() < vmax.z() + bufferLength) {
          REAL rhoGH = density * gravAccel * (vmax.z() - sph_pos.z());
          density *= pow(1 + rhoGH * gravScale / P0, 1.0 / gamma);
        }
        sph_particle->setType(SPHParticleType::BOUNDARY);
        sph_particle->setDensity(density);
        sph_particle->setLocalPosition(local_pos);
      } else {
        if (sph_pos.z() >= vmin.z()) {
          REAL rhoGH = density * gravAccel * (vmax.z() - sph_pos.z());
          density *= pow(1 + rhoGH * gravScale / P0, 1.0 / gamma);
        }
        sph_particle->setType(SPHParticleType::FREE);
        sph_particle->setDensity(density);
        sph_particle->setLocalPosition(local_pos);
      } // End domain if
    }   // End if not isGhost
  }

  // Remove the SPH particle that are inside the DEM particles
  // but not in the ghost layer of each particle
  removeRedundantSPHParticles(allSPHParticles);

  // Compute volume, pressure, viscosity of each particle
  for (auto& sph_particle : allSPHParticles) {
    sph_particle->initialize();
  }

  return allSPHParticles;

} // generateSPHParticleNoBottom

// For drainage problem 
// (middle layers only)
template <int dim>
SPHParticlePArray
SPHParticleCreator::generateSPHParticleMiddleLayers(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles)
{
  // Determine length parameter and SPH point spacing
  REAL xMin = spatialDomain.minCorner().x();
  REAL xMax = spatialDomain.maxCorner().x();
  auto sphLength = xMin - xMax;
  auto spaceInterval = util::getParam<REAL>("spaceInterval");
  auto numSPHPoint = sphLength/spaceInterval + 1;

  auto small_value = 0.01 * spaceInterval;
  //auto smoothLength = 0.0;
  auto kernelSize = 0.0;

  // The location of the water
  auto zMinWater = util::getParam<REAL>("waterZmin");
  auto zMaxWater = util::getParam<REAL>("waterZmax");

  // SPH parameters
  auto gravAccel = util::getParam<REAL>("gravAccel");
  auto gravScale = util::getParam<REAL>("gravScale");
  auto numLayers = util::getParam<int>("numLayers");
  auto gamma = util::getParam<REAL>("gamma");
  auto P0 = util::getParam<REAL>("P0");
  auto sphInitialDensity = util::getParam<REAL>("SPHInitialDensity");

  auto sphMass = computeMass<dim>(sphInitialDensity, sphLength, numSPHPoint);

  // get the dimensions of the sph domain
  Vec vmin = spatialDomain.minCorner();
  Vec vmax = spatialDomain.maxCorner();

  // Create an linearly spaced arrays of x/y/zcoords
  std::vector<REAL> xCoords, yCoords, zCoords;
  createCoords<dim>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                  zCoords);

  /*
  std::cout << "x:(" << xCoords.size() << ")";
  for (auto x: xCoords) { std::cout << x << " "; }
  std::cout << std::endl;
  std::cout << "y:(" << yCoords.size() << ")";
  for (auto y: yCoords) { std::cout << y << " "; }
  std::cout << std::endl;
  std::cout << "z:(" << zCoords.size() << ")";
  for (auto z: zCoords) { std::cout << z << " "; }
  std::cout << std::endl;
  */

  // Create the initial allSPHParticleVec
  SPHParticlePArray allSPHParticles =
    createParticleArray<dim>(sphMass, sphInitialDensity, xCoords, yCoords,
                         zCoords);

  //std::cout << "Created particle array: " << d_allSPHParticleVec.size() << "\n";

  // create and store SPHParticle objects into d_sphParticleVec
  for (auto& sph_particle : allSPHParticles) {
    dem::Vec local_pos = 0;
    dem::Vec sph_pos = sph_particle->getInitPosition();
    bool isInside = false;
    bool insideGhostLayer = false;

    for (auto& dem_particle : allDEMParticles) {

      if (sph_particle->isInsideDEMParticle<dim>(kernelSize, dem_particle,
                                               local_pos, insideGhostLayer)) {
        isInside = true;
        sph_particle->setType(SPHParticleType::GHOST);
        sph_particle->setLocalPosition(local_pos);
        sph_particle->setDEMParticle(dem_particle.get());
        break;
      }
    } // end dem particles

    // free/boundary particles
    if (!isInside) {
      REAL density = sphInitialDensity;
      REAL bufferLength = spaceInterval - small_value;
      if (sph_particle->isOutsideDomainWithoutZTop<dim>(bufferLength, vmin, vmax)) {
        sph_particle->setType(SPHParticleType::BOUNDARY);
        sph_particle->setDensity(density);
        sph_particle->setLocalPosition(local_pos);
      } else {
        if (sph_pos.z() >= zMinWater && sph_pos.z() <= zMaxWater) {
          REAL rhoGH = density * gravAccel * (vmax.z() - sph_pos.z());
          density *= pow(1 + rhoGH * gravScale / P0, 1.0 / gamma);
          sph_particle->setType(SPHParticleType::FREE);
          sph_particle->setDensity(density);
          sph_particle->setLocalPosition(local_pos);
        } else {
          sph_particle->setType(SPHParticleType::NONE);
          sph_particle->setDensity(density);
          sph_particle->setLocalPosition(local_pos);
        }
      } // End domain if
    }   // End if not isGhost
  }

  // Remove the SPH particle that are inside the DEM particles
  // but not in the ghost layer of each particle
  removeRedundantSPHParticles(allSPHParticles);

  // Compute volume, pressure, viscosity of each particle
  for (auto& sph_particle : allSPHParticles) {
    sph_particle->initialize();
  }

  return allSPHParticles;

} // generateSPHParticleMiddleLayers

template <>
REAL
SPHParticleCreator::computeMass<2>(const double& density, const double& length,
                                    const std::size_t& numPts) const
{
  auto L_over_N = length / numPts;
  return density * L_over_N * L_over_N;
}

template <>
REAL
SPHParticleCreator::computeMass<3>(const double& density, const double& length,
                                    const std::size_t& numPts) const
{
  auto L_over_N = length / numPts;
  return density * L_over_N * L_over_N * L_over_N;
}

template <>
void
SPHParticleCreator::createCoords<3>(const Vec& vmin, const Vec& vmax,
                                     const REAL& spaceInterval,
                                     const int& numLayers,
                                     std::vector<REAL>& xCoords,
                                     std::vector<REAL>& yCoords,
                                     std::vector<REAL>& zCoords) const
{
  REAL xmin = vmin.x();
  REAL ymin = vmin.y();
  REAL zmin = vmin.z();
  REAL xmax = vmax.x();
  REAL ymax = vmax.y();
  REAL zmax = vmax.z();

  // Create the domain buffer length
  REAL bufferLength = spaceInterval * numLayers;

  // Modify the domain size using the buffer length
  REAL xminBuffered = xmin - bufferLength;
  REAL xmaxBuffered = xmax + bufferLength;
  REAL yminBuffered = ymin - bufferLength;
  REAL ymaxBuffered = ymax + bufferLength;
  REAL zminBuffered = zmin - bufferLength;
  REAL zmaxBuffered = zmax + bufferLength;

  // Create an linearly spaced array of x/zcoords from x/zminBuffered to
  // x/zmaxBuffered
  xCoords =
    util::linspaceApprox<REAL>(xminBuffered, xmaxBuffered, spaceInterval);
  yCoords =
    util::linspaceApprox<REAL>(yminBuffered, ymaxBuffered, spaceInterval);
  zCoords =
    util::linspaceApprox<REAL>(zminBuffered, zmaxBuffered, spaceInterval);
}

template <>
void
SPHParticleCreator::createCoords<2>(const Vec& vmin, const Vec& vmax,
                                     const REAL& spaceInterval,
                                     const int& numLayers,
                                     std::vector<REAL>& xCoords,
                                     std::vector<REAL>& yCoords,
                                     std::vector<REAL>& zCoords) const
{
  createCoords<3>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                  zCoords);
  yCoords.clear();
}

template <>
SPHParticlePArray
SPHParticleCreator::createParticleArray<3>(const REAL& mass,
                                            const REAL& density,
                                            const std::vector<REAL>& xCoords,
                                            const std::vector<REAL>& yCoords,
                                            const std::vector<REAL>& zCoords)
{
  SPHParticlePArray allSPHParticles;

  ParticleID partID = 0;
  dem::Vec local_pos = 0;
  for (const auto& zCoord : zCoords) {
    for (const auto& yCoord : yCoords) {
      for (const auto& xCoord : xCoords) {
        allSPHParticles.push_back(std::make_shared<sph::SPHParticle>(
          partID, mass, density, xCoord, yCoord, zCoord, local_pos,
          SPHParticleType::NONE));
        partID++;
      }
    }
  }

  // std::cout << "Created " << d_allSPHParticleVec.size() << " 3D particles."
  // << std::endl;

  // Compiler should be able to interpret this as a move
  return allSPHParticles;
}

template <>
SPHParticlePArray
SPHParticleCreator::createParticleArray<2>(const REAL& mass,
                                            const REAL& density,
                                            const std::vector<REAL>& xCoords,
                                            const std::vector<REAL>& yCoords,
                                            const std::vector<REAL>& zCoords)
{
  SPHParticlePArray allSPHParticles;

  ParticleID partID = 0;
  dem::Vec local_pos = 0;
  for (const auto& zCoord : zCoords) {
    for (const auto& xCoord : xCoords) {
      allSPHParticles.push_back(std::make_shared<sph::SPHParticle>(
        partID, mass, density, xCoord, 0.0, zCoord, local_pos,
        SPHParticleType::NONE));
      partID++;
    }
  }

  // std::cout << "Created " << d_allSPHParticleVec.size() << " 2D particles."
  // << std::endl;

  // Compiler should be able to interpret this as a move
  return allSPHParticles;
}

void
SPHParticleCreator::removeRedundantSPHParticles(SPHParticlePArray& allSPHParticles)
{
  allSPHParticles.erase(
    std::remove_if(allSPHParticles.begin(), allSPHParticles.end(),
                   [](SPHParticleP particle) {
                     if (particle->getType() == SPHParticleType::NONE) {
                       return true;
                     }
                     return false;
                   }),
    allSPHParticles.end());
}


template SPHParticlePArray SPHParticleCreator::generateSPHParticleDam<2>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);
template SPHParticlePArray SPHParticleCreator::generateSPHParticleDam<3>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);
template SPHParticlePArray SPHParticleCreator::generateSPHParticleNoBottom<2>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);
template SPHParticlePArray SPHParticleCreator::generateSPHParticleNoBottom<3>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);
template SPHParticlePArray SPHParticleCreator::generateSPHParticleMiddleLayers<2>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);
template SPHParticlePArray SPHParticleCreator::generateSPHParticleMiddleLayers<3>(
  const dem::Box& spatialDomain, const dem::DEMParticlePArray& allDEMParticles);

} // end namespace sph