#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Geometry/Ellipsoid.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Util/Utility.h>
#include <Core/Math/Vec.h>

using namespace dem;

using DEMParticleShape = DEMParticle::DEMParticleShape;
using ParticleParameters = DEMParticleCreator::ParticleParameters;


DEMParticlePArray
DEMParticleCreator::generateDEMParticles(std::size_t layerFlag,
                                         DEMParticleShape particleShape,
                                         const Box& spatialDomain, 
                                         Gradation& gradation)
{
  if (layerFlag == 0) {
    return generateDEMParticles<0>(particleShape, spatialDomain, gradation);
  } else if (layerFlag == 1) {
    return generateDEMParticles<1>(particleShape, spatialDomain, gradation);
  } 
  return generateDEMParticles<2>(particleShape, spatialDomain, gradation);
}

/**
 * Create a set of DEM particles:
 *    layerFlag = 0 => A single DEM particle
 *    layerFlag = 1 => A single horizontal layer of DEM particles
 *    layerFlag = 2 => Multiple horizontal layers of DEM particles
 */
template <std::size_t layerFlag>
DEMParticlePArray
DEMParticleCreator::generateDEMParticles(DEMParticleShape particleShape,
                                         const Box& spatialDomain, 
                                         Gradation& gradation)
{
  // Get particle parameters
  ParticleParameters params = getParticleParameters(gradation);

  // Create the particles
  DEMParticlePArray particles = createParticles<layerFlag>(particleShape,
    spatialDomain, gradation, params);

  return particles;

} // generateDEMParticle


/**
 * This structure is used to pass parameters used in particle creation.
 * The spacing and offset is different for spheres.  It's just based 
 * on the maximum particle radius for artibitrary ellipsoids.
 */
ParticleParameters
DEMParticleCreator::getParticleParameters(const Gradation& gradation)
{
  ParticleParameters params;

  // Get number of layers
  params.numLayers = static_cast<int>(util::getParam<REAL>("particleLayers"));

  // Get material parameters
  params.youngModulus = util::getParam<REAL>("young");
  params.poissonRatio = util::getParam<REAL>("poisson");

  // Get max particle size
  params.maxDiameter = gradation.getPtclMaxRadius()*2.0;

  // Get randomization flags
  params.randomOrientation = 
    static_cast<bool>(util::getParam<REAL>("randomOrientation"));
  params.randomRadiusRatio = 
    static_cast<bool>(util::getParam<REAL>("randomRadiusRatio"));
  params.randomVelocity = 
    static_cast<bool>(util::getParam<REAL>("randomVelocity"));

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
 * Create a single particle at the center of the domain box
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<0>(DEMParticleShape particleShape,
                                       const Box& spatialDomain, 
                                       Gradation& gradation, 
                                       const ParticleParameters& params)
{
  DEMParticlePArray particles;
  std::size_t particleNum = 0;

  DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 
    particleShape, DEMParticle::DEMParticleType::FREE,
    spatialDomain.center(), gradation, 
    params.youngModulus, params.poissonRatio,
    params.randomOrientation, params.randomRadiusRatio, params.randomVelocity);

  particles.push_back(particle);

  return particles;
}

/** 
 * Create a single horizontal layer of particles at the z-location
 * of the center of the domain
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<1>(DEMParticleShape particleShape,
                                       const Box& spatialDomain, 
                                       Gradation& gradation, 
                                       const ParticleParameters& params)
{
  auto z_cen = spatialDomain.center().z();

  auto x_min = spatialDomain.minCorner().x() + params.edge;
  auto x_max = spatialDomain.maxCorner().x() - params.edge;
  auto xCoords = util::linspaceApprox<REAL>(x_min, x_max, params.maxDiameter);

  auto y_min = spatialDomain.minCorner().y() + params.edge;
  auto y_max = spatialDomain.maxCorner().y() - params.edge;
  auto yCoords = util::linspaceApprox<REAL>(y_min, y_max, params.maxDiameter);

  std::cout << params << "\n";

  DEMParticlePArray particles;
  std::size_t particleNum = 0;

  std::cout << "Num (x, y) = (" << xCoords.size() << "," << yCoords.size() << ")";
  for (auto x : xCoords) {
    std::cout << x << " ";
  }
  std::cout << std::endl;

  for (auto xCoord : xCoords) {
    for (auto yCoord : yCoords) {
      DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 
        particleShape, DEMParticle::DEMParticleType::FREE,
        Vec(xCoord, yCoord, z_cen), gradation, params.youngModulus,
        params.poissonRatio, params.randomOrientation, params.randomRadiusRatio, 
        params.randomVelocity);
      particles.push_back(particle);
    }
  }

  return particles;
}

/** 
 * Create multiple horizontal layers of particles from a starting 
 * z-position to an ending z-position inside the domain.
 * Particles in each layer are offset by a small amount in a zig-zag
 * manner.
 */
template <>
DEMParticlePArray
DEMParticleCreator::createParticles<2>(DEMParticleShape particleShape,
                                       const Box& spatialDomain, 
                                       Gradation& gradation, 
                                       const ParticleParameters& params)
{
  // Get the z-location parameters
  REAL z_min = 0;
  try {
    z_min = util::getParam<REAL>("floatMinZ");
  } catch (const std::out_of_range& err) {
    z_min = spatialDomain.minCorner().z();
  }

  REAL z_max = z_min + params.numLayers * params.maxDiameter;
  auto zCoords = util::linspace<REAL>(z_min, z_max, params.numLayers-1);

  /*
  REAL z_max = 0;
  try {
    z_max = util::getParam<REAL>("floatMaxZ");
  } catch (const std::out_of_range& err) {
    z_max = spatialDomain.maxCorner().z();
  }
  z_min += params.maxDiameter;
  z_max -= params.maxDiameter;
  auto zCoords = util::linspaceApprox<REAL>(z_min, z_max, params.maxDiameter);
  */

  auto x_min = spatialDomain.minCorner().x() + params.edge;
  auto x_max = spatialDomain.maxCorner().x() - params.edge;
  auto xCoords = util::linspaceApprox<REAL>(x_min, x_max, params.maxDiameter);

  /*
  std::cout << "maxDia = " << params.maxDiameter
            << " edge = " << params.edge
            << " xmin = " << x_min << " xmax = " << x_max
            << "x: ";
  for (auto x : xCoords) {
    std::cout << x << " ";
  }
  std::cout << std::endl;
  */

  auto y_min = spatialDomain.minCorner().y() + params.edge;
  auto y_max = spatialDomain.maxCorner().y() - params.edge;
  auto yCoords = util::linspaceApprox<REAL>(y_min, y_max, params.maxDiameter);

  DEMParticlePArray particles;
  std::size_t particleNum = 0;
  auto offset = params.offset;

  for (auto zCoord : zCoords) {
    for (auto xCoord : xCoords) {
      xCoord += offset;
      for (auto yCoord : yCoords) {
        yCoord += offset;
        DEMParticleP particle = std::make_shared<DEMParticle>(++particleNum, 
          particleShape, DEMParticle::DEMParticleType::FREE,
          Vec(xCoord, yCoord, zCoord), gradation, params.youngModulus,
          params.poissonRatio, params.randomOrientation, params.randomRadiusRatio,
          params.randomVelocity);
        particles.push_back(particle);
      }
    }
    offset *= -1;
  }

  return particles;
}

/** 
 * Non-obvious side effect: Converts particle type to BOUNDARY_PERIODIC if 
 * it is on the boundary 
 */
DEMParticlePArray
DEMParticleCreator::generatePeriodicDEMParticles(DEMParticlePArray& particles,
                                                 Box& spatialDomain,
                                                 REAL marginFactor, 
                                                 REAL faceShiftFactor) 
{
  // Find the largest radius ellipsoid in the input set and set
  // the margin to be marginFactor times  that value
  auto minmaxIter = std::minmax_element(particles.begin(), particles.end(),
   [](const DEMParticleP& p1, const DEMParticleP& p2){
     auto max_p1 = std::max({p1->radiusA(), p1->radiusB(), p1->radiusC()});
     auto max_p2 = std::max({p2->radiusA(), p2->radiusB(), p2->radiusC()});
     return max_p1 < max_p2;
   });
  auto minRadius = std::min({(*minmaxIter.first)->radiusA(),
                             (*minmaxIter.first)->radiusB(),
                             (*minmaxIter.first)->radiusC()});
  auto maxRadius = std::max({(*minmaxIter.second)->radiusA(),
                             (*minmaxIter.second)->radiusB(),
                             (*minmaxIter.second)->radiusC()});
  auto boundaryMargin = marginFactor*maxRadius;
  auto faceShift = faceShiftFactor*minRadius;

  // std::cout << "min rad = " << minRadius << " max rad = " << maxRadius << "\n";
  // std::cout << "Extra boundary margin = " << boundaryMargin 
  //           << " face shift = " << faceShift << "\n";
  
  // Create an oriented box from the spatial domain
  // and set up the faces
  Box shrunkDomain(
    spatialDomain.minCorner() + Vec(faceShift, faceShift, faceShift),
    spatialDomain.maxCorner() - Vec(faceShift, faceShift, faceShift));
  OrientedBox box(shrunkDomain);
  std::vector<Vec> vertices = box.vertices();
  constexpr std::array<std::array<int, 4>, 3> faceIndices = {{
      {{0, 4, 7, 3}} // x-
    , {{0, 1, 5, 4}} // y-
    , {{0, 3, 2, 1}} // z-
  }};
  std::vector<Face> faces;
  for (const auto& indices : faceIndices) {
    int i0 = indices[0]; int i1 = indices[1]; 
    int i2 = indices[2]; int i3 = indices[3];
    Vec  v0 = vertices[i0]; Vec  v1 = vertices[i1];
    Vec  v2 = vertices[i2]; Vec  v3 = vertices[i3];
    Face face(v0, v1, v2, v3);
    if (!face.isValid()) {
      std::cout << "**ERROR** The face " << face << " is invalid but"
                << " continuing anyway\n.";
    }
    faces.push_back(face);
  }

  // Compute starting ID for extra particles
  auto maxIter = std::max_element(particles.begin(), particles.end(),
   [](const DEMParticleP& p1, const DEMParticleP& p2){
     return p1->getId() < p2->getId();
   });
  auto particleID = (*maxIter)->getId();
  //auto particleID = particles.size();

  // Check intersections
  REAL widthX = shrunkDomain.dimX();
  REAL widthY = shrunkDomain.dimY();
  REAL widthZ = shrunkDomain.dimZ();
  DEMParticlePArray extraParticles;
  for (const auto& particle : particles) {

    auto position = particle->currentPosition();

    auto axis_a = particle->currentAxisA();
    auto axis_b = particle->currentAxisB();
    auto axis_c = particle->currentAxisC();

    auto radius_a = particle->radiusA();
    auto radius_b = particle->radiusB();
    auto radius_c = particle->radiusC();

    auto id = particle->getId();

    // Create an ellipsoid object for ellipsoid-face intersection tests
    // *TODO* Generalize to sphere and other particle shapes.
    Ellipsoid ellipsoid(id, position, axis_a, axis_b, axis_c, 
                        radius_a, radius_b, radius_c);

    int faceID = 1;
    for (const auto& face : faces) {
      auto status = ellipsoid.intersects(face);
      // std::cout << "Face = " << face << "\n";
      // std::cout << "status = " << std::boolalpha << status.first
      //           << " face " << static_cast<int>(status.second.first)
      //           << " , " << status.second.second << "\n";
      if (status.first) {
        std::vector<Vec> translations;

        switch (static_cast<Boundary::BoundaryID>(faceID)) {
          case Boundary::BoundaryID::NONE:
            break;
          case Boundary::BoundaryID::XMINUS:
          {
            Vec shift(widthX + boundaryMargin, 0, 0);
            //std::cout << " Loc: x-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, boundaryMargin, Vec(0, 1, 1),
                                 widthX, widthY, widthZ, 
                                 face, status, translations);
            break;
          }
          case Boundary::BoundaryID::XPLUS:
            break;
          case Boundary::BoundaryID::YMINUS:
          {
            Vec shift(0, widthY + boundaryMargin, 0);
            //std::cout << " Loc: y-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, boundaryMargin, Vec(1, 0, 1),
                                 widthX, widthY, widthZ, 
                                 face, status, translations);
            break;
          }
          case Boundary::BoundaryID::YPLUS:
            break;
          case Boundary::BoundaryID::ZMINUS:
          {
            Vec shift(0, 0, widthZ + boundaryMargin);
            //std::cout << " Loc: z-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, boundaryMargin, Vec(1, 1, 0),
                                 widthX, widthY, widthZ, 
                                 face, status, translations);
            break;
          }
          case Boundary::BoundaryID::ZPLUS:
            break;
        }

        // Create copies
        particle->setType(DEMParticle::DEMParticleType::BOUNDARY_PERIODIC);
        for (const auto& translation : translations) {
          DEMParticleP newParticle = std::make_shared<DEMParticle>(*particle);
          newParticle->setCurrentPosition(particle->currentPosition() + 
                                          translation);
          newParticle->setPreviousPosition(particle->previousPosition() + 
                                           translation);
          newParticle->setId(++particleID);
          newParticle->computeAndSetGlobalCoef();
          extraParticles.push_back(newParticle);
        }
      }
      faceID += 2;
    }
  }

  // Remove duplicates
  removeDuplicates(extraParticles);

  // Update the spatial domain
  Box expandedDomain(
    shrunkDomain.minCorner(),
    shrunkDomain.maxCorner() + Vec(boundaryMargin, boundaryMargin, boundaryMargin));
  spatialDomain = expandedDomain;

  return extraParticles;
}

void
DEMParticleCreator::addExtraTranslations(const Vec& shift, 
                                         REAL boundaryMargin,
                                         Vec inPlaneDiag,
                                         REAL widthX, REAL widthY, REAL widthZ,
                                         const Face& face,
                                         const IntersectionStatus status, 
                                         std::vector<Vec>& translations) 
{
  if (status.second.first == Face::Location::VERTEX) {
    int vertIndex = status.second.second;
    //std::cout << "is vertex " << vertIndex << "\n";
    if (vertIndex == 0) {
      for (int ii = 1; ii < 4; ++ii) {
        Vec inPlane = face.vertex[ii] - face.vertex[vertIndex];
        auto length = inPlane.length();
        inPlane.normalizeInPlace();
        if (ii == 2) {
          inPlane *= length;
          inPlane += (inPlaneDiag*boundaryMargin);
        } else {
          inPlane *= (length + boundaryMargin);
        }
        Vec outOfPlane = inPlane + shift;
        translations.push_back(inPlane);
        translations.push_back(outOfPlane);
      }
    } else if (vertIndex == 1 || vertIndex == 3) {
      Vec normal = Vec(1,1,1) - inPlaneDiag;
      Vec vec1 = normal * Vec(widthX + boundaryMargin, 
        widthY + boundaryMargin, widthZ + boundaryMargin);
      Vec vec2 = face.vertex[2] - face.vertex[vertIndex];
      auto length = vec2.length();
      vec2.normalizeInPlace();
      vec2 *= (length + boundaryMargin);
      //std::cout << "vec1 = " << vec1 << " vec2 = " << vec2 << "\n";
      translations.push_back(vec1 + vec2);
    }
  } 
  else if (status.second.first == Face::Location::EDGE) {
    int edgeIndex = status.second.second;
    //std::cout << "is edge " << edgeIndex << "\n";
    if (edgeIndex == 0 || edgeIndex == 3) {
      int oppIndex = (edgeIndex+3) % 4;
      Vec inPlane = face.vertex[oppIndex] - face.vertex[edgeIndex];
      auto length = inPlane.length() + boundaryMargin;
      inPlane.normalizeInPlace();
      inPlane *= length;
      Vec outOfPlane = inPlane + shift;
      translations.push_back(inPlane);
      translations.push_back(outOfPlane);
    }
  }
}

void
DEMParticleCreator::removeDuplicates(DEMParticlePArray& input)
{
  std::vector<Vec> seen;
  auto newEnd = std::remove_if(input.begin(), input.end(),
    [&seen](const DEMParticleP& particle)
    {
      Vec pos = particle->currentPosition();
      if (std::find_if(seen.begin(), seen.end(), 
                       [&pos](const Vec& seen_pos) {return pos == seen_pos;}) 
           != seen.end()) {
        //std::cout << "Found : " << pos << "\n";
        return true;
      }
      //std::cout << "Inserted : " << pos << "\n";
      seen.push_back(pos);
      return false;
    });
  input.erase(newEnd, input.end());
}

DEMParticlePArray
DEMParticleCreator::updatePeriodicDEMParticles(const OrientedBox& periodicDomain, 
                                               DEMParticlePArray& particles)
{
  auto e0 = periodicDomain.axis(0) * (periodicDomain.extent(0) * 2);
  auto e1 = periodicDomain.axis(1) * (periodicDomain.extent(1) * 2);
  auto e2 = periodicDomain.axis(2) * (periodicDomain.extent(2) * 2);

  std::vector<Vec> vertices = periodicDomain.vertices();
  constexpr std::array<std::array<int, 4>, 6> faceIndices = {{
      {{0, 4, 7, 3}} // x-
    , {{1, 2, 6, 5}} // x+
    , {{0, 1, 5, 4}} // y-
    , {{2, 3, 7, 6}} // y+
    , {{0, 3, 2, 1}} // z-
    , {{4, 5, 6, 7}}  // z+

  }};
  std::vector<Face> faces;
  for (const auto& indices : faceIndices) {
    int i0 = indices[0]; int i1 = indices[1]; 
    int i2 = indices[2]; int i3 = indices[3];
    Face face(vertices[i0], vertices[i1], vertices[i2], vertices[i3]);
    faces.push_back(face);
  }

  // Check intersections
  DEMParticlePArray extraParticles;
  for (const auto& particle : particles) {

    if (particle->getType() != DEMParticle::DEMParticleType::FREE) {
      continue;
    }

    auto position = particle->currentPosition();

    auto axis_a = particle->currentAxisA();
    auto axis_b = particle->currentAxisB();
    auto axis_c = particle->currentAxisC();

    auto radius_a = particle->radiusA();
    auto radius_b = particle->radiusB();
    auto radius_c = particle->radiusC();

    auto id = particle->getId();

    // Create an ellipsoid object for ellipsoid-face intersection tests
    Ellipsoid ellipsoid(id, position, axis_a, axis_b, axis_c, 
                        radius_a, radius_b, radius_c);

    int faceID = 1;
    for (const auto& face : faces) {
      auto status = ellipsoid.intersects(face);
      // std::cout << "Face = " << face << "\n";
      // std::cout << "status = " << std::boolalpha << status.first
      //           << " face " << static_cast<int>(status.second.first)
      //           << " , " << status.second.second << "\n";
      if (status.first) {
        std::vector<Vec> translations;

        switch (static_cast<Boundary::BoundaryID>(faceID)) {
          case Boundary::BoundaryID::NONE:
            break;
          case Boundary::BoundaryID::XMINUS:
          {
            Vec shift = e0;
            //std::cout << " Loc: x-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
          case Boundary::BoundaryID::XPLUS:
          {
            Vec shift = -e0;
            //std::cout << " Loc: x-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
          case Boundary::BoundaryID::YMINUS:
          {
            Vec shift = e1;
            //std::cout << " Loc: y-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
          case Boundary::BoundaryID::YPLUS:
          {
            Vec shift = -e1;
            //std::cout << " Loc: y-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
          case Boundary::BoundaryID::ZMINUS:
          {
            Vec shift = e2;
            //std::cout << " Loc: z-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
          case Boundary::BoundaryID::ZPLUS:
          {
            Vec shift = -e2;
            //std::cout << " Loc: z-: " << shift << "\n";
            translations.push_back(shift);
            addExtraTranslations(shift, face, status, translations);
            break;
          }
            break;
        }

        // Create copies
        particle->setType(DEMParticle::DEMParticleType::BOUNDARY_PERIODIC);
        for (const auto& translation : translations) {
          DEMParticleP newParticle = std::make_shared<DEMParticle>(*particle);
          newParticle->setCurrentPosition(particle->currentPosition() + 
                                          translation);
          newParticle->setPreviousPosition(particle->previousPosition() + 
                                           translation);
          newParticle->setId(particle->getId());
          newParticle->computeAndSetGlobalCoef();
          extraParticles.push_back(newParticle);
        }
      }
      faceID += 1;
    }
  }

  // Remove duplicates
  removeDuplicates(extraParticles);

  return extraParticles;
}

void
DEMParticleCreator::addExtraTranslations(const Vec& shift, 
                                         const Face& face,
                                         const IntersectionStatus status, 
                                         std::vector<Vec>& translations) 
{
  if (status.second.first == Face::Location::VERTEX) {
    int vertIndex = status.second.second;
    //std::cout << "is vertex " << vertIndex << "\n";
    for (int ii = 1; ii < 4; ++ii) {
      if (vertIndex != ii) {
        Vec inPlane = face.vertex[ii] - face.vertex[vertIndex];
        Vec outOfPlane = inPlane + shift;
        translations.push_back(inPlane);
        translations.push_back(outOfPlane);
      }
    }
  } 
  else if (status.second.first == Face::Location::EDGE) {
    int edgeIndex = status.second.second;
    //std::cout << "is edge " << edgeIndex << "\n";
    int oppIndex = (edgeIndex+3) % 4;
    Vec inPlane = face.vertex[oppIndex] - face.vertex[edgeIndex];
    Vec outOfPlane = inPlane + shift;
    translations.push_back(inPlane);
    translations.push_back(outOfPlane);
  }
}

namespace dem {
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<0>(
  DEMParticleShape shape, const Box& spatialDomain, Gradation& gradation);
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<1>(
  DEMParticleShape shape, const Box& spatialDomain, Gradation& gradation);
template DEMParticlePArray DEMParticleCreator::generateDEMParticles<2>(
  DEMParticleShape shape, const Box& spatialDomain, Gradation& gradation);
}
