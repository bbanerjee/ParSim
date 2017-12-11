#include <Simulations/DEM/PeriodicBCComputeStressStrain.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <InputOutput/DEMParticleFileReader.h>
#include <Boundary/BoundaryFileReader.h>
#include <Core/Util/Utility.h>

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"


#include <experimental/filesystem>

using namespace dem;
using util::combine;
namespace fs = std::experimental::filesystem;

void
PeriodicBCComputeStressStrain::execute(DiscreteElements* dem)
{
  // 1) Select the list of files containing data from a simulation
  std::string directory = util::getFilename("inputDataDirectory");
  std::string domain_filename = "oriented_domain";
  std::string particle_filename = "particle";
  std::vector<std::string> filenames;
  for (auto& file : fs::directory_iterator(directory)) {
    std::string filename = fs::path(file).filename();
    auto found = filename.find(".vtu");
    if (found != std::string::npos) {
      filenames.push_back(filename);
    }
  }
  std::sort(filenames.begin(), filenames.end());
  std::vector<std::string> domain_files;
  std::vector<std::string> particle_files;
  for (auto& filename : filenames) {
    auto found = filename.find(domain_filename);
    if (found != std::string::npos) {
      domain_files.push_back(directory + "/" + filename);
    }
    found = filename.find(particle_filename);
    if (found != std::string::npos) {
      particle_files.push_back(directory + "/" + filename);
    }
  }
  for (auto& filename : domain_files) {
    std::cout << filename << "\n";
  }
  for (auto& filename : particle_files) {
    std::cout << filename << "\n";
  }

  // 2) For each file
  //  a) Read the output boundary data from a periodic BC simulation
  //  b) Read the output particle data from a periodic BC simulation
  REAL youngModulus = util::getParam<REAL>("young");
  REAL poissonRatio = util::getParam<REAL>("poisson");

  BoundaryFileReader boundaryReader;
  DEMParticleFileReader particleReader;
  int iteration = 0;
  for (auto& particle_file : particle_files) {
    OrientedBox spatialDomain(Box(Vec(0,0,0), Vec(1,1,1)));
    boundaryReader.readVTK(domain_files[iteration], spatialDomain);
    std::cout << "Domain:" << spatialDomain << "\n";

    DEMParticlePArray particles;
    particleReader.readVTK(particle_file, youngModulus, poissonRatio, particles);
    //for (auto& particle :particles) {
    //  std::cout << "Particles:" << *particle << "\n";
    //}

    createTessellation(particles);

  }
}

void
PeriodicBCComputeStressStrain::createTessellation(
  const DEMParticlePArray& particles) 
{
  std::stringstream coords;
  coords << 3 << " " << particles.size() << " ";
  for (const auto& particle : particles) {
    coords << particle->currentPosition().x() << " "
           << particle->currentPosition().y() << " "
           << particle->currentPosition().z() << " ";
  }
  std::cout << "Coords: \n" << coords.str() << "\n";

  orgQhull::RboxPoints box;
  box.appendPoints(coords);

  orgQhull::Qhull hull;
  hull.runQhull(box, "d Qbb Qt i");

  std::stringstream elements;
  hull.setOutputStream(&elements);
  hull.outputQhull();
  std::cout << "Elements: \n" << elements.str() << "\n";
}