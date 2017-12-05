#include <Simulations/DEM/PeriodicBCAxisymmetricStrainDriven.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
PeriodicBCAxisymmetricStrainDriven::execute(DiscreteElements* dem)
{
  // Read the boundary conditions in all processes
  // The file should be small enough
  std::string bcFile = util::getFilename("demBoundaryConditionFilename");
  if (bcFile == "none") {
    std::cout << "**ERROR** Input boundary condition file " << bcFile
              << " has not been specified in the input file\n";
    return;
  }
  dem->readBoundaryConditions(bcFile);

  dem->allowPatchDomainResize(Boundary::BoundaryID::XMINUS);
  dem->allowPatchDomainResize(Boundary::BoundaryID::XPLUS);
  dem->allowPatchDomainResize(Boundary::BoundaryID::YMINUS);
  dem->allowPatchDomainResize(Boundary::BoundaryID::YPLUS);
  dem->allowPatchDomainResize(Boundary::BoundaryID::ZMINUS);
  dem->allowPatchDomainResize(Boundary::BoundaryID::ZPLUS);

  std::ofstream progressInf;
  std::ofstream balancedInf;

  std::string outputFolder(".");
  std::string outputParticleCSVFile;
  REAL distX = 0, distY = 0, distZ = 0;

  if (dem->getMPIRank() == 0) {
    std::string boundaryFile = util::getFilename("boundaryFilename");
    std::string particleFile = util::getFilename("particleFilename");
    dem->readBoundary(boundaryFile);
    dem->readParticles(particleFile);

    // Create the output writer in the master process
    // <outputFolder> filename.pe3d </outputFolder>
    auto folderName =  util::getFilename("outputFolder");
    outputFolder = util::createOutputFolder(folderName);
    dem->createOutputWriter(outputFolder, 0);

    dem->openProgressOutputFile(progressInf, folderName + ".progress");
    dem->openProgressOutputFile(balancedInf, folderName + ".balanced");

    dem->writeBoundaryToFile();
    dem->writeBoundaryToFile(OrientedBox(dem->getSpatialDomain()));
    dem->writePatchGridToFile();
    dem->writeParticlesToFile(0);
    outputParticleCSVFile = combine("output_particles_", 0, 3);
    dem->printParticlesCSV(outputFolder, outputParticleCSVFile, 0);
    dem->printBoundary();
    dem->printBoundaryContacts();
    dem->getStartDimension(distX, distY, distZ);
  }

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  // Scatter the particles
  dem->scatterParticles();

  // Create a modifiable oriented bounding box from the
  // spatial domain that will be modified by each process
  const auto spatialDomain = dem->getSpatialDomain();
  OrientedBox periodicDomain(spatialDomain);

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");

  auto netStep = endStep - startStep + 1;
  auto netSnap = endSnap - startSnap + 1;

  auto iterSnap = startSnap;
  auto iteration = startStep;

  auto startTime = util::getParam<REAL>("timeAccrued");
  auto timeStep = util::getParam<REAL>("timeStep");
  //auto endTime = util::getParam<REAL>("endTime");
  auto endTime = startTime + (endStep - startStep)*timeStep;

  auto curTime = startTime;

  std::size_t numOverlaps = 0;
  //int numBisections = 0;
  while (iteration <= endStep && curTime < endTime) {

    // Communicate ghost particles to patches
    dem->communicateGhostParticles(iteration);

    // Get the particles and domain
    auto& patchParticles = dem->getModifiableParticleVec();

    //if (numOverlaps > 0 && numBisections < 5) {
    //  timeStep *= 0.5;
    //  ++numBisections;
    //}

    //std::cout << "Before calc: timeStep = " << timeStep << "\n";
    timeStep = dem->calcTimeStep(timeStep); 

    // Substep until number of overlapping particles is zero
    dem->findContact(iteration);
    numOverlaps = dem->numOverlappingParticles();

    // Switch particle type of patch particles from periodic to fixed
    //dem->switchParticleType(DEMParticle::DEMParticleType::BOUNDARY_PERIODIC,
    //                        DEMParticle::DEMParticleType::FIXED,
    //                        patchParticles);

    /*
    if (dem->isBoundaryProcess()) {
      dem->findBoundaryContacts(iteration);
    }
    */

    dem->initializeForces();
    dem->internalForce(timeStep, iteration);

    /*
    if (dem->isBoundaryProcess()) {
      dem->boundaryForce(timeStep, iteration);
    }
    */

    dem->updateParticles(timeStep, iteration);


    // Switch back particle type of patch particles from fixed to periodic
    //dem->switchParticleType(DEMParticle::DEMParticleType::FIXED,
    //                        DEMParticle::DEMParticleType::BOUNDARY_PERIODIC,
    //                        patchParticles);

    // Apply the particle boundary conditions to each set of patch particles
    // and update the periodic boundary
    dem->applyParticleBC(curTime, spatialDomain, periodicDomain, patchParticles);

    // Identify free particles that have crossed the periodic domain
    // boundary
    DEMParticleCreator particleCreator;
    auto newParticles = 
      particleCreator.updatePeriodicDEMParticles(periodicDomain, patchParticles);

    //dem->gatherBoundaryContacts(); 
    dem->updatePatchBox();

    if (iteration % (netStep / netSnap) == 0) {

      dem->gatherParticles();
      dem->gatherEnergy();

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile();
        dem->writeBoundaryToFile(periodicDomain);
        dem->printBoundary();
        dem->writePatchGridToFile();
        dem->writeParticlesToFile(iterSnap);
        outputParticleCSVFile = combine("output_particles_", iterSnap, 3);
        dem->printParticlesCSV(outputFolder, outputParticleCSVFile, 0);
        //dem->printBoundaryContacts();
        dem->appendToProgressOutputFile(progressInf, iteration, timeStep, distX, distY, distZ);
        std::cout << "Iteration = " << iteration 
                  << " cur time = " << curTime
                  << " timeStep = " << timeStep << "\n";
      }

      dem->printContact(combine(outputFolder, "contact_", iterSnap, 3));
      ++iterSnap;
    }

    dem->releaseReceivedParticles();
    dem->migrateParticles(iteration);

    ++iteration;
    curTime += timeStep;
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(balancedInf);
  }
}
