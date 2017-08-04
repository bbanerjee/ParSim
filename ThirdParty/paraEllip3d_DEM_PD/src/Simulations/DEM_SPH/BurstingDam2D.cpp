#include <Simulations/DEM_SPH/BurstingDam2D.h>
#include <Core/Util/Utility.h>

void 
BurstingDam2D::execute(DiscreteElements* dem, SmoothParticleHydrodynamics* sph) 
{     
  if (dem->getMPIRank() == 0) {
    auto boundaryFile = InputParameters::get().datafile["boundaryFile"];
    auto particleFile = InputParameters::get().datafile["particleFile"];
    dem->readBoundary(boundaryFile);
    dem->readParticle(particleFile);
    sph->generateSPHParticle2D();
    dem->openDepositProg(progressInf, "deposit_progress");
    sph->openSPHTecplot(sphTecplotInf, "sph_results.dat");
  }

  // scatter particles only once; also updates grid for the first time
  dem->scatterDEMSPHParticle(); 
  pd->scatterSPHParticle(dem->getAllContainer()); 
  dem->calcNeighborRanks();

  // sph parameters
  auto gravAccel = util::getParam<REAL>("gravAccel");
  auto gravScale = util::getParam<REAL>("gravScale");
  auto waterLength = util::getParam<REAL>("waterLength");
  auto nSPHPoint = util::getParam<std::size_t>("nSPHPoint");
  
  REAL D = 5.0*gravAccel*gravScale*waterLength;
  REAL space_interval = waterLength/static_cast<REAL>(nSPHPoint-1);
  REAL smoothLength = 1.5*space_interval;
  REAL kernelSize = 3.0*smoothLength;
  REAL sphCellSize = dem->getGradation().getPtclMaxRadius() + kernelSize;
  REAL one_devide_h = 1.0/smoothLength;
  REAL factor_kernel = 7.0/(478.0*Pi*smoothLength*smoothLength);
  REAL factor_kernel_gradient = factor_kernel/smoothLength;
  REAL qmin = 0.7593;	// this is fixed for quintic kernel function
  REAL Wqmin = sph->kernelFunction(qmin);
  REAL p1 = 4;	// for the Lennard-Jones boundary forces 
  REAL p2 = 2;

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  auto timeStep = util::getParam<REAL>("timeStep");
  auto timeAccrued = util::getParam<REAL>("timeAccrued");

  auto netStep = endStep - startStep + 1;
  auto netSnap = endSnap - startSnap + 1;
  REAL timeCount = 0;
  REAL timeIncr  = timeStep * netStep;
  REAL timeTotal = timeAccrued + timeStep * netStep;

  REAL iteration = startStep;
  std::size_t iterSnap = startSnap;
  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    / Create the output writer in the master process
    // <outputFolder> rigidInc.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    dem->createOutputWriter(outputFolder, iterSnap-1);

    dem->writeBoundaryToFile();
    dem->writeGridToFile();
    dem->writeParticlesToFile(iterSnap);
    dem->printBdryContact();
    sph->printSPHTecplot(sphTecplotInf, iterSnap-1);
  }

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);
  std::cerr << "Proc = " << s_mpiRank << " outputFolder = " << outputFolder << "\n";

  // initialization SPH velocity based on equation(4.3) in 
  // http://www.artcompsci.org/vol_1/v1_web/node34.html
  sph->initialSPHVelocity2D();	

  while (timeAccrued < timeTotal) { 

    bool toCheckTime = (iteration + 1) % (netStep / netSnap) == 0;

    // update position and density of SPH particles based on equation (4.1)
    sph->updateSPHLeapFrogPositionDensity();
    sph->migrateSPHParticle();

    REAL commuT = migraT = gatherT = totalT = 0;  
    REAL time1 = time2 = time3 = 0;  
    REAL time0 = MPI_Wtime();

    dem->commuParticle();
    sph->commuSPHParticle();

    if (toCheckTime) {
      time1 = MPI_Wtime(); 
      commuT = time1 - time0;
    }

    // use values from last step, must call before findContact (which clears data)
    dem->calcTimeStep(); 
    dem->findContact();
    if (dem->isBdryProcess()) dem->findBdryContact();

    dem->clearContactForce();
    dem->internalForce();
    if (dem->isBdryProcess()) dem->boundaryForce();

    dem->dragForce();

    sph->calculateSPHDensityDotVelocityDotLinkedList2D();

    // For the 2D simulation
    // fix the y for all SPH points
    for (auto& particle : sph->getSPHParticleVec()) {
      switch (particle->getType()) {

      // free sph particle
      case 1: 
        particle->fixY();
        break;

      // ghost sph particle
      case 2: 
        break;

      // boundary sph particle
      case 3: 
        particle->fixXYZ();
        break;

      default:
        std::cout << "SPH particle type of pta should be 1, 2 or 3!" << std::endl;
        exit(-1);
      } // switch
    }

    dem->updateParticle();   

    // update velocity of SPH particles based on equation (4.2)
    sph->updateSPHLeapFrogVelocity();	

    // updateGridMaxZ();	// do not update grid in DEM-SPH coupling model
    // updateGridMaxZ() is for deposition or explosion. If they go out of side walls, particles are discarded.
    // updateGrid() updates all six directions, thus side walls may "disappear" if particles go far out of side walls 
    // and cause some grids to extrude out of side walls.

    timeCount += timeStep;
    timeAccrued += timeStep;
    if (timeCount >= timeIncr/netSnap) { 

      if (toCheckTime) {
        time2 = MPI_Wtime();
      }

      dem->gatherParticle();
      sph->gatherSPHParticle();
      dem->gatherBdryContact();
      dem->gatherEnergy(); 
      
      if (toCheckTime) {
        time3 = MPI_Wtime(); 
        gatherT = time3 - time2;
      }

      if (dem->getMPIRank() == 0) {
        dem->writeBoundaryToFile();
        dem->writeGridToFile();
        dem->writeParticlesToFile(iterSnap);
        dem->printBdryContact();
        sph->printSPHTecplot(sphTecplotInf, iterSnap);
        printDepositProg(progressInf);
      }
      printContact(combine(outputFolder, "bursting_contact_", iterSnap, 3));
    
      timeCount = 0;
      ++iterSnap;
    }

    // late release because printContact refers to received particles
    dem->releaseRecvParticle(); 
    // late release since calculateSPHDensityDotVelocityDotLinkedList3D 
    // will use these communicated sph particles
    sph->releaseRecvSPHParticle(); 

    if (toCheckTime) {
      time1 = MPI_Wtime();
    }
    dem->migrateParticle();

    if (toCheckTime) {
      time2 = MPI_Wtime(); 
      migraT = time2 - time1; 
      totalT = time2 - time0;
    }

    if (dem->getMPIRank() == 0 && toCheckTime) {
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT 
               << std::setw(OWID) << migraT 
               << std::setw(OWID) << totalT - commuT - migraT 
               << std::setw(OWID) << totalT 
               << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;
    }

    ++iteration;
  } 

  if (dem->getMPIRank() == 0) {
    dem->closeProg(progressInf);
  }
}

