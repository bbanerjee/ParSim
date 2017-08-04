  void Assembly::drainageProblem() 
  {     
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      generateSPHParticleNoBottom3D();
      openDepositProg(progressInf, "deposit_progress");
      openSPHTecplot(sphTecplotInf, "sph_results.dat");
    }
    scatterDEMSPHParticle(); // scatter particles only once; also updates grid for the first time
    calcNeighborRanks();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1;
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    // sph parameters
    D = 5*(dem::Parameter::getSingleton().parameter["gravAccel"])*(dem::Parameter::getSingleton().parameter["gravScale"])*(allContainer.getMaxCorner().getZ()-allContainer.getMinCorner().getZ());
    space_interval = dem::Parameter::getSingleton().parameter["spaceInterval"];
    smoothLength = 1.5*space_interval;
    kernelSize = 3*smoothLength;
    one_devide_h = 1.0/smoothLength;
    factor_kernel = 1.0/(120.0*Pi*smoothLength*smoothLength*smoothLength);
    factor_kernel_gradient = 1.0/(120.0*Pi*smoothLength*smoothLength*smoothLength*smoothLength);
    REAL qmin = 0.7593;	// this is fixed for quintic kernel function
    Wqmin = kernelFunction(qmin);
    p1 = 4;	// for the Lennard-Jones boundary forces 
    p2 = 2;
    sphCellSize = gradation.getPtclMaxRadius() + kernelSize;

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    /**/REAL timeCount = 0;
    /**/timeAccrued = static_cast<REAL> (dem::Parameter::getSingleton().parameter["timeAccrued"]);
    /**/REAL timeIncr  = timeStep * netStep;
    /**/REAL timeTotal = timeAccrued + timeStep * netStep;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "bursting_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "bursting_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "bursting_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "bursting_bdrycntc_", iterSnap -1, 3));
      printSPHTecplot(sphTecplotInf, iterSnap-1);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "compuT" << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
	
//printSPHParticle(strcat(combineString(cstr0, "SPH_particle_", mpiRank, 3), ".dat"));

      initialSPHVelocity3D();
    /**/while (timeAccrued < timeTotal) { 
      bool toCheckTime = (iteration + 1) % (netStep / netSnap) == 0;

      // update position and density of SPH particles based on equation (4.1)
      updateSPHLeapFrogPositionDensity();
      migrateSPHParticle();

      commuT = migraT = gatherT = totalT = 0;  time0 = MPI_Wtime();
      commuParticle();
      commuSPHParticle();

      if (toCheckTime) time2 = MPI_Wtime(); commuT = time2 - time0;

      /**/calcTimeStep(); // use values from last step, must call before findContact (which clears data)
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      dragForce();

      calculateSPHDensityDotVelocityDotLinkedList3D();

      updateParticle();   
      updateSPHLeapFrogVelocity();	// update velocity of SPH particles based on equation (4.2)

//      updateGridMaxZ();	// do not update grid in DEM-SPH coupling model
      // updateGridMaxZ() is for deposition or explosion. If they go out of side walls, particles are discarded.
      // updateGrid() updates all six directions, thus side walls may "disappear" if particles go far out of side walls 
      // and cause some grids to extrude out of side walls.

      /**/timeCount += timeStep;
      /**/timeAccrued += timeStep;
      /**/if (timeCount >= timeIncr/netSnap) { 
	if (toCheckTime) time1 = MPI_Wtime();
	gatherParticle();
	gatherSPHParticle();	// particles that are out of container have been deleted
	gatherBdryContact();
	gatherEnergy(); if (toCheckTime) time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "bursting_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "bursting_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "bursting_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "bursting_bdrycntc_", iterSnap, 3));
	  printDepositProg(progressInf);
	  printSPHTecplot(sphTecplotInf, iterSnap);
	}
	printContact(combineString(cstr, "bursting_contact_", iterSnap, 3));
      
	/**/timeCount = 0;
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      releaseRecvSPHParticle(); // late release since calculateSPHDensityDotVelocityDotLinkedList3D will use these communicated sph particles

      if (toCheckTime) time1 = MPI_Wtime();
      migrateParticle();

      if (toCheckTime) time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && toCheckTime) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT - commuT - migraT << std::setw(OWID) << totalT 
		 <<std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;
      ++iteration;
    } 
  
    if (mpiRank == 0) closeProg(progressInf);
  }

