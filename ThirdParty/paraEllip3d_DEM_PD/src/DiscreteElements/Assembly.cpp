//
//                                    x1(1)
//                          ---------------------
//                         /                    /|
//                        /                    / |
//                       /                    /  |
//                      /        z2(6)       /   |
//                     /                    /    |
//                    /                    /     |                    z
//                   /                    /      |                    |
//                  |---------------------       |                    |
//            y1(3) |                    | y2(4) |                    |____ y
//                  |                    |       /                   /
//                  |                    |      /                   /
//                  |        x2(2)       |     /                   x
//                  |                    |    /
//                  |                    |   /
//                  |                    |  /
//                  |                    | /
//                  |                    |/
//                  ----------------------
//                             z1(5)
//
//
//    It is preferable to use the description of surface x1, x2, y1, y2, z1, z2,
//    where x1 < x2, y1 < y2, z1 < z2. We also use surface 1, 2, 3, 4, 5, 6 accordingly.
//

#include <DiscreteElements/Assembly.h>
#include <Core/Const/const.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <cassert>
#include <utility>
#include <algorithm>
#include <sys/time.h>
#include <omp.h>

//#define BINNING
//#define TIME_PROFILE

//static time_t timeStamp; // for file timestamping
//static struct timeval time_w1, time_w2; // for wall-clock time record
//static struct timeval time_p1, time_p2; // for internal wall-clock time profiling, can be used on any piece of code
//static struct timeval time_r1, time_r2; // for internal wall-clock time profiling for contact resolution only (excluding space search)

namespace dem {

  struct timeval timediff(const struct timeval &time1, const struct timeval &time2) {
    struct timeval diff;
    diff.tv_sec = time2.tv_sec - time1.tv_sec; 
    if( ( diff.tv_usec = time2.tv_usec - time1.tv_usec) < 0) {
      diff.tv_usec += 1000000;
      --diff.tv_sec;
    }
    return(diff); 
  }


  long int timediffmsec(const struct timeval &time1, const struct timeval &time2) {
    struct timeval diff = timediff(time1, time2);
    return(diff.tv_sec * 1000000 + diff.tv_usec); 
  }


  REAL timediffsec(const struct timeval &time1, const struct timeval &time2) {
    return( (REAL) timediffmsec(time1,time2) / 1.0e+6);
  }


  char *combineString(char *cstr, const char *str, std::size_t num, std::size_t width) {
    std::string obj(str);
    std::stringstream ss;
    ss << std::setw(width) << std::setfill('0') << std::right << num;
    obj += ss.str();
    return strcpy( cstr, obj.c_str() );
  }

  // input:   number percentage smaller from data file
  // output:  mass percentage smaller to disk file debugInf
  // purpose: let mass percentage smaller satisfy particle size distribution curve
  // method:  use trial and error method on number percentage until mass percentage is satisfied
  void Assembly::tuneMassPercentage() 
  {
    if (mpiRank == 0) {
      REAL minX = dem::Parameter::getSingleton().parameter["minX"];
      REAL minY = dem::Parameter::getSingleton().parameter["minY"];
      REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
      REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
      REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
      REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
      std::size_t particleLayers = dem::Parameter::getSingleton().parameter["particleLayers"];

      setContainer(Rectangle(minX, minY, minZ, maxX, maxY, maxZ));

      buildBoundary(5, "deposit_boundary_ini");
    
      std::size_t sieveNum = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["sieveNum"]);
      std::vector<REAL> percent(sieveNum), size(sieveNum);
      std::vector<std::pair<REAL, REAL> > &grada = dem::Parameter::getSingleton().gradation;
      assert(grada.size() == sieveNum);
      for (std::size_t i = 0; i < sieveNum; ++i) {
	percent[i] = grada[i].first;
	size[i] = grada[i].second;
      }
      REAL ratioBA = dem::Parameter::getSingleton().parameter["ratioBA"];
      REAL ratioCA = dem::Parameter::getSingleton().parameter["ratioCA"];
      setGradation(Gradation(sieveNum, percent, size, ratioBA, ratioCA));

      generateParticle(particleLayers, "float_particle_ini"); 

      // statistics of mass distribution
      Gradation massGrad = gradation;
      std::vector<REAL> &massPercent = massGrad.getPercent();
      std::vector<REAL> &massSize    = massGrad.getSize();
      for (std::size_t i = 0; i < massPercent.size(); ++i) massPercent[i] = 0;

      for (auto itr = allParticleVec.cbegin(); itr != allParticleVec.cend(); ++itr)
	for (int i = massPercent.size()-1; i >= 0 ; --i) { // do not use size_t for descending series
	  if ( (*itr)->getA() <= massSize[i] )
	    massPercent[i] += (*itr)->getMass();
	}
      REAL totalMass = massPercent[0];
      for (std::size_t i = 0; i < massPercent.size(); ++i) massPercent[i] /= totalMass;
      debugInf << std::endl << "mass percentage of particles:" << std::endl
	       << std::setw(OWID) << massPercent.size() << std::endl;
      for (std::size_t i = 0; i < massPercent.size(); ++i)
	debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID) << massSize[i] << std::endl;
    }

  }


  void Assembly::depositIntoContainer() 
  {
    if (mpiRank == 0) {
      REAL minX = dem::Parameter::getSingleton().parameter["minX"];
      REAL minY = dem::Parameter::getSingleton().parameter["minY"];
      REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
      REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
      REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
      REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
      std::size_t particleLayers = dem::Parameter::getSingleton().parameter["particleLayers"];

      setContainer(Rectangle(minX, minY, minZ, maxX, maxY, maxZ));

      buildBoundary(5, "deposit_boundary_ini");
    
      std::size_t sieveNum = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["sieveNum"]);
      std::vector<REAL> percent(sieveNum), size(sieveNum);
      std::vector<std::pair<REAL, REAL> > &grada = dem::Parameter::getSingleton().gradation;
      assert(grada.size() == sieveNum);
      for (std::size_t i = 0; i < sieveNum; ++i) {
	percent[i] = grada[i].first;
	size[i] = grada[i].second;
      }
      REAL ratioBA = dem::Parameter::getSingleton().parameter["ratioBA"];
      REAL ratioCA = dem::Parameter::getSingleton().parameter["ratioCA"];
      setGradation(Gradation(sieveNum, percent, size, ratioBA, ratioCA));

      generateParticle(particleLayers, "float_particle_ini"); 
    }

    deposit("deposit_boundary_ini",
	    "float_particle_ini");

    if (mpiRank == 0) {
      setContainer(Rectangle(allContainer.getMinCorner().getX(),
			     allContainer.getMinCorner().getY(),
			     allContainer.getMinCorner().getZ(),
			     allContainer.getMaxCorner().getX(),
			     allContainer.getMaxCorner().getY(),
			     dem::Parameter::getSingleton().parameter["trimHeight"]));
      buildBoundary(6, "trim_boundary_ini");
      char cstr[50];
      std::size_t endSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
      trim(false,
	   combineString(cstr, "deposit_particle_", endSnap, 3),
	   "trim_particle_ini");
    }
  }


  void Assembly::resumeDepositIntoContainer() 
  {
    deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
	    dem::Parameter::getSingleton().datafile["particleFile"].c_str());
  
    if (mpiRank == 0) {
      setContainer(Rectangle(allContainer.getMinCorner().getX(),
			     allContainer.getMinCorner().getY(),
			     allContainer.getMinCorner().getZ(),
			     allContainer.getMaxCorner().getX(),
			     allContainer.getMaxCorner().getY(),
			     dem::Parameter::getSingleton().parameter["trimHeight"]));
      buildBoundary(6, "trim_boundary_ini");
      char cstr[50];
      std::size_t endSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
      trim(false,
	   combineString(cstr, "deposit_particle_", endSnap, 3),
	   "trim_particle_ini");
    }
    
  }


  void Assembly::proceedFromPreset() 
  {
    deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
	    dem::Parameter::getSingleton().datafile["particleFile"].c_str());
  }


  void Assembly::deposit(const char *inputBoundary,
			 const char *inputParticle) 
  {     
    if (mpiRank == 0) {
      readBoundary(inputBoundary); 
      readParticle(inputParticle);
      openDepositProg(progressInf, "deposit_progress");
    }
    scatterParticle(); // scatter particles only once; also updates grid for the first time

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1;
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    /**/REAL timeCount = 0;
    /**/timeAccrued = static_cast<REAL> (dem::Parameter::getSingleton().parameter["timeAccrued"]);
    /**/REAL timeTotal = timeAccrued + timeStep * netStep;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "deposit_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "deposit_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "deposit_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "deposit_bdrycntc_", iterSnap -1, 3));
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "compuT" << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    /**/while (timeAccrued < timeTotal) { 
      //while (iteration <= endStep) {
      bool toCheckTime = (iteration + 1) % (netStep / netSnap) == 0;

      commuT = migraT = gatherT = totalT = 0;  time0 = MPI_Wtime();
      commuParticle(); if (toCheckTime) time2 = MPI_Wtime(); commuT = time2 - time0;

      /**/calcTimeStep(); // use values from last step, must call before findConact (which clears data)
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();   
      updateGrid(); // universal; updateGridMaxZ() for deposition only

      /**/timeCount += timeStep;
      /**/timeAccrued += timeStep;
      /**/if (timeCount >= timeTotal/netSnap) { 
	//if (iteration % (netStep / netSnap) == 0) {
	if (toCheckTime) time1 = MPI_Wtime();
	gatherParticle();
	gatherBdryContact();
	gatherEnergy(); if (toCheckTime) time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "deposit_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "deposit_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "deposit_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "deposit_bdrycntc_", iterSnap, 3));
	  printDepositProg(progressInf);
	}
	printContact(combineString(cstr, "deposit_contact_", iterSnap, 3));
      
	/**/timeCount = 0;
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      if (toCheckTime) time1 = MPI_Wtime();
      migrateParticle(); if (toCheckTime) time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && toCheckTime) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT - commuT - migraT << std::setw(OWID) << totalT 
		 <<std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;
      ++iteration;
    } 
  
    if (mpiRank == 0) closeProg(progressInf);
  }


  void Assembly::expandCavityParticle() 
  {
    if (mpiRank == 0) {
      const char *inputParticle = dem::Parameter::getSingleton().datafile["particleFile"].c_str();
      REAL percent = dem::Parameter::getSingleton().parameter["expandPercent"];
      readParticle(inputParticle); 
    
      REAL x1 = dem::Parameter::getSingleton().parameter["cavityMinX"];
      REAL y1 = dem::Parameter::getSingleton().parameter["cavityMinY"];
      REAL z1 = dem::Parameter::getSingleton().parameter["cavityMinZ"];
      REAL x2 = dem::Parameter::getSingleton().parameter["cavityMaxX"];
      REAL y2 = dem::Parameter::getSingleton().parameter["cavityMaxY"];
      REAL z2 = dem::Parameter::getSingleton().parameter["cavityMaxZ"];
    
      ParticlePArray cavityParticleVec;
      Vec center;
    
      for (auto it = allParticleVec.cbegin(); it != allParticleVec.cend(); ++it ){
	center=(*it)->getCurrPos();
	if(center.getX() > x1 && center.getX() < x2 &&
	   center.getY() > y1 && center.getY() < y2 &&
	   center.getZ() > z1 && center.getZ() < z2 ) {
	  (*it)->expand(percent);
	  cavityParticleVec.push_back(*it);
	}
      }
    
      printParticle("cavity_particle_ini", cavityParticleVec);
      printParticle("expand_particle_ini");
    }
  
    deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
	    "expand_particle_ini");
  }


  void Assembly::resumeExpandCavityParticle() 
  {
    deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
	    dem::Parameter::getSingleton().datafile["particleFile"].c_str());
  }


  void Assembly::isotropic() 
  {
    std::size_t isotropicType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["isotropicType"]);
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openCompressProg(progressInf, "isotropic_progress");
      openCompressProg(balancedInf, "isotropic_balanced");
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1;
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL sigmaEnd, sigmaInc, sigmaVar;
    std::size_t sigmaDiv;
  
    sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
    sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];
    std::vector<REAL> &sigmaPath = dem::Parameter::getSingleton().sigmaPath;
    std::size_t sigma_i = 0;

    if (isotropicType == 1) 
      sigmaVar = sigmaEnd;
    else if (isotropicType == 2) {
      REAL sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
      sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
      sigmaVar = sigmaStart;
    } else if (isotropicType == 3) {
      sigmaVar = sigmaPath[sigma_i];
      sigmaInc = (sigmaPath[sigma_i+1] -  sigmaPath[sigma_i])/sigmaDiv;
      sigmaEnd = sigmaPath[sigmaPath.size()-1];
    }

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "isotropic_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "isotropic_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "isotropic_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "isotropic_bdrycntc_", iterSnap -1, 3));
      printBoundary(combineString(cstr0, "isotropic_boundary_", iterSnap - 1, 3));
      getStartDimension(distX, distY, distZ);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    while (iteration <= endStep) {
      commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
      updateBoundary(sigmaVar, "isotropic");
      updateGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	time1 = MPI_Wtime();
	gatherParticle();
	gatherEnergy(); time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "isotropic_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "isotropic_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "isotropic_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "isotropic_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "isotropic_boundary_", iterSnap, 3));
	  printCompressProg(progressInf, distX, distY, distZ);
	}
	printContact(combineString(cstr, "isotropic_contact_", iterSnap, 3));      
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      time1 = MPI_Wtime();
      migrateParticle(); time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (isotropicType == 1) {
	if (tractionErrorTol(sigmaVar, "isotropic")) {
	  if (mpiRank == 0) {
	    printParticle("isotropic_particle_end");
	    printBdryContact("isotropic_bdrycntc_end");
	    printBoundary("isotropic_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      } else if (isotropicType == 2) {
	if (tractionErrorTol(sigmaVar, "isotropic")) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVar += sigmaInc;
	}
	if (tractionErrorTol(sigmaEnd, "isotropic")) {
	  if (mpiRank == 0) {
	    printParticle("isotropic_particle_end");
	    printBdryContact("isotropic_bdrycntc_end");
	    printBoundary("isotropic_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      }

      if (isotropicType == 3) {
	if (tractionErrorTol(sigmaVar, "isotropic")) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVar += sigmaInc;
	  if (sigmaVar == sigmaPath[sigma_i+1]) {
	    sigmaVar = sigmaPath[++sigma_i];
	    sigmaInc = (sigmaPath[sigma_i+1] -  sigmaPath[sigma_i])/sigmaDiv;
	  }
	}
	if (tractionErrorTol(sigmaEnd, "isotropic")) {
	  if (mpiRank == 0) {
	    printParticle("isotropic_particle_end");
	    printBdryContact("isotropic_bdrycntc_end");
	    printBoundary("isotropic_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      }

      ++iteration;
    } 
  
    if (mpiRank == 0) {
      closeProg(progressInf);
      closeProg(balancedInf);
    }
  }


  void Assembly::odometer() 
  {
    std::size_t odometerType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["odometerType"]);
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openCompressProg(progressInf, "odometer_progress");
      openCompressProg(balancedInf, "odometer_balanced");
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1;
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL sigmaEnd, sigmaInc, sigmaVar;
    std::size_t sigmaDiv;
  
    sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
    sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];
    std::vector<REAL> &sigmaPath = dem::Parameter::getSingleton().sigmaPath;
    std::size_t sigma_i = 0;

    if (odometerType == 1) {
      REAL sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
      sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
      sigmaVar = sigmaStart;
    } else if (odometerType == 2) {
      sigmaVar = sigmaPath[sigma_i];
      sigmaInc = (sigmaPath[sigma_i+1] -  sigmaPath[sigma_i])/sigmaDiv;
      sigmaEnd = sigmaPath[sigmaPath.size()-1];
    }

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "odometer_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "odometer_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "odometer_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "odometer_bdrycntc_", iterSnap -1, 3));
      printBoundary(combineString(cstr0, "odometer_boundary_", iterSnap - 1, 3));
      getStartDimension(distX, distY, distZ);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    while (iteration <= endStep) {
      commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
      updateBoundary(sigmaVar, "odometer");
      updateGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	time1 = MPI_Wtime();
	gatherParticle();
	gatherEnergy(); time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "odometer_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "odometer_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "odometer_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "odometer_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "odometer_boundary_", iterSnap, 3));
	  printCompressProg(progressInf, distX, distY, distZ);
	}
	printContact(combineString(cstr, "odometer_contact_", iterSnap, 3));      
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      time1 = MPI_Wtime();
      migrateParticle(); time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (odometerType == 1) {
	if (tractionErrorTol(sigmaVar, "odometer")) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVar += sigmaInc;
	}
	if (tractionErrorTol(sigmaEnd, "odometer")) {
	  if (mpiRank == 0) {
	    printParticle("odometer_particle_end");
	    printBdryContact("odometer_bdrycntc_end");
	    printBoundary("odometer_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      } else if (odometerType == 2) {
	if (tractionErrorTol(sigmaVar, "odometer")) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVar += sigmaInc;
	  if (sigmaVar == sigmaPath[sigma_i+1]) {
	    sigmaVar = sigmaPath[++sigma_i];
	    sigmaInc = (sigmaPath[sigma_i+1] -  sigmaPath[sigma_i])/sigmaDiv;
	  }
	}
	if (tractionErrorTol(sigmaEnd, "odometer")) {
	  if (mpiRank == 0) {
	    printParticle("odometer_particle_end");
	    printBdryContact("odometer_bdrycntc_end");
	    printBoundary("odometer_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      }

      ++iteration;
    } 
  
    if (mpiRank == 0) {
      closeProg(progressInf);
      closeProg(balancedInf);
    }
  }


  void Assembly::triaxial() 
  {
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openCompressProg(progressInf, "triaxial_progress");
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1; 
    REAL sigmaConf = dem::Parameter::getSingleton().parameter["sigmaConf"];
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "triaxial_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "triaxial_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "triaxial_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "triaxial_bdrycntc_", iterSnap -1, 3));
      printBoundary(combineString(cstr0, "triaxial_boundary_", iterSnap - 1, 3));
      getStartDimension(distX, distY, distZ);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    while (iteration <= endStep) {
      commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      // displacement control relies on constant time step, so do not call calcTimeStep().
      //calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
      updateBoundary(sigmaConf, "triaxial");
      updateGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	time1 = MPI_Wtime();
	gatherParticle();
	gatherEnergy(); time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "triaxial_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "triaxial_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "triaxial_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "triaxial_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "triaxial_boundary_", iterSnap, 3));
	  //printCompressProg(progressInf, distX, distY, distZ); // redundant
	}
	printContact(combineString(cstr, "triaxial_contact_", iterSnap, 3));      
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      time1 = MPI_Wtime();
      migrateParticle(); time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (mpiRank == 0 && iteration % 10 == 0)
	printCompressProg(progressInf, distX, distY, distZ);

      // no break condition, just through top/bottom displacement control
      ++iteration;
    } 

    if (mpiRank == 0) {
      printParticle("triaxial_particle_end");
      printBdryContact("triaxial_bdrycntc_end");
      printBoundary("triaxial_boundary_end");
      printCompressProg(progressInf, distX, distY, distZ);
    }
  
    if (mpiRank == 0) closeProg(progressInf);
  }


  void Assembly::planeStrain() 
  {
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openCompressProg(progressInf, "plnstrn_progress");
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1; 
    REAL sigmaConf = dem::Parameter::getSingleton().parameter["sigmaConf"];
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "plnstrn_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "plnstrn_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "plnstrn_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "plnstrn_bdrycntc_", iterSnap -1, 3));
      printBoundary(combineString(cstr0, "plnstrn_boundary_", iterSnap - 1, 3));
      getStartDimension(distX, distY, distZ);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    while (iteration <= endStep) {
      commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      // displacement control relies on constant time step, so do not call calcTimeStep().
      //calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
      updateBoundary(sigmaConf, "plnstrn");
      updateGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	time1 = MPI_Wtime();
	gatherParticle();
	gatherEnergy(); time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "plnstrn_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "plnstrn_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "plnstrn_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "plnstrn_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "plnstrn_boundary_", iterSnap, 3));
	  //printCompressProg(progressInf, distX, distY, distZ); // redundant
	}
	printContact(combineString(cstr, "plnstrn_contact_", iterSnap, 3));      
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      time1 = MPI_Wtime();
      migrateParticle(); time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (mpiRank == 0 && iteration % 10 == 0)
	printCompressProg(progressInf, distX, distY, distZ);

      // no break condition, just through top/bottom displacement control
      ++iteration;
    } 

    if (mpiRank == 0) {
      printParticle("plnstrn_particle_end");
      printBdryContact("plnstrn_bdrycntc_end");
      printBoundary("plnstrn_boundary_end");
      printCompressProg(progressInf, distX, distY, distZ);
    }
  
    if (mpiRank == 0) closeProg(progressInf);
  }


  void Assembly::trueTriaxial() 
  {
    std::size_t trueTriaxialType = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["trueTriaxialType"]);
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openCompressProg(progressInf, "trueTriaxial_progress");
      openCompressProg(balancedInf, "trueTriaxial_balanced");
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1; 
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    REAL sigmaStart, sigmaEndZ, sigmaEndX, sigmaEndY;
    REAL sigmaDiv, sigmaIncZ, sigmaIncX, sigmaIncY, sigmaVarZ, sigmaVarX, sigmaVarY;
    REAL sigmaInit[3], sigmaEnd, sigmaInc, sigmaVar;
    std::size_t changeDirc;
    sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];

    if (trueTriaxialType == 1) {
      sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
      sigmaEndZ  = dem::Parameter::getSingleton().parameter["sigmaEndZ"];
      sigmaEndX  = dem::Parameter::getSingleton().parameter["sigmaEndX"];
      sigmaEndY  = dem::Parameter::getSingleton().parameter["sigmaEndY"];
      sigmaIncZ  = (sigmaEndZ - sigmaStart) / sigmaDiv;
      sigmaIncX  = (sigmaEndX - sigmaStart) / sigmaDiv;
      sigmaIncY  = (sigmaEndY - sigmaStart) / sigmaDiv;
      sigmaVarZ  = sigmaStart;
      sigmaVarX  = sigmaStart;
      sigmaVarY  = sigmaStart;
    } else if (trueTriaxialType == 2) {
      sigmaInit[0] = dem::Parameter::getSingleton().parameter["sigmaStartX"];
      sigmaInit[1] = dem::Parameter::getSingleton().parameter["sigmaStartY"];
      sigmaInit[2] = dem::Parameter::getSingleton().parameter["sigmaStartZ"];
      sigmaEnd     = dem::Parameter::getSingleton().parameter["sigmaEnd"];
      changeDirc   = dem::Parameter::getSingleton().parameter["changeDirc"];
      sigmaInc     = (sigmaEnd - sigmaInit[changeDirc]) / sigmaDiv;
      sigmaVar     = sigmaInit[changeDirc];
    }

    REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "trueTriaxial_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "trueTriaxial_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "trueTriaxial_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "trueTriaxial_bdrycntc_", iterSnap -1, 3));
      printBoundary(combineString(cstr0, "trueTriaxial_boundary_", iterSnap - 1, 3));
      getStartDimension(distX, distY, distZ);
    }
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    while (iteration <= endStep) {
      commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary

      if (trueTriaxialType == 1)
	updateBoundary(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY);
      else if (trueTriaxialType == 2) {
	REAL sigmaX = 0.0, sigmaY = 0.0, sigmaZ = 0.0;
	if (changeDirc == 0) {
	  sigmaX = sigmaVar;     sigmaY = sigmaInit[1]; sigmaZ = sigmaInit[2];
	} else if (changeDirc == 1) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaVar;     sigmaZ = sigmaInit[2];
	} else if (changeDirc == 2) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaInit[1]; sigmaZ = sigmaVar;
	}
	updateBoundary(sigmaZ, "trueTriaxial", sigmaX, sigmaY);
      }

      updateGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	time1 = MPI_Wtime();
	gatherParticle();
	gatherEnergy(); time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "trueTriaxial_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "trueTriaxial_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "trueTriaxial_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "trueTriaxial_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "trueTriaxial_boundary_", iterSnap, 3));
	  printCompressProg(progressInf, distX, distY, distZ);
	}
	printContact(combineString(cstr, "trueTriaxial_contact_", iterSnap, 3));      
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      time1 = MPI_Wtime();
      migrateParticle(); time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (trueTriaxialType == 1) {
	if (tractionErrorTol(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY)) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVarZ += sigmaIncZ;
	  sigmaVarX += sigmaIncX;
	  sigmaVarY += sigmaIncY;
	}
	if (tractionErrorTol(sigmaEndZ, "trueTriaxial", sigmaEndX, sigmaEndY)) {
	  if (mpiRank == 0) {
	    printParticle("trueTriaxial_particle_end");
	    printBdryContact("trueTriaxial_bdrycntc_end");
	    printBoundary("trueTriaxial_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      } else if (trueTriaxialType == 2) {
	REAL sigmaX = 0.0, sigmaY = 0.0, sigmaZ = 0.0;
	if (changeDirc == 0) {
	  sigmaX = sigmaVar;     sigmaY = sigmaInit[1]; sigmaZ = sigmaInit[2];
	} else if (changeDirc == 1) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaVar;     sigmaZ = sigmaInit[2];
	} else if (changeDirc == 2) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaInit[1]; sigmaZ = sigmaVar;
	}
	if (tractionErrorTol(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
	  if (mpiRank == 0) printCompressProg(balancedInf, distX, distY, distZ);
	  sigmaVar += sigmaInc;
	} 
   
	if (changeDirc == 0) {
	  sigmaX = sigmaEnd;     sigmaY = sigmaInit[1]; sigmaZ = sigmaInit[2];
	} else if (changeDirc == 1) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaEnd;     sigmaZ = sigmaInit[2];
	} else if (changeDirc == 2) {
	  sigmaX = sigmaInit[0]; sigmaY = sigmaInit[1]; sigmaZ = sigmaEnd;
	}
	if (tractionErrorTol(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
	  if (mpiRank == 0) {
	    printParticle("trueTriaxial_particle_end");
	    printBdryContact("trueTriaxial_bdrycntc_end");
	    printBoundary("trueTriaxial_boundary_end");
	    printCompressProg(balancedInf, distX, distY, distZ);
	  }
	  break;
	}
      }

      ++iteration;
    } 

    if (mpiRank == 0) {
      printParticle("trueTriaxial_particle_end");
      printBdryContact("trueTriaxial_bdrycntc_end");
      printBoundary("trueTriaxial_boundary_end");
      printCompressProg(progressInf, distX, distY, distZ);
    }

    if (mpiRank == 0) {
      closeProg(progressInf);
      closeProg(balancedInf);
    }
  }


  bool Assembly::tractionErrorTol(REAL sigma, std::string type, REAL sigmaX, REAL sigmaY) {
    // sigma implies sigmaZ
    REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];

    std::map<std::string, REAL> normalForce;
    REAL x1, x2, y1, y2, z1, z2;
    // do not use mergeBoundaryVec because each process calls this function.
    for(auto it = boundaryVec.cbegin(); it != boundaryVec.cend(); ++it) {
      std::size_t id = (*it)->getId();
      Vec normal = (*it)->getNormalForce();
      Vec point  = (*it)->getPoint();
      switch (id) {
      case 1: 
	normalForce["x1"] = fabs(normal.getX());
	x1 = point.getX();
	break;
      case 2:
	normalForce["x2"] = normal.getX();
	x2 = point.getX();
	break;
      case 3:
	normalForce["y1"] = fabs(normal.getY());
	y1 = point.getY();
	break;
      case 4:
	normalForce["y2"] = normal.getY();
	y2 = point.getY();
	break;
      case 5:
	normalForce["z1"] = fabs(normal.getZ());
	z1 = point.getZ();
	break;
      case 6:
	normalForce["z2"] = normal.getZ();
	z2 = point.getZ();
	break;
      }
    }
    REAL areaX = (y2 - y1) * (z2 - z1);
    REAL areaY = (z2 - z1) * (x2 - x1);
    REAL areaZ = (x2 - x1) * (y2 - y1);

    if (type.compare("isotropic") == 0)
      return (    fabs(normalForce["x1"]/areaX - sigma) / sigma <= tol
		  && fabs(normalForce["x2"]/areaX - sigma) / sigma <= tol
		  && fabs(normalForce["y1"]/areaY - sigma) / sigma <= tol
		  && fabs(normalForce["y2"]/areaY - sigma) / sigma <= tol
		  && fabs(normalForce["z1"]/areaZ - sigma) / sigma <= tol
		  && fabs(normalForce["z2"]/areaZ - sigma) / sigma <= tol );

    else if (type.compare("odometer") == 0)
      return ( fabs(normalForce["z1"]/areaZ - sigma) / sigma <= tol
	       && fabs(normalForce["z2"]/areaZ - sigma) / sigma <= tol );

    else if (type.compare("triaxial") == 0)
      return true; // always near equilibrium

    else if (type.compare("trueTriaxial") == 0)
      return (    fabs(normalForce["x1"]/areaX - sigmaX) / sigmaX <= tol
		  && fabs(normalForce["x2"]/areaX - sigmaX) / sigmaX <= tol
		  && fabs(normalForce["y1"]/areaY - sigmaY) / sigmaY <= tol
		  && fabs(normalForce["y2"]/areaY - sigmaY) / sigmaY <= tol
		  && fabs(normalForce["z1"]/areaZ - sigma)  / sigma  <= tol
		  && fabs(normalForce["z2"]/areaZ - sigma)  / sigma  <= tol );

    return false;
  }


  void Assembly::coupleWithSonicFluid() 
  {
    std::ofstream particleInf;
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      openDepositProg(progressInf, "couple_progress");
      openParticleProg(particleInf, "particle_progress");
      /*1*/ fluid.initParameter(allContainer, gradation);
      /*2*/ fluid.initialize();
    }
    scatterParticle();

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1;
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    REAL timeCount = 0;
    timeAccrued = static_cast<REAL> (dem::Parameter::getSingleton().parameter["timeAccrued"]);
    REAL timeTotal = timeAccrued + timeStep * netStep;
    if (mpiRank == 0) {
      plotBoundary(strcat(combineString(cstr0, "couple_bdryplot_", iterSnap - 1, 3), ".dat"));
      plotGrid(strcat(combineString(cstr0, "couple_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "couple_particle_", iterSnap - 1, 3));
      printBdryContact(combineString(cstr0, "couple_bdrycntc_", iterSnap -1, 3));
      /*3*/ fluid.plot(strcat(combineString(cstr0, "couple_fluidplot_", iterSnap -1, 3), ".dat")); 
    }
    /*
    if (mpiRank == 0)
      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;
    */
    while (timeAccrued < timeTotal) { 
      //REAL time0 = MPI_Wtime();
      commuParticle(); 
      //REAL time2 = MPI_Wtime(); 
      //REAL commuT = time2 - time0;

      calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();

      /*4*/ fluid.getParticleInfo(particleVec); // not allParticleVec
      /*5*/ fluid.runOneStep();
      /*6*/ fluid.calcParticleForce(particleVec, particleInf); // not allParticleVec
      /*7*/ fluid.penalize();

      internalForce();
      if (isBdryProcess()) boundaryForce();

      updateParticle();
      updateGridMaxZ();

      timeCount += timeStep;
      timeAccrued += timeStep;
      if (timeCount >= timeTotal/netSnap) { 
	//REAL time1 = MPI_Wtime();
	gatherParticle();
	gatherBdryContact();
	gatherEnergy(); 
        //time2 = MPI_Wtime(); 
        //REAL gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "couple_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "couple_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "couple_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "couple_bdrycntc_", iterSnap, 3));
	  printDepositProg(progressInf);
	  /*8*/ fluid.plot(strcat(combineString(cstr, "couple_fluidplot_", iterSnap, 3), ".dat"));
	}
	printContact(combineString(cstr, "couple_contact_", iterSnap, 3));
      
	timeCount = 0;
	++iterSnap;
      }

      releaseRecvParticle(); // late release because printContact refers to received particles
      //REAL time1 = MPI_Wtime();
      migrateParticle(); 
      //time2 = MPI_Wtime(); 
      //REAL migraT = time2 - time1; 
      //REAL totalT = time2 - time0;
      /*
      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;
      */

      ++iteration;
    } 
  
    if (mpiRank == 0) {
      closeProg(progressInf);
      closeProg(particleInf);
    }
  }


  // particleLayers:
  // 0 - one free particle
  // 1 - a horizontal layer of free particles
  // 2 - multiple layers of free particles
  void Assembly::generateParticle(std::size_t particleLayers,
				  const char *genParticle)
  {
    REAL young = dem::Parameter::getSingleton().parameter["young"];
    REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

    REAL x,y,z;
    std::size_t particleNum = 0;
    REAL diameter = gradation.getPtclMaxRadius()*2.0;

    REAL offset   = 0;
    REAL edge     = diameter;
    if (gradation.getSize().size() == 1 &&
	gradation.getPtclRatioBA() == 1.0 && 
	gradation.getPtclRatioCA() == 1.0) {
      edge   = diameter*2.0;
      offset = diameter*0.25;
    }
  
    REAL x1 = allContainer.getMinCorner().getX() + edge;
    REAL y1 = allContainer.getMinCorner().getY() + edge;
    REAL z1 = allContainer.getMinCorner().getZ() + diameter;
    REAL x2 = allContainer.getMaxCorner().getX() - edge;
    REAL y2 = allContainer.getMaxCorner().getY() - edge;
    //REAL z2 = allContainer.getMaxCorner().getZ() - diameter;
    REAL z2 = dem::Parameter::getSingleton().parameter["floatMaxZ"] - diameter;
    REAL x0 = allContainer.getCenter().getX();
    REAL y0 = allContainer.getCenter().getY();
    REAL z0 = allContainer.getCenter().getZ();

    if (particleLayers == 0) {      // just one free particle
      ParticleP newptcl = std::make_shared<Particle>(particleNum+1, 0, Vec(x0,y0,z0), 
                                                     gradation, young, poisson);
      allParticleVec.push_back(newptcl);
      particleNum++;
    }
    else if (particleLayers == 1) { // a horizontal layer of free particles
      for (x = x1; x - x2 < EPS; x += diameter)
	for (y = y1; y - y2 < EPS; y += diameter) {
	  ParticleP newptcl = std::make_shared<Particle>(particleNum+1, 0, Vec(x,y,z0), 
                                                         gradation, young, poisson);
	  allParticleVec.push_back(newptcl);
	  particleNum++;
	}
    }
    else if (particleLayers == 2) { // multiple layers of free particles
      for (z = z1; z - z2 < EPS; z += diameter) {
	for (x = x1 + offset; x - x2 < EPS; x += diameter)
	  for (y = y1 + offset; y - y2 < EPS; y += diameter) {
	    ParticleP newptcl = std::make_shared<Particle>(particleNum+1, 0, Vec(x,y,z), 
                                                           gradation, young, poisson);
	    allParticleVec.push_back(newptcl);
	    particleNum++;
	  }	
	offset *= -1;
      }
    } 

    printParticle(genParticle); 
  }
  

  void Assembly::trimOnly() {
    if (mpiRank == 0) {
      readBoundary(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
      trim(true,
	   dem::Parameter::getSingleton().datafile["particleFile"].c_str(),
	   "trim_particle_end");
    }
  }


  void Assembly::trim(bool toRebuild,
		      const char *inputParticle,
		      const char *trmParticle)
  {
    if (toRebuild) readParticle(inputParticle);
    trimHistoryNum = allParticleVec.size();

    Vec  v1 = allContainer.getMinCorner();
    Vec  v2 = allContainer.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
    REAL maxR = gradation.getPtclMaxRadius();
 
    // BB: Feb 2, 2017:
    // Not an efficient operation
    // Better approach may be to use a list if random access of vector 
    // members is not needed
    allParticleVec.erase(
      std::remove_if(allParticleVec.begin(), allParticleVec.end(),
                    [&x1, &y1, &z1, &x2, &y2, &z2, &maxR](ParticleP particle) {
                      Vec center = particle->getCurrPos();
                      if (center.getX() < x1 || center.getX() > x2 ||
	                  center.getY() < y1 || center.getY() > y2 ||
	                  center.getZ() < z1 || center.getZ() + maxR > z2) {
                        return true;
                      }
                      return false;
                    }),
      allParticleVec.end()
    );

    /*
    ParticlePArray::iterator itr;
    Vec center;
    for (auto itr = allParticleVec.begin(); itr != allParticleVec.end(); ) {
      center=(*itr)->getCurrPos();
      if(center.getX() < x1 || center.getX() > x2 ||
	 center.getY() < y1 || center.getY() > y2 ||
	 center.getZ() < z1 || center.getZ() + maxR > z2)
	{
	  delete (*itr); // release memory
	  itr = allParticleVec.erase(itr); 
	}
      else
	++itr;
    }
    */
  
    printParticle(trmParticle);
  }


  void Assembly::findParticleInRectangle(const Rectangle &container,
			  		 const ParticlePArray &inputParticle,
			  		 ParticlePArray &foundParticle) {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
    for (std::size_t pt = 0; pt < inputParticle.size(); ++pt) {
      Vec center = inputParticle[pt]->getCurrPos();
      // it is critical to use EPS
      if (center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	  center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	  center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS)
	foundParticle.push_back(inputParticle[pt]);
    }
  }


  void Assembly::findPeriParticleInRectangle(const Rectangle &container,
			  		     const std::vector<periDynamics::PeriParticle*> &inputParticle,
			  		     std::vector<periDynamics::PeriParticle*> &foundParticle) {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
    for (std::size_t pt = 0; pt < inputParticle.size(); ++pt) {
      Vec center = inputParticle[pt]->getCurrPosition();
      // it is critical to use EPS
      if (center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	  center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	  center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS)
	foundParticle.push_back(inputParticle[pt]);
    }
  }

  void Assembly::removeParticleOutRectangle() {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();

    // BB: Feb 2, 2017:
    // Not an efficient operation
    // Better approach may be to use a list if random access of vector 
    // members is not needed
    REAL epsilon = EPS;
    particleVec.erase(
      std::remove_if(particleVec.begin(), particleVec.end(),
                    [&x1, &y1, &z1, &x2, &y2, &z2, &epsilon](ParticleP particle) {
                      Vec center = particle->getCurrPos();
                      if ( !(center.getX() - x1 >= -epsilon && 
                             center.getX() - x2 < -epsilon &&
	                     center.getY() - y1 >= -epsilon && 
                             center.getY() - y2 < -epsilon &&
	                     center.getZ() - z1 >= -epsilon && 
                             center.getZ() - z2 < -epsilon) ) {
                        return true;
                      }
                      return false;
                    }),
      particleVec.end()
    );

    /*
    ParticlePArray::iterator itr;
    Vec center;
    //std::size_t flag = 0;

    for (itr = particleVec.begin(); itr != particleVec.end(); ) {
      center=(*itr)->getCurrPos();
      // it is critical to use EPS
      if ( !(center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	     center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	     center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS) )
	{
	  // debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
	  // << " removed=" << std::setw(3) << (*itr)->getId();	
	  // flag = 1;
	  delete (*itr); // release memory
	  itr = particleVec.erase(itr); 
	}
      else
	++itr;
    }
    */

    /*
      if (flag == 1) {
      debugInf << " now " << particleVec.size() << ": ";
      for (ParticlePArray::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

  }


  void Assembly::removePeriParticleOutRectangle() {
    Vec  v1 = container.getMinCorner();
    Vec  v2 = container.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();

    std::vector<periDynamics::PeriParticle*>::iterator itr;
    Vec center;
    //std::size_t flag = 0;

    for (itr = periParticleVec.begin(); itr != periParticleVec.end(); ) {
      center=(*itr)->getCurrPosition();
      // it is critical to use EPS
      if ( !(center.getX() - x1 >= -EPS && center.getX() - x2 < -EPS &&
	     center.getY() - y1 >= -EPS && center.getY() - y2 < -EPS &&
	     center.getZ() - z1 >= -EPS && center.getZ() - z2 < -EPS) )
	{
	  /*
	    debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank
	    << " removed=" << std::setw(3) << (*itr)->getId();	
	    flag = 1;
	  */
	  delete (*itr); // release memory
	  itr = periParticleVec.erase(itr); 
	}
      else
	++itr;
    }
    /*
      if (flag == 1) {
      debugInf << " now " << particleVec.size() << ": ";
      for (ParticlePArray::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

  }

  REAL Assembly::getPtclMaxX(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return -1/EPS;

    auto it = inputParticle.cbegin();
    REAL x0 = (*it)->getCurrPos().getX();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getX() > x0 )
	x0 = (*it)->getCurrPos().getX();
    }
    return x0;
  }


  REAL Assembly::getPtclMinX(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return 1/EPS;

    auto it = inputParticle.cbegin();
    REAL x0 = (*it)->getCurrPos().getX();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getX() < x0 )
	x0 = (*it)->getCurrPos().getX();
    }
    return x0;
  }


  REAL Assembly::getPtclMaxY(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return -1/EPS;

    auto it = inputParticle.cbegin();
    REAL y0 = (*it)->getCurrPos().getY();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getY() > y0 )
	y0 = (*it)->getCurrPos().getY();
    }
    return y0;
  }


  REAL Assembly::getPtclMinY(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return 1/EPS;

    auto it = inputParticle.cbegin();
    REAL y0 = (*it)->getCurrPos().getY();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getY() < y0 )
	y0 = (*it)->getCurrPos().getY();
    }
    return y0;
  }


  REAL Assembly::getPtclMaxZ(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return -1/EPS;

    auto it = inputParticle.cbegin();
    REAL z0 = (*it)->getCurrPos().getZ();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getZ() > z0 )
	z0 = (*it)->getCurrPos().getZ();
    }
    return z0;
  }


  REAL Assembly::getPtclMinZ(const ParticlePArray &inputParticle) const {
    if (inputParticle.size() == 0)
      return 1/EPS;

    auto it = inputParticle.cbegin();
    REAL z0 = (*it)->getCurrPos().getZ();
    for (; it != inputParticle.cend(); ++it) {
      if ( (*it)->getCurrPos().getZ() < z0 )
	z0 = (*it)->getCurrPos().getZ();
    }
    return z0;
  }


  void Assembly::setCommunicator(boost::mpi::communicator &comm) {
    boostWorld = comm;
    mpiWorld = MPI_Comm(comm);
    mpiProcX = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcX"]);
    mpiProcY = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcY"]);
    mpiProcZ = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcZ"]);
  
    // create Cartesian virtual topology (unavailable in boost.mpi) 
    int ndim = 3;
    int dims[3] = {mpiProcX, mpiProcY, mpiProcZ};
    int periods[3] = {0, 0, 0};
    int reorder = 0; // mpiRank not reordered
    MPI_Cart_create(mpiWorld, ndim, dims, periods, reorder, &cartComm);
    MPI_Comm_rank(cartComm, &mpiRank); 
    MPI_Comm_size(cartComm, &mpiSize);
    MPI_Cart_coords(cartComm, mpiRank, ndim, mpiCoords);
    mpiTag = 0;
    assert(mpiRank == boostWorld.rank());
    //debugInf << mpiRank << " " << mpiCoords[0] << " " << mpiCoords[1] << " " << mpiCoords[2] << std::endl;

    for (int iRank = 0; iRank < mpiSize; ++iRank) {
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      if (coords[0] == 0 || coords[0] == mpiProcX - 1 ||
	  coords[1] == 0 || coords[1] == mpiProcY - 1 ||
	  coords[2] == 0 || coords[2] == mpiProcZ - 1)
	bdryProcess.push_back(iRank);
    }
  }


  void Assembly::scatterParticle() {
    // partition particles and send to each process
    if (mpiRank == 0) { // process 0
      setGrid(Rectangle(grid.getMinCorner().getX(),
			grid.getMinCorner().getY(),
			grid.getMinCorner().getZ(),
			grid.getMaxCorner().getX(),
			grid.getMaxCorner().getY(),
			getPtclMaxZ(allParticleVec) + gradation.getPtclMaxRadius() ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      ParticlePArray tmpParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
	tmpParticleVec.clear(); // do not release memory!
	int ndim = 3;
	int coords[3];
	MPI_Cart_coords(cartComm, iRank, ndim, coords);
	Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
			    v1.getY() + vspan.getY() / mpiProcY * coords[1],
			    v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
			    v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
			    v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
			    v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
	findParticleInRectangle(container, allParticleVec, tmpParticleVec);
	if (iRank != 0)
	  reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
	if (iRank == 0) {
	  particleVec.resize(tmpParticleVec.size());
	  for (auto i = 0u; i < particleVec.size(); ++i)
	    particleVec[i] = std::make_shared<Particle>(*tmpParticleVec[i]); // default synthesized copy constructor
	} // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, particleVec);
    }

    // content of allParticleVec may need to be printed, so do not clear it. 
    //if (mpiRank == 0) releaseGatheredParticle();

    // broadcast necessary info
    broadcast(boostWorld, gradation, 0);
    broadcast(boostWorld, boundaryVec, 0);
    broadcast(boostWorld, allContainer, 0);
    broadcast(boostWorld, grid, 0);
  }


  // this is to scatter the dem and sph particle
  // two point: (1) Partition the peri-points before any bonds constructed, all peri-points will be treated as free peri-points. 
  //		    This is also what the coupling model shows, i.e. even the peri-points that are bonded to DEM particles still have properties of free peri-points
  //            (2) Before partition there are no any bonds in PeriParticle, after partition, constructNeighbor() is called, then these peri-bonds, boundary-bonds,
  //		    and peri-DEM-bonds will be generated
  void Assembly::scatterDEMPeriParticle() {
    // partition particles and send to each process
    if (mpiRank == 0) { // process 0
      setGrid(Rectangle(allContainer.getMinCorner().getX() - point_interval*0.2,
			allContainer.getMinCorner().getY() - point_interval*0.2,
			allContainer.getMinCorner().getZ() - point_interval*0.2,
			allContainer.getMaxCorner().getX() + point_interval*0.2,
			allContainer.getMaxCorner().getY() + point_interval*0.2,
			allContainer.getMaxCorner().getZ() + point_interval*0.2 ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;

      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      ParticlePArray tmpParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
	tmpParticleVec.clear(); // do not release memory!
	int ndim = 3;
	int coords[3];
	MPI_Cart_coords(cartComm, iRank, ndim, coords);
	Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
			    v1.getY() + vspan.getY() / mpiProcY * coords[1],
			    v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
			    v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
			    v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
			    v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
	findParticleInRectangle(container, allParticleVec, tmpParticleVec);
	if (iRank != 0)
	  reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
	  // before send, the SPHParticle.demParticle == NULL, since NULL is assigned when SPHParticle is created
	if (iRank == 0) {
	  particleVec.resize(tmpParticleVec.size());
	  for (auto i = 0u; i < particleVec.size(); ++i)
	    particleVec[i] = std::make_shared<Particle>(*tmpParticleVec[i]);  // default synthesized copy constructor
	} // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, particleVec);
    }

    // content of allParticleVec may need to be printed, so do not clear it. 
    //if (mpiRank == 0) releaseGatheredParticle();

    ///////////////////////////////////////////////////////////////////////////////
    // partition peri-points and send to each process
    if (mpiRank == 0) { // process 0
      setGrid(Rectangle(allContainer.getMinCorner().getX() - point_interval*0.2,
			allContainer.getMinCorner().getY() - point_interval*0.2,
			allContainer.getMinCorner().getZ() - point_interval*0.2,
			allContainer.getMaxCorner().getX() + point_interval*0.2,
			allContainer.getMaxCorner().getY() + point_interval*0.2,
			allContainer.getMaxCorner().getZ() + point_interval*0.2 ));
    
      Vec v1 = grid.getMinCorner();
      Vec v2 = grid.getMaxCorner();
      Vec vspan = v2 - v1;
      boost::mpi::request *reqs = new boost::mpi::request [mpiSize - 1];
      std::vector<periDynamics::PeriParticle*> tmpPeriParticleVec;
      for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
	tmpPeriParticleVec.clear(); // do not release memory!
	int ndim = 3;
	int coords[3];
	MPI_Cart_coords(cartComm, iRank, ndim, coords);
	Rectangle container(v1.getX() + vspan.getX() / mpiProcX * coords[0],
			    v1.getY() + vspan.getY() / mpiProcY * coords[1],
			    v1.getZ() + vspan.getZ() / mpiProcZ * coords[2],
			    v1.getX() + vspan.getX() / mpiProcX * (coords[0] + 1),
			    v1.getY() + vspan.getY() / mpiProcY * (coords[1] + 1),
			    v1.getZ() + vspan.getZ() / mpiProcZ * (coords[2] + 1));
	findPeriParticleInRectangle(container, allPeriParticleVec, tmpPeriParticleVec);
	if (iRank != 0)
	  reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, tmpPeriParticleVec); // non-blocking send
	if (iRank == 0) {
	  periParticleVec.resize(tmpPeriParticleVec.size());
	  for (auto i = 0u; i < periParticleVec.size(); ++i)
	    periParticleVec[i] = new periDynamics::PeriParticle(*tmpPeriParticleVec[i]); // default synthesized copy constructor
	} // now particleVec do not share memeory with allParticleVec
      }
      boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
      delete [] reqs;

    } else { // other processes except 0
      boostWorld.recv(0, mpiTag, periParticleVec);
    }

    // broadcast necessary info
    broadcast(boostWorld, gradation, 0);
    broadcast(boostWorld, boundaryVec, 0);
    broadcast(boostWorld, allContainer, 0);
    broadcast(boostWorld, grid, 0);
    broadcast(boostWorld, point_interval, 0);
    broadcast(boostWorld, maxHorizonSize, 0);
  } // scatterDEMPeriParticle



  bool Assembly::isBdryProcess() {
    return (mpiCoords[0] == 0 || mpiCoords[0] == mpiProcX - 1 ||
	    mpiCoords[1] == 0 || mpiCoords[1] == mpiProcY - 1 ||
	    mpiCoords[2] == 0 || mpiCoords[2] == mpiProcZ - 1);
  }


  void Assembly::commuParticle() 
  {
    // determine container of each process
    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;
    container = Rectangle(v1.getX() + vspan.getX() / mpiProcX * mpiCoords[0],
			  v1.getY() + vspan.getY() / mpiProcY * mpiCoords[1],
			  v1.getZ() + vspan.getZ() / mpiProcZ * mpiCoords[2],
			  v1.getX() + vspan.getX() / mpiProcX * (mpiCoords[0] + 1),
			  v1.getY() + vspan.getY() / mpiProcY * (mpiCoords[1] + 1),
			  v1.getZ() + vspan.getZ() / mpiProcZ * (mpiCoords[2] + 1));

    // find neighboring blocks
    rankX1 = -1; rankX2 = -1; rankY1 = -1; rankY2 = -1; rankZ1 = -1; rankZ2 = -1;
    rankX1Y1 = -1; rankX1Y2 = -1; rankX1Z1 = -1; rankX1Z2 = -1; 
    rankX2Y1 = -1; rankX2Y2 = -1; rankX2Z1 = -1; rankX2Z2 = -1; 
    rankY1Z1 = -1; rankY1Z2 = -1; rankY2Z1 = -1; rankY2Z2 = -1; 
    rankX1Y1Z1 = -1; rankX1Y1Z2 = -1; rankX1Y2Z1 = -1; rankX1Y2Z2 = -1; 
    rankX2Y1Z1 = -1; rankX2Y1Z2 = -1; rankX2Y2Z1 = -1; rankX2Y2Z2 = -1;
    // x1: -x direction
    int neighborCoords[3] = {mpiCoords[0], mpiCoords[1], mpiCoords[2]};
    --neighborCoords[0];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1);
    // x2: +x direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2);
    // y1: -y direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1);
    // y2: +y direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2);
    // z1: -z direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankZ1);
    // z2: +z direction
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankZ2);
    // x1y1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1);
    // x1y2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2);
    // x1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z1);
    // x1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z2);
    // x2y1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1);
    // x2y2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2);
    // x2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z1);
    // x2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z2);
    // y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z1);
    // y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z2);
    // y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z1);
    // y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z2);
    // x1y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z1);
    // x1y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z2);
    // x1y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z1);
    // x1y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    --neighborCoords[0]; ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z2);
    // x2y1z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z1);
    // x2y1z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; --neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z2);
    // x2y2z1
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1]; --neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z1);
    // x2y2z2
    neighborCoords[0] = mpiCoords[0];
    neighborCoords[1] = mpiCoords[1];
    neighborCoords[2] = mpiCoords[2];
    ++neighborCoords[0]; ++neighborCoords[1]; ++neighborCoords[2];
    MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z2);

    // if found, communicate with neighboring blocks
    ParticlePArray particleX1, particleX2;
    ParticlePArray particleY1, particleY2;
    ParticlePArray particleZ1, particleZ2;
    ParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2; 
    ParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2; 
    ParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2; 
    ParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2; 
    ParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
    v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
    v2 = container.getMaxCorner();   
    //debugInf << "rank=" << mpiRank << ' ' << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '  << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
    REAL cellSize = gradation.getPtclMaxRadius() * 2;
    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX(), v1.getY(), v1.getZ(), 
			    v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1, particleVec, particleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			    v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2, particleVec, particleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY(), v1.getZ(), 
			    v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerY1, particleVec, particleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			    v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerY2, particleVec, particleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ(),
			    v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerZ1, particleVec, particleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			    v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerZ2, particleVec, particleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX(), v1.getY(), v1.getZ(),
			      v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			      v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX(), v1.getY(), v1.getZ(),
			      v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			      v1.getX() + cellSize, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			      v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
			      v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			      v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
			      v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY(), v1.getZ(),
			      v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			      v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			      v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
			      v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX(), v1.getY(), v1.getZ(),
				v1.getX() + cellSize, v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
				v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
				v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
				v1.getX() + cellSize, v2.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
				v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
				v2.getX(), v1.getY() + cellSize, v2.getZ());
      findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
				v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX() - cellSize, v2.getY() - cellSize, v2.getZ() - cellSize,
				v2.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  particleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
    }

    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // merge: particles inside container (at front) + particles from neighoring blocks (at end)
    recvParticleVec.clear();
    // 6 surfaces
    if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
    if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
    if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
    if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
    if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
    if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

    mergeParticleVec.clear();
    mergeParticleVec = particleVec; // duplicate pointers, pointing to the same memory
    mergeParticleVec.insert(mergeParticleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

    /*
      ParticlePArray testParticleVec;
      testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
      debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4) << mpiRank 
      << " ptclNum=" << std::setw(4) << particleVec.size() 
      << " surface="
      << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
      << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
      << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()  
      << " recv="
      << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
      << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
      << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size() 
      << " rNum="    
      << std::setw(4) << recvParticleVec.size() << ": ";   

      for (ParticlePArray::const_iterator it = testParticleVec.begin(); it != testParticleVec.end();++it)
      debugInf << (*it)->getId() << ' ';
      debugInf << std::endl;
      testParticleVec.clear();
    */
  }


  // before the transfer to peri-points, their bondVec should be cleared
  void Assembly::commuPeriParticle() 
  {
    // determine container of each process
    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;
    container = Rectangle(v1.getX() + vspan.getX() / mpiProcX * mpiCoords[0],
			  v1.getY() + vspan.getY() / mpiProcY * mpiCoords[1],
			  v1.getZ() + vspan.getZ() / mpiProcZ * mpiCoords[2],
			  v1.getX() + vspan.getX() / mpiProcX * (mpiCoords[0] + 1),
			  v1.getY() + vspan.getY() / mpiProcY * (mpiCoords[1] + 1),
			  v1.getZ() + vspan.getZ() / mpiProcZ * (mpiCoords[2] + 1));

    // if found, communicate with neighboring blocks
    std::vector<periDynamics::PeriParticle*> periParticleX1, periParticleX2;
    std::vector<periDynamics::PeriParticle*> periParticleY1, periParticleY2;
    std::vector<periDynamics::PeriParticle*> periParticleZ1, periParticleZ2;
    std::vector<periDynamics::PeriParticle*> periParticleX1Y1, periParticleX1Y2, periParticleX1Z1, periParticleX1Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX2Y1, periParticleX2Y2, periParticleX2Z1, periParticleX2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleY1Z1, periParticleY1Z2, periParticleY2Z1, periParticleY2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX1Y1Z1, periParticleX1Y1Z2, periParticleX1Y2Z1, periParticleX1Y2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX2Y1Z1, periParticleX2Y1Z2, periParticleX2Y2Z1, periParticleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
    v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
    v2 = container.getMaxCorner();   
    //debugInf << "rank=" << mpiRank << ' ' << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '  << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
    REAL cellSize = std::max(4*maxHorizonSize, gradation.getPtclMaxRadius()+3*point_interval);	// constructNeighbor() is based on 2*horizonSize, 
    // refer to PeriParitcle.calcAcceleration(), in this function the deformationGradient and sigma of the other peri-points in the bond are needed,
    // this means that the deformationGradient, sigma and Kinv of the other peri-point should also be calculated even if it is in recvParticleVec, thus
    // we need to transfer 2*cellSize peri-points, the peri-points in the outer cell are used to calculate the deformationGradient, sigma and Kinv
    // of the peri-points in inner cell
    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX(), v1.getY(), v1.getZ(), 
			    v1.getX() + cellSize, v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX1, periParticleVec, periParticleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  periParticleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rperiParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			    v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2, periParticleVec, periParticleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  periParticleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rperiParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY(), v1.getZ(), 
			    v2.getX(), v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerY1, periParticleVec, periParticleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  periParticleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rperiParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			    v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerY2, periParticleVec, periParticleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  periParticleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rperiParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ(),
			    v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerZ1, periParticleVec, periParticleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  periParticleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rperiParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			    v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerZ2, periParticleVec, periParticleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  periParticleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rperiParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX(), v1.getY(), v1.getZ(),
			      v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerX1Y1, periParticleVec, periParticleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  periParticleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rperiParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			      v1.getX() + cellSize, v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX1Y2, periParticleVec, periParticleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  periParticleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rperiParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX(), v1.getY(), v1.getZ(),
			      v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX1Z1, periParticleVec, periParticleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  periParticleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rperiParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			      v1.getX() + cellSize, v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX1Z2, periParticleVec, periParticleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  periParticleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rperiParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			      v2.getX(), v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerX2Y1, periParticleVec, periParticleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  periParticleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rperiParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
			      v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2Y2, periParticleVec, periParticleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  periParticleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rperiParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
			      v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX2Z1, periParticleVec, periParticleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  periParticleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rperiParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
			      v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2Z2, periParticleVec, periParticleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  periParticleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rperiParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY(), v1.getZ(),
			      v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerY1Z1, periParticleVec, periParticleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  periParticleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rperiParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
			      v2.getX(), v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerY1Z2, periParticleVec, periParticleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  periParticleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rperiParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
			      v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerY2Z1, periParticleVec, periParticleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  periParticleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rperiParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
			      v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerY2Z2, periParticleVec, periParticleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  periParticleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rperiParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX(), v1.getY(), v1.getZ(),
				v1.getX() + cellSize, v1.getY() + cellSize, v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX1Y1Z1, periParticleVec, periParticleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  periParticleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rperiParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX(), v1.getY(), v2.getZ() - cellSize,
				v1.getX() + cellSize, v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerX1Y1Z2, periParticleVec, periParticleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  periParticleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rperiParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX(), v2.getY() - cellSize, v1.getZ(),
				v1.getX() + cellSize, v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX1Y2Z1, periParticleVec, periParticleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  periParticleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rperiParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX(), v2.getY() - cellSize, v2.getZ() - cellSize,
				v1.getX() + cellSize, v2.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerX1Y2Z2, periParticleVec, periParticleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  periParticleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rperiParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX() - cellSize, v1.getY(), v1.getZ(),
				v2.getX(), v1.getY() + cellSize, v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX2Y1Z1, periParticleVec, periParticleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  periParticleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rperiParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX() - cellSize, v1.getY(), v2.getZ() - cellSize,
				v2.getX(), v1.getY() + cellSize, v2.getZ());
      findPeriParticleInRectangle(containerX2Y1Z2, periParticleVec, periParticleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  periParticleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rperiParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX() - cellSize, v2.getY() - cellSize, v1.getZ(),
				v2.getX(), v2.getY(), v1.getZ() + cellSize);
      findPeriParticleInRectangle(containerX2Y2Z1, periParticleVec, periParticleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  periParticleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rperiParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX() - cellSize, v2.getY() - cellSize, v2.getZ() - cellSize,
				v2.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2Y2Z2, periParticleVec, periParticleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  periParticleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rperiParticleX2Y2Z2);
    }

    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // merge: periParticles inside container (at front) + periParticles from neighoring blocks (at end)
    recvPeriParticleVec.clear();
    // 6 surfaces
    if (rankX1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1.begin(), rperiParticleX1.end());
    if (rankX2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2.begin(), rperiParticleX2.end());
    if (rankY1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1.begin(), rperiParticleY1.end());
    if (rankY2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2.begin(), rperiParticleY2.end());
    if (rankZ1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleZ1.begin(), rperiParticleZ1.end());
    if (rankZ2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleZ2.begin(), rperiParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1.begin(), rperiParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2.begin(), rperiParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Z1.begin(), rperiParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Z2.begin(), rperiParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1.begin(), rperiParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2.begin(), rperiParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Z1.begin(), rperiParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Z2.begin(), rperiParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1Z1.begin(), rperiParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1Z2.begin(), rperiParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2Z1.begin(), rperiParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2Z2.begin(), rperiParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1Z1.begin(), rperiParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1Z2.begin(), rperiParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2Z1.begin(), rperiParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2Z2.begin(), rperiParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1Z1.begin(), rperiParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1Z2.begin(), rperiParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2Z1.begin(), rperiParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2Z2.begin(), rperiParticleX2Y2Z2.end());

    mergePeriParticleVec.clear();
    mergePeriParticleVec = periParticleVec; // duplicate pointers, pointing to the same memory
    mergePeriParticleVec.insert(mergePeriParticleVec.end(), recvPeriParticleVec.begin(), recvPeriParticleVec.end());

    /*
      ParticlePArray testParticleVec;
      testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
      testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
      debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4) << mpiRank 
      << " ptclNum=" << std::setw(4) << particleVec.size() 
      << " surface="
      << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
      << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
      << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()  
      << " recv="
      << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
      << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
      << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size() 
      << " rNum="    
      << std::setw(4) << recvParticleVec.size() << ": ";   

      for (ParticlePArray::const_iterator it = testParticleVec.begin(); it != testParticleVec.end();++it)
      debugInf << (*it)->getId() << ' ';
      debugInf << std::endl;
      testParticleVec.clear();
    */
  }

  void Assembly::releaseRecvParticle() {
    // release memory of received particles
    /*
    for (ParticlePArray::iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      delete (*it);
    */
    recvParticleVec.clear();
    // 6 surfaces
    rParticleX1.clear();
    rParticleX2.clear();
    rParticleY1.clear();
    rParticleY2.clear();
    rParticleZ1.clear();
    rParticleZ2.clear();
    // 12 edges
    rParticleX1Y1.clear();
    rParticleX1Y2.clear();
    rParticleX1Z1.clear();
    rParticleX1Z2.clear();
    rParticleX2Y1.clear();
    rParticleX2Y2.clear();
    rParticleX2Z1.clear();
    rParticleX2Z2.clear();
    rParticleY1Z1.clear();
    rParticleY1Z2.clear();
    rParticleY2Z1.clear();
    rParticleY2Z2.clear();
    // 8 vertices
    rParticleX1Y1Z1.clear();
    rParticleX1Y1Z2.clear();
    rParticleX1Y2Z1.clear();
    rParticleX1Y2Z2.clear();
    rParticleX2Y1Z1.clear();
    rParticleX2Y1Z2.clear();
    rParticleX2Y2Z1.clear();
    rParticleX2Y2Z2.clear();
  }

  void Assembly::releaseRecvPeriParticle() {
//    // delete those periBonds between recvPeriParticle
//    for (std::vector<periDynamics::PeriParticle*>::iterator it = periParticleVec.begin(); it != periParticleVec.end(); ++it){
//	(*it)->eraseRecvPeriBonds();
//    }
    for(std::vector<periDynamics::PeriBond*>::iterator pt=recvPeriBondVec.begin(); pt!=recvPeriBondVec.end(); pt++){
	delete (*pt);
    }
    recvPeriBondVec.clear();
    std::vector<periDynamics::PeriBond*>().swap(recvPeriBondVec); // actual memory release

    // release memory of received particles
    for (std::vector<periDynamics::PeriParticle*>::iterator it = recvPeriParticleVec.begin(); it != recvPeriParticleVec.end(); ++it){
      delete (*it);
    }
    recvPeriParticleVec.clear();
    std::vector<periDynamics::PeriParticle*>().swap(recvPeriParticleVec); // actual memory release
    // 6 surfaces
    rperiParticleX1.clear();
    rperiParticleX2.clear();
    rperiParticleY1.clear();
    rperiParticleY2.clear();
    rperiParticleZ1.clear();
    rperiParticleZ2.clear();
    // 12 edges
    rperiParticleX1Y1.clear();
    rperiParticleX1Y2.clear();
    rperiParticleX1Z1.clear();
    rperiParticleX1Z2.clear();
    rperiParticleX2Y1.clear();
    rperiParticleX2Y2.clear();
    rperiParticleX2Z1.clear();
    rperiParticleX2Z2.clear();
    rperiParticleY1Z1.clear();
    rperiParticleY1Z2.clear();
    rperiParticleY2Z1.clear();
    rperiParticleY2Z2.clear();
    // 8 vertices
    rperiParticleX1Y1Z1.clear();
    rperiParticleX1Y1Z2.clear();
    rperiParticleX1Y2Z1.clear();
    rperiParticleX1Y2Z2.clear();
    rperiParticleX2Y1Z1.clear();
    rperiParticleX2Y1Z2.clear();
    rperiParticleX2Y2Z1.clear();
    rperiParticleX2Y2Z2.clear();
  }

  void Assembly::releasePeriBondVec(){
    for(std::vector<periDynamics::PeriBond*>::iterator pt=periBondVec.begin(); pt!=periBondVec.end(); pt++){
	delete (*pt);
    }
    periBondVec.clear();
    std::vector<periDynamics::PeriBond*>().swap(periBondVec); // actual memory release
  }

  void Assembly::updateGrid() {
    updateGridMinX();
    updateGridMaxX();
    updateGridMinY();
    updateGridMaxY();
    updateGridMinZ();
    updateGridMaxZ();
  }


  void Assembly::updateGridMinX() {
    REAL pMinX = getPtclMinX(particleVec);
    REAL minX = 0;
    MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  
    setGrid(Rectangle(minX - gradation.getPtclMaxRadius(),
		      grid.getMinCorner().getY(),
		      grid.getMinCorner().getZ(),
		      grid.getMaxCorner().getX(),
		      grid.getMaxCorner().getY(),
		      grid.getMaxCorner().getZ() ));
  }


  void Assembly::updateGridMaxX() {
    REAL pMaxX = getPtclMaxX(particleVec);
    REAL maxX = 0;
    MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
  
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      grid.getMinCorner().getY(),
		      grid.getMinCorner().getZ(),
		      maxX + gradation.getPtclMaxRadius(),
		      grid.getMaxCorner().getY(),
		      grid.getMaxCorner().getZ() ));
  }


  void Assembly::updateGridMinY() {
    REAL pMinY = getPtclMinY(particleVec);
    REAL minY = 0;
    MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      minY - gradation.getPtclMaxRadius(),
		      grid.getMinCorner().getZ(),
		      grid.getMaxCorner().getX(),
		      grid.getMaxCorner().getY(),
		      grid.getMaxCorner().getZ() ));
  }


  void Assembly::updateGridMaxY() {
    REAL pMaxY = getPtclMaxY(particleVec);
    REAL maxY = 0;
    MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
  
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      grid.getMinCorner().getY(),
		      grid.getMinCorner().getZ(),
		      grid.getMaxCorner().getX(),
		      maxY + gradation.getPtclMaxRadius(),
		      grid.getMaxCorner().getZ() ));
  }


  void Assembly::updateGridMinZ() {
    REAL pMinZ = getPtclMinZ(particleVec);
    REAL minZ = 0;
    MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      grid.getMinCorner().getY(),
		      minZ - gradation.getPtclMaxRadius(),
		      grid.getMaxCorner().getX(),
		      grid.getMaxCorner().getY(),
		      grid.getMaxCorner().getZ() ));
  }


  void Assembly::updateGridMaxZ() {
    // update compute grids adaptively due to particle motion
    REAL pMaxZ = getPtclMaxZ(particleVec);
    REAL maxZ = 0;
    MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

    // no need to broadcast grid as it is updated in each process
    setGrid(Rectangle(grid.getMinCorner().getX(),
		      grid.getMinCorner().getY(),
		      grid.getMinCorner().getZ(),
		      grid.getMaxCorner().getX(),
		      grid.getMaxCorner().getY(),
		      maxZ + gradation.getPtclMaxRadius() ));
  }

  void Assembly::updatePeriGrid(){
    std::vector<periDynamics::PeriParticle*>::const_iterator pit = periParticleVec.begin();
    REAL pMaxX=(*pit)->getCurrPosition().getX(); 
    REAL pMinX=(*pit)->getCurrPosition().getX(); 
    REAL pMaxY=(*pit)->getCurrPosition().getY(); 
    REAL pMinY=(*pit)->getCurrPosition().getY(); 
    REAL pMaxZ=(*pit)->getCurrPosition().getZ(); 
    REAL pMinZ=(*pit)->getCurrPosition().getZ();

    REAL tmpx, tmpy, tmpz;
    dem::Vec tmp_xyz;
    for(pit=periParticleVec.begin(); pit!=periParticleVec.end(); pit++){
	tmp_xyz = (*pit)->getCurrPosition();
	tmpx=tmp_xyz.getX(); tmpy=tmp_xyz.getY(); tmpz=tmp_xyz.getZ();
	if( pMaxX<tmpx )
	    pMaxX=tmpx;
   	if( pMinX>tmpx )
	    pMinX=tmpx;
	if( pMaxY<tmpy )
	    pMaxY=tmpy;
   	if( pMinY>tmpy )
	    pMinY=tmpy;
	if( pMaxZ<tmpz )
	    pMaxZ=tmpz;
   	if( pMinZ>tmpz )
	    pMinZ=tmpz;  
    }
    REAL maxX=0; REAL maxY=0; REAL maxZ=0;
    REAL minX=0; REAL minY=0; REAL minZ=0;
    MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

    // no need to broadcast grid as it is updated in each process
    setGrid(Rectangle(minX - point_interval*0.2,
		      minY - point_interval*0.2,
		      minZ - point_interval*0.2,
		      maxX + point_interval*0.2,
		      maxY + point_interval*0.2,
		      maxZ + point_interval*0.2 ));
  }

  void Assembly::migrateParticle() 
  {
    Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
    REAL segX = vspan.getX() / mpiProcX;
    REAL segY = vspan.getY() / mpiProcY;
    REAL segZ = vspan.getZ() / mpiProcZ;
    Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
    Vec v2 = container.getMaxCorner();  

    // if a neighbor exists, transfer particles crossing the boundary in between.
    ParticlePArray particleX1, particleX2;
    ParticlePArray particleY1, particleY2;
    ParticlePArray particleZ1, particleZ2;
    ParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2; 
    ParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2; 
    ParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2; 
    ParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2; 
    ParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX() - segX, v1.getY(), v1.getZ(), 
			    v1.getX(), v2.getY(), v2.getZ());
      findParticleInRectangle(containerX1, particleVec, particleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  particleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX(), v1.getY(), v1.getZ(),
			    v2.getX() + segX, v2.getY(), v2.getZ());
      findParticleInRectangle(containerX2, particleVec, particleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  particleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY() - segY, v1.getZ(), 
			    v2.getX(), v1.getY(), v2.getZ());
      findParticleInRectangle(containerY1, particleVec, particleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  particleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY(), v1.getZ(),
			    v2.getX(), v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerY2, particleVec, particleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  particleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ() - segZ,
			    v2.getX(), v2.getY(), v1.getZ());
      findParticleInRectangle(containerZ1, particleVec, particleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  particleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ(),
			    v2.getX(), v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerZ2, particleVec, particleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  particleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX() - segX, v1.getY() - segY, v1.getZ(),
			      v1.getX(), v1.getY(), v2.getZ());
      findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  particleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX() - segX, v2.getY(), v1.getZ(),
			      v1.getX(), v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  particleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX() - segX, v1.getY(), v1.getZ() -segZ,
			      v1.getX(), v2.getY(), v1.getZ());
      findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  particleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX() - segX, v1.getY(), v2.getZ(),
			      v1.getX(), v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  particleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX(), v1.getY() - segY, v1.getZ(),
			      v2.getX() + segX, v1.getY(), v2.getZ());
      findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  particleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX(), v2.getY(), v1.getZ(),
			      v2.getX() + segX, v2.getY() + segY, v2.getZ());
      findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  particleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX(), v1.getY(), v1.getZ() - segZ,
			      v2.getX() + segX, v2.getY(), v1.getZ());
      findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  particleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX(), v1.getY(), v2.getZ(),
			      v2.getX() + segX, v2.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  particleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY() - segY, v1.getZ() - segZ,
			      v2.getX(), v1.getY(), v1.getZ());
      findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  particleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY() - segY, v2.getZ(),
			      v2.getX(), v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  particleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY(), v1.getZ() - segZ,
			      v2.getX(), v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  particleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY(), v2.getZ(),
			      v2.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  particleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX() - segX, v1.getY() - segY, v1.getZ() - segZ,
				v1.getX(), v1.getY(), v1.getZ());
      findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  particleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX() - segX, v1.getY() - segY, v2.getZ(),
				v1.getX(), v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  particleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX() - segX, v2.getY(), v1.getZ() - segZ,
				v1.getX(), v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  particleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX() - segX, v2.getY(), v2.getZ(),
				v1.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  particleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX(), v1.getY() - segY, v1.getZ() - segZ,
				v2.getX() + segX, v1.getY(), v1.getZ());
      findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  particleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX(), v1.getY() - segY, v2.getZ(),
				v2.getX() + segX, v1.getY(), v2.getZ() + segZ);
      findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  particleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX(), v2.getY(), v1.getZ() - segZ,
				v2.getX() + segX, v2.getY() + segY, v1.getZ());
      findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  particleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX(), v2.getY(), v2.getZ(),
				v2.getX() + segX, v2.getY() + segY, v2.getZ() + segZ);
      findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  particleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
    }
    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // delete outgoing particles
    removeParticleOutRectangle();

    // add incoming particles
    recvParticleVec.clear(); // new use of recvParticleVec
    // 6 surfaces
    if (rankX1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(), rParticleX1.end());
    if (rankX2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(), rParticleX2.end());
    if (rankY1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(), rParticleY1.end());
    if (rankY2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(), rParticleY2.end());
    if (rankZ1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(), rParticleZ1.end());
    if (rankZ2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(), rParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(), rParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(), rParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(), rParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(), rParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(), rParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(), rParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(), rParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(), rParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(), rParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(), rParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(), rParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(), rParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(), rParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(), rParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(), rParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(), rParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(), rParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(), rParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(), rParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(), rParticleX2Y2Z2.end());

    particleVec.insert(particleVec.end(), recvParticleVec.begin(), recvParticleVec.end());

    /*
      if (recvParticleVec.size() > 0) {    
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank 
      << "   added=";
      for (ParticlePArray::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << " now " << particleVec.size() << ": ";
      for (ParticlePArray::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

    // do not release memory of received particles because they are part of and managed by particleVec
    // 6 surfaces
    rParticleX1.clear();
    rParticleX2.clear();
    rParticleY1.clear();
    rParticleY2.clear();
    rParticleZ1.clear();
    rParticleZ2.clear();
    // 12 edges
    rParticleX1Y1.clear();
    rParticleX1Y2.clear();
    rParticleX1Z1.clear();
    rParticleX1Z2.clear();
    rParticleX2Y1.clear();
    rParticleX2Y2.clear();
    rParticleX2Z1.clear();
    rParticleX2Z2.clear();
    rParticleY1Z1.clear();
    rParticleY1Z2.clear();
    rParticleY2Z1.clear();
    rParticleY2Z2.clear();
    // 8 vertices
    rParticleX1Y1Z1.clear();
    rParticleX1Y1Z2.clear();
    rParticleX1Y2Z1.clear();
    rParticleX1Y2Z2.clear();
    rParticleX2Y1Z1.clear();
    rParticleX2Y1Z2.clear();
    rParticleX2Y2Z1.clear();
    rParticleX2Y2Z2.clear();

    recvParticleVec.clear();
  }


  void Assembly::migratePeriParticle() 
  {
    Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
    REAL segX = vspan.getX() / mpiProcX;
    REAL segY = vspan.getY() / mpiProcY;
    REAL segZ = vspan.getZ() / mpiProcZ;
    Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
    Vec v2 = container.getMaxCorner();  

    // if a neighbor exists, transfer particles crossing the boundary in between.
    std::vector<periDynamics::PeriParticle*> periParticleX1, periParticleX2;
    std::vector<periDynamics::PeriParticle*> periParticleY1, periParticleY2;
    std::vector<periDynamics::PeriParticle*> periParticleZ1, periParticleZ2;
    std::vector<periDynamics::PeriParticle*> periParticleX1Y1, periParticleX1Y2, periParticleX1Z1, periParticleX1Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX2Y1, periParticleX2Y2, periParticleX2Z1, periParticleX2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleY1Z1, periParticleY1Z2, periParticleY2Z1, periParticleY2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX1Y1Z1, periParticleX1Y1Z2, periParticleX1Y2Z1, periParticleX1Y2Z2; 
    std::vector<periDynamics::PeriParticle*> periParticleX2Y1Z1, periParticleX2Y1Z2, periParticleX2Y2Z1, periParticleX2Y2Z2; 
    boost::mpi::request reqX1[2], reqX2[2];
    boost::mpi::request reqY1[2], reqY2[2];
    boost::mpi::request reqZ1[2], reqZ2[2];
    boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
    boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
    boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
    boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
    boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

    // 6 surfaces
    if (rankX1 >= 0) { // surface x1
      Rectangle containerX1(v1.getX() - segX, v1.getY(), v1.getZ(), 
			    v1.getX(), v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX1, periParticleVec, periParticleX1);
      reqX1[0] = boostWorld.isend(rankX1, mpiTag,  periParticleX1);
      reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rperiParticleX1);
    }
    if (rankX2 >= 0) { // surface x2
      Rectangle containerX2(v2.getX(), v1.getY(), v1.getZ(),
			    v2.getX() + segX, v2.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2, periParticleVec, periParticleX2);
      reqX2[0] = boostWorld.isend(rankX2, mpiTag,  periParticleX2);
      reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rperiParticleX2);
    }
    if (rankY1 >= 0) {  // surface y1
      Rectangle containerY1(v1.getX(), v1.getY() - segY, v1.getZ(), 
			    v2.getX(), v1.getY(), v2.getZ());
      findPeriParticleInRectangle(containerY1, periParticleVec, periParticleY1);
      reqY1[0] = boostWorld.isend(rankY1, mpiTag,  periParticleY1);
      reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rperiParticleY1);
    }
    if (rankY2 >= 0) {  // surface y2
      Rectangle containerY2(v1.getX(), v2.getY(), v1.getZ(),
			    v2.getX(), v2.getY() + segY, v2.getZ());
      findPeriParticleInRectangle(containerY2, periParticleVec, periParticleY2);
      reqY2[0] = boostWorld.isend(rankY2, mpiTag,  periParticleY2);
      reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rperiParticleY2);
    }
    if (rankZ1 >= 0) {  // surface z1
      Rectangle containerZ1(v1.getX(), v1.getY(), v1.getZ() - segZ,
			    v2.getX(), v2.getY(), v1.getZ());
      findPeriParticleInRectangle(containerZ1, periParticleVec, periParticleZ1);
      reqZ1[0] = boostWorld.isend(rankZ1, mpiTag,  periParticleZ1);
      reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rperiParticleZ1);
    }
    if (rankZ2 >= 0) {  // surface z2
      Rectangle containerZ2(v1.getX(), v1.getY(), v2.getZ(),
			    v2.getX(), v2.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerZ2, periParticleVec, periParticleZ2);
      reqZ2[0] = boostWorld.isend(rankZ2, mpiTag,  periParticleZ2);
      reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rperiParticleZ2);
    }
    // 12 edges
    if (rankX1Y1 >= 0) { // edge x1y1
      Rectangle containerX1Y1(v1.getX() - segX, v1.getY() - segY, v1.getZ(),
			      v1.getX(), v1.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX1Y1, periParticleVec, periParticleX1Y1);
      reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag,  periParticleX1Y1);
      reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rperiParticleX1Y1);
    }
    if (rankX1Y2 >= 0) { // edge x1y2
      Rectangle containerX1Y2(v1.getX() - segX, v2.getY(), v1.getZ(),
			      v1.getX(), v2.getY() + segY, v2.getZ());
      findPeriParticleInRectangle(containerX1Y2, periParticleVec, periParticleX1Y2);
      reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag,  periParticleX1Y2);
      reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rperiParticleX1Y2);
    }
    if (rankX1Z1 >= 0) { // edge x1z1
      Rectangle containerX1Z1(v1.getX() - segX, v1.getY(), v1.getZ() -segZ,
			      v1.getX(), v2.getY(), v1.getZ());
      findPeriParticleInRectangle(containerX1Z1, periParticleVec, periParticleX1Z1);
      reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag,  periParticleX1Z1);
      reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rperiParticleX1Z1);
    }
    if (rankX1Z2 >= 0) { // edge x1z2
      Rectangle containerX1Z2(v1.getX() - segX, v1.getY(), v2.getZ(),
			      v1.getX(), v2.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX1Z2, periParticleVec, periParticleX1Z2);
      reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag,  periParticleX1Z2);
      reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rperiParticleX1Z2);
    }
    if (rankX2Y1 >= 0) { // edge x2y1
      Rectangle containerX2Y1(v2.getX(), v1.getY() - segY, v1.getZ(),
			      v2.getX() + segX, v1.getY(), v2.getZ());
      findPeriParticleInRectangle(containerX2Y1, periParticleVec, periParticleX2Y1);
      reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag,  periParticleX2Y1);
      reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rperiParticleX2Y1);
    }
    if (rankX2Y2 >= 0) { // edge x2y2
      Rectangle containerX2Y2(v2.getX(), v2.getY(), v1.getZ(),
			      v2.getX() + segX, v2.getY() + segY, v2.getZ());
      findPeriParticleInRectangle(containerX2Y2, periParticleVec, periParticleX2Y2);
      reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag,  periParticleX2Y2);
      reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rperiParticleX2Y2);
    }
    if (rankX2Z1 >= 0) { // edge x2z1
      Rectangle containerX2Z1(v2.getX(), v1.getY(), v1.getZ() - segZ,
			      v2.getX() + segX, v2.getY(), v1.getZ());
      findPeriParticleInRectangle(containerX2Z1, periParticleVec, periParticleX2Z1);
      reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag,  periParticleX2Z1);
      reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rperiParticleX2Z1);
    }
    if (rankX2Z2 >= 0) { // edge x2z2
      Rectangle containerX2Z2(v2.getX(), v1.getY(), v2.getZ(),
			      v2.getX() + segX, v2.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX2Z2, periParticleVec, periParticleX2Z2);
      reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag,  periParticleX2Z2);
      reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rperiParticleX2Z2);
    }
    if (rankY1Z1 >= 0) { // edge y1z1
      Rectangle containerY1Z1(v1.getX(), v1.getY() - segY, v1.getZ() - segZ,
			      v2.getX(), v1.getY(), v1.getZ());
      findPeriParticleInRectangle(containerY1Z1, periParticleVec, periParticleY1Z1);
      reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag,  periParticleY1Z1);
      reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rperiParticleY1Z1);
    }
    if (rankY1Z2 >= 0) { // edge y1z2
      Rectangle containerY1Z2(v1.getX(), v1.getY() - segY, v2.getZ(),
			      v2.getX(), v1.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerY1Z2, periParticleVec, periParticleY1Z2);
      reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag,  periParticleY1Z2);
      reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rperiParticleY1Z2);
    }
    if (rankY2Z1 >= 0) { // edge y2z1
      Rectangle containerY2Z1(v1.getX(), v2.getY(), v1.getZ() - segZ,
			      v2.getX(), v2.getY() + segY, v1.getZ());
      findPeriParticleInRectangle(containerY2Z1, periParticleVec, periParticleY2Z1);
      reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag,  periParticleY2Z1);
      reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rperiParticleY2Z1);
    }
    if (rankY2Z2 >= 0) { // edge y2z2
      Rectangle containerY2Z2(v1.getX(), v2.getY(), v2.getZ(),
			      v2.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findPeriParticleInRectangle(containerY2Z2, periParticleVec, periParticleY2Z2);
      reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag,  periParticleY2Z2);
      reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rperiParticleY2Z2);
    }
    // 8 vertices
    if (rankX1Y1Z1 >= 0) { // edge x1y1z1
      Rectangle containerX1Y1Z1(v1.getX() - segX, v1.getY() - segY, v1.getZ() - segZ,
				v1.getX(), v1.getY(), v1.getZ());
      findPeriParticleInRectangle(containerX1Y1Z1, periParticleVec, periParticleX1Y1Z1);
      reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag,  periParticleX1Y1Z1);
      reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rperiParticleX1Y1Z1);
    }
    if (rankX1Y1Z2 >= 0) { // edge x1y1z2
      Rectangle containerX1Y1Z2(v1.getX() - segX, v1.getY() - segY, v2.getZ(),
				v1.getX(), v1.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX1Y1Z2, periParticleVec, periParticleX1Y1Z2);
      reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag,  periParticleX1Y1Z2);
      reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rperiParticleX1Y1Z2);
    }
    if (rankX1Y2Z1 >= 0) { // edge x1y2z1
      Rectangle containerX1Y2Z1(v1.getX() - segX, v2.getY(), v1.getZ() - segZ,
				v1.getX(), v2.getY() + segY, v1.getZ());
      findPeriParticleInRectangle(containerX1Y2Z1, periParticleVec, periParticleX1Y2Z1);
      reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag,  periParticleX1Y2Z1);
      reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rperiParticleX1Y2Z1);
    }
    if (rankX1Y2Z2 >= 0) { // edge x1y2z2
      Rectangle containerX1Y2Z2(v1.getX() - segX, v2.getY(), v2.getZ(),
				v1.getX(), v2.getY() + segY, v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX1Y2Z2, periParticleVec, periParticleX1Y2Z2);
      reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag,  periParticleX1Y2Z2);
      reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rperiParticleX1Y2Z2);
    }
    if (rankX2Y1Z1 >= 0) { // edge x2y1z1
      Rectangle containerX2Y1Z1(v2.getX(), v1.getY() - segY, v1.getZ() - segZ,
				v2.getX() + segX, v1.getY(), v1.getZ());
      findPeriParticleInRectangle(containerX2Y1Z1, periParticleVec, periParticleX2Y1Z1);
      reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag,  periParticleX2Y1Z1);
      reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rperiParticleX2Y1Z1);
    }
    if (rankX2Y1Z2 >= 0) { // edge x2y1z2
      Rectangle containerX2Y1Z2(v2.getX(), v1.getY() - segY, v2.getZ(),
				v2.getX() + segX, v1.getY(), v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX2Y1Z2, periParticleVec, periParticleX2Y1Z2);
      reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag,  periParticleX2Y1Z2);
      reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rperiParticleX2Y1Z2);
    }
    if (rankX2Y2Z1 >= 0) { // edge x2y2z1
      Rectangle containerX2Y2Z1(v2.getX(), v2.getY(), v1.getZ() - segZ,
				v2.getX() + segX, v2.getY() + segY, v1.getZ());
      findPeriParticleInRectangle(containerX2Y2Z1, periParticleVec, periParticleX2Y2Z1);
      reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag,  periParticleX2Y2Z1);
      reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rperiParticleX2Y2Z1);
    }
    if (rankX2Y2Z2 >= 0) { // edge x2y2z2
      Rectangle containerX2Y2Z2(v2.getX(), v2.getY(), v2.getZ(),
				v2.getX() + segX, v2.getY() + segY, v2.getZ() + segZ);
      findPeriParticleInRectangle(containerX2Y2Z2, periParticleVec, periParticleX2Y2Z2);
      reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag,  periParticleX2Y2Z2);
      reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rperiParticleX2Y2Z2);
    }
    // 6 surfaces
    if (rankX1 >= 0) boost::mpi::wait_all(reqX1, reqX1 + 2);
    if (rankX2 >= 0) boost::mpi::wait_all(reqX2, reqX2 + 2);
    if (rankY1 >= 0) boost::mpi::wait_all(reqY1, reqY1 + 2);
    if (rankY2 >= 0) boost::mpi::wait_all(reqY2, reqY2 + 2);
    if (rankZ1 >= 0) boost::mpi::wait_all(reqZ1, reqZ1 + 2);
    if (rankZ2 >= 0) boost::mpi::wait_all(reqZ2, reqZ2 + 2);
    // 12 edges
    if (rankX1Y1 >= 0) boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
    if (rankX1Y2 >= 0) boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);  
    if (rankX1Z1 >= 0) boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
    if (rankX1Z2 >= 0) boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
    if (rankX2Y1 >= 0) boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
    if (rankX2Y2 >= 0) boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);  
    if (rankX2Z1 >= 0) boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
    if (rankX2Z2 >= 0) boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2); 
    if (rankY1Z1 >= 0) boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
    if (rankY1Z2 >= 0) boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
    if (rankY2Z1 >= 0) boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
    if (rankY2Z2 >= 0) boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2); 
    // 8 vertices
    if (rankX1Y1Z1 >= 0) boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
    if (rankX1Y1Z2 >= 0) boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
    if (rankX1Y2Z1 >= 0) boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
    if (rankX1Y2Z2 >= 0) boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
    if (rankX2Y1Z1 >= 0) boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
    if (rankX2Y1Z2 >= 0) boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
    if (rankX2Y2Z1 >= 0) boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
    if (rankX2Y2Z2 >= 0) boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);  

    // delete outgoing particles
    removePeriParticleOutRectangle();

    // add incoming particles
    recvPeriParticleVec.clear();
    // 6 surfaces
    if (rankX1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1.begin(), rperiParticleX1.end());
    if (rankX2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2.begin(), rperiParticleX2.end());
    if (rankY1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1.begin(), rperiParticleY1.end());
    if (rankY2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2.begin(), rperiParticleY2.end());
    if (rankZ1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleZ1.begin(), rperiParticleZ1.end());
    if (rankZ2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleZ2.begin(), rperiParticleZ2.end());
    // 12 edges
    if (rankX1Y1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1.begin(), rperiParticleX1Y1.end());
    if (rankX1Y2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2.begin(), rperiParticleX1Y2.end());
    if (rankX1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Z1.begin(), rperiParticleX1Z1.end());
    if (rankX1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Z2.begin(), rperiParticleX1Z2.end());
    if (rankX2Y1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1.begin(), rperiParticleX2Y1.end());
    if (rankX2Y2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2.begin(), rperiParticleX2Y2.end());
    if (rankX2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Z1.begin(), rperiParticleX2Z1.end());
    if (rankX2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Z2.begin(), rperiParticleX2Z2.end());
    if (rankY1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1Z1.begin(), rperiParticleY1Z1.end());
    if (rankY1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY1Z2.begin(), rperiParticleY1Z2.end());
    if (rankY2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2Z1.begin(), rperiParticleY2Z1.end());
    if (rankY2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleY2Z2.begin(), rperiParticleY2Z2.end());
    // 8 vertices
    if (rankX1Y1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1Z1.begin(), rperiParticleX1Y1Z1.end());
    if (rankX1Y1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y1Z2.begin(), rperiParticleX1Y1Z2.end());
    if (rankX1Y2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2Z1.begin(), rperiParticleX1Y2Z1.end());
    if (rankX1Y2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX1Y2Z2.begin(), rperiParticleX1Y2Z2.end());
    if (rankX2Y1Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1Z1.begin(), rperiParticleX2Y1Z1.end());
    if (rankX2Y1Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y1Z2.begin(), rperiParticleX2Y1Z2.end());
    if (rankX2Y2Z1 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2Z1.begin(), rperiParticleX2Y2Z1.end());
    if (rankX2Y2Z2 >= 0) recvPeriParticleVec.insert(recvPeriParticleVec.end(), rperiParticleX2Y2Z2.begin(), rperiParticleX2Y2Z2.end());

    for(std::vector<periDynamics::PeriParticle*>::iterator pit=recvPeriParticleVec.begin(); pit!=recvPeriParticleVec.end(); pit++)
	(*pit)->constructMatrixMember();

    periParticleVec.insert(periParticleVec.end(), recvPeriParticleVec.begin(), recvPeriParticleVec.end());

    /*
      if (recvParticleVec.size() > 0) {    
      debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2) << mpiRank 
      << "   added=";
      for (ParticlePArray::const_iterator it = recvParticleVec.begin(); it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << " now " << particleVec.size() << ": ";
      for (ParticlePArray::const_iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
      debugInf << std::endl;
      }
    */

    // do not release memory of received particles because they are part of and managed by particleVec
    // 6 surfaces
    rParticleX1.clear();
    rParticleX2.clear();
    rParticleY1.clear();
    rParticleY2.clear();
    rParticleZ1.clear();
    rParticleZ2.clear();
    // 12 edges
    rParticleX1Y1.clear();
    rParticleX1Y2.clear();
    rParticleX1Z1.clear();
    rParticleX1Z2.clear();
    rParticleX2Y1.clear();
    rParticleX2Y2.clear();
    rParticleX2Z1.clear();
    rParticleX2Z2.clear();
    rParticleY1Z1.clear();
    rParticleY1Z2.clear();
    rParticleY2Z1.clear();
    rParticleY2Z2.clear();
    // 8 vertices
    rParticleX1Y1Z1.clear();
    rParticleX1Y1Z2.clear();
    rParticleX1Y2Z1.clear();
    rParticleX1Y2Z2.clear();
    rParticleX2Y1Z1.clear();
    rParticleX2Y1Z2.clear();
    rParticleX2Y2Z1.clear();
    rParticleX2Y2Z2.clear();

    recvPeriParticleVec.clear();
  }

  void Assembly::gatherParticle() {
    // update allParticleVec: process 0 collects all updated particles from each other process  
    if (mpiRank != 0) {// each process except 0
      boostWorld.send(0, mpiTag, particleVec);
    }
    else { // process 0
      // allParticleVec is cleared before filling with new data
      releaseGatheredParticle();

      // duplicate particleVec so that it is not destroyed by allParticleVec in next iteration,
      // otherwise it causes memory error.
      ParticlePArray dupParticleVec(particleVec.size());
      for (std::size_t i = 0; i < dupParticleVec.size(); ++i)
	dupParticleVec[i] = std::make_shared<Particle>(*particleVec[i]);

      // fill allParticleVec with dupParticleVec and received particles
      allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(), dupParticleVec.end());

      ParticlePArray tmpParticleVec;
      long gatherRam = 0;
      for (int iRank = 1; iRank < mpiSize; ++iRank) {

	tmpParticleVec.clear();// do not destroy particles!
	boostWorld.recv(iRank, mpiTag, tmpParticleVec);
	allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(), tmpParticleVec.end());
	gatherRam += tmpParticleVec.size();

      }
      //debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = " << gatherRam * sizeof(Particle) << std::endl;
    }
  }


  void Assembly::releaseGatheredParticle() {
    // clear allParticleVec, avoid long time memory footprint.
    /*
    for (ParticlePArray::iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
      delete (*it);
    */
    allParticleVec.clear();
    ParticlePArray().swap(allParticleVec); // actual memory release
  }

  void Assembly::gatherPeriParticle() {
    // update allPeriParticleVec: process 0 collects all updated particles from each other process  
    for(std::vector<periDynamics::PeriParticle*>::iterator it=periParticleVec.begin(); it!=periParticleVec.end(); it++){
      (*it)->assignSigma();	// store Matrix sigma to each value
    }
    if (mpiRank != 0) {// each process except 0
      boostWorld.send(0, mpiTag, periParticleVec);
    }
    else { // process 0
      // allPeriParticleVec is cleared before filling with new data
      releaseGatheredPeriParticle();

      // duplicate PeriParticleVec so that it is not destroyed by allPeriParticleVec in next iteration,
      // otherwise it causes memory error.
      std::vector<periDynamics::PeriParticle*> dupPeriParticleVec(periParticleVec.size());
      for (std::size_t i = 0; i < dupPeriParticleVec.size(); ++i){
	dupPeriParticleVec[i] = new periDynamics::PeriParticle(*periParticleVec[i]);
	dupPeriParticleVec[i]->releaseBondVec();
      }

      // fill allParticleVec with dupParticleVec and received particles
      allPeriParticleVec.insert(allPeriParticleVec.end(), dupPeriParticleVec.begin(), dupPeriParticleVec.end());
      std::vector<periDynamics::PeriParticle*> tmpPeriParticleVec;
      long gatherRam = 0;
      for (int iRank = 1; iRank < mpiSize; ++iRank) {
	tmpPeriParticleVec.clear();// do not destroy particles!
	boostWorld.recv(iRank, mpiTag, tmpPeriParticleVec);
	allPeriParticleVec.insert(allPeriParticleVec.end(), tmpPeriParticleVec.begin(), tmpPeriParticleVec.end());
	gatherRam += tmpPeriParticleVec.size();

      }
      //debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = " << gatherRam * sizeof(Particle) << std::endl;
    }
  }

  void Assembly::releaseGatheredPeriParticle() {
    // clear allPeriParticleVec, avoid long time memory footprint.
    for (std::vector<periDynamics::PeriParticle*>::iterator it = allPeriParticleVec.begin(); it != allPeriParticleVec.end(); ++it)
      delete (*it);
    allPeriParticleVec.clear();
    std::vector<periDynamics::PeriParticle*>().swap(allPeriParticleVec); // actual memory release
  }

  void Assembly::gatherBdryContact() {
    if (isBdryProcess()) {
      if (mpiRank != 0)
	boostWorld.send(0, mpiTag, boundaryVec);
    }

    if (mpiRank == 0) {
      mergeBoundaryVec.clear();
      BoundaryPArray().swap(mergeBoundaryVec); // actual memory release
      mergeBoundaryVec = boundaryVec; 

      BoundaryPArray tmpBoundaryVec;   
      for (std::size_t it = 0; it < bdryProcess.size(); ++it) {
	if (bdryProcess[it] != 0) {// not root process
	  tmpBoundaryVec.clear();  // do not destroy particles!
	  boostWorld.recv(bdryProcess[it], mpiTag, tmpBoundaryVec);
	  // merge tmpBoundaryVec into mergeBoundaryVec
	  assert(tmpBoundaryVec.size() == mergeBoundaryVec.size());
	  for (std::size_t jt = 0; jt < tmpBoundaryVec.size(); ++jt)
	    mergeBoundaryVec[jt]->getContactInfo().insert(   \
							  mergeBoundaryVec[jt]->getContactInfo().end(), \
							  tmpBoundaryVec[jt]->getContactInfo().begin(), \
							  tmpBoundaryVec[jt]->getContactInfo().end() );
	}    
      }

      // must update after collecting all boundary contact info
      for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	(*it)->updateStatForce();
    }
  }
  

  void Assembly::printBdryContact(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: printBdryContact" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
  
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      (*it)->printContactInfo(ofs);
    }
  
    ofs.close();
  }


  void Assembly::gatherEnergy() {
    calcTransEnergy();
    calcRotatEnergy();
    calcKinetEnergy();
    calcGraviEnergy(allContainer.getMinCorner().getZ());
    calcMechaEnergy();
  }


  void Assembly::closeProg(std::ofstream &ofs) {
    ofs.close();
  }


  void  Assembly::openDepositProg(std::ofstream &ofs, const char *str) {
    ofs.open(str);
    if(!ofs) { debugInf << "stream error: openDepositProg" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << std::setw(OWID) << "iteration"
	<< std::setw(OWID) << "normal_x1"
	<< std::setw(OWID) << "normal_x2"
	<< std::setw(OWID) << "normal_y1"
	<< std::setw(OWID) << "normal_y2"
	<< std::setw(OWID) << "normal_z1"
	<< std::setw(OWID) << "normal_z2"
    
	<< std::setw(OWID) << "contact_x1"
	<< std::setw(OWID) << "contact_x2"
	<< std::setw(OWID) << "contact_y1"
	<< std::setw(OWID) << "contact_y2"
	<< std::setw(OWID) << "contact_z1"
	<< std::setw(OWID) << "contact_z2"
	<< std::setw(OWID) << "contact_inside"
    
	<< std::setw(OWID) << "penetr_x1"
	<< std::setw(OWID) << "penetr_x2"
	<< std::setw(OWID) << "penetr_y1"
	<< std::setw(OWID) << "penetr_y2"
	<< std::setw(OWID) << "penetr_z1"
	<< std::setw(OWID) << "penetr_z2"

	<< std::setw(OWID) << "avgNormal"
	<< std::setw(OWID) << "avgShear"
	<< std::setw(OWID) << "avgPenetr"
    
	<< std::setw(OWID) << "transEnergy"
	<< std::setw(OWID) << "rotatEnergy"
	<< std::setw(OWID) << "kinetEnergy"
	<< std::setw(OWID) << "graviEnergy"
	<< std::setw(OWID) << "mechaEnergy"

	<< std::setw(OWID) << "vibra_est_dt"
	<< std::setw(OWID) << "impact_est_dt"
	<< std::setw(OWID) << "actual_dt"
	<< std::setw(OWID) << "accruedTime"
    
	<< std::endl;
  }


  void Assembly::printDepositProg(std::ofstream &ofs) {
    REAL var[6];
  
    // normalForce
    for (std::size_t i = 0; i < 6; ++i)
      var[i] = 0;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      Vec normal = (*it)->getNormalForce();
      switch (id) {
      case 1: 
	var[0] = fabs(normal.getX());
	break;
      case 2:
	var[1] = normal.getX();
	break;
      case 3:
	var[2] = fabs(normal.getY());
	break;
      case 4:
	var[3] = normal.getY();
	break;
      case 5:
	var[4] = fabs(normal.getZ());
	break;
      case 6:
	var[5] = normal.getZ();
	break;
      }
    }
    ofs << std::setw(OWID) << iteration;
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << var[i];
  
    // contactNum
    for (std::size_t i = 0; i < 6; ++i)
      var[i] = 0;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      var[id - 1] = (*it)->getContactNum();
    }
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << static_cast<std::size_t> (var[i]);
    ofs << std::setw(OWID) << allContactNum;
  
    // avgPenetr
    for (std::size_t i = 0; i < 6; ++i)
      var[i] = 0;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      var[id - 1] = (*it)->getAvgPenetr();
    }
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << var[i];
  
    // average data
    ofs << std::setw(OWID) << avgNormal
	<< std::setw(OWID) << avgShear
	<< std::setw(OWID) << avgPenetr;

    // energy
    ofs << std::setw(OWID) << transEnergy
	<< std::setw(OWID) << rotatEnergy
	<< std::setw(OWID) << kinetEnergy
	<< std::setw(OWID) << graviEnergy
	<< std::setw(OWID) << mechaEnergy;

    // time
    ofs << std::setw(OWID) << vibraTimeStep
	<< std::setw(OWID) << impactTimeStep
	<< std::setw(OWID) << timeStep
	<< std::setw(OWID) << timeAccrued;

    ofs << std::endl;
  }


  void Assembly::getStartDimension(REAL &distX, REAL &distY, REAL &distZ) {
    REAL x1, x2, y1, y2, z1, z2;
    // use boundaryVec
    for(BoundaryPArray::const_iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it) {
      switch ((*it)->getId()) {
      case 1: 
	x1 = (*it)->getPoint().getX();
	break;
      case 2:
	x2 = (*it)->getPoint().getX();
	break;
      case 3:
	y1 = (*it)->getPoint().getY();
	break;
      case 4:
	y2 = (*it)->getPoint().getY();
	break;
      case 5:
	z1 = (*it)->getPoint().getZ();
	break;
      case 6:
	z2 = (*it)->getPoint().getZ();
	break;
      }
    }
    distX = x2 - x1;
    distY = y2 - y1;
    distZ = z2 - z1;
  }


  void  Assembly::openCompressProg(std::ofstream &ofs, const char *str) {
    ofs.open(str);
    if(!ofs) { debugInf << "stream error: openCompressProg" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << std::setw(OWID) << "iteration"
	<< std::setw(OWID) << "traction_x1"
	<< std::setw(OWID) << "traction_x2"
	<< std::setw(OWID) << "traction_y1"
	<< std::setw(OWID) << "traction_y2"
	<< std::setw(OWID) << "traction_z1"
	<< std::setw(OWID) << "traction_z2"
	<< std::setw(OWID) << "mean_stress"

	<< std::setw(OWID) << "bulk_volume"
	<< std::setw(OWID) << "density"
	<< std::setw(OWID) << "epsilon_x"
	<< std::setw(OWID) << "epsilon_y"
	<< std::setw(OWID) << "epsilon_z"
	<< std::setw(OWID) << "epsilon_v"
	<< std::setw(OWID) << "void_ratio"
	<< std::setw(OWID) << "porosity"

	<< std::setw(OWID) << "velocity_x1"
	<< std::setw(OWID) << "velocity_x2"
	<< std::setw(OWID) << "velocity_y1"
	<< std::setw(OWID) << "velocity_y2"
	<< std::setw(OWID) << "velocity_z1"
	<< std::setw(OWID) << "velocity_z2"
    
	<< std::setw(OWID) << "contact_x1"
	<< std::setw(OWID) << "contact_x2"
	<< std::setw(OWID) << "contact_y1"
	<< std::setw(OWID) << "contact_y2"
	<< std::setw(OWID) << "contact_z1"
	<< std::setw(OWID) << "contact_z2"
	<< std::setw(OWID) << "contact_inside"
    
	<< std::setw(OWID) << "penetr_x1"
	<< std::setw(OWID) << "penetr_x2"
	<< std::setw(OWID) << "penetr_y1"
	<< std::setw(OWID) << "penetr_y2"
	<< std::setw(OWID) << "penetr_z1"
	<< std::setw(OWID) << "penetr_z2"
    
	<< std::setw(OWID) << "avgNormal"
	<< std::setw(OWID) << "avgShear"
	<< std::setw(OWID) << "avgPenetr"

	<< std::setw(OWID) << "transEnergy"
	<< std::setw(OWID) << "rotatEnergy"
	<< std::setw(OWID) << "kinetEnergy"
	<< std::setw(OWID) << "graviEnergy"
	<< std::setw(OWID) << "mechaEnergy"

	<< std::setw(OWID) << "vibra_est_dt"
	<< std::setw(OWID) << "impact_est_dt"
	<< std::setw(OWID) << "actual_dt"
    
	<< std::endl;
  }


  void Assembly::printCompressProg(std::ofstream &ofs, REAL distX, REAL distY, REAL distZ) {
    REAL x1, x2, y1, y2, z1, z2;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      switch ((*it)->getId()) {
      case 1: 
	x1 = (*it)->getPoint().getX();
	break;
      case 2:
	x2 = (*it)->getPoint().getX();
	break;
      case 3:
	y1 = (*it)->getPoint().getY();
	break;
      case 4:
	y2 = (*it)->getPoint().getY();
	break;
      case 5:
	z1 = (*it)->getPoint().getZ();
	break;
      case 6:
	z2 = (*it)->getPoint().getZ();
	break;
      }
    }
    REAL areaX = (y2 - y1) * (z2 - z1);
    REAL areaY = (z2 - z1) * (x2 - x1);
    REAL areaZ = (x2 - x1) * (y2 - y1);
    REAL bulkVolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
    REAL voidRatio = bulkVolume / getParticleVolume() - 1;

    REAL var[6], vel[6];
    // normalForce
    for (std::size_t i = 0; i < 6; ++i) {
      var[i] = 0;
      vel[i] = 0;
    }
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      Vec normal = (*it)->getNormalForce();
      Vec veloc  = (*it)->getVeloc();
      switch (id) {
      case 1: 
	var[0] = fabs(normal.getX()) / areaX;
	vel[0] = veloc.getX();
	break;
      case 2:
	var[1] = normal.getX() / areaX;
	vel[1] = veloc.getX();
	break;
      case 3:
	var[2] = fabs(normal.getY()) / areaY;
	vel[2] = veloc.getY();
	break;
      case 4:
	var[3] = normal.getY() / areaY;
	vel[3] = veloc.getY();
	break;
      case 5:
	var[4] = fabs(normal.getZ()) / areaZ;
	vel[4] = veloc.getZ();
	break;
      case 6:
	var[5] = normal.getZ() / areaZ;
	vel[5] = veloc.getZ();
	break;
      }
    }
    ofs << std::setw(OWID) << iteration;
    REAL avg = 0;
    for (std::size_t i = 0; i < 6; ++i) {
      ofs << std::setw(OWID) << var[i];
      avg += var[i];
    }
    ofs << std::setw(OWID) << avg/6;  

    // volume
    ofs << std::setw(OWID) << bulkVolume
	<< std::setw(OWID) << getMass() / bulkVolume
	<< std::setw(OWID) << 1-(x2-x1)/distX
	<< std::setw(OWID) << 1-(y2-y1)/distY
	<< std::setw(OWID) << 1-(z2-z1)/distZ
	<< std::setw(OWID) << 3-(x2-x1)/distX - (y2-y1)/distY - (z2-z1)/distZ
	<< std::setw(OWID) << voidRatio
	<< std::setw(OWID) << voidRatio / (1 + voidRatio);

    // velocity
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << vel[i];  

    // contactNum
    for (std::size_t i = 0; i < 6; ++i)
      var[i] = 0;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      var[id - 1] = (*it)->getContactNum();
    }
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << static_cast<std::size_t> (var[i]);
    ofs << std::setw(OWID) << allContactNum;
  
    // avgPenetr
    for (std::size_t i = 0; i < 6; ++i)
      var[i] = 0;
    for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
      std::size_t id = (*it)->getId();
      var[id - 1] = (*it)->getAvgPenetr();
    }
    for (std::size_t i = 0; i < 6; ++i)
      ofs << std::setw(OWID) << var[i];
  
    // average data
    ofs << std::setw(OWID) << avgNormal
	<< std::setw(OWID) << avgShear
	<< std::setw(OWID) << avgPenetr;

    // energy
    ofs << std::setw(OWID) << transEnergy
	<< std::setw(OWID) << rotatEnergy
	<< std::setw(OWID) << kinetEnergy
	<< std::setw(OWID) << graviEnergy
	<< std::setw(OWID) << mechaEnergy;

    // time
    ofs << std::setw(OWID) << vibraTimeStep
	<< std::setw(OWID) << impactTimeStep
	<< std::setw(OWID) << timeStep;

    ofs << std::endl;
  }

  void  Assembly::openPeriProgress(std::ofstream &ofs, const char *str) {
    ofs.open(str);
    if(!ofs) { debugInf << "stream error: openPeriProgress" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << "Title = \"Particle Information\"" << std::endl;
    ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"KE\" \"P\" \"Mises\" \"Volume\" \"horizonSize\"" << std::endl;
  }

  void Assembly::printPeriProgress(std::ofstream &ofs, const int iframe) const {
	ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	REAL sigma11, sigma12; // sigma13;
	REAL /*sigma21,*/ sigma22, sigma23;
	REAL sigma31, /*sigma32,*/ sigma33;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = allPeriParticleVec.begin(); pt!= allPeriParticleVec.end(); pt++) {
//	    if( (*pt)->getInitPosition().getX()> 0.5*(dem::Parameter::getSingleton().parameter["Xmin"]+dem::Parameter::getSingleton().parameter["Xmax"]) )
//		continue;
	    sigma11=(*pt)->getSigma11(); sigma12=(*pt)->getSigma12(); 
            // sigma13=(*pt)->getSigma13();
	    // sigma21=(*pt)->getSigma21(); 
            sigma22=(*pt)->getSigma22(); sigma23=(*pt)->getSigma23();
	    sigma31=(*pt)->getSigma31(); 
            // sigma32=(*pt)->getSigma32(); 
            sigma33=(*pt)->getSigma33();
	    pressure = sigma11 + sigma22 + sigma33;
	    vonMisesStress = sqrt(0.5*( (sigma11-sigma22)*(sigma11-sigma22)
			         + (sigma22-sigma33)*(sigma22-sigma33)
				 + (sigma11-sigma33)*(sigma11-sigma33) )
			     + 3*( sigma12*sigma12				 + sigma23*sigma23
				 + sigma31*sigma31 ));
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::setw(20) << (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getDisplacement().getY() 
		<< std::setw(20) << (*pt)->getDisplacement().getZ() 
		<< std::setw(20) << (*pt)->getVelocity().getX()
		<< std::setw(20) << (*pt)->getVelocity().getY() 
		<< std::setw(20) << (*pt)->getVelocity().getZ() 
		<< std::setw(20) << vfabs((*pt)->getVelocity())
		<< std::setw(20) << pressure
		<< std::setw(20) << vonMisesStress
		<< std::setw(20) << (*pt)->getParticleVolume()
		<< std::setw(20) << (*pt)->getHorizonSize()
		<< std::endl;
	    ofs.flush();
	}
  }

  void Assembly::printPeriProgressHalf(std::ofstream &ofs, const int iframe) const {
	ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	REAL sigma11, sigma12; // sigma13;
	REAL /*sigma21,*/ sigma22, sigma23;
	REAL sigma31, /*sigma32,*/ sigma33;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = allPeriParticleVec.begin(); pt!= allPeriParticleVec.end(); pt++) {
	    if( (*pt)->getInitPosition().getX()> 0.5*(dem::Parameter::getSingleton().parameter["Xmin"]+dem::Parameter::getSingleton().parameter["Xmax"]) )
		continue;
	    sigma11=(*pt)->getSigma11(); sigma12=(*pt)->getSigma12(); 
            // sigma13=(*pt)->getSigma13();
	    // sigma21=(*pt)->getSigma21(); 
            sigma22=(*pt)->getSigma22(); sigma23=(*pt)->getSigma23();
	    sigma31=(*pt)->getSigma31(); 
            // sigma32=(*pt)->getSigma32(); 
            sigma33=(*pt)->getSigma33();
	    pressure = sigma11 + sigma22 + sigma33;
	    vonMisesStress = sqrt(0.5*( (sigma11-sigma22)*(sigma11-sigma22)
			         + (sigma22-sigma33)*(sigma22-sigma33)
				 + (sigma11-sigma33)*(sigma11-sigma33) )
			     + 3*( sigma12*sigma12				 + sigma23*sigma23
				 + sigma31*sigma31 ));
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::setw(20) << (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getDisplacement().getY() 
		<< std::setw(20) << (*pt)->getDisplacement().getZ() 
		<< std::setw(20) << (*pt)->getVelocity().getX()
		<< std::setw(20) << (*pt)->getVelocity().getY() 
		<< std::setw(20) << (*pt)->getVelocity().getZ() 
		<< std::setw(20) << vfabs((*pt)->getVelocity())
		<< std::setw(20) << pressure
		<< std::setw(20) << vonMisesStress
		<< std::setw(20) << (*pt)->getParticleVolume()
		<< std::setw(20) << (*pt)->getHorizonSize()
		<< std::endl;
	    ofs.flush();
	}
  }

  void  Assembly::openParticleProg(std::ofstream &ofs, const char *str) {
    ofs.open(str);
    if(!ofs) { debugInf << "stream error: openParticleProg" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    ofs << std::setw(OWID) << "iteration"
	<< std::setw(OWID) << "accruedTime"
	<< std::setw(OWID) << "penalFx"
	<< std::setw(OWID) << "penalFy"
	<< std::setw(OWID) << "penalFz"
	<< std::setw(OWID) << "pressureFx"
	<< std::setw(OWID) << "pressureFy"
	<< std::setw(OWID) << "pressureFz"
	<< std::setw(OWID) << "penalMx"
	<< std::setw(OWID) << "penalMy"
	<< std::setw(OWID) << "penalMz"
	<< std::setw(OWID) << "pressureMx"
	<< std::setw(OWID) << "pressureMy"
	<< std::setw(OWID) << "pressureMz"
	<< std::setw(OWID) << "accelX"
	<< std::setw(OWID) << "accelY"
	<< std::setw(OWID) << "accelZ"
	<< std::setw(OWID) << "velocX"
	<< std::setw(OWID) << "velocY"
	<< std::setw(OWID) << "velocZ"
	<< std::endl;
  }


  void Assembly::readParticle(const char *inputParticle) {

    REAL young = dem::Parameter::getSingleton().parameter["young"];
    REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

    std::ifstream ifs(inputParticle);
    if(!ifs) { debugInf << "stream error: readParticle" << std::endl; exit(-1); }
    std::size_t particleNum;
    ifs >> particleNum;
    std::string str;
    ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str
	>> str >> str >> str >> str >> str >> str >> str >> str >> str >> str
	>> str >> str >> str >> str >> str >> str >> str >> str >> str;
  
    /*
    ParticlePArray::iterator it;
    for(it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
      delete (*it);
    */
    allParticleVec.clear();

    std::size_t id, type;
    REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
    REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
    for (std::size_t i = 0; i < particleNum; ++i){
      ifs >> id >> type >> a >> b >> c >> px >> py >> pz 
          >> dax >> day >> daz >> dbx >> dby >> dbz >> dcx >> dcy >> dcz
	  >> vx >> vy >> vz >> omx >> omy >> omz >> fx >> fy >> fz >> mx >> my >> mz;
      ParticleP pt= std::make_shared<Particle>(id, type, Vec(a,b,c), Vec(px,py,pz), 
                                               Vec(dax,day,daz), Vec(dbx,dby,dbz), 
                                               Vec(dcx,dcy,dcz), young, poisson);
    
      // optional settings for a particle's initial status
      if ( (static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
	pt->setPrevVeloc(Vec(vx,vy,vz));
	pt->setCurrVeloc(Vec(vx,vy,vz));
	pt->setPrevOmga(Vec(omx,omy,omz));
	pt->setCurrOmga(Vec(omx,omy,omz));
	pt->setForce(Vec(fx,fy,fz));  // initial force
	pt->setMoment(Vec(mx,my,mz)); // initial moment
      }
    
      //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial force
      //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial moment

      allParticleVec.push_back(pt);
    }
  
    std::size_t sieveNum;
    ifs >> sieveNum;
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    for (std::size_t i = 0; i < sieveNum; ++i)
      ifs >> percent[i] >> size[i];
    REAL ratio_ba, ratio_ca;
    ifs >> ratio_ba >> ratio_ca;
    setGradation(Gradation(sieveNum, percent, size, ratio_ba, ratio_ca));
    ifs.close();
  }


  void Assembly::printParticle(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: printParticle" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << allParticleVec.size() << std::endl;
    ofs << std::setw(OWID) << "id"
	<< std::setw(OWID) << "type"
	<< std::setw(OWID) << "radius_a"
	<< std::setw(OWID) << "radius_b"
	<< std::setw(OWID) << "radius_c"
	<< std::setw(OWID) << "position_x"
	<< std::setw(OWID) << "position_y"
	<< std::setw(OWID) << "position_z"
	<< std::setw(OWID) << "axle_a_x"
	<< std::setw(OWID) << "axle_a_y"
	<< std::setw(OWID) << "axle_a_z"
	<< std::setw(OWID) << "axle_b_x"
	<< std::setw(OWID) << "axle_b_y"
	<< std::setw(OWID) << "axle_b_z"
	<< std::setw(OWID) << "axle_c_x"
	<< std::setw(OWID) << "axle_c_y"
	<< std::setw(OWID) << "axle_c_z"
	<< std::setw(OWID) << "velocity_x"
	<< std::setw(OWID) << "velocity_y"
	<< std::setw(OWID) << "velocity_z"
	<< std::setw(OWID) << "omga_x"
	<< std::setw(OWID) << "omga_y"
	<< std::setw(OWID) << "omga_z"
	<< std::setw(OWID) << "force_x"
	<< std::setw(OWID) << "force_y"
	<< std::setw(OWID) << "force_z"
	<< std::setw(OWID) << "moment_x"
	<< std::setw(OWID) << "moment_y"
	<< std::setw(OWID) << "moment_z"
	<< std::endl;
  
    Vec vObj;
    ParticlePArray::const_iterator  it;
    for (it = allParticleVec.begin(); it != allParticleVec.end();++it)  {
      ofs << std::setw(OWID) << (*it)->getId()
	  << std::setw(OWID) << (*it)->getType()
	  << std::setw(OWID) << (*it)->getA()
	  << std::setw(OWID) << (*it)->getB()
	  << std::setw(OWID) << (*it)->getC();
    
      vObj=(*it)->getCurrPos();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecA();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecB();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecC();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrVeloc();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrOmga();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getForce();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getMoment();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ() << std::endl;
    }
  
    std::size_t sieveNum = gradation.getSieveNum();
    std::vector<REAL> percent = gradation.getPercent();
    std::vector<REAL> size    = gradation.getSize();
    ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
    for (std::size_t i = 0; i < sieveNum; ++i)
      ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i] << std::endl;
    ofs << std::endl << std::setw(OWID) << gradation.getPtclRatioBA() << std::setw(OWID) << gradation.getPtclRatioCA() << std::endl;

    ofs.close();
  }


  void Assembly::printParticle(const char *str, ParticlePArray  &particleVec) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: printParticle" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << particleVec.size() << std::endl;
    ofs << std::setw(OWID) << "id"
	<< std::setw(OWID) << "type"
	<< std::setw(OWID) << "radius_a"
	<< std::setw(OWID) << "radius_b"
	<< std::setw(OWID) << "radius_c"
	<< std::setw(OWID) << "position_x"
	<< std::setw(OWID) << "position_y"
	<< std::setw(OWID) << "position_z"
	<< std::setw(OWID) << "axle_a_x"
	<< std::setw(OWID) << "axle_a_y"
	<< std::setw(OWID) << "axle_a_z"
	<< std::setw(OWID) << "axle_b_x"
	<< std::setw(OWID) << "axle_b_y"
	<< std::setw(OWID) << "axle_b_z"
	<< std::setw(OWID) << "axle_c_x"
	<< std::setw(OWID) << "axle_c_y"
	<< std::setw(OWID) << "axle_c_z"
	<< std::setw(OWID) << "velocity_x"
	<< std::setw(OWID) << "velocity_y"
	<< std::setw(OWID) << "velocity_z"
	<< std::setw(OWID) << "omga_x"
	<< std::setw(OWID) << "omga_y"
	<< std::setw(OWID) << "omga_z"
	<< std::setw(OWID) << "force_x"
	<< std::setw(OWID) << "force_y"
	<< std::setw(OWID) << "force_z"
	<< std::setw(OWID) << "moment_x"
	<< std::setw(OWID) << "moment_y"
	<< std::setw(OWID) << "moment_z"
	<< std::endl;
  
    Vec vObj;
    for (auto it = particleVec.cbegin(); it != particleVec.cend();++it)  {
      ofs << std::setw(OWID) << (*it)->getId()
	  << std::setw(OWID) << (*it)->getType()
	  << std::setw(OWID) << (*it)->getA()
	  << std::setw(OWID) << (*it)->getB()
	  << std::setw(OWID) << (*it)->getC();
    
      vObj=(*it)->getCurrPos();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecA();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecB();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrDirecC();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrVeloc();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getCurrOmga();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getForce();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ();
    
      vObj=(*it)->getMoment();
      ofs << std::setw(OWID) << vObj.getX()
	  << std::setw(OWID) << vObj.getY()
	  << std::setw(OWID) << vObj.getZ() << std::endl;
    }
  
    ofs.close();
  }


  void Assembly::readBoundary(const char *str) {
    std::ifstream ifs(str);
    if(!ifs) { debugInf << "stream error: readBoundary" << std::endl; exit(-1); }

    REAL x1, y1, z1, x2, y2, z2;
    ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
    setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
    // compute grid assumed to be the same as container, change in scatterParticle() if necessary.
    setGrid(Rectangle(x1, y1, z1, x2, y2, z2)); 

    boundaryVec.clear();
    BoundaryP bptr;
    std::size_t boundaryNum;
    std::size_t type;
    ifs >> boundaryNum;
    for(std::size_t i = 0; i < boundaryNum; ++i) {
      ifs >> type;
      if(type == 1) // plane boundary
	bptr = std::make_shared<planeBoundary>(type, ifs);
      else if(type == 2) // cylindrical boundary
	bptr = std::make_shared<cylinderBoundary>(type, ifs);

      boundaryVec.push_back(bptr);
    }

    ifs.close();
  }


  void Assembly::printBoundary(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: printBoundary" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
  
    Vec v1 = allContainer.getMinCorner();
    Vec v2 = allContainer.getMaxCorner();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
  
    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
	<< std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl << std::endl
	<< std::setw(OWID) << boundaryVec.size() << std::endl;

    for(BoundaryPArray::const_iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
      (*it)->print(ofs);
  
    ofs.close();
  }


  void Assembly::plotBoundary(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: plotBoundary" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    REAL x1, y1, z1, x2, y2, z2;
    x1 = allContainer.getMinCorner().getX();
    y1 = allContainer.getMinCorner().getY();
    z1 = allContainer.getMinCorner().getZ();
    x2 = allContainer.getMaxCorner().getX();
    y2 = allContainer.getMaxCorner().getY();
    z2 = allContainer.getMaxCorner().getZ();

    ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
    ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
    ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
    ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
    ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
    ofs << "1 2 3 4 5 6 7 8" << std::endl;

    ofs.close();
  }


  void Assembly::plotGrid(const char *str) const {
    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: plotGrid" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;

    ofs << "ZONE N=" << (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1)
	<< ", E=" << mpiProcX * mpiProcY * mpiProcZ << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;

    std::vector<Vec> coords((mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1));
    std::size_t index = 0;
    for (auto i = 0; i < mpiProcX + 1; ++i)
      for (auto j = 0; j < mpiProcY + 1; ++j)
	for (auto k = 0; k < mpiProcZ + 1; ++k)
	  coords[index++] = Vec(v1.getX() + vspan.getX() / mpiProcX * i,
				v1.getY() + vspan.getY() / mpiProcY * j,
				v1.getZ() + vspan.getZ() / mpiProcZ * k);

    for (auto i = 0; i < (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1); ++i)
      ofs << std::setw(OWID) << coords[i].getX() 
	  << std::setw(OWID) << coords[i].getY() 
	  << std::setw(OWID) << coords[i].getZ() << std::endl;

    for (int iRank = 0; iRank < mpiSize; ++iRank) {
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, 3, coords);

      int id4 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + coords[2];
      int id1 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + coords[2];
      int id3 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + coords[2];
      int id2 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + coords[2];

      int id8 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + (coords[2]+1);
      int id5 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + coords[1]*(mpiProcZ+1) + (coords[2]+1);
      int id7 = 1 + coords[0]*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + (coords[2]+1);
      int id6 = 1 + (coords[0]+1)*(mpiProcZ+1)*(mpiProcY+1) + (coords[1]+1)*(mpiProcZ+1) + (coords[2]+1);

      ofs << std::setw(8) << id1 << std::setw(8) << id2 << std::setw(8) << id3 << std::setw(8) << id4 
	  << std::setw(8) << id5 << std::setw(8) << id6 << std::setw(8) << id7 << std::setw(8) << id8 << std::endl;
    }

    ofs.close();
  }


  void Assembly::findContact() { // various implementations
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];

    if (ompThreads == 1) { // non-openmp single-thread version, time complexity bigO(n x n), n is the number of particles.  
      contactVec.clear();
    
#ifdef TIME_PROFILE
      REAL time_r = 0; // time consumed in contact resolution, i.e., tmpContact.isOverlapped()
      gettimeofday(&time_p1, NULL); 
#endif
    
      std::size_t num1 = particleVec.size();      // particles inside container
      std::size_t num2 = mergeParticleVec.size(); // paticles inside container (at front) + particles from neighboring blocks (at end)
      for (std::size_t i = 0; i < num1; ++i) {    // NOT (num1 - 1), in parallel situation where one particle could contact received particles!
	Vec u = particleVec[i]->getCurrPos();
	for (std::size_t j = i + 1; j < num2; ++j){
	  Vec v = mergeParticleVec[j]->getCurrPos();
	  if ( ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA())
	       && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	    Contact tmpContact(particleVec[i].get(), mergeParticleVec[j].get()); // a local and temparory object
#ifdef TIME_PROFILE
	    gettimeofday(&time_r1, NULL); 
#endif
	    if(tmpContact.isOverlapped())
	      contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
#ifdef TIME_PROFILE
	    gettimeofday(&time_r2, NULL); 
	    time_r += timediffsec(time_r1, time_r2);
#endif
	  }
	}
      }	
    
#ifdef TIME_PROFILE
      gettimeofday(&time_p2, NULL);
      debugInf<< std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
#endif

    }

    else if (ompThreads > 1) { // openmp implementation: various loop scheduling - (static), (static,1), (dynamic), (dynamic,1)
      contactVec.clear();
    
#ifdef TIME_PROFILE
      gettimeofday(&time_p1, NULL); 
#endif
    
      std::size_t i, j;
      Vec u, v;
      std::size_t num1 = particleVec.size();
      std::size_t num2 = mergeParticleVec.size();  
      int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
    
#pragma omp parallel for num_threads(ompThreads) private(i, j, u, v) shared(num1, num2) schedule(dynamic)
      for (i = 0; i < num1; ++i) { 
	u = particleVec[i]->getCurrPos();
	for (j = i + 1; j < num2; ++j) {
	  v = mergeParticleVec[j]->getCurrPos();
	  if ( ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA() )
	       && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
	       && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
	       && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
	    Contact tmpContact(particleVec[i].get(), mergeParticleVec[j].get()); // a local and temparory object
	    if(tmpContact.isOverlapped())
#pragma omp critical
	      contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
	  }
	}
      }
    
#ifdef TIME_PROFILE
      gettimeofday(&time_p2, NULL);
      debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
#endif

    } // end of openmp implementation

  }


  void Assembly::internalForce(){
    REAL pAvg[3], sum[3];
    for (std::size_t i = 0; i < 3; ++i) {
      pAvg[i] = 0;
      sum[i]  = 0;
    }

    if(contactVec.size() > 0) {
      for (ContactArray::iterator it = contactVec.begin(); it != contactVec.end(); ++it)
	it->checkinPrevTgt(contactTgtVec); // checkin previous tangential force and displacment    
    
#ifdef TIME_PROFILE
      gettimeofday(&time_p1,NULL); 
#endif 

      contactTgtVec.clear(); // contactTgtVec must be cleared before filling in new values.
      for (ContactArray::iterator it = contactVec.begin(); it != contactVec.end(); ++ it){
	it->contactForce();             // cannot be parallelized as it may change a particle's force simultaneously.
	it->checkoutTgt(contactTgtVec); // checkout current tangential force and displacment
	pAvg[0] += it->getNormalForce();
	pAvg[1] += it->getTgtForce();
	pAvg[2] += it->getPenetration();
      }
      for (std::size_t i = 0; i < 3; ++i)
	pAvg[i] /= contactVec.size();
    
#ifdef TIME_PROFILE
      gettimeofday(&time_p2,NULL);
      debugInf << std::setw(OWID) << "internalForce=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::endl; 
#endif   
    }

    MPI_Reduce(pAvg, sum, 3, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
    avgNormal = sum[0]/mpiSize;
    avgShear  = sum[1]/mpiSize;
    avgPenetr = sum[2]/mpiSize;
  }


  void Assembly::updateParticle() {
    for(ParticlePArray::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      (*it)->update();
  }


  void Assembly::updateBoundary(REAL sigma, std::string type, REAL sigmaX, REAL sigmaY) {
    if (mpiRank == 0) {
      REAL x1, x2, y1, y2, z1, z2;
      for(BoundaryPArray::const_iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it) {
	switch ((*it)->getId()) {
	case 1: 
	  x1 = (*it)->getPoint().getX();
	  break;
	case 2:
	  x2 = (*it)->getPoint().getX();
	  break;
	case 3:
	  y1 = (*it)->getPoint().getY();
	  break;
	case 4:
	  y2 = (*it)->getPoint().getY();
	  break;
	case 5:
	  z1 = (*it)->getPoint().getZ();
	  break;
	case 6:
	  z2 = (*it)->getPoint().getZ();
	  break;
	}
      }
      REAL areaX = (y2 - y1) * (z2 - z1);
      REAL areaY = (z2 - z1) * (x2 - x1);
      REAL areaZ = (x2 - x1) * (y2 - y1);

      if (type.compare("isotropic") == 0) {
	for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	  (*it)->updateIsotropic(sigma, areaX, areaY, areaZ);
      } else if (type.compare("odometer") == 0) {
	for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	  (*it)->updateOdometer(sigma, areaX, areaY, areaZ);
      } else if (type.compare("triaxial") == 0) {
	for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	  (*it)->updateTriaxial(sigma, areaX, areaY, areaZ);
      } else if (type.compare("plnstrn") == 0) {
	for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	  (*it)->updatePlaneStrain(sigma, areaX, areaY, areaZ);
      } else if (type.compare("trueTriaxial") == 0) {
	for(BoundaryPArray::iterator it = mergeBoundaryVec.begin(); it != mergeBoundaryVec.end(); ++it)
	  (*it)->updateTrueTriaxial(sigma, areaX, areaY, areaZ, sigmaX, sigmaY);
      }

      // update boundaryVec from mergeBoundaryVec and remove contactInfo to reduce MPI transmission
      boundaryVec = mergeBoundaryVec; 
      for(BoundaryPArray::iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
	(*it)->clearContactInfo();

      // update allContainer
      for(BoundaryPArray::const_iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it) {
	switch ((*it)->getId()) {
	case 1: 
	  x1 = (*it)->getPoint().getX();
	  break;
	case 2:
	  x2 = (*it)->getPoint().getX();
	  break;
	case 3:
	  y1 = (*it)->getPoint().getY();
	  break;
	case 4:
	  y2 = (*it)->getPoint().getY();
	  break;
	case 5:
	  z1 = (*it)->getPoint().getZ();
	  break;
	case 6:
	  z2 = (*it)->getPoint().getZ();
	  break;
	}
      }
      setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
    }

    broadcast(boostWorld, boundaryVec, 0);

  }


  void Assembly::clearContactForce() {
    for(ParticlePArray::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
      (*it)->clearContactForce();
  }


  void Assembly::findBdryContact() {
    for(BoundaryPArray::iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
      (*it)->findBdryContact(particleVec);
  }


  void Assembly::boundaryForce() {
    for(BoundaryPArray::iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
      (*it)->boundaryForce(boundaryTgtMap);
  }


  void Assembly::printContact(char *str) const
  {
    // There are two implementions of printContact
    // implementation 1: parallel IO, each process prints to a data file using a shared pointer.
    //                   and use post-processing tool to remove redundant info.
    MPI_Status status;
    MPI_File contactFile;
    MPI_File_open(mpiWorld, str, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &contactFile);
    if(boostWorld.rank() == 0 && !contactFile) { std::cout << "stream error: printContact" << std::endl; exit(-1);}

    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);

    for (ContactArray::const_iterator it = contactVec.begin(); it != contactVec.end(); ++it)
      inf << std::setw(OWID) << it->getP1()->getId()
	  << std::setw(OWID) << it->getP2()->getId()
	  << std::setw(OWID) << it->getPoint1().getX()
	  << std::setw(OWID) << it->getPoint1().getY()
	  << std::setw(OWID) << it->getPoint1().getZ()
	  << std::setw(OWID) << it->getPoint2().getX()
	  << std::setw(OWID) << it->getPoint2().getY()
	  << std::setw(OWID) << it->getPoint2().getZ()
	  << std::setw(OWID) << it->getRadius1()
	  << std::setw(OWID) << it->getRadius2()
	  << std::setw(OWID) << it->getPenetration()
	  << std::setw(OWID) << it->getTgtDisp()
	  << std::setw(OWID) << it->getContactRadius()
	  << std::setw(OWID) << it->getR0()
	  << std::setw(OWID) << it->getE0()
	  << std::setw(OWID) << it->getNormalForce()
	  << std::setw(OWID) << it->getTgtForce()
	  << std::setw(OWID) << ( it->getPoint1().getX() + it->getPoint2().getX() )/2
	  << std::setw(OWID) << ( it->getPoint1().getY() + it->getPoint2().getY() )/2
	  << std::setw(OWID) << ( it->getPoint1().getZ() + it->getPoint2().getZ() )/2
	  << std::setw(OWID) << it->normalForceVec().getX()
	  << std::setw(OWID) << it->normalForceVec().getY()
	  << std::setw(OWID) << it->normalForceVec().getZ()
	  << std::setw(OWID) << it->tgtForceVec().getX()
	  << std::setw(OWID) << it->tgtForceVec().getY()
	  << std::setw(OWID) << it->tgtForceVec().getZ()
	  << std::setw(OWID) << it->getVibraTimeStep()
	  << std::setw(OWID) << it->getImpactTimeStep()
	  << std::endl;

    int length = (OWID*28 + 1) *contactVec.size();
    // write a file at a location specified by a shared file pointer (blocking, collective)
    // note MPI_File_write_shared is non-collective
    MPI_File_write_ordered(contactFile, const_cast<char*> (inf.str().c_str()), length, MPI_CHAR, &status);
    MPI_File_close(&contactFile);

    // implementation 2: each process prints to an individual file.
    //                   use post-processing tool to merge files and remove redundance.
    /*
    char csuf[10];
    combineString(csuf, ".p", mpiRank, 5);
    strcat(str, csuf);

    std::ofstream ofs(str);
    if(!ofs) { debugInf << "stream error: printContact" << std::endl; exit(-1); }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
  
    ofs << std::setw(OWID) << contactVec.size() << std::endl;
    ofs << std::setw(OWID) << "ptcl_1"
	<< std::setw(OWID) << "ptcl_2"
	<< std::setw(OWID) << "point1_x"
	<< std::setw(OWID) << "point1_y"
	<< std::setw(OWID) << "point1_z"
	<< std::setw(OWID) << "point2_x"
	<< std::setw(OWID) << "point2_y"
	<< std::setw(OWID) << "point2_z"
	<< std::setw(OWID) << "radius_1"
	<< std::setw(OWID) << "radius_2"
	<< std::setw(OWID) << "penetration"
	<< std::setw(OWID) << "tangt_disp"
	<< std::setw(OWID) << "contact_radius"
	<< std::setw(OWID) << "R0"
	<< std::setw(OWID) << "E0"
	<< std::setw(OWID) << "normal_force"
	<< std::setw(OWID) << "tangt_force"
	<< std::setw(OWID) << "contact_x"
	<< std::setw(OWID) << "contact_y"
	<< std::setw(OWID) << "contact_z"
	<< std::setw(OWID) << "normal_x"
	<< std::setw(OWID) << "normal_y"
	<< std::setw(OWID) << "normal_z"
	<< std::setw(OWID) << "tangt_x"
	<< std::setw(OWID) << "tangt_y"
	<< std::setw(OWID) << "tangt_z"
	<< std::setw(OWID) << "vibra_t_step"
	<< std::setw(OWID) << "impact_t_step"
	<< std::endl;
  
    ContactArray::const_iterator it;
    for (it = contactVec.begin(); it != contactVec.end(); ++it)
      ofs << std::setw(OWID) << it->getP1()->getId()
	  << std::setw(OWID) << it->getP2()->getId()
	  << std::setw(OWID) << it->getPoint1().getX()
	  << std::setw(OWID) << it->getPoint1().getY()
	  << std::setw(OWID) << it->getPoint1().getZ()
	  << std::setw(OWID) << it->getPoint2().getX()
	  << std::setw(OWID) << it->getPoint2().getY()
	  << std::setw(OWID) << it->getPoint2().getZ()
	  << std::setw(OWID) << it->getRadius1()
	  << std::setw(OWID) << it->getRadius2()
	  << std::setw(OWID) << it->getPenetration()
	  << std::setw(OWID) << it->getTgtDisp()
	  << std::setw(OWID) << it->getContactRadius()
	  << std::setw(OWID) << it->getR0()
	  << std::setw(OWID) << it->getE0()
	  << std::setw(OWID) << it->getNormalForce()
	  << std::setw(OWID) << it->getTgtForce()
	  << std::setw(OWID) << ( it->getPoint1().getX() + it->getPoint2().getX() )/2
	  << std::setw(OWID) << ( it->getPoint1().getY() + it->getPoint2().getY() )/2
	  << std::setw(OWID) << ( it->getPoint1().getZ() + it->getPoint2().getZ() )/2
	  << std::setw(OWID) << it->normalForceVec().getX()
	  << std::setw(OWID) << it->normalForceVec().getY()
	  << std::setw(OWID) << it->normalForceVec().getZ()
	  << std::setw(OWID) << it->tgtForceVec().getX()
	  << std::setw(OWID) << it->tgtForceVec().getY()
	  << std::setw(OWID) << it->tgtForceVec().getZ()
	  << std::setw(OWID) << it->getVibraTimeStep()
	  << std::setw(OWID) << it->getImpactTimeStep()
	  << std::endl;
    ofs.close();
    */
  }


  void Assembly::calcTransEnergy() {
    REAL pEngy = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin(); it != particleVec.end(); ++it) {
      if ((*it)->getType() == 0)
	pEngy += (*it)->getTransEnergy();
    }
    MPI_Reduce(&pEngy, &transEnergy, 1, MPI_DOUBLE, MPI_SUM, 0,  mpiWorld);
  }


  void Assembly::calcRotatEnergy() {
    REAL pEngy = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin(); it != particleVec.end(); ++it) {
      if ((*it)->getType() == 0)
	pEngy += (*it)->getRotatEnergy();
    }
    MPI_Reduce(&pEngy, &rotatEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
  }


  void Assembly::calcKinetEnergy() {
    REAL pEngy = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin();it != particleVec.end(); ++it) {
      if ((*it)->getType() == 0)
	pEngy += (*it)->getKinetEnergy();
    }
    MPI_Reduce(&pEngy, &kinetEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
  }


  void Assembly::calcGraviEnergy(REAL ref) {
    REAL pEngy = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin();it != particleVec.end(); ++it) {
      if ((*it)->getType() == 0)
	pEngy += (*it)->getPotenEnergy(ref);
    }
    MPI_Reduce(&pEngy, &graviEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
  }


  void Assembly::calcMechaEnergy() {
    mechaEnergy = kinetEnergy +  graviEnergy;
  }


  REAL Assembly::getMass() const {
    REAL var = 0;
    for (ParticlePArray::const_iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
      var += (*it)->getMass();
    return var;
  }


  REAL Assembly::getParticleVolume() const {
    REAL var = 0;
    for (ParticlePArray::const_iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
      if ((*it)->getType() == 0)
	var += (*it)->getVolume();
    return var;
  }


  void Assembly::calcTimeStep() {
    calcVibraTimeStep();
    calcImpactTimeStep();
    calcContactNum();

    REAL CFL = 0.5;
    std::valarray<REAL> dt(3);
    dt[0] = dem::Parameter::getSingleton().parameter["timeStep"];
    dt[1] = CFL * vibraTimeStep;
    dt[2] = CFL * impactTimeStep;

    timeStep = dt.min();
  }


  void Assembly::calcContactNum() {
    std::size_t pContactNum = contactVec.size();
    MPI_Reduce(&pContactNum, &allContactNum, 1, MPI_INT, MPI_SUM, 0, mpiWorld);
  }


  void Assembly::calcVibraTimeStep() {
    REAL pTimeStep = 1/EPS;
    if (contactVec.size() == 0)
      pTimeStep = 1/EPS;
    else {
      ContactArray::const_iterator it = contactVec.begin();
      pTimeStep = it->getVibraTimeStep();
      for (++it; it != contactVec.end(); ++it) {
	REAL val = it->getVibraTimeStep(); 
	pTimeStep = val < pTimeStep ? val : pTimeStep;
      }
    }

    MPI_Allreduce(&pTimeStep, &vibraTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  }


  void Assembly::calcImpactTimeStep() {
    REAL pTimeStep = 1/EPS;
    if (contactVec.size() == 0)
      pTimeStep = 1/EPS;
    else {
      ContactArray::const_iterator it = contactVec.begin();
      pTimeStep = it->getImpactTimeStep();
      for (++it; it != contactVec.end(); ++it) {
	REAL val = it->getImpactTimeStep(); 
	pTimeStep = val < pTimeStep ? val : pTimeStep;
      }
    }

    MPI_Allreduce(&pTimeStep, &impactTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  }


  REAL Assembly::getAvgTransVelocity() const {
    REAL avgv = 0;
    std::size_t count = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin(); it != particleVec.end(); ++it)
      if ((*it)->getType() == 0) {
	avgv += vfabs((*it)->getCurrVeloc());
	++count;
      }
    return avgv /= count;
  }


  REAL Assembly::getAvgRotatVelocity() const {
    REAL avgv = 0;
    std::size_t count = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin(); it != particleVec.end(); ++it)
      if ((*it)->getType() == 0) {
	avgv += vfabs((*it)->getCurrOmga());
	++count;
      }
    return avgv /= count;
  }


  REAL Assembly::getAvgForce() const {
    REAL avgv = 0;
    std::size_t count = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin(); it != particleVec.end(); ++it)
      if ((*it)->getType() == 0) {
	avgv += vfabs((*it)->getForce());
	++count;
      }
    return avgv/count;
  }


  REAL Assembly::getAvgMoment() const {
    REAL avgv = 0;
    std::size_t count = 0;
    ParticlePArray::const_iterator it;
    for (it = particleVec.begin();it != particleVec.end(); ++it)
      if ((*it)->getType() == 0) {
	avgv += vfabs((*it)->getMoment());
	++count;
      }
    return avgv /= count;
  }


  void Assembly::buildBoundary(std::size_t boundaryNum,
			       const char *boundaryFile)
  {
    std::ofstream ofs(boundaryFile);
    if(!ofs) { debugInf << "stream error: buildBoundary" << std::endl; exit(-1);}
    ofs.setf(std::ios::scientific, std::ios::floatfield);

    Vec  v1 = allContainer.getMinCorner();
    Vec  v2 = allContainer.getMaxCorner();
    Vec  v0 = allContainer.getCenter();
    REAL x1 = v1.getX();
    REAL y1 = v1.getY();
    REAL z1 = v1.getZ();
    REAL x2 = v2.getX();
    REAL y2 = v2.getY();
    REAL z2 = v2.getZ();
    REAL x0 = v0.getX();
    REAL y0 = v0.getY();
    REAL z0 = v0.getZ();

    ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
	<< std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl << std::endl
	<< std::setw(OWID) << boundaryNum << std::endl << std::endl;
  
    if (boundaryNum == 1) {   // only a bottom boundary, i.e., boundary 5
      ofs << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl
      
	  << std::setw(OWID) << 5
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z1 << std::endl << std::endl;
    
    }
    else if (boundaryNum == 5) { // no top boundary, i.e., no boundary 6
      // boundary 1
      ofs << std::setw(OWID) << 1
	  << std::setw(OWID) << 1 << std::endl
      
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x1
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z0 << std::endl

	  << std::setw(OWID) << " "
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 1 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z2 << std::endl << std::endl
      
	// boundary 2
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 1 << std::endl

	  << std::setw(OWID) << 2
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x2     
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z0 << std::endl

	  << std::setw(OWID) << " "
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 1 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z2 << std::endl << std::endl
      
	// boundary 3
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 1 << std::endl

	  << std::setw(OWID) << 3
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y1
	  << std::setw(OWID) << z0 << std::endl

	  << std::setw(OWID) << " "
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 1 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z2 << std::endl << std::endl
      
	// boundary 4
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 1 << std::endl

	  << std::setw(OWID) << 4
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x0      
	  << std::setw(OWID) << y2
	  << std::setw(OWID) << z0 << std::endl

	  << std::setw(OWID) << " "
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 1 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z2 << std::endl << std::endl
      
	// boundary 5
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 5
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z1 << std::endl << std::endl;
    }
    else if (boundaryNum == 6) { // all 6 boundaries
      // boundary 1
      ofs << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 1      
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x1
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z0 << std::endl << std::endl

	// boundary 2
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 2
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x2     
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z0 << std::endl << std::endl
      
	// boundary 3
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 3
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y1
	  << std::setw(OWID) << z0 << std::endl << std::endl
      
	// boundary 4
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 4
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << x0      
	  << std::setw(OWID) << y2
	  << std::setw(OWID) << z0 << std::endl << std::endl
      
	// boundary 5
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 5
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << -1
	  << std::setw(OWID) << x0
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z1 << std::endl << std::endl
      
	// boundary 6
	  << std::setw(OWID) << 1
	  << std::setw(OWID) << 0 << std::endl

	  << std::setw(OWID) << 6
	  << std::setw(OWID) << 0 
	  << std::setw(OWID) << 0
	  << std::setw(OWID) << 1    
	  << std::setw(OWID) << x0      
	  << std::setw(OWID) << y0
	  << std::setw(OWID) << z2 << std::endl << std::endl;
    }
  
    ofs.close();
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// periDynamics part //////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
    void Assembly::initialPeriDynamics(const char* inputFile){
/* // not used in DEM-PD coupling code, i.e. rigidInclusion()
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	std::cout << "Problem Initilization " << std::endl;
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	std::cout << "Read data file ..." << std::endl;
	readPeriDynamicsData(inputFile);
	std::cout << "Calculate particle volume ..." << std::endl;
	calcParticleVolume();
	// writeMeshCheckVolume("checkv.dat"); exit(1);
	std::cout << "Calculate horizon size ..." << std::endl;
	calcHorizonSize();
	std::cout << "Construct neighor list ..." << std::endl;
	constructNeighbor();
	std::cout << "Calculate Kinv ..." << std::endl;
	calcParticleKinv();
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	    (*pt)->initial();
	}

	// prescrible the essential boundary condition
	std::cout << "Prescribe the boundary condition ..." << std::endl;
	prescribeEssentialBoundaryCondition(0);

	// calculate the stress at each particle
	std::cout << "Calculate particle stress ..." << std::endl;
	calcParticleStress();

	// calculate the acceleration for each particle
	std::cout << "Calculate particle acceleration ..." << std::endl;
	calcParticleAcceleration();

	// traction boundary 
	ApplyExternalForce(0);

	// apply initial velocity boundary condition, in this case give the cubicTopBoundaryVec particles initial velocities
	//for(std::vector<PeriParticle*>::iterator pt=cubicTopBoundaryVec.begin(); pt!=cubicTopBoundaryVec.end(); pt++){
	//	(*pt)->setInitVelocity(Vec(0.0, 0.0, 1.0));
	//}
*/
    } // end initialPeriDynamics


    void Assembly::prescribeEssentialBoundaryCondition(const int istep){

//	for(std::vector<periDynamics::PeriParticle*>::iterator pt=bottomBoundaryVec.begin(); pt!=bottomBoundaryVec.end(); pt++){
//		(*pt)->prescribeBottomDisplacement(0.0);	// fix z displacement in the z direction
//	}

	//REAL dispz;
	//if(istep <= 200) {
	//	dispz = 1.8*0.05*double(istep)/200.;
	//}else {
	//	dispz = 1.8*0.05;
	//}

	//for(std::vector<periDynamics::PeriParticle*>::iterator pt=topBoundaryVec.begin(); pt!=topBoundaryVec.end(); pt++){
	//	(*pt)->prescribeTopDisplacement(dispz);	// fix z displacement in the z direction
	//}

    } // end prescribeEssentialBoundaryCondition()


    void Assembly::solve(const char* outputFile){
/* // not used in the DEM-PD coupling code, i.e. rigidInclusion()
	// open the tecplot file for output
	std::ofstream ofs(outputFile);
	int iframe = 0;
	writeParticleTecplot(ofs,iframe);
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	std::cout << "Start of the time loop " << std::endl;
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	std::ofstream datafile("uxyz.dat");
	datafile.setf(std::ios::scientific, std::ios::floatfield);
	datafile.precision(10);
	datafile << "VARIABLES = \"Time step\", \"UX\", \"UY\", \"UZ\"" << std::endl;

	int nsteps = 10000;	// is not true
	for(int istep = 1; istep <= nsteps; istep++) {
			
	    runFirstHalfStep();

	    prescribeEssentialBoundaryCondition(istep);

	    checkBondParticleAlive();

	    calcParticleStress();

	    calcParticleAcceleration();

	    ApplyExternalForce(istep);

	    runSecondHalfStep();

//	    if( istep % printInterval == 0) {
//		std::cout << "*** current time step is	" << istep << std::endl;
//		iframe++;
//		writeParticleTecplot(ofs,iframe);
//		datafile << istep
//		<< std::setw(20) << periParticleVec[568]->getDisplacement().getX()
//		<< std::setw(20) << periParticleVec[568]->getDisplacement().getY()
//		<< std::setw(20) << periParticleVec[568]->getDisplacement().getZ() << std::endl;
//	    }
//	    if( istep % 200 == 0) {
//		writeDisplacementData("ux.dat","uy.dat","uz.dat");
//	    }

	} // time loop
	ofs.close();
	datafile.close();
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	std::cout << "Simulation Finished !" << std::endl;
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	writeDisplacementData("ux.dat","uy.dat","uz.dat");
*/
    } // end solve()


    void Assembly::writeDisplacementData(const char *outputFilex, const char *outputFiley, const char *outputFilez) {
/* // not used in the DEM-PD coupling code, i.e. rigidInclusion()
	// displacment along the x axis
	std::ofstream ofs(outputFilex);
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "VARIABLES = \"X\", \"UX\"" << std::endl;
//	for(int index = 0; index < 5; index++){
//	    int node = Uxindex[index];
//	    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getX()
//		<< std::setw(20) << periParticleVec[node]->getDisplacement().getX() << std::endl;
//	}
	ofs.flush();
	ofs.close();

	//dispalcement along the y axis
	ofs.open(outputFiley);
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "VARIABLES = \"Y\", \"UY\"" << std::endl;
//	for(int index = 0; index < 5; index++){
//	    int node = Uyindex[index];
//	    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getY()
//		<< std::setw(20) << periParticleVec[node]->getDisplacement().getY() << std::endl;
//	}
	ofs.flush();
	ofs.close();

	//dispalcement along the z axis
	ofs.open(outputFilez);
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "VARIABLES = \"Z\", \"UZ\"" << std::endl;
//	for(int index = 0; index < 39; index++){
//	    int node = Uzindex[index];
//	    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getZ()
//		<< std::setw(20) << periParticleVec[node]->getDisplacement().getZ() << std::endl;
//	}
	ofs.flush();
	ofs.close();
*/
    } // end writeDisplacementData


    void Assembly::runFirstHalfStep(){ 
	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int i;
	num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num) schedule(dynamic)
	for(i=0; i<num; i++){
	    periParticleVec[i]->updateDisplacement();
	}
    } // end runFirstHalfStep()


    void Assembly::runSecondHalfStep(){ 
	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int i;
	num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num) schedule(dynamic)
	for(i=0; i<num; i++){
	    periParticleVec[i]->updateVelocity();
	}
    } // end runSecondHalfStep()


    void Assembly::constructNeighbor(){	// this function should be called after scattering in each cpu to construct peri-bonds
	// neighbor - searches and constructs the particle neighborlists
	// construct neighborlist for all particles ...
	// compute the weighting function for all particles ...
	if(periParticleVec.empty())
	    return;
	for(std::vector<periDynamics::PeriParticle*>::iterator i_nt=periParticleVec.begin(); i_nt!=periParticleVec.end(); i_nt++){
	    (*i_nt)->clearPeriBonds();	// bondVec should be empty at this time	    
	}
	periBondVec.clear();
        std::size_t num = periParticleVec.size(); 
        int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
        std::size_t i_nt, j_nt;
#pragma omp parallel for num_threads(ompThreads) private(i_nt, j_nt) shared(num) schedule(dynamic)
	for(i_nt=0; i_nt<num-1; i_nt++){

	    Vec coord0_i = periParticleVec[i_nt] ->getInitPosition();
	    REAL horizonSize_i = periParticleVec[i_nt]->getHorizonSize();
	    for(j_nt=i_nt+1; j_nt<num; j_nt++){
		Vec coord0_j = periParticleVec[j_nt]->getInitPosition();
		REAL tmp_length = vfabs(coord0_i-coord0_j);

		REAL horizonSize_j = periParticleVec[j_nt]->getHorizonSize();
		REAL horizonSize_ij = (horizonSize_i+horizonSize_j)*0.5;//This will lead to the fact that horizion is not a sphere!!! 

		REAL ratio = tmp_length/horizonSize_ij;

		// establish the neighbor list
		if(ratio <= 2.0){ 
		    // create bond
		    periDynamics::PeriBond* bond_pt = new periDynamics::PeriBond(tmp_length, periParticleVec[i_nt], periParticleVec[j_nt]);
#pragma omp critical
{
		    periParticleVec[i_nt]->pushBackBondVec(bond_pt);
		    periParticleVec[j_nt]->pushBackBondVec(bond_pt);
		    periBondVec.push_back(bond_pt);
}
		    REAL factor = 3.0/(2.0*Pi*horizonSize_ij*horizonSize_ij*horizonSize_ij); // for the factor of 3d window function
				
		    // weighting function (influence function)
		    if(ratio < 1.0){
			bond_pt->setWeight( factor*(2.0/3.0-ratio*ratio+0.5*ratio*ratio*ratio) );
		    }
		    else{
			bond_pt->setWeight( factor*(2.0-ratio)*(2.0-ratio)*(2.0-ratio)/6.0 );
		    }		     
		} // if(ratio<2.0)

	    } // end j_nt
	} // end i_nt

    } // end constNeighbor()


    void Assembly::findRecvPeriBonds(){	// this function should be called after commuPeriParticle() in each cpu to construct peri-bonds
	// neighbor - searches and constructs the particle neighborlists
	// construct neighborlist for all particles ...
	// compute the weighting function for all particles ...
	if(recvPeriParticleVec.empty())
	    return;
	for(std::vector<periDynamics::PeriParticle*>::iterator i_nt=recvPeriParticleVec.begin(); i_nt<recvPeriParticleVec.end(); i_nt++){
	    (*i_nt)->clearPeriBonds();	// bondVec should be empty at this time	    
	}
	recvPeriBondVec.clear();
	for(std::vector<periDynamics::PeriParticle*>::iterator i_nt=recvPeriParticleVec.begin(); i_nt<recvPeriParticleVec.end(); i_nt++){
	    Vec coord0_i = (*i_nt) ->getInitPosition();
	    REAL horizonSize_i = (*i_nt)->getHorizonSize();
	    for(std::vector<periDynamics::PeriParticle*>::iterator j_nt=periParticleVec.begin(); j_nt<periParticleVec.end(); j_nt++){
		Vec coord0_j = (*j_nt)->getInitPosition();					
		REAL tmp_length = vfabs(coord0_i-coord0_j);					

		REAL horizonSize_j = (*j_nt)->getHorizonSize();
		REAL horizonSize_ij = (horizonSize_i+horizonSize_j)*0.5;//This will lead to the fact that horizion is not a sphere!!! 

		REAL ratio = tmp_length/horizonSize_ij;

		// establish the neighbor list
		if(ratio <= 2.0){ 
		    // create bond
		    periDynamics::PeriBond* bond_pt = new periDynamics::PeriBond(tmp_length, *i_nt, *j_nt);
		    (*j_nt)->pushBackBondVec(bond_pt);
		    (*i_nt)->pushBackBondVec(bond_pt);	// i_nt is in recvPeriParticleVec, this is to calculate the deformationGradient, sigma and Kinv
							// for the peri-points in inner cell of recvPeriParticleVec, refer to commuPeriParticle()
//		    bond_pt->setIsRecv();	// bond_pt->isRecv=true;
		    recvPeriBondVec.push_back(bond_pt);

		    REAL factor = 3.0/(2.0*Pi*horizonSize_ij*horizonSize_ij*horizonSize_ij); // for the factor of 3d window function
				
		    // weighting function (influence function)
		    if(ratio < 1.0){
			bond_pt->setWeight( factor*(2.0/3.0-ratio*ratio+0.5*ratio*ratio*ratio) );
		    }
		    else{
			bond_pt->setWeight( factor*(2.0-ratio)*(2.0-ratio)*(2.0-ratio)/6.0 );
		    }		     
		} // if(ratio<2.0)

	    } // end j_nt
	} // end i_nt

	// since the calculation of PeriParticle.calcAcceleration() needs to know the deformationGradient, sigma, and Kinv of the peri-points in peri-bonds
	// even these peri-points are in recvPeriParticleVec, thus commuPeriParticle() communicate two cellSize peri-points, the deformationGradient, sigma and Kinv
	// of the inner peri-points needs to be calculated exactly, thus the peri-bonds between the recvPeriParticles should also be constructed
	for(std::vector<periDynamics::PeriParticle*>::iterator i_nt=recvPeriParticleVec.begin(); i_nt<recvPeriParticleVec.end()-1; i_nt++){
	    Vec coord0_i = (*i_nt) ->getInitPosition();
	    REAL horizonSize_i = (*i_nt)->getHorizonSize();
	    for(std::vector<periDynamics::PeriParticle*>::iterator j_nt=i_nt+1; j_nt<recvPeriParticleVec.end(); j_nt++){
		Vec coord0_j = (*j_nt)->getInitPosition();			// need to use "<" since if recvPeriParticleVec is empty, then j_nt 
		REAL tmp_length = vfabs(coord0_i-coord0_j);			// will exceed the limit of recvPeriParticleVec, then segmentational fault

		REAL horizonSize_j = (*j_nt)->getHorizonSize();
		REAL horizonSize_ij = (horizonSize_i+horizonSize_j)*0.5;//This will lead to the fact that horizion is not a sphere!!! 

		REAL ratio = tmp_length/horizonSize_ij;
		// establish the neighbor list
		if(ratio <= 2.0){ 
		    // create bond
		    periDynamics::PeriBond* bond_pt = new periDynamics::PeriBond(tmp_length, *i_nt, *j_nt);
		    (*i_nt)->pushBackBondVec(bond_pt);
		    (*j_nt)->pushBackBondVec(bond_pt);
//		    bond_pt->setIsRecv();	// bond_pt->isRecv=true;
		    recvPeriBondVec.push_back(bond_pt);

		    REAL factor = 3.0/(2.0*Pi*horizonSize_ij*horizonSize_ij*horizonSize_ij); // for the factor of 3d window function
				
		    // weighting function (influence function)
		    if(ratio < 1.0){
			bond_pt->setWeight( factor*(2.0/3.0-ratio*ratio+0.5*ratio*ratio*ratio) );
		    }
		    else{
			bond_pt->setWeight( factor*(2.0-ratio)*(2.0-ratio)*(2.0-ratio)/6.0 );
		    }		     
		} // if(ratio<2.0)

	    } // end j_nt
	} // end i_nt

    } // end findRecvPeriBonds()


    void Assembly::readPeriDynamicsData(const char *InputFile){	// should only be called by master cpu
	// readData - reads controlling parameters, particle positions and mesh connectivities
	// @param char * - reference of the input file name
    	
	std::ifstream ifs(InputFile);
	if(!ifs) {
	    std::cout << "stream error!" << std::endl; exit(-1);
	}
		
	ifs >> ndim >> nPeriParticle >> nele;
	// read particle information, create and store PeriParticle objects into periParticleVec
	for(int ip=0; ip<nPeriParticle; ip++){
	    REAL tmp_x, tmp_y, tmp_z;
	    int tmp_int;
	    ifs >> tmp_int >> tmp_x >> tmp_y >> tmp_z;
	    periDynamics::PeriParticle* tmp_pt = new periDynamics::PeriParticle(tmp_x, tmp_y, tmp_z);
	    allPeriParticleVecInitial.push_back(tmp_pt);

	}
		
	// read the connectivity information
	connectivity = new int*[nele];
	for(int iel=0; iel<nele; iel++){
		connectivity[iel] = new int[8];
	}
		
	for(int iel=0; iel<nele; iel++){
	    int tmp_int;
	    ifs >> tmp_int;
	    for(int node=0; node<8; node++){
	   	ifs >> connectivity[iel][node]; 
	    }
	}
/* // this fixedPeriParticleVec needs to be found in another way, such as by locations
// read the inside fixed peri-points, new added part
int nFixed, inode;
ifs >> nFixed;
fixedPeriParticleVec.clear();
for(int ipt=0; ipt<nFixed; ipt++){
    ifs >> inode;
    fixedPeriParticleVec.push_back(periParticleVecInitial[inode-1]);
}
*/


/*
  	std::ofstream ofs("interfacePeriPoints.dat");
  	if(!ofs) {
    	    cout << "stream error!" << endl; exit(-1);
  	}

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = interfacePeriParticleVec.begin(); pt!= interfacePeriParticleVec.end(); pt++) {

	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs.flush();
	}
	ofs.close();
*/

	point_interval = vfabs( allPeriParticleVecInitial[connectivity[0][0]]->getInitPosition() 
			      - allPeriParticleVecInitial[connectivity[0][1]]->getInitPosition() );

//	// read particle indices for linear elasticity verfication, can be deleted
//	for(int ip=0; ip < 39; ip++) {ifs >> Uzindex[ip];}
//	for(int ip=0; ip < 5; ip++) {ifs >> Uyindex[ip];}
//	for(int ip=0; ip < 5; ip++) {ifs >> Uxindex[ip];}

    } // readPeriDynamicsData()


    void Assembly::writeMesh(const char *outputFile){
/* // not used in DEM-PD coupling code, i.e. rigidInclustion()
	std::ofstream ofs(outputFile);
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "Title = \"Mesh Checking\"" << std::endl;
	ofs << "VARIABLES = \"X\", \"Y\",\"Z\"" << std::endl;
	ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT ET = BRICK" << std::endl;
	for(int node = 0; node < nPeriParticle; node++){
	    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getX()
		<< std::setw(20) << periParticleVec[node]->getInitPosition().getY() 
		<< std::setw(20) << periParticleVec[node]->getInitPosition().getZ() << std::endl;
	}
	for(int iel = 0; iel < nele; iel++){
	    for(int node = 0; node < 8; node++){
		ofs << std::setw(10) << connectivity[iel][node]; 
	    }
	    ofs << std::endl;
	}
	ofs.close();
*/
    } // end writeMesh


    void Assembly::writeMeshCheckVolume(const char *outputFile){
/* // not used in DEM-PD coupling code, i.e. rigidInclusion()
	std::ofstream ofs(outputFile);
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "Title = \"Volume Checking\"" << std::endl;
	ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"V\"" << std::endl;
	ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT ET = BRICK" << std::endl;
	// Output the coordinates and the array information
	for(std::vector<periDynamics::PeriParticle*>::iterator pt = periParticleVecInitial.begin(); pt!= periParticleVecInitial.end(); pt++) {
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() 
		<< std::setw(20) << (*pt)->getInitPosition().getZ() 
		<< std::setw(20) << (*pt)->getParticleVolume()
		<< std::endl;
	}
	for(int iel = 0; iel < nele; iel++){
	    for(int node = 0; node < 8; node++){
		ofs << std::setw(10) << connectivity[iel][node]; 
	    }
	    ofs << std::endl;
	}
	ofs.close();
*/

    } // end writeMeshCheckVolume


    void Assembly::writeParticleTecplot(std::ofstream &ofs, const int iframe) const{
	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	if(iframe == 0) {
	    ofs << "Title = \"Particle Information\"" << std::endl;
	    ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"KE\" \"P\" \"Mises\"" << std::endl;
	}
	ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	Matrix sigma;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = allPeriParticleVec.begin(); pt!= allPeriParticleVec.end(); pt++) {
	    sigma = (*pt)->getSigma();
	    pressure = sigma(1,1) + sigma(2,2) + sigma(3,3);
	    vonMisesStress = sqrt(0.5*( (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))
			         + (sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))
				 + (sigma(1,1)-sigma(3,3))*(sigma(1,1)-sigma(3,3)) )
			     + 3*( sigma(1,2)*sigma(1,2)				 + sigma(2,3)*sigma(2,3)
				 + sigma(3,1)*sigma(3,1) ));
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::setw(20) << (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getDisplacement().getY() 
		<< std::setw(20) << (*pt)->getDisplacement().getZ() 
		<< std::setw(20) << (*pt)->getVelocity().getX()
		<< std::setw(20) << (*pt)->getVelocity().getY() 
		<< std::setw(20) << (*pt)->getVelocity().getZ() 
		<< std::setw(20) << vfabs((*pt)->getVelocity())
		<< std::setw(20) << pressure
		<< std::setw(20) << vonMisesStress
		<< std::endl;
	    ofs.flush();
	}

    } // end writeparticleTecplot

// this function has some problems, allPeriParticleVecInitial does not exist anymore after the first gather
    void Assembly::printPeriDomain(const char* str) const{
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "Title = \"Particle Information\"" << std::endl;
	ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" \"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\"" << std::endl;
//	ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" << nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" << std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	Matrix sigma;
	Matrix Kinv_tmp;
	Matrix deformationG;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = periParticleVec.begin(); pt!= periParticleVec.end(); pt++) {
	    sigma = (*pt)->getSigma();
	    pressure = sigma(1,1) + sigma(2,2) + sigma(3,3);
	    vonMisesStress = sqrt(0.5*( (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))
			         + (sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))
				 + (sigma(1,1)-sigma(3,3))*(sigma(1,1)-sigma(3,3)) )
			     + 3*( sigma(1,2)*sigma(1,2)				 + sigma(2,3)*sigma(2,3)
				 + sigma(3,1)*sigma(3,1) ));
	    Kinv_tmp = (*pt)->getParticleKinv();
	    deformationG = (*pt)->getDeformationGradient();
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::setw(20) << (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getDisplacement().getY() 
		<< std::setw(20) << (*pt)->getDisplacement().getZ() 
		<< std::setw(20) << (*pt)->getVelocity().getX()
		<< std::setw(20) << (*pt)->getVelocity().getY() 
		<< std::setw(20) << (*pt)->getVelocity().getZ() 
		<< std::setw(20) << (*pt)->getAcceleration().getX()
		<< std::setw(20) << (*pt)->getAcceleration().getY() 
		<< std::setw(20) << (*pt)->getAcceleration().getZ()
		<< std::setw(20) << vfabs((*pt)->getVelocity())
		<< std::setw(20) << pressure
		<< std::setw(20) << vonMisesStress
		<< std::setw(20) << (*pt)->getParticleVolume()
		<< std::setw(20) << (*pt)->getHorizonSize()
		<< std::setw(20) << (*pt)->getBondsNumber()
		<< std::setw(20) << dem::det(Kinv_tmp)
		<< std::setw(20) << dem::det(deformationG)
		<< std::endl;
	    ofs.flush();
	}
//	for(int iel = 0; iel < nele; iel++){
//	    for(int node = 0; node < 8; node++){
//		ofs << std::setw(10) << connectivity[iel][node]; 
//	    }
//	    ofs << std::endl;
//	}

	ofs.close();

    } // end printPeriDomain


// this function has some problems, allPeriParticleVecInitial does not exist anymore after the first gather
    void Assembly::printRecvPeriDomain(const char* str) const{
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "Title = \"Particle Information\"" << std::endl;
	ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" \"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\"" << std::endl;
//	ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" << nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" << std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	Matrix sigma;
	Matrix Kinv_tmp;
	Matrix deformationG;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = recvPeriParticleVec.begin(); pt!= recvPeriParticleVec.end(); pt++) {
	    sigma = (*pt)->getSigma();
	    pressure = sigma(1,1) + sigma(2,2) + sigma(3,3);
	    vonMisesStress = sqrt(0.5*( (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))
			         + (sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))
				 + (sigma(1,1)-sigma(3,3))*(sigma(1,1)-sigma(3,3)) )
			     + 3*( sigma(1,2)*sigma(1,2)				 + sigma(2,3)*sigma(2,3)
				 + sigma(3,1)*sigma(3,1) ));
	    Kinv_tmp = (*pt)->getParticleKinv();
	    deformationG = (*pt)->getDeformationGradient();
	    ofs << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::setw(20) << (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getDisplacement().getY() 
		<< std::setw(20) << (*pt)->getDisplacement().getZ() 
		<< std::setw(20) << (*pt)->getVelocity().getX()
		<< std::setw(20) << (*pt)->getVelocity().getY() 
		<< std::setw(20) << (*pt)->getVelocity().getZ() 
		<< std::setw(20) << (*pt)->getAcceleration().getX()
		<< std::setw(20) << (*pt)->getAcceleration().getY() 
		<< std::setw(20) << (*pt)->getAcceleration().getZ()
		<< std::setw(20) << vfabs((*pt)->getVelocity())
		<< std::setw(20) << pressure
		<< std::setw(20) << vonMisesStress
		<< std::setw(20) << (*pt)->getParticleVolume()
		<< std::setw(20) << (*pt)->getHorizonSize()
		<< std::setw(20) << (*pt)->getBondsNumber()
		<< std::setw(20) << dem::det(Kinv_tmp)
		<< std::setw(20) << dem::det(deformationG)
		<< std::endl;
	    ofs.flush();
	}
//	for(int iel = 0; iel < nele; iel++){
//	    for(int node = 0; node < 8; node++){
//		ofs << std::setw(10) << connectivity[iel][node]; 
//	    }
//	    ofs << std::endl;
//	}

	ofs.close();

    } // end printRecvPeriDomain

void Assembly::printPeriParticle(const char* str) const{
  std::ofstream ofs(str);
  if(!ofs) {
    std::cout << "stream error!" << std::endl; exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  
  ofs << std::setw(OWID) << "Node ID"
      << std::setw(OWID) << "U1"
      << std::setw(OWID) << "U2"
      << std::setw(OWID) << "U3"
      << std::setw(OWID) << "S.Mises"
      << std::setw(OWID) << "S11"
      << std::setw(OWID) << "S22"
      << std::setw(OWID) << "S33"
      << std::setw(OWID) << "S12"
      << std::setw(OWID) << "S13"
      << std::setw(OWID) << "S23"
      << std::setw(OWID) << "X"
      << std::setw(OWID) << "Y"
      << std::setw(OWID) << "Z"
      << std::endl;
  
	int id = 0;
	Matrix sigma;
	REAL vonMisesStress;
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = allPeriParticleVec.begin(); pt!= allPeriParticleVec.end(); pt++) {
	    sigma = (*pt)->getSigma();
	    id++;
	    vonMisesStress = sqrt(0.5*( (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))
			         + (sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))
				 + (sigma(1,1)-sigma(3,3))*(sigma(1,1)-sigma(3,3)) )
			     + 3*( sigma(1,2)*sigma(1,2)				 + sigma(2,3)*sigma(2,3)
				 + sigma(3,1)*sigma(3,1) ));
	    ofs << std::setw(OWID) << id
		<< std::setw(OWID) << (*pt)->getDisplacement().getX()
		<< std::setw(OWID) << (*pt)->getDisplacement().getY()
		<< std::setw(OWID) << (*pt)->getDisplacement().getZ()
		<< std::setw(OWID) << vonMisesStress
		<< std::setw(OWID) << sigma(1,1)
		<< std::setw(OWID) << sigma(2,2)
		<< std::setw(OWID) << sigma(3,3)
		<< std::setw(OWID) << sigma(1,2)
		<< std::setw(OWID) << sigma(1,3)
		<< std::setw(OWID) << sigma(2,3)
		<< std::setw(OWID) << (*pt)->getInitPosition().getX()
		<< std::setw(OWID) << (*pt)->getInitPosition().getY()
		<< std::setw(OWID) << (*pt)->getInitPosition().getZ()
		<< std::endl;
	}
  
  ofs.close();
}



    void Assembly::printPeriDomainSphere(const char* str) const{
/* // not used in DEM-PD coupling code, i.e. rigidInclusion()
  	std::ofstream ofs(str);
  	if(!ofs) {
    	    std::cout << "stream error!" << endl; exit(-1);
  	}

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(10);
	ofs << "Title = \"Particle Information\"" << std::endl;
	ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Sigma_rr\" \"Sigma_tt\" \"Tao_rt\" " << std::endl;
	ofs << "ZONE T = \"periDomain\" "<< std::endl;
	// Output the coordinates and the array information
	REAL pressure, vonMisesStress;
	Matrix sigma;
	// this is not to print all the peri-points, it just print the peri-points in the interface between the sphere and the box
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = interfacePeriParticleVec.begin(); pt!= interfacePeriParticleVec.end(); pt++) {
	    sigma = (*pt)->getSigma();
	    Vec currPosition = (*pt)->getCurrPos();
	    REAL R = vfabs(currPosition);
	    REAL theta = acos(currPosition.getZ()/R);

	    // for phi should notice that x can be zero
	    REAL phi;
	    if( abs(currPosition.getX())<EPS ){
		if( currPosition.getY()>0 ){
		    phi = Pi*0.5;
		}
		else{
		    phi = Pi*1.5;
		}
	    }
	    else {
	    	phi = atan(currPosition.getY()/currPosition.getX());	
	    }

	    Matrix Qtrans(3,3);	// the transfor tensor from cartesian to spherical coordinate
				// refer to http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm 
	    Qtrans(1,1) = sin(theta)*cos(phi); Qtrans(1,2) = sin(theta)*sin(phi); Qtrans(1,3) = cos(theta);
	    Qtrans(2,1) = cos(theta)*cos(phi); Qtrans(2,2) = cos(theta)*sin(phi); Qtrans(2,3) = -sin(theta);
	    Qtrans(3,1) = -sin(phi);	       Qtrans(3,2) = cos(phi);	 	  Qtrans(3,3) = 0;

	    Matrix sigma_sphere = Qtrans*sigma*trans(Qtrans);

	    ofs << std::setw(20) << currPosition.getX()
		<< std::setw(20) << currPosition.getY()
		<< std::setw(20) << currPosition.getZ()
		<< std::setw(20) << sigma_sphere(1,1)
		<< std::setw(20) << sigma_sphere(2,2)
		<< std::setw(20) << sigma_sphere(1,2)
		<< std::endl;
	    ofs.flush();
	}

	ofs.close();
*/
    } // end printPeriDomainSphere



    void Assembly::calcDeformationGradient(){	
		// calcDeformationGradient - calculates the deformation gradient for all peri-particles
    }

    void Assembly::calcHorizonSize(){
	maxHorizonSize = 0;
	for(int iel = 0; iel < nele; iel++){
	    // get initial positions of the particles in this connectivity
 	    int n1 = connectivity[iel][0]; 	// index of first node
	    int n2 = connectivity[iel][1];
	    int n4 = connectivity[iel][3];
	    int n5 = connectivity[iel][4];

	    Vec coord1 = allPeriParticleVecInitial[n1-1]->getInitPosition();	// the index of input file is starting from 1
	    Vec coord2 = allPeriParticleVecInitial[n2-1]->getInitPosition();
	    Vec coord4 = allPeriParticleVecInitial[n4-1]->getInitPosition();
	    Vec coord5 = allPeriParticleVecInitial[n5-1]->getInitPosition();

	    REAL tmp1, tmp2, tmp3;
	    tmp1 = vfabs(coord2-coord1);
	    tmp2 = vfabs(coord4-coord1);
	    tmp3 = vfabs(coord5-coord1);

	    REAL tmpmax;	// max number in tmp1, tmp2 and tmp3
	    tmpmax = std::max( std::max(tmp1, tmp2), std::max(tmp2, tmp3) );
	    if(1.5075*tmpmax>maxHorizonSize) maxHorizonSize=1.5075*tmpmax;	    

	    for(int node = 0; node < 8; node++){
	        allPeriParticleVecInitial[connectivity[iel][node]-1]->replaceHorizonSizeIfLarger(1.5075*tmpmax);
	    }
	}
	std::cout << "maxHorizonSize: " << maxHorizonSize << std::endl;

    } // end calcHorizonSize()


    void Assembly::setInitIsv(){

    	REAL isv_tmp;
	if(dem::Parameter::getSingleton().parameter["typeConstitutive"] == 1){ // 1---implicit, 2---explicit
	    isv_tmp = dem::Parameter::getSingleton().parameter["Chi"];
	}
	else{
	    isv_tmp = dem::Parameter::getSingleton().parameter["c"];
	}

	for(std::vector<periDynamics::PeriParticle*>::iterator pt=allPeriParticleVecInitial.begin(); pt!=allPeriParticleVecInitial.end(); pt++){
	    (*pt)->setInitIsv(isv_tmp);
	}

    } // setInitIsv()


    void Assembly::calcParticleVolume() {
		
	int *numofpieces;
	REAL *particleVolume;
	REAL xi, eta, zeta;
	numofpieces = new int[nPeriParticle];
	for(int node = 0; node < nPeriParticle; numofpieces[node] = 0, node++);	// initialize
	int nip = 2;
	Matrix gp_loc3D;
	Matrix gp_weight3D;
	periDynamics::gauss3D( nip, gp_loc3D, gp_weight3D);
	particleVolume = new REAL[nPeriParticle];
	for(int node = 0; node < nPeriParticle; particleVolume[node] = 0.0, node++); // initialize
	Matrix xl(3,8);
	Matrix shp;

	for(int iel = 0; iel < nele; iel++) {

/*	    // check if this element contains the peri-point that is in the sphere (fixedPeriParticleVec)
	    int isContain = 0;	// if contains

	    // check the first node
	    for(std::vector<periDynamics::PeriParticle*>::iterator pt=fixedPeriParticleVec.begin(); pt!=fixedPeriParticleVec.end(); pt++){
		if(periParticleVecInitial[connectivity[iel][0]-1] == (*pt) || periParticleVecInitial[connectivity[iel][1]-1] == (*pt)
		|| periParticleVecInitial[connectivity[iel][2]-1] == (*pt) || periParticleVecInitial[connectivity[iel][3]-1] == (*pt)
		|| periParticleVecInitial[connectivity[iel][4]-1] == (*pt) || periParticleVecInitial[connectivity[iel][5]-1] == (*pt)
		|| periParticleVecInitial[connectivity[iel][6]-1] == (*pt) || periParticleVecInitial[connectivity[iel][7]-1] == (*pt) ){
		    isContain = 1;
		    break;
		}
	    }
	    if(isContain == 1){
		continue;
	    }
*/

	    for(int node = 0; node < 8; node++) {
		int nposition = connectivity[iel][node]-1;
		xl(1,node+1) = allPeriParticleVecInitial[nposition]->getInitPosition().getX();
		xl(2,node+1) = allPeriParticleVecInitial[nposition]->getInitPosition().getY();
		xl(3,node+1) = allPeriParticleVecInitial[nposition]->getInitPosition().getZ();
		numofpieces[nposition] += 1;
	    }
	    for(int ik = 0; ik < 8; ik++) {
		xi =   gp_loc3D(1,ik+1);
		eta =  gp_loc3D(2,ik+1);
		zeta = gp_loc3D(3,ik+1);
		// call fem function to get the volume
		REAL xsj = 0.0;
		periDynamics::shp3d(xi, eta, zeta, xl, shp, xsj);
		for(int node = 0; node < 8; node++) {
		    int nposition = connectivity[iel][node] - 1;
		    particleVolume[nposition] = particleVolume[nposition] + gp_weight3D(1,ik+1)*xsj/8.0;
		}
	    }
	}

	//Commented to compare result with the fortran code
	for(int node = 0; node < nPeriParticle; node++) {
	    particleVolume[node] = particleVolume[node]*8.0/(REAL(numofpieces[node]));
	}

	// store the particle volume into the object PeriParticle
	for(int node = 0; node < nPeriParticle; node++) {
	    allPeriParticleVecInitial[node]->setParticleVolume(particleVolume[node]);
	}

	delete[] numofpieces;
	delete[] particleVolume;
    } // end calcParticleVolume


/*    void Assembly::checkBondParticleAlive(){
	// compute the bond length and check for whether a bond or particle is alive or not
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	    if( (*pt)->getIsAlive() ){ // particle not alive, then go to next particle
		(*pt)->checkParticleAlive(dem::Parameter::getSingleton().parameter["bondStretchLimit"]);
	    }

	} // end particle

    } // end checkBondParticleAlive()
*/

    void Assembly::checkBondParticleAlive(){
	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int pt;
	num = periBondVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num) schedule(dynamic)
	for(pt=0; pt<num; pt++){
	    periBondVec[pt]->checkIfAlive();
	}

	// compute the bond length and check for whether a bond or particle is alive or not
	num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num) schedule(dynamic)
	for(pt=0; pt<num; pt++){
	    if( periParticleVec[pt]->getIsAlive() ){ // particle not alive, then go to next particle
		periParticleVec[pt]->checkParticleAlive();
	    }
	} // end particle

    } // end checkBondParticleAlive()

    void Assembly::calcParticleKinv(){ // this should be called after commuDEMPeriParticle(), and find peri-bonds between communicated periParticles
	// Compute the inverse of the shape tensor K
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++) {
	    (*pt)->calcParticleKinv();
	} 
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++) {
	    (*pt)->calcParticleKinv();
	} 

    } // end calcParticleKinv()	


    void Assembly::calcParticleStress(){
	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int i;
	num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num) schedule(dynamic)
	for(i=0; i<num; i++){
	    periParticleVec[i]->calcParticleStress();
	} 
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++){
	    (*pt)->calcParticleStress();
    	} 
    } // calcParticleStress()


    void Assembly::calcParticleAcceleration(){
	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int i;
	num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num) schedule(dynamic)
	for(i=0; i<num; i++){
	    periParticleVec[i]->calcParticleAcceleration();
	}
    } // calcParticleAcceleration()


    void Assembly::calcRecvParticleKinv(){ // this should be called after commuDEMPeriParticle(), and find peri-bonds between communicated periParticles
	// Compute the inverse of the shape tensor K
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++) {
	    (*pt)->calcParticleKinv();
	} 

    } // end calcParticleKinv()	

    void Assembly::ApplyExternalForce(int istep) {
	// deal with the external force, applied at the top of the boundary
	REAL factor = 0.0;
	REAL rampStep = dem::Parameter::getSingleton().parameter["rampStep"];
	if(istep <= rampStep) {
	    factor = REAL(istep)/REAL(rampStep);
	}
	else {
	    factor = REAL(1.0);
	}
//	for(std::vector<periDynamics::PeriParticle*>::iterator pt=topBoundaryVec.begin(); pt!=topBoundaryVec.end(); pt++){
//	    (*pt)->setAcceleration(Vec(0.0,0.0,factor*dem::Parameter::getSingleton().parameter["bodyDensity"]));
//	}
    } // end ApplyExternalForce


    // delete those peri-points that are inside sand particles
    // actually what are deleted are the peri-points that are inside
    // a enlarged sand particle to avoid peri-points go to inside the 
    // sand after a few steps. Here the enlarged sand particle is a little
    // larger than the sand particle, the delta = 0.5*point_interval and 
    // delta is the difference of principle lengths
    void Assembly::removeInsidePeriParticles(){

	bool is_inside;
	allPeriParticleVec.clear();
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=allPeriParticleVecInitial.begin(); pt!=allPeriParticleVecInitial.end(); pt++){
	    Vec coord_pt = (*pt)->getCurrPosition();
	    is_inside = false;	// if this peri-point is inside sand particles
//	    if(coord_pt.getX()==0 && coord_pt.getY()==0 && coord_pt.getZ()==1.5){
//		is_inside = true;
//	    }

	    // remove the inside peri-points that are in the box mesh
	    for(ParticlePArray::iterator dem_pt=allParticleVec.begin(); dem_pt!=allParticleVec.end(); dem_pt++){
		REAL a = (*dem_pt)->getA()+0.5*point_interval;	// enlarged sand particle
		REAL b = (*dem_pt)->getB()+0.5*point_interval;
		REAL c = (*dem_pt)->getC()+0.5*point_interval;

		Vec coord_pt_tmp = (*dem_pt)->globalToLocal(coord_pt-(*dem_pt)->getCurrPos());	// this is important to get local coordinate

		REAL x_pt = coord_pt_tmp.getX();
		REAL y_pt = coord_pt_tmp.getY();
		REAL z_pt = coord_pt_tmp.getZ();

		if( (x_pt/a)*(x_pt/a) + (y_pt/b)*(y_pt/b) + (z_pt/c)*(z_pt/c) < 1 ){ // this peri-points is inside sand particle
		    is_inside = true;
		    break;	// do not need to check other sand particles
		}
	    }

	    if(is_inside == false){ // this peri-point is not within any sand particles
		allPeriParticleVec.push_back(*pt);
	    }
	}



/*
// remove the inside peri-points that are in the box mesh
//	bool is_inside;
//	periParticleVec.clear();
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVecInitial.begin(); pt!=periParticleVecInitial.end(); pt++){
	    Vec coord_pt = (*pt)->getCurrPos();
	    is_inside = false;	// if this peri-point is inside sand particles
	    for(std::vector<particle*>::iterator dem_pt=ParticleVec.begin(); dem_pt!=ParticleVec.end(); dem_pt++){
//		REAL a = (*dem_pt)->getA()+0.5*point_interval;	// enlarged sand particle
//		REAL b = (*dem_pt)->getB()+0.5*point_interval;
//		REAL c = (*dem_pt)->getC()+0.5*point_interval;
		
		// since we know the model, the sphere's radius is 0.1m, then we can delete points within radius as 0.095m
		REAL a = 0.099;	// enlarged sand particle
		REAL b = 0.099;
		REAL c = 0.099;

		coord_pt = (*dem_pt)->localVec(coord_pt-(*dem_pt)->getCurrPos());	// this is important to get local coordinate

		REAL x_pt = coord_pt.getX();
		REAL y_pt = coord_pt.getY();
		REAL z_pt = coord_pt.getZ();

		if( (x_pt/a)*(x_pt/a) + (y_pt/b)*(y_pt/b) + (z_pt/c)*(z_pt/c) < 1 ){ // this peri-points is inside sand particle
		    is_inside = true;
		    break;	// do not need to check other sand particles
		}
	    }


	    if(is_inside == false){ // this peri-point is not within any sand particles
		periParticleVec.push_back(*pt);
	    }
	}
*/


    } // end removeInsidePeriParticles()

   void Assembly::clearPeriDEMBonds(){
	for(std::vector<PeriDEMBond*>::iterator bt=periDEMBondVec.begin(); bt!=periDEMBondVec.end(); bt++){
	    delete (*bt);
	}
	periDEMBondVec.clear();
        std::vector<PeriDEMBond*>().swap(periDEMBondVec); // actual memory release
	plan_gravity = dem::Parameter::getSingleton().parameter["periDensity"]*point_interval*point_interval*9.8;	// rho*l^2*g, used in PeriDEMBond.cpp
   }

   void Assembly::eraseBrokenPeriDEMBonds(){
	for(std::vector<PeriDEMBond*>::iterator bt=periDEMBondVec.begin(); bt!=periDEMBondVec.end(); ){
	    if( (*bt)->getIsAlive() ){	// still alive
		bt++;
	    }
	    else{
	      delete (*bt);
	      bt = periDEMBondVec.erase(bt);	// release memory
	    }
	}
   }

   // if a peri-point is first time (i.e. previous step peri-point is not bonded to dem particle) entering the influence domain of a DEM particle, 
   // then construct peri-DEM bond between this peri-point and this DEM particle
   void Assembly::findPeriDEMBonds(){

	eraseBrokenPeriDEMBonds();
	// construct sand peri-bonds
	// here mergePeriParticleVec is used, since if we have a DEM particle that is crossing the boundary between two cpus, then sand-peri bonds between 
	// the communicated peri-points can provide force on the DEM particle from the peri-points in neighboring cpus
	// mergeParticleVec is used, since if a DEM particle is crossing the boundary between two cpus, then these sand-peri bonds between the communicated 
	// DEM particles can provide force on peri-points from the DEM particles in neighboring cpus
	REAL delta = point_interval*3.2;	// 3 layers of peri-points

	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int peri_pt;
	num = mergePeriParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num, delta) schedule(dynamic)
	for(peri_pt=0; peri_pt<num; peri_pt++){
	    Vec xyz_peri = mergePeriParticleVec[peri_pt]->getCurrPosition();

	    for(auto dem_pt=mergeParticleVec.begin(); dem_pt!=mergeParticleVec.end(); dem_pt++){
		// check and construct the periDEMBondVec in this particle
		REAL ra = (*dem_pt)->getA();
		REAL rb = (*dem_pt)->getB();
		REAL rc = (*dem_pt)->getC();
	    	Vec xyz_peri_tmp = (*dem_pt)->globalToLocal(xyz_peri-(*dem_pt)->getCurrPos());	// this is very important, since all calculations below for ellipsoid
	    	REAL x_peri = xyz_peri_tmp.getX();						// are based on the local coordinate of the ellipsoid
	    	REAL y_peri = xyz_peri_tmp.getY();
	    	REAL z_peri = xyz_peri_tmp.getZ();
	    	REAL x_peri_de = x_peri/(ra+delta);						// are based on the local coordinate of the ellipsoid
	    	REAL y_peri_de = y_peri/(rb+delta);
	    	REAL z_peri_de = z_peri/(rc+delta);
		if( x_peri_de*x_peri_de + y_peri_de*y_peri_de + z_peri_de*z_peri_de < 1 ){
		    // check if (*peri_pt) is first time entering the influence domain of (*dem_pt)
		    if( mergePeriParticleVec[peri_pt]->isBonded((*dem_pt)->getId()) ) continue;	// already bonded		

	    	    // calcualte the local coordinates of the projector on the surface of the ellipsoid particle
	    	    Vec projector_local;

	    	    // get the projector of this peri-point, *peri_pt, on the surface of this ellipsoid, *dem_pt

		    // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in the same octant
		    REAL kd = sqrt(ra*ra*rb*rb*rc*rc/(x_peri*x_peri*rb*rb*rc*rc+y_peri*y_peri*ra*ra*rc*rc+z_peri*z_peri*ra*ra*rb*rb));
		    REAL dd = sqrt(x_peri*x_peri+y_peri*y_peri+z_peri*z_peri)*(1-kd);

		    REAL n1 = 2*x_peri/(ra+dd)/(ra+dd);
		    REAL n2 = 2*y_peri/(rb+dd)/(rb+dd);
		    REAL n3 = 2*z_peri/(rc+dd)/(rc+dd);		

		    REAL aa = n1*n1*rb*rb*rc*rc+n2*n2*ra*ra*rc*rc+n3*n3*ra*ra*rb*rb;
		    REAL bb = 2*(n1*x_peri*rb*rb*rc*rc+n2*y_peri*ra*ra*rc*rc+n3*z_peri*ra*ra*rb*rb);
		    REAL cc = rb*rb*rc*rc*x_peri*x_peri+ra*ra*rc*rc*y_peri*y_peri+ra*ra*rb*rb*z_peri*z_peri-ra*ra*rb*rb*rc*rc;

		    REAL kn = (-bb+sqrt(bb*bb-4*aa*cc))*0.5/aa;

//std::cout << "kn: " << kn << std::endl;


		    projector_local.setX(x_peri+kn*n1);
		    projector_local.setY(y_peri+kn*n2);
		    projector_local.setZ(z_peri+kn*n3);


//		    // this is to use the cross point instead of the normal point as the projector point
//		    projector_local.setX(x_peri*kd);
//		    projector_local.setY(y_peri*kd);
//		    projector_local.setZ(z_peri*kd);

		    PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, (*dem_pt).get(), mergePeriParticleVec[peri_pt]);

//		    // this is used to test the coupled force model, October 10, 2014
//		    // in this test model, the sand-peri-points will move along the dem-particle
//		    PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, *dem_pt, *peri_pt);

#pragma omp critical
		    periDEMBondVec.push_back(bond_tmp);

		}
	    }
	}

   } // findPeriDEMBonds()


    // construct boundary bonds, sand bonds, July 15, 2014
    void Assembly::constructBoundarySandPeriBonds(){
/*
	// delete previous bonds vector
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=topBoundaryBondVec.begin(); pt!=topBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=bottomBoundaryBondVec.begin(); pt!=bottomBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=leftBoundaryBondVec.begin(); pt!=leftBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=rightBoundaryBondVec.begin(); pt!=rightBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=frontBoundaryBondVec.begin(); pt!=frontBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}
    	for(std::vector<dem::PeriBoundaryBond*>::iterator pt=backBoundaryBondVec.begin(); pt!=backBoundaryBondVec.end(); pt++){
	    delete (*pt);
    	}

    	topBoundaryBondVec.clear();
    	bottomBoundaryBondVec.clear();
    	leftBoundaryBondVec.clear();
    	rightBoundaryBondVec.clear();
    	frontBoundaryBondVec.clear();
    	backBoundaryBondVec.clear();

    	for(std::vector<dem::PeriDEMBond*>::iterator pt=periDEMBondVec.begin(); pt!=periDEMBondVec.end(); pt++){
	    delete (*pt);
    	}
    	periDEMBondVec.clear();


	REAL delta = point_interval*3;	// 3 layers of peri-points

	// boundary coordinates
	REAL x_min = getApt(3).getX();
 	REAL x_max = getApt(1).getX();
	REAL y_min = getApt(4).getY();
	REAL y_max = getApt(2).getY();
	REAL z_min = getApt(6).getZ();
	REAL z_max = getApt(5).getZ();

// no boundaries in pile penetration
//	// construct boundary peri-bonds
//	PeriBoundaryBond* bond_tmp;
//	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){ // overlap over peri-points
//	    Vec xyz_pt = (*pt)->getCurrPos();
//	    REAL x_pt = xyz_pt.getX();
//	    REAL y_pt = xyz_pt.getY();
//	    REAL z_pt = xyz_pt.getZ();
//
//	    if( x_pt-x_min<delta ){	// back boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_min,y_pt,z_pt), (*pt));
//		backBoundaryBondVec.push_back(bond_tmp);
//	    }
//	    if( x_max-x_pt<delta ){	// front boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_max,y_pt,z_pt), (*pt));
//		frontBoundaryBondVec.push_back(bond_tmp);
//	    }
//	    if( y_pt-y_min<delta ){	// left boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_min,z_pt), (*pt));
//		leftBoundaryBondVec.push_back(bond_tmp);
//	    }
//	    if( y_max-y_pt<delta ){	// right boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_max,z_pt), (*pt));
//		rightBoundaryBondVec.push_back(bond_tmp);
//	    }
//	    if( z_pt-z_min<delta ){	// bottom boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_min), (*pt));
//		bottomBoundaryBondVec.push_back(bond_tmp);
//	    }
//	    if( z_max-z_pt<delta ){	// top boundary points
//		bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_max), (*pt));
//		topBoundaryBondVec.push_back(bond_tmp);
//	    }
//
//	} 
//

	// construct sand peri-bonds
	for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=interfacePeriParticleVec.begin(); peri_pt!=interfacePeriParticleVec.end(); peri_pt++){
	    Vec xyz_peri = (*peri_pt)->getCurrPos();

	    for(ParticlePArray::iterator dem_pt=ParticleVec.begin(); dem_pt!=ParticleVec.end(); dem_pt++){
		// check and construct the periDEMBondVec in this particle
		REAL ra = (*dem_pt)->getA();
		REAL rb = (*dem_pt)->getB();
		REAL rc = (*dem_pt)->getC();
	    	xyz_peri = (*dem_pt)->localVec(xyz_peri-(*dem_pt)->getCurrPos());	// this is very important, since all calculations below for ellipsoid
	    	REAL x_peri = xyz_peri.getX();						// are based on the local coordinate of the ellipsoid
	    	REAL y_peri = xyz_peri.getY();
	    	REAL z_peri = xyz_peri.getZ();
//		if( (x_peri/(ra+delta))*(x_peri/(ra+delta))+(y_peri/(rb+delta))*(y_peri/(rb+delta))+(z_peri/(rc+delta))*(z_peri/(rc+delta)) < 1 ){
	    	    // calcualte the local coordinates of the projector on the surface of the ellipsoid particle
	    	    Vec projector_local;

	    	    // get the projector of this peri-point, *peri_pt, on the surface of this ellipsoid, *dem_pt

		    // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in the same octant
		    REAL kd = sqrt(ra*ra*rb*rb*rc*rc/(x_peri*x_peri*rb*rb*rc*rc+y_peri*y_peri*ra*ra*rc*rc+z_peri*z_peri*ra*ra*rb*rb));
//		    REAL dd = sqrt(x_peri*x_peri+y_peri*y_peri+z_peri*z_peri)*(1-kd);
//
//		    REAL n1 = 2*x_peri/(ra+dd)/(ra+dd);
//		    REAL n2 = 2*y_peri/(rb+dd)/(rb+dd);
//		    REAL n3 = 2*z_peri/(rc+dd)/(rc+dd);		
//
//		    REAL aa = n1*n1*rb*rb*rc*rc+n2*n2*ra*ra*rc*rc+n3*n3*ra*ra*rb*rb;
//		    REAL bb = 2*(n1*x_peri*rb*rb*rc*rc+n2*y_peri*ra*ra*rc*rc+n3*z_peri*ra*ra*rb*rb);
//		    REAL cc = rb*rb*rc*rc*x_peri*x_peri+ra*ra*rc*rc*y_peri*y_peri+ra*ra*rb*rb*z_peri*z_peri-ra*ra*rb*rb*rc*rc;
//
//		    REAL kn = (-bb+sqrt(bb*bb-4*aa*cc))*0.5/aa;
//
//std::cout << "kn: " << kn << std::endl;
//
//
//		    projector_local.setX(x_peri+kn*n1);
//		    projector_local.setY(y_peri+kn*n2);
//		    projector_local.setZ(z_peri+kn*n3);
//
//
		    // this is to use the cross point instead of the normal point as the projector point
		    projector_local.setX(x_peri*kd);
		    projector_local.setY(y_peri*kd);
		    projector_local.setZ(z_peri*kd);

		    PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, *dem_pt, *peri_pt);

//		    // this is used to test the coupled force model, October 10, 2014
//		    // in this test model, the sand-peri-points will move along the dem-particle
//		    PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, *dem_pt, *peri_pt);

		    periDEMBondVec.push_back(bond_tmp);

//		}
	    }
	}
	
*/
    } // end constructBoundarySandPeriBonds


    void Assembly::applyCoupledForces(){

/* no boundaries in the pile penetration
	// calculate the coupled forces between boundaries
	// calculate current length of the bonds
	// boundary coordinates
	REAL x_min = getApt(3).getX();
 	REAL x_max = getApt(1).getX();
	REAL y_min = getApt(4).getY();
	REAL y_max = getApt(2).getY();
	REAL z_min = getApt(6).getZ();
	REAL z_max = getApt(5).getZ();
	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=bottomBoundaryBondVec.begin(); bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(z_min, 6);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=topBoundaryBondVec.begin(); bond_pt!=topBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(z_max, 5);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=backBoundaryBondVec.begin(); bond_pt!=backBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(x_min, 3);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=frontBoundaryBondVec.begin(); bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(x_max, 1);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=leftBoundaryBondVec.begin(); bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(y_min, 4);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=rightBoundaryBondVec.begin(); bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(y_max, 2);	
		// pass the boundary coordinate and the number of the boundary
	}
*/
	// calculate coupled forces between sand particles and peri-points
	for(std::vector<PeriDEMBond*>::iterator bond_pt=periDEMBondVec.begin(); bond_pt!=periDEMBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce();
	}
    } // end applyCoupledForces()



// this is used to test the coupling model. In this test model, the sand-peri-points will move along the 
// dem-particles. And the peri-points have no influence on the dem-particles. October 10, 2014
    void Assembly::applyCoupledBoundary(){

/* no boundaries in the pile penetration
	// calculate the coupled forces between boundaries
	// calculate current length of the bonds
	// boundary coordinates
	REAL x_min = getApt(3).getX();
 	REAL x_max = getApt(1).getX();
	REAL y_min = getApt(4).getY();
	REAL y_max = getApt(2).getY();
	REAL z_min = getApt(6).getZ();
	REAL z_max = getApt(5).getZ();
	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=bottomBoundaryBondVec.begin(); bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(z_min, 6);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=topBoundaryBondVec.begin(); bond_pt!=topBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(z_max, 5);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=backBoundaryBondVec.begin(); bond_pt!=backBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(x_min, 3);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=frontBoundaryBondVec.begin(); bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(x_max, 1);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=leftBoundaryBondVec.begin(); bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(y_min, 4);	
		// pass the boundary coordinate and the number of the boundary
	}

	for(std::vector<PeriBoundaryBond*>::iterator bond_pt=rightBoundaryBondVec.begin(); bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondForce(y_max, 2);	
		// pass the boundary coordinate and the number of the boundary
	}
*/

	// calculate coupled forces between sand particles and peri-points
	for(std::vector<PeriDEMBond*>::iterator bond_pt=periDEMBondVec.begin(); bond_pt!=periDEMBondVec.end(); bond_pt++){
	    (*bond_pt)->applyBondBoundary();
	}


    } // end applyCoupledBoundary()


  void Assembly::constructPeriMatrix(){	// construct the Matrix members in periParticleVec, 
					// since currently the transfer of pointer array in Matrix is not implemented well
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	    (*pt)->constructMatrixMember();	// construct these Matrix members
	}

  } // constructPeriMatrix()


  void Assembly::constructRecvPeriMatrix(){	// construct the Matrix members in periParticleVec, 
					// since currently the transfer of pointer array in Matrix is not implemented well
	for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++){
	    (*pt)->constructMatrixMember();	// construct these Matrix members
	}

  } // constructRecvPeriMatrix()


  void Assembly::findBoundaryPeriParticles(){	// it does not matter since this will be only called once
      bottomBoundaryVec.clear();
      frontBoundaryVec.clear();
      leftBoundaryVec.clear();
      topBoundaryInnerVec.clear();
      topBoundaryEdgeVec.clear();
      topBoundaryCornerVec.clear();

      REAL x1 = dem::Parameter::getSingleton().parameter["Xmin"]; REAL x2 = dem::Parameter::getSingleton().parameter["Xmax"];
      REAL y1 = dem::Parameter::getSingleton().parameter["Ymin"]; REAL y2 = dem::Parameter::getSingleton().parameter["Ymax"];
      REAL z1 = dem::Parameter::getSingleton().parameter["Zmin"]; REAL z2 = dem::Parameter::getSingleton().parameter["Zmax"];
    for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=periParticleVec.begin(); peri_pt!=periParticleVec.end(); peri_pt++){
	REAL x_peri = (*peri_pt)->getInitPosition().getX();
	REAL y_peri = (*peri_pt)->getInitPosition().getY();
	REAL z_peri = (*peri_pt)->getInitPosition().getZ();
	if(z_peri==z2){	// top boundary
//	    if( (x_peri==x2&&y_peri==y2) || x_peri==x1&&y_peri==y2 || x_peri==x2&&y_peri==y1 || x_peri==x1&&y_peri==y1){	// corner points
//		topBoundaryCornerVec.push_back(*peri_pt);
//	    }
//	    else if(x_peri==x1 || y_peri==y1 || x_peri==x2 || y_peri==y2){	// edge points
//		topBoundaryEdgeVec.push_back(*peri_pt);
//	    }
//	    else {	// inner points
//	    	topBoundaryInnerVec.push_back(*peri_pt);
//	    }
	}
	else if(z_peri<z1+4.5*point_interval){	// bottom boundary
	    bottomBoundaryVec.push_back(*peri_pt);
	}

	if(x_peri<x1+4.5*point_interval || x_peri>x2-4.5*point_interval){
	    leftBoundaryVec.push_back(*peri_pt);
	}
	if(y_peri<y1+4.5*point_interval || y_peri>y2-4.5*point_interval){
	    frontBoundaryVec.push_back(*peri_pt);
	}
    }
	
  } // findBoundaryPeriParticles()


  void Assembly::findFixedPeriParticles(){	// it does not matter since this will be only called once
    fixedPeriParticleVec.clear();
    REAL radius = dem::Parameter::getSingleton().parameter["fixRadius"];
    radius = radius+point_interval*0.2;	// in order to find all points on the spherical surface
    REAL x0 = dem::Parameter::getSingleton().parameter["periFixCentroidX"];
    REAL y0 = dem::Parameter::getSingleton().parameter["periFixCentroidY"];
    REAL z0 = dem::Parameter::getSingleton().parameter["periFixCentroidZ"];
    for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=periParticleVec.begin(); peri_pt!=periParticleVec.end(); peri_pt++){
	REAL x_peri = (*peri_pt)->getInitPosition().getX();
	REAL y_peri = (*peri_pt)->getInitPosition().getY();
	REAL z_peri = (*peri_pt)->getInitPosition().getZ();
	if( (x_peri-x0)*(x_peri-x0)+(y_peri-y0)*(y_peri-y0)+(z_peri-z0)*(z_peri-z0)<radius*radius+EPS ){
	    fixedPeriParticleVec.push_back(*peri_pt);
	}
    }
	
  } // findFixedPeriParticles()

  void Assembly::applyPeriBoundaryCondition(){
    dem::Vec zero_vec = dem::Vec(0,0,0);
//    for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=fixedPeriParticleVec.begin(); peri_pt!=fixedPeriParticleVec.end(); peri_pt++){
//        (*peri_pt)->prescribeDisplacement(zero_vec);
//    }

	int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];
	int num;	// number of peri-points
	int peri_pt;
	num = bottomBoundaryVec.size();
//#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num) schedule(dynamic)
//    for(peri_pt=0; peri_pt<num; peri_pt++){
//        bottomBoundaryVec[peri_pt]->prescribeBottomDisplacement(0);
//    }
	num = leftBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num) schedule(dynamic)
    for(peri_pt=0; peri_pt<num; peri_pt++){
        leftBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
    }
	num = frontBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num) schedule(dynamic)
    for(peri_pt=0; peri_pt<num; peri_pt++){
        frontBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
    }

  } // applyPeriBoundaryCondition()

  void Assembly::applyTractionBoundary(int g_iteration){
	// set traction boundary
	REAL force = dem::Parameter::getSingleton().parameter["periForce"];
	REAL rampStep = dem::Parameter::getSingleton().parameter["rampStep"];
	REAL framp = 0;
	if(g_iteration <= rampStep){
	    framp = g_iteration/rampStep;
	}
	else {
	    framp = 1;
	}
	framp = framp*force;
 	dem::Vec tmp_vec = dem::Vec(0,0,framp);
	for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=topBoundaryInnerVec.begin(); peri_pt!=topBoundaryInnerVec.end(); peri_pt++){
	    (*peri_pt)->addAccelerationByForce(tmp_vec);
	}
	for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=topBoundaryEdgeVec.begin(); peri_pt!=topBoundaryEdgeVec.end(); peri_pt++){
	    (*peri_pt)->addAccelerationByForce(tmp_vec*0.5);
	}
	for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=topBoundaryCornerVec.begin(); peri_pt!=topBoundaryCornerVec.end(); peri_pt++){
	    (*peri_pt)->addAccelerationByForce(tmp_vec*0.25);
	}
	for(std::vector<periDynamics::PeriParticle*>::iterator peri_pt=bottomBoundaryVec.begin(); peri_pt!=bottomBoundaryVec.end(); peri_pt++){
	    (*peri_pt)->addAccelerationByForce(-tmp_vec);
	}

  } // applyTractionBoundary()

  void Assembly::rigidInclusion() 
  { 
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      distX = dem::Parameter::getSingleton().parameter["Xmax"] - dem::Parameter::getSingleton().parameter["Xmin"];
      distY = dem::Parameter::getSingleton().parameter["Ymax"] - dem::Parameter::getSingleton().parameter["Ymin"];
      distZ = dem::Parameter::getSingleton().parameter["Zmax"] - dem::Parameter::getSingleton().parameter["Zmin"];
      REAL x1 = dem::Parameter::getSingleton().parameter["Xmin"];  REAL x2 = dem::Parameter::getSingleton().parameter["Xmax"];
      REAL y1 = dem::Parameter::getSingleton().parameter["Ymin"];  REAL y2 = dem::Parameter::getSingleton().parameter["Ymax"];
      REAL z1 = dem::Parameter::getSingleton().parameter["Zmin"];  REAL z2 = dem::Parameter::getSingleton().parameter["Zmax"];
      setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
      setGrid(Rectangle(x1, y1, z1, x2, y2, z2));	// compute grid assumed to be the same as container, change in scatterParticle() if necessary.
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      readPeriDynamicsData(dem::Parameter::getSingleton().datafile["periFile"].c_str());
      openCompressProg(progressInf, "rigidInc_progress");
      openPeriProgress(periProgInf, "rigidInc_peri.dat");
      calcParticleVolume();		// volume and horizon size are related to the mesh, July 14, 2014
      calcHorizonSize();	 	// so these two should be calculated before peri-points that are within sand particles are deleted
      removeInsidePeriParticles();	// delete those peri-points that are inside sand particles
    }
    scatterDEMPeriParticle();
    constructPeriMatrix();	// construct the Matrix members in periParticleVec
    constructNeighbor();

//if(mpiRank==0){
//printPeriDomain("peridomain_rank_0.dat");
//printRecvPeriDomain("peridomain_recv_rank_0.dat");	// bondVec is not empty
//}
//if(mpiRank==1){
//printPeriDomain("peridomain_rank_1.dat");
//printRecvPeriDomain("peridomain_recv_rank_1.dat");	// bondVec is empty
//}

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1; 
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    //REAL time0, time1, time2, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    if (mpiRank == 0) {
      plotGrid(strcat(combineString(cstr0, "rigidInc_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "rigidInc_particle_", iterSnap - 1, 3));
      printPeriProgress(periProgInf, 0); 
    }
//    if (mpiRank == 0)
//      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
//	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;

    commuParticle();		// the commuPeriParticle() have to be called after commuParticle(), since 
    commuPeriParticle();	// rankX1... are calculated in commuParticle()
    constructRecvPeriMatrix();	// construct the Matrix members in recvPeriParticle

//if(mpiRank==0){
//printPeriDomain("peridomain_rank_0.dat");
//printRecvPeriDomain("peridomain_recv_rank_0.dat");	// bondVec is not empty
//}
//if(mpiRank==1){
//printPeriDomain("peridomain_rank_1.dat");
//printRecvPeriDomain("peridomain_recv_rank_1.dat");	// bondVec is empty
//}

    findRecvPeriBonds();	// find the peri-bonds between periParticleVec and recvPeriParticle
    calcParticleKinv();
    for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	(*pt)->initial();
    }
    for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++){
	(*pt)->initial();
    }
    calcParticleStress();
    calcParticleAcceleration();
    releaseRecvParticle();
    releaseRecvPeriParticle();
    releasePeriBondVec();	// free the bondVec in periParticleVec and bondVec.clear() before constructNeighbor() again

/*
  	std::ofstream ofs_top("topBoundaryInnerPeriPoints.dat");
  	if(!ofs_top) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs_top.setf(std::ios::scientific, std::ios::floatfield);
	ofs_top.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = topBoundaryInnerVec.begin(); pt!= topBoundaryInnerVec.end(); pt++) {

	    ofs_top << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs_top.flush();
	}
	ofs_top.close();

  	std::ofstream ofs_bot("topBoundaryEdgePeriPoints.dat");
  	if(!ofs_bot) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs_bot.setf(std::ios::scientific, std::ios::floatfield);
	ofs_bot.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = topBoundaryEdgeVec.begin(); pt!= topBoundaryEdgeVec.end(); pt++) {

	    ofs_bot << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs_bot.flush();
	}
	ofs_bot.close();


 	std::ofstream ofs_fix("topBoundaryCornerPeriPoints.dat");
  	if(!ofs_fix) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs_fix.setf(std::ios::scientific, std::ios::floatfield);
	ofs_fix.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = topBoundaryCornerVec.begin(); pt!= topBoundaryCornerVec.end(); pt++) {

	    ofs_fix << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs_fix.flush();
	}
	ofs_fix.close();


 	std::ofstream ofs2("bottomBoundaryPoints.dat");
  	if(!ofs2) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs2.setf(std::ios::scientific, std::ios::floatfield);
	ofs2.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = bottomBoundaryVec.begin(); pt!= bottomBoundaryVec.end(); pt++) {

	    ofs2<< std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs2.flush();
	}
	ofs2.close();


 	std::ofstream ofs3("leftxBoundaryPoints.dat");
  	if(!ofs3) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs3.setf(std::ios::scientific, std::ios::floatfield);
	ofs3.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = leftBoundaryVec.begin(); pt!= leftBoundaryVec.end(); pt++) {

	    ofs3<< std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs3.flush();
	}
	ofs3.close();


 	std::ofstream ofs4("frontyBoundaryPoints.dat");
  	if(!ofs4) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs4.setf(std::ios::scientific, std::ios::floatfield);
	ofs4.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = frontBoundaryVec.begin(); pt!= frontBoundaryVec.end(); pt++) {

	    ofs4<< std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs4.flush();
	}
	ofs4.close();

if(mpiRank==0){
 	std::ofstream ofs5("fixedPeriPoints_rank_0.dat");
  	if(!ofs5) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs5.setf(std::ios::scientific, std::ios::floatfield);
	ofs5.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = fixedPeriParticleVec.begin(); pt!= fixedPeriParticleVec.end(); pt++) {

	    ofs5<< std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs5.flush();
	}
	ofs5.close();
}
if(mpiRank==1){
 	std::ofstream ofs5("fixedPeriPoints_rank_1.dat");
  	if(!ofs5) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs5.setf(std::ios::scientific, std::ios::floatfield);
	ofs5.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = fixedPeriParticleVec.begin(); pt!= fixedPeriParticleVec.end(); pt++) {

	    ofs5<< std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs5.flush();
	}
	ofs5.close();
}

  	std::ofstream ofs6("periParticleVec.dat");
  	if(!ofs6) {
    	    std::cout << "stream error!" << std::endl; exit(-1);
  	}

	ofs6.setf(std::ios::scientific, std::ios::floatfield);
	ofs6.precision(10);
	for(std::vector<periDynamics::PeriParticle*>::const_iterator pt = periParticleVec.begin(); pt!= periParticleVec.end(); pt++) {

	    ofs6 << std::setw(20) << (*pt)->getInitPosition().getX() + (*pt)->getDisplacement().getX()
		<< std::setw(20) << (*pt)->getInitPosition().getY() + (*pt)->getDisplacement().getY()
		<< std::setw(20) << (*pt)->getInitPosition().getZ() + (*pt)->getDisplacement().getZ()
		<< std::endl;
	    ofs6.flush();
	}
	ofs6.close();
*/
    while (iteration <= endStep) {
      findBoundaryPeriParticles();	// boundary peri-points are determined by their initial positions, however they are still needed to be determined
      findFixedPeriParticles();		// every step since some peri-points can pass from one cpu to another
      //commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
      commuParticle(); 
      //time2 = MPI_Wtime(); commuT = time2 - time0;

      // displacement control relies on constant time step, so do not call calcTimeStep().
      //calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();
      runFirstHalfStep();

      applyPeriBoundaryCondition();
      commuPeriParticle();
      constructRecvPeriMatrix();
      constructNeighbor();
      findRecvPeriBonds();

      calcParticleStress();
      calcParticleAcceleration();	// peri-point acceleration is set zero at the beginning

      applyTractionBoundary(iteration-startStep);
      runSecondHalfStep();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
//      updateBoundary(sigmaConf, "triaxial");
//      updateGrid();
      updatePeriGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	//time1 = MPI_Wtime();
	gatherParticle();
	gatherPeriParticle();
	gatherEnergy(); 
        //time2 = MPI_Wtime(); gatherT = time2 - time1;

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "rigidInc_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "rigidInc_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "rigidInc_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "rigidInc_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "rigidInc_boundary_", iterSnap, 3));
	  //printCompressProg(progressInf, distX, distY, distZ); // redundant
	  printPeriProgress(periProgInf, iterSnap);
	}
	printContact(combineString(cstr, "rigidInc_contact_", iterSnap, 3));      
	++iterSnap;
      }
      releaseRecvParticle(); // late release because printContact refers to received particles
      releaseRecvPeriParticle();	// release recvPeriParticleVec and also remove peri-bonds between recvPeriParticle
      releasePeriBondVec();
      //time1 = MPI_Wtime();
      migrateParticle();
      migratePeriParticle();
      //time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
//      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
//	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
//		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (mpiRank == 0 && iteration % 10 == 0)
	printCompressProg(progressInf, distX, distY, distZ);

      // no break condition, just through top/bottom displacement control
      ++iteration;
    } 

    if (mpiRank == 0) {
      printParticle("rigidInc_particle_end");
      printBdryContact("rigidInc_bdrycntc_end");
      printBoundary("rigidInc_boundary_end");
      printCompressProg(progressInf, distX, distY, distZ);
//      printPeriProgress(periProgInf, iterSnap);
    }
    if (mpiRank == 0){
      closeProg(progressInf);
      closeProg(periProgInf);
    }

  } // rigidInclusion


  void Assembly::pullOutPeri() 
  { 
    REAL distX, distY, distZ;
    if (mpiRank == 0) {
      distX = dem::Parameter::getSingleton().parameter["Xmax"] - dem::Parameter::getSingleton().parameter["Xmin"];
      distY = dem::Parameter::getSingleton().parameter["Ymax"] - dem::Parameter::getSingleton().parameter["Ymin"];
      distZ = dem::Parameter::getSingleton().parameter["Zmax"] - dem::Parameter::getSingleton().parameter["Zmin"];
      REAL x1 = dem::Parameter::getSingleton().parameter["Xmin"];  REAL x2 = dem::Parameter::getSingleton().parameter["Xmax"];
      REAL y1 = dem::Parameter::getSingleton().parameter["Ymin"];  REAL y2 = dem::Parameter::getSingleton().parameter["Ymax"];
      REAL z1 = dem::Parameter::getSingleton().parameter["Zmin"];  REAL z2 = dem::Parameter::getSingleton().parameter["Zmax"];
      setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
      setGrid(Rectangle(x1, y1, z1, x2, y2, z2));	// compute grid assumed to be the same as container, change in scatterParticle() if necessary.
      readParticle(dem::Parameter::getSingleton().datafile["particleFile"].c_str());
      readPeriDynamicsData(dem::Parameter::getSingleton().datafile["periFile"].c_str());
      openCompressProg(progressInf, "rigidInc_progress");
      openPeriProgress(periProgInf, "rigidInc_peri.dat");
      openPeriProgress(periProgInfHalf, "rigidInc_peri_half.dat");
      calcParticleVolume();		// volume and horizon size are related to the mesh, July 14, 2014
      calcHorizonSize();	 	// so these two should be calculated before peri-points that are within sand particles are deleted
      removeInsidePeriParticles();	// delete those peri-points that are inside sand particles
    }
    scatterDEMPeriParticle();
    constructPeriMatrix();	// construct the Matrix members in periParticleVec
    constructNeighbor();

//if(mpiRank==0){
//printPeriDomain("peridomain_rank_0.dat");
//printRecvPeriDomain("peridomain_recv_rank_0.dat");	// bondVec is not empty
//}
//if(mpiRank==1){
//printPeriDomain("peridomain_rank_1.dat");
//printRecvPeriDomain("peridomain_recv_rank_1.dat");	// bondVec is empty
//}

    std::size_t startStep = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startStep"]);
    std::size_t endStep   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endStep"]);
    std::size_t startSnap = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["startSnap"]);
    std::size_t endSnap   = static_cast<std::size_t> (dem::Parameter::getSingleton().parameter["endSnap"]);
    std::size_t netStep   = endStep - startStep + 1;
    std::size_t netSnap   = endSnap - startSnap + 1; 
    timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

    //REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
    iteration = startStep;
    std::size_t iterSnap = startSnap;
    char cstr0[50];
    if (mpiRank == 0) {
      plotGrid(strcat(combineString(cstr0, "rigidInc_gridplot_", iterSnap - 1, 3), ".dat"));
      printParticle(combineString(cstr0, "rigidInc_particle_", iterSnap - 1, 3));
      printPeriProgress(periProgInf, 0); 
      printPeriProgressHalf(periProgInfHalf, 0); 
    }
//    if (mpiRank == 0)
//      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" << std::setw(OWID) << "migraT"
//	       << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" << std::endl;

    commuParticle();		// the commuPeriParticle() have to be called after commuParticle(), since 
    commuPeriParticle();	// rankX1... are calculated in commuParticle()
    constructRecvPeriMatrix();	// construct the Matrix members in recvPeriParticle

    findRecvPeriBonds();	// find the peri-bonds between periParticleVec and recvPeriParticle
    calcParticleKinv();
    for(std::vector<periDynamics::PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	(*pt)->initial();
    }
    for(std::vector<periDynamics::PeriParticle*>::iterator pt=recvPeriParticleVec.begin(); pt!=recvPeriParticleVec.end(); pt++){
	(*pt)->initial();
    }
    calcParticleStress();
    calcParticleAcceleration();
    releaseRecvParticle();
//    releaseRecvPeriParticle();
//    releasePeriBondVec();	// free the bondVec in periParticleVec and bondVec.clear() before constructNeighbor() again
    findBoundaryPeriParticles();
    clearPeriDEMBonds();
    findPeriDEMBonds();
    while (iteration <= endStep) {
//      findBoundaryPeriParticles();	// boundary peri-points are determined by their initial positions, however they are still needed to be determined
//      findFixedPeriParticles();		// every step since some peri-points can pass from one cpu to another
      //commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
//      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

      // displacement control relies on constant time step, so do not call calcTimeStep().
      //calcTimeStep(); // use values from last step, must call before findConact
      findContact();
      if (isBdryProcess()) findBdryContact();

      clearContactForce();
      internalForce();
      if (isBdryProcess()) boundaryForce();
      runFirstHalfStep();

      applyPeriBoundaryCondition();
//      commuPeriParticle();
//      constructRecvPeriMatrix();
//      constructNeighbor();
//      checkBondParticleAlive();
      findRecvPeriBonds();

      calcParticleStress();
      calcParticleAcceleration();	// peri-point acceleration is set zero at the beginning

      findPeriDEMBonds();	// based on current states of DEM particles and peri-points
      applyCoupledForces();
      runSecondHalfStep();

      updateParticle();
      gatherBdryContact(); // must call before updateBoundary
//      updateBoundary(sigmaConf, "triaxial");
//      updateGrid();
//      updatePeriGrid();
   
      if (iteration % (netStep / netSnap) == 0) {
	//time1 = MPI_Wtime();
	gatherParticle();
	gatherPeriParticle();
	gatherEnergy(); 
        //time2 = MPI_Wtime(); gatherT = time2 - time1;

//	checkBondParticleAlive();
//	findPeriDEMBonds();	// update peri-dem bonds

	char cstr[50];
	if (mpiRank == 0) {
	  plotBoundary(strcat(combineString(cstr, "rigidInc_bdryplot_", iterSnap, 3), ".dat"));
	  plotGrid(strcat(combineString(cstr, "rigidInc_gridplot_", iterSnap, 3), ".dat"));
	  printParticle(combineString(cstr, "rigidInc_particle_", iterSnap, 3));
	  printBdryContact(combineString(cstr, "rigidInc_bdrycntc_", iterSnap, 3));
	  printBoundary(combineString(cstr, "rigidInc_boundary_", iterSnap, 3));
	  //printCompressProg(progressInf, distX, distY, distZ); // redundant
	  printPeriProgress(periProgInf, iterSnap);
	  printPeriProgressHalf(periProgInfHalf, iterSnap);
	}
	printContact(combineString(cstr, "rigidInc_contact_", iterSnap, 3));      
	++iterSnap;
      }
//      releaseRecvParticle(); // late release because printContact refers to received particles
//      releaseRecvPeriParticle();	// release recvPeriParticleVec and also remove peri-bonds between recvPeriParticle
//      releasePeriBondVec();
      //time1 = MPI_Wtime();
//      migrateParticle();
//      migratePeriParticle();
      //time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
//      if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore gather and print time at this step
//	debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT << std::setw(OWID) << migraT
//		 << std::setw(OWID) << totalT << std::setw(OWID) << (commuT + migraT)/totalT*100 << std::endl;

      if (mpiRank == 0 && iteration % 10 == 0)
	printCompressProg(progressInf, distX, distY, distZ);

      // no break condition, just through top/bottom displacement control
      ++iteration;
    } 

    if (mpiRank == 0) {
      printParticle("rigidInc_particle_end");
      printBdryContact("rigidInc_bdrycntc_end");
      printBoundary("rigidInc_boundary_end");
      printCompressProg(progressInf, distX, distY, distZ);
//      printPeriProgress(periProgInf, iterSnap);
    }
    if (mpiRank == 0){
      closeProg(progressInf);
      closeProg(periProgInf);
      closeProg(periProgInfHalf);
    }

  } // pullOutPeri


} // namespace dem ends

/*
// create a specimen from discreate particles through floating and then gravitation,
// boundaries are composed of fixed particles.
void Assembly::deposit_PtclBdry(gradation& grad,
int   particleLayers,
REAL rsize,
int   totalSteps,  
int   snapNum,
int   interval,
const char *iniptclfile,   
const char *ParticleFile, 
const char *contactfile,
const char *progressfile, 
const char *debugfile)
{
if (grad.rorc == 1) {
RORC = grad.rorc;
container.setCenter(Vec(0,0,0));
container.setDimx(grad.dimn);
container.setDimy(grad.dimn);
container.setDimz(grad.dimn);
	
generate_p(grad, iniptclfile, particleLayers, rsize, 4.0);
deposit_p(totalSteps,        // totalSteps
snapNum,          // number of snapNum
interval,           // print interval
grad.dimn,          // dimension of particle-composed-boundary
rsize,              // relative container size
iniptclfile,        // input file, initial particles
ParticleFile,       // output file, resulted particles, including snapNum 
contactfile,        // output file, resulted contacts, including snapNum 
progressfile,       // output file, statistical info
debugfile);         // output file, debug info
}
}
*/

 /*
 // particleLayers:
 // 0 - one free particle
 // 1 - a horizontal layer of free particles
 // 2 - multiple layers of free particles
 // ht- how many times of size would be the floating height
 void Assembly::generate_p(gradation&  grad,
 const char *ParticleFile,
 int particleLayers,
 REAL rsize,
 REAL ht)
 {
 REAL x,y,z;
 Particle* newptcl;
 particleNum = 0;
 REAL wall=2.2; // wall - wall height; ht - free particle height
 REAL est =1.02;
 int grid=static_cast<int> (nearbyint(rsize*10)-1);  

 // grid: dimension of free particle array.
 // 7 - small dimn container
 // 9 - medium dimn container 
 // 11- large dimn container 

 REAL dimn=grad.dimn;
 // particle boundary 1
 x=dimn/2*(grid+1)/10;
 for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
 for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }

 // particle boundary 2
 y=dimn/2*(grid+1)/10;
 for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
 for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }

 // particle boundary 3
 x=-dimn/2*(grid+1)/10;
 for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
 for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }

 // particle boundary 4
 y=-dimn/2*(grid+1)/10;
 for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
 for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }

 // particle boundary 6
 z=-dimn/2;
 for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
 for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }

 if (particleLayers == 0) {      // just one free particle
 newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }
 else if (particleLayers == 1) { // a horizontal layer of free particles
 z=dimn/2;
 for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
 for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }
 }
 else if (particleLayers == 2) { // multiple layers of free particles
 for (z=dimn/2; z<dimn/2 + dimn*ht; z+=dimn/2/5)
 for (x=-dimn/2*(grid-1)/10; x<dimn/2*(grid-1)/10*est; x+=dimn/2/5)
 for (y=-dimn/2*(grid-1)/10; y<dimn/2*(grid-1)/10*est; y+=dimn/2/5) {
 newptcl = new Particle(particleNum+1, 0, Vec(x,y,z), grad, young, poisson);
 particleVec.push_back(newptcl);
 particleNum++;
 }	
 }
    
 printParticle(ParticleFile);
    
 }
 */

 // rule out
 /*
   void Assembly::plotCavity(const char *str) const {
   std::ofstream ofs(str);
   if(!ofs) { debugInf << "stream error: plotCavity" << std::endl; exit(-1); }
   ofs.setf(std::ios::scientific, std::ios::floatfield);
   ofs.precision(OPREC);

   REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
   x1 = cavity.getMinCorner().getX();
   y1 = cavity.getMinCorner().getY();
   z1 = cavity.getMinCorner().getZ();
   x2 = cavity.getMaxCorner().getX();
   y2 = cavity.getMaxCorner().getY();
   z2 = cavity.getMaxCorner().getZ();

   ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
   ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
   ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
   ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1 << std::endl;
   ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1 << std::endl;
   ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
   ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
   ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2 << std::endl;
   ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2 << std::endl;
   ofs << "1 2 3 4 5 6 7 8" << std::endl;

   ofs.close();
   }

   void Assembly::plotSpring(const char *str) const {
   std::ofstream ofs(str);
   if(!ofs) { debugInf << "stream error: plotSpring" << std::endl; exit(-1); }
   ofs.setf(std::ios::scientific, std::ios::floatfield);
   ofs.precision(OPREC);

   std::size_t totalMemParticle = 0;
   for (std::size_t i = 0; i < memBoundary.size(); ++i) 
   for (std::size_t j = 0; j < memBoundary[i].size(); ++j) 
   for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k) 
   ++totalMemParticle;
   std::size_t totalSpring = springVec.size();
   ofs << "ZONE N=" << totalMemParticle << ", E=" << totalSpring << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << std::endl;
   Particle *pt = NULL;
   Vec vt;
   for (std::size_t i = 0; i < memBoundary.size(); ++i) 
   for (std::size_t j = 0; j < memBoundary[i].size(); ++j) 
   for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k) {
   pt = memBoundary[i][j][k]; 
   vt = pt->getCurrPos();
   ofs << std::setw(OWID) << vt.getX() << std::setw(OWID) << vt.getY() << std::setw(OWID) << vt.getZ() << std::endl;
   }
   for (std::size_t i = 0; i < springVec.size(); ++i) {
   ofs << std::setw(OWID) << springVec[i]->getParticleId1() - trimHistoryNum  << std::setw(OWID) << springVec[i]->getParticleId2() - trimHistoryNum << std::endl;
   }

   ofs.close();
   }  

   void Assembly::printMemParticle(const char *str) const  {
   std::ofstream ofs(str);
   if(!ofs) { debugInf << "stream error: printMemParticle" << std::endl; exit(-1); }
   ofs.setf(std::ios::scientific, std::ios::floatfield);
   ofs.precision(OPREC);
  
   std::size_t totalMemParticle = 0;
   for (std::size_t i = 0; i < memBoundary.size(); ++i) 
   for (std::size_t j = 0; j < memBoundary[i].size(); ++j) 
   for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k) 
   ++totalMemParticle;
  
   ofs << std::setw(OWID) << totalMemParticle << std::setw(OWID) << 1 << std::endl;
   ofs << std::setw(OWID) << container.getCenter().getX()
   << std::setw(OWID) << container.getCenter().getY()
   << std::setw(OWID) << container.getCenter().getZ()
   << std::setw(OWID) << container.getDimx()
   << std::setw(OWID) << container.getDimy()
   << std::setw(OWID) << container.getDimz() << std::endl;
  
   ofs << std::setw(OWID) << "ID"
   << std::setw(OWID) << "type"
   << std::setw(OWID) << "radius_a"
   << std::setw(OWID) << "radius_b"
   << std::setw(OWID) << "radius_c"
   << std::setw(OWID) << "position_x"
   << std::setw(OWID) << "position_y"
   << std::setw(OWID) << "position_z"
   << std::setw(OWID) << "axle_a_x"
   << std::setw(OWID) << "axle_a_y"
   << std::setw(OWID) << "axle_a_z"
   << std::setw(OWID) << "axle_b_x"
   << std::setw(OWID) << "axle_b_y"
   << std::setw(OWID) << "axle_b_z"
   << std::setw(OWID) << "axle_c_x"
   << std::setw(OWID) << "axle_c_y"
   << std::setw(OWID) << "axle_c_z"
   << std::setw(OWID) << "velocity_x"
   << std::setw(OWID) << "velocity_y"
   << std::setw(OWID) << "velocity_z"
   << std::setw(OWID) << "omga_x"
   << std::setw(OWID) << "omga_y"
   << std::setw(OWID) << "omga_z"
   << std::setw(OWID) << "force_x"
   << std::setw(OWID) << "force_y"
   << std::setw(OWID) << "force_z"
   << std::setw(OWID) << "moment_x"
   << std::setw(OWID) << "moment_y"
   << std::setw(OWID) << "moment_z"
   << std::endl;
  
   Particle *it = NULL;
   Vec vObj;
   for (std::size_t i = 0; i < memBoundary.size(); ++i) 
   for (std::size_t j = 0; j < memBoundary[i].size(); ++j) 
   for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k) {
   it = memBoundary[i][j][k];
   ofs << std::setw(OWID) << it->getId()
   << std::setw(OWID) << it->getType()
   << std::setw(OWID) << it->getA()
   << std::setw(OWID) << it->getB()
   << std::setw(OWID) << it->getC();
	
   vObj=it->getCurrPos();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getCurrDirecA();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getCurrDirecB();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getCurrDirecC();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getCurrVeloc();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getCurrOmga();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getForce();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ();
	
   vObj=it->getMoment();
   ofs << std::setw(OWID) << vObj.getX()
   << std::setw(OWID) << vObj.getY()
   << std::setw(OWID) << vObj.getZ() << std::endl;
   }
   ofs.close();  
   }
 
   // vector elements are in the order of:
   // x1: inner, outer
   // x2: inner, outer
   // y1: inner, outer
   // y2: inner, outer
   // z1: inner, outer
   // z2: inner, outer
   void Assembly::checkMembrane(vector<REAL> &vx ) const {
   ParticlePArray vec1d;  // 1-dimension
   std::vector< ParticlePArray  > vec2d; // 2-dimension
   REAL in, out, tmp;
   REAL x1_in, x1_out, x2_in, x2_out;
   REAL y1_in, y1_out, y2_in, y2_out;
   REAL z1_in, z1_out, z2_in, z2_out;

   // surface x1
   vec2d = memBoundary[0];
   in = vec2d[0][0]->getCurrPos().getX();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getX();
   if (tmp < out) out = tmp;
   if (tmp > in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   x1_in  = in;
   x1_out = out;

   // surface x2
   vec2d.clear();
   vec2d = memBoundary[1];
   in = vec2d[0][0]->getCurrPos().getX();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getX();
   if (tmp > out) out = tmp;
   if (tmp < in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   x2_in  = in;
   x2_out = out;

   // surface y1
   vec2d.clear();
   vec2d = memBoundary[2];
   in = vec2d[0][0]->getCurrPos().getY();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getY();
   if (tmp < out) out = tmp;
   if (tmp > in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   y1_in  = in;
   y1_out = out;

   // surface y2
   vec2d.clear();
   vec2d = memBoundary[3];
   in = vec2d[0][0]->getCurrPos().getY();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getY();
   if (tmp > out) out = tmp;
   if (tmp < in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   y2_in  = in;
   y2_out = out;
  
   // surface z1
   vec2d.clear();
   vec2d = memBoundary[4];
   in = vec2d[0][0]->getCurrPos().getZ();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getZ();
   if (tmp < out) out = tmp;
   if (tmp > in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   z1_in  = in;
   z1_out = out;

   // surface z2
   vec2d.clear();
   vec2d = memBoundary[5];
   in = vec2d[0][0]->getCurrPos().getZ();
   out= in;
   for (std::size_t i = 0; i < vec2d.size(); ++i)
   for (std::size_t j = 0; j < vec2d[i].size(); ++j) {
   tmp = vec2d[i][j]->getCurrPos().getZ();
   if (tmp > out) out = tmp;
   if (tmp < in ) in  = tmp;
   }
   vx.push_back(in);
   vx.push_back(out);
   z2_in  = in;
   z2_out = out;

   }
  
   //  1. it is important and helpful to mark a member function as const
   //     if it does NOT change member data.
   //  2. when a constant member function traverses member data, it can
   //     NOT change the data.
   //  3. then if it traverses a member data of a list, it should use a
   //     const_iterator, otherwise compiler will give errors.
   //  4. a const_iterator such as it also guarantees that (*it) will NOT
   //     change any data. if (*it) call a modification function, the 
   //     compiler will give errors.

   // OPENMP_IMPL: 
   // 0: implementation 0, ts partitions, based on linked list
   // 1: implementation 1, ts partitions, based on vector
   // 2: implementation 2, no partition, each thread leaps by ts until completed
   // 3: implementation 3, no partition, each thread leaps by ts until num/2 and handles two particles.
   // 4: implementation 4, no partition, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)

   //start of def OPENMP 
   #ifdef OPENMP	

   #if OPENMP_IMPL == 0
   // OpenMP implementation 0: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
   // implementation is based on linked list, also works for vector but not efficient.
   void Assembly::findContact() { 
   contactVec.clear();
   int possContact = 0;
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p1,NULL); 
   #endif

   int tid;   // thread id
   int ts;    // number of threads
   int num;   // number of particles
   int tnum;  // number of particles per thread
   int i, j;
   Vec u, v;
   num = particleVec.size();
   ot = particleVec.begin();
   ParticlePArray::iterator ot, it, pt;
  
   #pragma omp parallel num_threads(nThreads) private(tid, ts, tnum, it, pt, i, j, u, v) shared(num) reduction(+: possContact)
   {
   tid = omp_get_thread_num();
   ts  = omp_get_num_threads();
   tnum = num / ts;  // divide itso ts partitions
   rnum = num % ts;  // remainder of the division
   it = ot;          // start particle of each thread
    
   // determine starting point and extend of each partition
   // this algorithm applies to both list and vector
   if (rnum == 0) {
   for (i = 0; i < tid * tnum; ++i)
   ++it;         // starting point of each partition
   }
   else {
   if (tid < rnum) {
   tnum += 1;    // tnum changed
   for (i = 0; i < tid * tnum ; ++i)
   ++it;
   }
   else {
   for (i = 0; i < rnum * (tnum + 1) + (tid - rnum) * tnum; ++ i)
   ++it;
   }
   }
    
   // explore each partition
   for (j = 0 ; j < tnum; ++j, ++it) { 
   u=(*it)->getCurrPos();
   for (pt = it, ++pt; pt != particleVec.end(); ++pt) {
   v=(*pt)->getCurrPos();
   if (   ( vfabs(v-u) < (*it)->getA() + (*pt)->getA())
   && ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are fixed particles
   && ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are free boundary particles
   && ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are ghost particles
   contact<Particle> tmpContact(*it, *pt); // a local and temparory object
   ++possContact;
   if(tmpContact.isOverlapped())
   #pragma omp critical
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   }
   }
   }
   }
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
   #endif
   possContactNum   = possContact;
   actualContactNum = contactVec.size();
   } // end of OpenMP implementation 0

   #elif OPENMP_IMPL == 1
   // OpenMP implementation 1: ts partitions, each thread handles a partition, max diff = n*n*(1-1/ts)/ts
   // implementation is based on vector index.
   void Assembly::findContact() { 
   contactVec.clear();
   int possContact = 0;
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p1,NULL); 
   #endif
   int tid;   // thread id
   int ts;    // number of threads
   int num;   // number of particles
   int start; // start particle index of each thread
   int end;   // last particle index of each thread
   int i, j;
   Vec u, v;
   num = particleVec.size();
  
   #pragma omp parallel num_threads(nThreads) private(tid, ts, start, end, i, j, u, v) shared(num) reduction(+: possContact)
   {
   tid = omp_get_thread_num();
   ts  = omp_get_num_threads();
   start = tid * num / ts;
   end   = (tid + 1) * num / ts - 1;
    
   // explore each partition
   for (i = start; i <= end; ++i) { 
   u = particleVec[i]->getCurrPos();
   for (j = i + 1; j < num; ++j) {
   v = particleVec[j]->getCurrPos();
   if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
   && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
   && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
   && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
   contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
   ++possContact;
   if(tmpContact.isOverlapped())
   #pragma omp critical
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   }
   }
   }
   }
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf <<  std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
   #endif
   possContactNum   = possContact;  
   actualContactNum = contactVec.size();
   } // end of OpenMP implementation 1

   #elif OPENMP_IMPL == 2
   // OpenMP implementation 2: no partitions, each thread leaps by ts until completed, max diff = n*(ts-1)/ts
   void Assembly::findContact() { 
   contactVec.clear();
   int possContact = 0;
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p1,NULL); 
   #endif
   int tid;   // thread id
   int ts;    // number of threads
   int num;   // number of particles
   int i, j;
   Vec u, v;
   num = particleVec.size();
  
   #pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, u, v) shared(num) reduction(+: possContact)
   {
   tid = omp_get_thread_num();
   ts  = omp_get_num_threads();
    
   // explore each partition
   for (i = tid; i < num; i += ts) { 
   u = particleVec[i]->getCurrPos();
   for (j = i + 1; j < num; ++j) {
   v = particleVec[j]->getCurrPos();
   if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
   && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
   && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
   && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
   contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
   ++possContact;
   if(tmpContact.isOverlapped())
   #pragma omp critical
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   }
   }
   }
   }
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
   #endif
   possContactNum   = possContact;  
   actualContactNum = contactVec.size();
   } // end of OpenMP implementation 2

   #elif OPENMP_IMPL == 3
   // OpenMP implementation 3: no partitions, each thread leaps by ts until num/2 and handles two particles, max diff = 0
   void Assembly::findContact() { 
   contactVec.clear();
   int possContact = 0;
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p1,NULL); 
   #endif
   int tid;   // thread id
   int ts;    // number of threads
   int num;   // number of particles
   int i, j, k;
   Vec u, v;
   num = particleVec.size();
  
   #pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, k, u, v) shared(num) reduction(+: possContact)
   {
   tid = omp_get_thread_num();
   ts  = omp_get_num_threads();
    
   // explore each partition, works whether num is odd or even
   for (i = tid; i <= num / 2; i += ts) {
   int inc = num - 1 - 2*i;
   if (inc == 0) inc = 1; // avoid infinite loop when num is odd
   for (k = i; k <= num - 1 - i; k += inc ) {
   u = particleVec[k]->getCurrPos();
   for (j = k + 1; j < num; ++j) {
   v = particleVec[j]->getCurrPos();
   if (   ( vfabs(v-u) < particleVec[k]->getA() + particleVec[j]->getA() )
   && ( particleVec[k]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
   && ( particleVec[k]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
   && ( particleVec[k]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
   contact<Particle> tmpContact(particleVec[k], particleVec[j]); // a local and temparory object
   ++possContact;
   if(tmpContact.isOverlapped())
   #pragma omp critical
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   }
   }
   }
   }
   }
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
   #endif
   possContactNum   = possContact;  
   actualContactNum = contactVec.size();
   } // end of OpenMP implementation 3

   #elif OPENMP_IMPL == 4
   // OpenMP implementation 4: no partitions, parallel for, various loop scheduling: (static), (static,1), (dynamic), (dynamic,1)
   void Assembly::findContact() { 
   contactVec.clear();
   int possContact = 0;
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p1,NULL); 
   #endif

   int num;   // number of particles
   int i, j;
   Vec u, v;
   num = particleVec.size();
  
   #pragma omp parallel for num_threads(nThreads) private(i, j, u, v) shared(num) reduction(+: possContact) schedule(dynamic)
   for (i = 0; i < num - 1; ++i) { 
   u = particleVec[i]->getCurrPos();
   for (j = i + 1; j < num; ++j) {
   v = particleVec[j]->getCurrPos();
   if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
   && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )      // not both are fixed particles
   && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )      // not both are free boundary particles
   && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) { // not both are ghost particles
   contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and temparory object
   ++possContact;
   if(tmpContact.isOverlapped())
   #pragma omp critical
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   }
   }
   }
  
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2); 
   #endif
   possContactNum   = possContact;  
   actualContactNum = contactVec.size();
   } // end of OpenMP implementation 4

   #endif

   #else //else of def OPENMP, i.e., serial versions start here:

   //start of ndef BINNING
   #ifndef BINNING
   void Assembly::findContact() { // serial version, O(n x n), n is the number of particles.
   contactVec.clear();
   possContactNum = 0;

   #ifdef TIME_PROFILE
   REAL time_r = 0; // time consumed in contact resolution, i.e., tmpContact.isOverlapped()
   gettimeofday(&time_p1,NULL); 
   #endif
    
   int num1 = particleVec.size();  // particles inside container
   int num2 = mergeParticleVec.size(); // particles inside container (at front) + particles from neighboring blocks (at end)
   for (int i = 0; i < num1 - 1; ++i) {
   Vec u = particleVec[i]->getCurrPos();
   for (int j = i + 1; j < num2; ++j) {
   Vec v = mergeParticleVec[j]->getCurrPos();
   if (   ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA())
   && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )      // not both are fixed particles
   && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )      // not both are free boundary particles
   && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )  ) { // not both are ghost particles
   Contact tmpContact(particleVec[i], mergeParticleVec[j]); // a local and temparory object
   ++possContactNum;
   #ifdef TIME_PROFILE
   gettimeofday(&time_r1,NULL); 
   #endif
   if(tmpContact.isOverlapped())
   contactVec.push_back(tmpContact);    // containers use value semantics, so a "copy" is pushed back.
   #ifdef TIME_PROFILE
   gettimeofday(&time_r2,NULL); 
   time_r += timediffsec(time_r1, time_r2);
   #endif
   }
   }
   }	
    
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
   #endif
    
   actualContactNum = contactVec.size();
   }

   //else of ndef BINNING
   #else
   void Assembly::findContact() { // serial version, binning methods, cell slightly larger than maximum particle
   contactVec.clear();
   possContactNum = 0;
  
   #ifdef TIME_PROFILE
   REAL time_r = 0;
   gettimeofday(&time_p1,NULL); 
   #endif
   REAL maxDiameter = gradation.getPtclMaxRadius() * 2;
   int  nx = floor (container.getDimx() / maxDiameter);
   int  ny = floor (container.getDimy() / maxDiameter);
   int  nz = floor (container.getDimz() *1.5 / maxDiameter);
   REAL dx = container.getDimx() / nx;
   REAL dy = container.getDimx() / ny;
   REAL dz = container.getDimx() *1.5 / nz;
   Vec  minCorner= container.getMinCorner();
   REAL x0 = minCorner.getX();
   REAL y0 = minCorner.getY();
   REAL z0 = minCorner.getZ();
  
   // 26 neighbors of each cell
   int neighbor[26][3];
   int count = 0;
   for (int i = -1; i < 2; ++i)
   for (int j = -1; j < 2; ++j)
   for (int k = -1; k < 2; ++k) {
   if (! (i == 0 && j == 0 && k==0 ) ) {
   neighbor[count][0] = i;
   neighbor[count][1] = j;
   neighbor[count][2] = k;
   ++count;
   }
   }
 
   // 4-dimensional array of cellVec
   typedef std::pair<bool, ParticlePArray > cellT;
   std::vector< std::vector< std::vector < cellT > > > cellVec;
   cellVec.resize(nx);
   for (int i = 0; i < cellVec.size(); ++i) {
   cellVec[i].resize(ny);
   for (int j = 0; j < cellVec[i].size(); ++j)
   cellVec[i][j].resize(nz);
   }
   // mark each cell as not searched
   for (int i = 0; i < nx; ++i)
   for (int j = 0; j < ny; ++j)
   for (int k = 0; k < nz; ++k)
   cellVec[i][j][k].first = false; // has not ever been searched

   // find particles in each cell
   Vec center;
   REAL x1, x2, y1, y2, z1, z2;
   for (int i = 0; i < nx; ++i)
   for (int j = 0; j < ny; ++j)
   for (int k = 0; k < nz; ++k) {
   x1 = x0 + dx * i;
   x2 = x0 + dx * (i + 1);
   y1 = y0 + dy * j;
   y2 = y0 + dy * (j + 1);
   z1 = z0 + dz * k;
   z2 = z0 + dz * (k + 1);
   for (int pt = 0; pt < particleVec.size(); ++pt) {
   center = particleVec[pt]->getCurrPos();
   if (center.getX() >= x1 && center.getX() < x2 &&
   center.getY() >= y1 && center.getY() < y2 &&
   center.getZ() >= z1 && center.getZ() < z2)
   cellVec[i][j][k].second.push_back( particleVec[pt] );
   }
   }
  
   // for each cell:
   Particle *it, *pt;
   Vec u, v;
   for (int i = 0; i < nx; ++i)
   for (int j = 0; j < ny; ++j)
   for (int k = 0; k < nz; ++k) {
   // for particles inside the cell	  
   for (int m = 0; m < cellVec[i][j][k].second.size(); ++m) {
   it = cellVec[i][j][k].second[m];
   u  = it->getCurrPos();
	  
   // for particles inside the cell itself   
   for (int n = m + 1; n < cellVec[i][j][k].second.size(); ++n) {
   //debugInf <<  i << " " << j << " " << k << " " << "m n size=" << m << " " << n << " " <<  cellVec[i][j][k].size() << std::endl;
   pt = cellVec[i][j][k].second[n];
   v  = pt->getCurrPos();
   if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
   ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
   ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
   ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
   contact<Particle> tmpContact(it, pt); // a local and temparory object
   ++possContactNum;
   #ifdef TIME_PROFILE
   gettimeofday(&time_r1,NULL); 
   #endif
   if(tmpContact.isOverlapped())
   contactVec.push_back(tmpContact);   // containers use value semantics, so a "copy" is pushed back.
   #ifdef TIME_PROFILE
   gettimeofday(&time_r2,NULL); 
   time_r += timediffsec(time_r1, time_r2);
   #endif
   }
   }
	  
   // for 26 neighboring cells
   for (int ncell = 0; ncell < 26; ++ncell ) {
   int ci = i + neighbor[ncell][0];
   int cj = j + neighbor[ncell][1];
   int ck = k + neighbor[ncell][2];
   if (ci > -1 && ci < nx && cj > -1 && cj < ny && ck > -1 && ck < nz && cellVec[ci][cj][ck].first == false ) {
   //debugInf << "i j k m ncell ci cj ck size contacts= " << i << " " << j << " " << k << " " << m  << " " << ncell << " " << ci << " " << cj << " " << ck << " " << cellVec[ci][cj][ck].second.size() << " "  << contactVec.size() << std::endl;
   ParticlePArray vt = cellVec[ci][cj][ck].second;
   for (int n = 0; n < vt.size(); ++n) {
   pt = vt[n];
   v  = pt->getCurrPos();
   if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
   ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed particles
   ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free boundary particles
   ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost particles
   contact<Particle> tmpContact(it, pt); // a local and temparory object
   ++possContactNum;
   #ifdef TIME_PROFILE
   gettimeofday(&time_r1,NULL); 
   #endif
   if(tmpContact.isOverlapped())
   contactVec.push_back(tmpContact);   // containers use value semantics, so a "copy" is pushed back.
   #ifdef TIME_PROFILE
   gettimeofday(&time_r2,NULL); 
   time_r += timediffsec(time_r1, time_r2);
   #endif
		  
   }
   }
   }
   }
   }
   cellVec[i][j][k].first = true; // searched, will not be searched again
	
   }
  
   #ifdef TIME_PROFILE
   gettimeofday(&time_p2,NULL);
   debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) << timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" << std::setw(OWID) << time_r; 
   #endif
  
   actualContactNum = contactVec.size();
   }

   //end of ndef BINNING
   #endif

   //end of def OPENMP 
   #endif


   Vec Assembly::getTopFreeParticlePosition() const {
   ParticlePArray::const_iterator it,jt,kt;
   it=particleVec.begin();
   while (it!=particleVec.end() && (*it)->getType()!=0)   // find the 1st free particle
   ++it;

   if (it==particleVec.end())    // no free particles
   return 0;

   jt=it; 
   kt=it;
    
   // two cases:
   // 1: 1st particle is not free
   // 2: 1st particle is free
   if (++kt!=particleVec.end()) { // case1: more than 2 particles; case 2: more than 1 particle
   for(++it;it!=particleVec.end();++it) {
   if ((*it)->getType()==0)
   if ((*it)->getCurrPos().getZ() > (*jt)->getCurrPos().getZ())
   jt=it;
   }
   return (*jt)->getCurrPos();
   }
   else {
   if ((*it)->getType()==0)  // case1: only 2 particles, the 2nd one is free; case2: only 1 particle
   return (*it)->getCurrPos();
   else
   return 0;
   }

   }



   REAL Assembly::ellipPileForce() {
   REAL val=0;
   for(ParticlePArray::iterator it=particleVec.begin();it!=particleVec.end();++it)
   if ((*it)->getType()==3) {
   val = (*it)->getForce().getZ();
   break;
   }
   return val;
   }

   Vec Assembly::ellipPileDimn() {
   Vec val;
   for(ParticlePArray::iterator it=particleVec.begin();it!=particleVec.end();++it)
   if ((*it)->getType()==3) {
   val = Vec((*it)->getA(), (*it)->getB(), (*it)->getC());
   break;
   }
   return val;
   }

   REAL Assembly::ellipPileTipZ() {
   REAL val=0;
   for(ParticlePArray::iterator it=particleVec.begin();it!=particleVec.end();++it)
   if ((*it)->getType()==3) {
   val = (*it)->getCurrPos().getZ()-(*it)->getA();
   break;
   }
   return val;
   }

   REAL Assembly::ellipPilePeneVol() {
   REAL val=0;
   if (getTopFreeParticlePosition().getZ()-ellipPileTipZ() <= 0)
   val=0;
   else{
   // low: a signed number as lower limit for volumetric integration
   REAL low=ellipPileTipZ() + ellipPileDimn().getX() - getTopFreeParticlePosition().getZ(); 
   REAL lowint=low-pow(low,3)/3.0/pow(ellipPileDimn().getX(),2);
   val = Pi * ellipPileDimn().getY() * ellipPileDimn().getZ()
   *(2.0/3*ellipPileDimn().getX()-lowint);
   }
   return val;
   }

   void Assembly::ellipPileUpdate() {
   for(ParticlePArray::iterator it=particleVec.begin();it!=particleVec.end();++it) {
   if ((*it)->getType()==3) {
   (*it)->setCurrVeloc(Vec(0, 0, -pileRate));
   (*it)->setCurrPos( (*it)->getPrevPos() + (*it)->getCurrVeloc() * timeStep);
   }
   }
   }





   void Assembly::springForce() {
   for (vector<Spring*>::iterator it = springVec.begin(); it != springVec.end(); ++it)
   (*it)->applyForce();
   }

   void Assembly::readCavityBoundary(const char *str) {
   std::ifstream ifs(str);
   if(!ifs) { debugInf << "stream error: readCavityBoundary" << std::endl; exit(-1); }  

   Boundary<Particle>* rbptr;
   int type;
   cavityBoundaryVec.clear();
   int boundaryNum;
   ifs >> boundaryNum;
   for(int i = 0; i < boundaryNum; i++) {
   ifs >> type;
   if(type == 1) // plane boundary
   rbptr = new planeBoundary<Particle>(ifs);
   cavityBoundaryVec.push_back(rbptr);
   }

   ifs.close();
   }


   void Assembly::printCavityBoundary(const char *str) const {
   std::ofstream ofs(str);
   if(!ofs) { debugInf << "stream error: printCavityBoundary" << std::endl; exit(-1); }
   ofs.setf(std::ios::scientific, std::ios::floatfield);
  
   ofs << std::setw(OWID) << cavityBoundaryVec.size() << std::endl;
   BoundaryPArray::const_iterator rt;
   for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
   (*rt)->display(ofs);
   ofs << std::endl;
  
   ofs.close();
   }



   void Assembly::findCavityContact() {
   BoundaryPArray::iterator rt;
   for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
   (*rt)->findBdryContact(allParticleVec);
   }


   void Assembly::cavityBoundaryForce() {
   BoundaryPArray::iterator rt;
   for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
   (*rt)->boundaryForce(boundaryTgtMap);
   }

   Vec Assembly::getNormalForce(int bdry) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==bdry)
   return (*it)->getNormalForce();
   }
   return 0;
   }

   Vec Assembly::getShearForce(int bdry) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==bdry)
   return (*it)->getShearForce();
   }
   return 0;
   }

   REAL Assembly::getAvgNormal(int bdry) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==bdry)
   return (*it)->getAvgNormal();
   }
   return 0;
   }

   Vec Assembly::getApt(int bdry) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==bdry)
   return (*it)->getApt();
   }
   return 0;
   }


   Vec Assembly::getDirc(int bdry) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==bdry)
   return (*it)->getDirc();
   }
   return 0;
   }

   REAL Assembly::getArea(int n) const {
   BoundaryPArray::const_iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==n)
   return (*it)->area;
   }
   return 0;
   }

   void Assembly::setArea(int n, REAL a) {
   BoundaryPArray::iterator it;
   for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
   if((*it)->getBdryID()==n)
   (*it)->area=a;
   }
   }

   REAL Assembly::getAvgPressure() const {
   BoundaryPArray::const_iterator rt;
   REAL avgpres=0;
   for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
   avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
   return avgpres/=boundaryVec.size();
   }

   // only update CoefOfLimits[0] for specified boundaries
   void Assembly::updateBoundary(int bn[], UPDATECTL rbctl[], int num) {
   for(int i=0;i<num;i++) {
   for(BoundaryPArray::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt) {
   if((*rt)->getBdryID()==bn[i]) {
   (*rt)->update(rbctl[i]);
   break;
   }
   }
   }
   }

   // update CoefOfLimits[1,2,3,4] for all 6 boundaries
   void Assembly::updateBoundary6() {
   for(BoundaryPArray::iterator rt=boundaryVec.begin();rt!=boundaryVec.end();++rt) {
   if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3) {
   for(BoundaryPArray::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
   if((*lt)->getBdryID()==4)
   (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==2)
   (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==5)
   (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==6)
   (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
   }
   }
   else if((*rt)->getBdryID()==2 || (*rt)->getBdryID()==4) {
   for(BoundaryPArray::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
   if((*lt)->getBdryID()==1)
   (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==3)
   (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==5)
   (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==6)
   (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
   }

   }
   else if((*rt)->getBdryID()==5 || (*rt)->getBdryID()==6) {
   for(BoundaryPArray::iterator lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
   if((*lt)->getBdryID()==1)
   (*rt)->CoefOfLimits[1].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==3)
   (*rt)->CoefOfLimits[2].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==2)
   (*rt)->CoefOfLimits[3].apt=(*lt)->CoefOfLimits[0].apt;
   else if((*lt)->getBdryID()==4)
   (*rt)->CoefOfLimits[4].apt=(*lt)->CoefOfLimits[0].apt;
   }

   }
	
   }
   }


   void Assembly::deposit_repose(int   interval,
   const char *inibdryfile,
   const char *ParticleFile, 
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   this->container = container;
  
   buildBoundary(5, inibdryfile); // container unchanged
  
   angleOfRepose(interval,           // print interval
   inibdryfile,        // input file, initial boundaries
   ParticleFile,       // output file, resulted particles, including snapNum 
   contactfile,        // output file, resulted contacts, including snapNum 
   progressfile,       // output file, statistical info
   debugfile);         // output file, debug info
   }

   void Assembly::angleOfRepose(int   interval,
   const char *inibdryfile,
   const char *ParticleFile, 
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   // pre_1: open streams for output.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: angleOfRepose" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf.precision(OPREC);
   progressInf << std::setw(OWID) << "iteration"
   << std::setw(OWID) << "poss_contact"
   << std::setw(OWID) << "actual_contact"
   << std::setw(OWID) << "penetration"
   << std::setw(OWID) << "avg_normal"
   << std::setw(OWID) << "avg_tangt"
   << std::setw(OWID) << "avg_velocity"
   << std::setw(OWID) << "avg_omga"
   << std::setw(OWID) << "avg_force"
   << std::setw(OWID) << "avg_moment"
   << std::setw(OWID) << "trans_energy"
   << std::setw(OWID) << "rotat_energy"
   << std::setw(OWID) << "kinet_energy"
   << std::setw(OWID) << "poten_energy"
   << std::setw(OWID) << "total_energy"
   << std::setw(OWID) << "void_ratio"
   << std::setw(OWID) << "porosity"
   << std::setw(OWID) << "coord_number"
   << std::setw(OWID) << "density"
   << std::setw(OWID) << "sigma_y1"
   << std::setw(OWID) << "sigma_y2"
   << std::setw(OWID) << "sigma_x1"
   << std::setw(OWID) << "sigma_x2"
   << std::setw(OWID) << "sigma_z1"
   << std::setw(OWID) << "sigma_z2"
   << std::setw(OWID) << "mean_stress"
   << std::setw(OWID) << "dimx"
   << std::setw(OWID) << "dimy"
   << std::setw(OWID) << "dimz"
   << std::setw(OWID) << "volume"
   << std::setw(OWID) << "epsilon_x"
   << std::setw(OWID) << "epsilon_y"
   << std::setw(OWID) << "epsilon_z"
   << std::setw(OWID) << "epsilon_v"
   << std::setw(OWID) << "vibra_t_step"
   << std::setw(OWID) << "impact_t_step"
   << std::setw(OWID) << "wall_time" << std::endl;
  
   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: angleOfRepose" << std::endl; exit(-1); }
   debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
   // pre_2. create boundaries from existing files.
   readBoundary(inibdryfile);

   // pre_3: define variables used in iterations.
   REAL avgNormal=0;
   REAL avgTangt=0;
   int  stepsnum=0;
   char stepsstr[4];
   bool toSnapshot=false;
   char stepsfp[50];
   REAL void_ratio=0;
   REAL bdry_penetr[7] = {0,0,0,0,0,0,0};
   int  bdry_cntnum[7] = {0,0,0,0,0,0,0};

   REAL maxRadius = gradation.getPtclMaxRadius();
   REAL maxDiameter = maxRadius * 2.0;
   REAL z0 = container.getMinCorner().getZ();
   ParticlePArray lastPtcls;
   Particle *newPtcl = NULL;
   int layers = 1; // how many layers of new particles to generate each time

   iteration = 0; 
   int particleNum = 0;
   REAL zCurr;
   gettimeofday(&time_w1,NULL);
   // iterations starting ...
   do
   {
   // 1. add particle
   if ( particleNum == 0 ) {
   zCurr = z0 + maxRadius;

   for ( int i = 0; i != layers; ++i) {
   newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   lastPtcls.push_back(newPtcl);
	  
   newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	  
   newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	  
   newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	  
   newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
   }
   toSnapshot = true;

   }
   else {
   vector<Particle*>::iterator it;
   bool allInContact = false;
   for ( it = lastPtcls.begin(); it != lastPtcls.end(); ++it) {
   if ( (*it)->isInContact() ) 
   allInContact = true;
   else {
   allInContact = false;
   break;
   }
   }

   if ( allInContact ) {

   lastPtcls.clear(); // do not delete those pointers to release memory; particleVec will do it.
   zCurr = getPtclMaxZ(allParticleVec) + maxDiameter;

   for ( int i = 0; i != layers; ++i) {
   newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   lastPtcls.push_back(newPtcl);
	    
   newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	    
   newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	    
   newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);
	    
   newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr + maxDiameter * i), gradation, young, poisson);
   particleVec.push_back(newPtcl);
   ++particleNum;
   //lastPtcls.push_back(newPtcl);	
   }
   toSnapshot = true;
   }
   }
   // 2. create possible boundary particles and contacts between particles.
   findContact();
   findBdryContact();
      
   // 3. set particle forces/moments as zero before each re-calculation,
   clearContactForce();	
      
   // 4. calculate contact forces/moments and apply them to particles.
   internalForce(avgNormal, avgTangt);
      
   // 5. calculate boundary forces/moments and apply them to particles.
   boundaryForce(bdry_penetr, bdry_cntnum);
      
   // 6. update particles' velocity/omga/position/orientation based on force/moment.
   updateParticle();
      
   // 7. (1) output particles and contacts information as snapNum.
   if (toSnapshot) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   time(&timeStamp);
   g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
   ++stepsnum;
   toSnapshot = false;
   }
      
   // 8. (2) output stress and strain info.
   if (iteration % interval == 0) {
   gettimeofday(&time_w2,NULL);
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*(getActualContactNum()
   +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
   +bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << "0"
   << std::setw(OWID) << getVibraTimeStep()
   << std::setw(OWID) << getImpactTimeStep()
   << std::setw(OWID) << timediffsec(time_w1,time_w2)
   << std::endl;
   }
      
   // 7. loop break conditions.
   ++iteration;
      
   } while (particleNum < 2000); //( zCurr < container.getMaxCorner().getZ() );  //(++iteration < totalSteps);
    
   // post_1. store the final snapshot of particles & contacts.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
   g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;

   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   void Assembly::scale_PtclBdry(int   totalSteps,  
   int   snapNum,
   int   interval,
   REAL dimn,
   REAL rsize,
   const char *iniptclfile,   
   const char *ParticleFile, 
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   deposit_p(totalSteps,        // totalSteps
   snapNum,          // number of snapNum
   interval,           // print interval
   dimn,               // dimension of particle-composed-boundary
   rsize,              // relative container size
   iniptclfile,        // input file, initial particles
   ParticleFile,       // output file, resulted particles, including snapNum 
   contactfile,        // output file, resulted contacts, including snapNum 
   progressfile,       // output file, statistical info
   debugfile);         // output file, debug info
   }


   // collapse a deposited specimen through gravitation
   void Assembly::collapse(int   totalSteps,  
   int   snapNum,
   int   interval,
   const char *iniptclfile,
   const char *initboundary,
   const char *ParticleFile,
   const char *contactfile,
   const char *progressfile,
   const char *debugfile)
   {
   buildBoundary(1,              // 1-only bottom boundary; 5-no top boundary;6-boxed 6 boundaries
   initboundary);  // output file, containing boundaries info
  
   deposit(totalSteps,        // number of iterations
   snapNum,          // number of snapNum
   interval,           // print interval
   iniptclfile,        // input file, initial particles
   initboundary,       // input file, boundaries
   ParticleFile,       // output file, resulted particles, including snapNum 
   contactfile,        // output file, resulted contacts, including snapNum 
   progressfile,       // output file, statistical info
   debugfile);         // output file, debug info
   }

  


   // make a cavity inside the sample and remove particles in the cavity
   void Assembly::trimCavity(bool toRebuild,
   const char *ParticleFile,
   const char *cavParticleFile)
   {
   if (toRebuild) readParticle(ParticleFile);
   trimHistoryNum = allParticleVec.size();

   REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
   x1 = cavity.getMinCorner().getX();
   y1 = cavity.getMinCorner().getY();
   z1 = cavity.getMinCorner().getZ();
   x2 = cavity.getMaxCorner().getX();
   y2 = cavity.getMaxCorner().getY();
   z2 = cavity.getMaxCorner().getZ();
   x0 = cavity.getCenter().getX();
   y0 = cavity.getCenter().getY();
   z0 = cavity.getCenter().getZ();
 
   ParticlePArray::iterator itr;
   Vec center;
   REAL delta = gradation.getPtclMaxRadius();

   for (itr = particleVec.begin(); itr != particleVec.end(); ) {
   center=(*itr)->getCurrPos();
   if(center.getX() + delta  >= x1 && center.getX() - delta <= x2 &&
   center.getY() + delta  >= y1 && center.getY() - delta <= y2 &&
   center.getZ() + delta  >= z1 && center.getZ() - delta <= z2 )
   {
   delete (*itr); // release memory
   itr = particleVec.erase(itr); 
   }
   else
   ++itr;
   }
  
   printParticle(cavParticleFile);
   }


   // expand partcile size by some percentage for particles inside cavity
   void Assembly::expandCavityParticle(bool toRebuild,
   REAL percent,
   const char *cavityptclfile,
   const char *ParticleFile,
   const char *newptclfile)
   {
   if (toRebuild) readParticle(ParticleFile);
   trimHistoryNum = allParticleVec.size();

   REAL x1,x2,y1,y2,z1,z2;
   x1 = cavity.getMinCorner().getX();
   y1 = cavity.getMinCorner().getY();
   z1 = cavity.getMinCorner().getZ();
   x2 = cavity.getMaxCorner().getX();
   y2 = cavity.getMaxCorner().getY();
   z2 = cavity.getMaxCorner().getZ();
 
   ParticlePArray::iterator itr;
   Vec center;

   int cavityPtclNum = 0;
   for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ) {
   center=(*itr)->getCurrPos();
   if(center.getX() > x1 && center.getX() < x2 &&
   center.getY() > y1 && center.getY() < y2 &&
   center.getZ() > z1 && center.getZ() < z2 )
   ++cavityPtclNum;
   }

   printCavityParticle(cavityPtclNum, cavityptclfile);

   for (itr = particleVec.begin(); itr != particleVec.end(); ++itr ) {
   center=(*itr)->getCurrPos();
   if(center.getX() > x1 && center.getX() < x2 &&
   center.getY() > y1 && center.getY() < y2 &&
   center.getZ() > z1 && center.getZ() < z2 )
   (*itr)->expand(percent);
   }

   printParticle(newptclfile);
   }


   void Assembly::printCavityParticle(int total, const char *str) const {
   std::ofstream ofs(str);
   if(!ofs) { debugInf << "stream error: printCavityParticle" << std::endl; exit(-1); }
   ofs.setf(std::ios::scientific, std::ios::floatfield);
   ofs.precision(OPREC);
   ofs << std::setw(OWID) << total << std::setw(OWID) << 1 << std::endl;
   ofs << std::setw(OWID) << cavity.getCenter().getX()
   << std::setw(OWID) << cavity.getCenter().getY()
   << std::setw(OWID) << cavity.getCenter().getZ()
   << std::setw(OWID) << cavity.getDimx()
   << std::setw(OWID) << cavity.getDimy()
   << std::setw(OWID) << cavity.getDimz() << std::endl;
  
   ofs << std::setw(OWID) << "ID"
   << std::setw(OWID) << "type"
   << std::setw(OWID) << "radius_a"
   << std::setw(OWID) << "radius_b"
   << std::setw(OWID) << "radius_c"
   << std::setw(OWID) << "position_x"
   << std::setw(OWID) << "position_y"
   << std::setw(OWID) << "position_z"
   << std::setw(OWID) << "axle_a_x"
   << std::setw(OWID) << "axle_a_y"
   << std::setw(OWID) << "axle_a_z"
   << std::setw(OWID) << "axle_b_x"
   << std::setw(OWID) << "axle_b_y"
   << std::setw(OWID) << "axle_b_z"
   << std::setw(OWID) << "axle_c_x"
   << std::setw(OWID) << "axle_c_y"
   << std::setw(OWID) << "axle_c_z"
   << std::setw(OWID) << "velocity_x"
   << std::setw(OWID) << "velocity_y"
   << std::setw(OWID) << "velocity_z"
   << std::setw(OWID) << "omga_x"
   << std::setw(OWID) << "omga_y"
   << std::setw(OWID) << "omga_z"
   << std::setw(OWID) << "force_x"
   << std::setw(OWID) << "force_y"
   << std::setw(OWID) << "force_z"
   << std::setw(OWID) << "moment_x"
   << std::setw(OWID) << "moment_y"
   << std::setw(OWID) << "moment_z"
   << std::endl;

   REAL x1,x2,y1,y2,z1,z2;
   x1 = cavity.getMinCorner().getX();
   y1 = cavity.getMinCorner().getY();
   z1 = cavity.getMinCorner().getZ();
   x2 = cavity.getMaxCorner().getX();
   y2 = cavity.getMaxCorner().getY();
   z2 = cavity.getMaxCorner().getZ();
  
   Vec tmp;
   ParticlePArray::const_iterator  it;
   for (it=particleVec.begin();it!=particleVec.end();++it)  {
   Vec center=(*it)->getCurrPos();
   if(center.getX() > x1 && center.getX() < x2 &&
   center.getY() > y1 && center.getY() < y2 &&
   center.getZ() > z1 && center.getZ() < z2 ) {

   ofs << std::setw(OWID) << (*it)->getId()
   << std::setw(OWID) << (*it)->getType()
   << std::setw(OWID) << (*it)->getA()
   << std::setw(OWID) << (*it)->getB()
   << std::setw(OWID) << (*it)->getC();
    
   tmp=(*it)->getCurrPos();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getCurrDirecA();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getCurrDirecB();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getCurrDirecC();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getCurrVeloc();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getCurrOmga();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getForce();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ();
    
   tmp=(*it)->getMoment();
   ofs << std::setw(OWID) << tmp.getX()
   << std::setw(OWID) << tmp.getY()
   << std::setw(OWID) << tmp.getZ() << std::endl;
   }
   }
  
   ofs.close();
   }


   // bdrymum = 6 by default
   // the variable existMaxID is important because cavity and container
   // use the same boundaryTgtMap.
   void Assembly::buildCavityBoundary(int existMaxId, const char *boundaryFile)
   {
   std::ofstream ofs(boundaryFile);
   if(!ofs) { debugInf << "stream error: buildCavityBoundary" << std::endl; exit(-1); }

   REAL x1,x2,y1,y2,z1,z2,x0,y0,z0;
   x1 = cavity.getMinCorner().getX();
   y1 = cavity.getMinCorner().getY();
   z1 = cavity.getMinCorner().getZ();
   x2 = cavity.getMaxCorner().getX();
   y2 = cavity.getMaxCorner().getY();
   z2 = cavity.getMaxCorner().getZ();
   x0 = cavity.getCenter().getX();
   y0 = cavity.getCenter().getY();
   z0 = cavity.getCenter().getZ();

   int boundaryNum = 6;

   ofs.setf(std::ios::scientific, std::ios::floatfield);
   ofs << std::setw(OWID) << 0
   << std::setw(OWID) << boundaryNum << std::endl << std::endl;

   // boundary 1
   ofs << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 1
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0     
   << std::setw(OWID) << x2
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0    
   << std::setw(OWID) << y1
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y2
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << x0     
   << std::setw(OWID) << y0    
   << std::setw(OWID) << z2 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 
   << std::setw(OWID) << -1
   << std::setw(OWID) << x0     
   << std::setw(OWID) << y0     
   << std::setw(OWID) << z1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl
    
   // boundary 2
   << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 2
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0     
   << std::setw(OWID) << x0    
   << std::setw(OWID) << y2
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << x2
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << x1 
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0     
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << x0     
   << std::setw(OWID) << y0    
   << std::setw(OWID) << z2 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 
   << std::setw(OWID) << -1
   << std::setw(OWID) << x0     
   << std::setw(OWID) << y0      
   << std::setw(OWID) << z1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl
    
   // boundary 3
   << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 3
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0     
   << std::setw(OWID) << x1
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0     
   << std::setw(OWID) << y1
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0  
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0       
   << std::setw(OWID) << y2
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0     
   << std::setw(OWID) << z2 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 
   << std::setw(OWID) << -1
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0      
   << std::setw(OWID) << z1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl
    
   // boundary 4
   << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 4
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0     
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y1
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << x2
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << -1 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << x1 
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0     
   << std::setw(OWID) << z2 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 
   << std::setw(OWID) << -1
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0      
   << std::setw(OWID) << z1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl
    
   // boundary 5
   << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 5
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1     
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0
   << std::setw(OWID) << z2 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << x2
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << -1 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << x1 
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y2
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y1
   << std::setw(OWID) << z0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl
    
   // boundary 6
   << std::setw(OWID) << 1 << std::endl
   << std::setw(OWID) << existMaxId + 6
   << std::setw(OWID) << 5
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1    
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y0
   << std::setw(OWID) << z1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << 0
   << std::setw(OWID) << x2
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << -1 
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << x1 
   << std::setw(OWID) << y0
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y2
   << std::setw(OWID) << z0      
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl
    
   << std::setw(OWID) << 1
   << std::setw(OWID) << 0
   << std::setw(OWID) << -1
   << std::setw(OWID) << 0 
   << std::setw(OWID) << x0      
   << std::setw(OWID) << y1
   << std::setw(OWID) << z0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::endl << std::endl; 

   ofs.close();
   }


   // create boundary particles and springs connecting those boundary particles
   void Assembly::createMemParticle(REAL rRadius,
   bool toRebuild,
   const char *ParticleFile,
   const char *allParticle)
   {
   if (toRebuild) readParticle(ParticleFile);

   REAL radius = gradation.getMinPtclRadius();
   if (gradation.getSize().size() == 1 &&
   gradation.getPtclRatioBA() == 1.0 && 
   gradation.getPtclRatioCA() == 1.0)
   radius *= rRadius; // determine how tiny the boundary particles are
   REAL diameter = radius*2;
   Vec v1 = allContainer.getMinCorner();
   Vec v2 = allContainer.getMaxCorner();
   Vec v0 = allContainer.getCenter();
   REAL x1 = v1.getX();
   REAL y1 = v1.getY();
   REAL z1 = v1.getZ();
   REAL x2 = v2.getX();
   REAL y2 = v2.getY();
   REAL z2 = v2.getZ();
   REAL x0 = v0.getX();
   REAL y0 = v0.getY();
   REAL z0 = v0.getZ();

   Particle* newptcl = NULL;
   REAL x, y, z;
  
   ParticlePArray vec1d;  // 1-dimension
   std::vector< ParticlePArray  > vec2d; // 2-dimension
   Spring* newSpring = NULL;
   int memPtclIndex = trimHistoryNum;
   // process in the order of surfaces: x1 x2 y1 y2 z1 z2
   // surface x1
   x = x1 - radius;
   for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
   vec1d.clear();
   for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }

   // surface x2     
   vec2d.clear();
   x = x2 + radius;
   for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
   vec1d.clear();
   for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }
 
   // surface y1
   vec2d.clear();
   y = y1 - radius;
   for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
   vec1d.clear();
   for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }

   // surface y2
   vec2d.clear();
   y = y2 + radius;
   for (z = z1 + radius; z <= z2 - radius + EPS; z += diameter) {
   vec1d.clear();
   for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }

   // surface z1
   vec2d.clear();
   z = z1 - radius;
   for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
   vec1d.clear();
   for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }

   // surface z2
   vec2d.clear();
   z = z2 + radius;
   for (y = y1 + radius; y <= y2 - radius + EPS; y += diameter) {
   vec1d.clear();
   for (x = x1 + radius; x <= x2 - radius + EPS; x += diameter) {
   newptcl = new Particle(++memPtclIndex, 5, Vec(x,y,z), radius, young, poisson);
   vec1d.push_back(newptcl);
   particleVec.push_back(newptcl);
   }
   vec2d.push_back(vec1d);
   }
   memBoundary.push_back(vec2d);
   for (int i = 0; i != vec2d.size() ; ++i)
   for (int j = 0; j != vec2d[i].size() ; ++ j) {
   if (j + 1 < vec2d[i].size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i][j+1], memYoung);
   springVec.push_back(newSpring);
   }
   if (i + 1 < vec2d.size() ) {
   newSpring = new Spring(*vec2d[i][j], *vec2d[i+1][j], memYoung);
   springVec.push_back(newSpring);
   }
   }

   // membrane particles at the edges of each surface, for example,
   // x1y1 means particles on surface x1 connecting to particles on surface y1
   ParticlePArray x1y1;
   ParticlePArray x1y2;
   ParticlePArray x1z1;
   ParticlePArray x1z2;

   ParticlePArray x2y1;
   ParticlePArray x2y2;
   ParticlePArray x2z1;
   ParticlePArray x2z2;

   ParticlePArray y1x1;
   ParticlePArray y1x2;
   ParticlePArray y1z1;
   ParticlePArray y1z2;

   ParticlePArray y2x1;
   ParticlePArray y2x2;
   ParticlePArray y2z1;
   ParticlePArray y2z2;

   ParticlePArray z1x1;
   ParticlePArray z1x2;
   ParticlePArray z1y1;
   ParticlePArray z1y2;

   ParticlePArray z2x1;
   ParticlePArray z2x2;
   ParticlePArray z2y1;
   ParticlePArray z2y2;

   // find edge particles for each surface
   // memBoundary[0, 1, 2, 3, 4, 5] correspond to 
   // surface     x1 x2 y1 y2 z1 z2 respectively
   // surface x1
   vec2d.clear();
   vec2d = memBoundary[0];
   x1z1  = vec2d[0];
   x1z2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   x1y1.push_back(vec2d[i][0]);
   x1y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }
   // surface x2
   vec2d.clear();
   vec2d = memBoundary[1];
   x2z1  = vec2d[0];
   x2z2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   x2y1.push_back(vec2d[i][0]);
   x2y2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }
   // surface y1
   vec2d.clear();
   vec2d = memBoundary[2];
   y1z1  = vec2d[0];
   y1z2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   y1x1.push_back(vec2d[i][0]);
   y1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }
   // surface y2
   vec2d.clear();
   vec2d = memBoundary[3];
   y2z1  = vec2d[0];
   y2z2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   y2x1.push_back(vec2d[i][0]);
   y2x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }
   // surface z1
   vec2d.clear();
   vec2d = memBoundary[4];
   z1y1  = vec2d[0];
   z1y2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   z1x1.push_back(vec2d[i][0]);
   z1x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }
   // surface z2
   vec2d.clear();
   vec2d = memBoundary[5];
   z2y1  = vec2d[0];
   z2y2  = vec2d[vec2d.size() - 1];
   for (int i = 0; i < vec2d.size(); ++i) {
   z2x1.push_back(vec2d[i][0]);
   z2x2.push_back(vec2d[i][ vec2d[i].size() - 1 ]);  
   }

   // create springs connecting 12 edges of a cube
   // 4 edges on surface x1
   assert(x1y1.size() == y1x1.size());
   for (int i = 0; i < x1y1.size(); ++i) {
   newSpring = new Spring(*x1y1[i], *y1x1[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x1y2.size() == y2x1.size());
   for (int i = 0; i < x1y2.size(); ++i) {
   newSpring = new Spring(*x1y2[i], *y2x1[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x1z1.size() == z1x1.size());
   for (int i = 0; i < x1z1.size(); ++i) {
   newSpring = new Spring(*x1z1[i], *z1x1[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x1z2.size() == z2x1.size());
   for (int i = 0; i < x1z2.size(); ++i) {
   newSpring = new Spring(*x1z2[i], *z2x1[i], memYoung);
   springVec.push_back(newSpring);
   }
   // 4 edges on surface x2  
   assert(x2y1.size() == y1x2.size());
   for (int i = 0; i < x2y1.size(); ++i) {
   newSpring = new Spring(*x2y1[i], *y1x2[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x2y2.size() == y2x2.size());
   for (int i = 0; i < x2y2.size(); ++i) {
   newSpring = new Spring(*x2y2[i], *y2x2[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x2z1.size() == z1x2.size());
   for (int i = 0; i < x2z1.size(); ++i) {
   newSpring = new Spring(*x2z1[i], *z1x2[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(x2z2.size() == z2x2.size());
   for (int i = 0; i < x2z2.size(); ++i) {
   newSpring = new Spring(*x2z2[i], *z2x2[i], memYoung);
   springVec.push_back(newSpring);
   }
   // 2 edges on surface y1 
   assert(y1z1.size() == z1y1.size());
   for (int i = 0; i < y1z1.size(); ++i) {
   newSpring = new Spring(*y1z1[i], *z1y1[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(y1z2.size() == z2y1.size());
   for (int i = 0; i < y1z2.size(); ++i) {
   newSpring = new Spring(*y1z2[i], *z2y1[i], memYoung);
   springVec.push_back(newSpring);
   }
   // 2 edges on surface y2
   assert(y2z1.size() == z1y2.size());
   for (int i = 0; i < y2z1.size(); ++i) {
   newSpring = new Spring(*y2z1[i], *z1y2[i], memYoung);
   springVec.push_back(newSpring);
   }
   assert(y2z2.size() == z2y2.size());
   for (int i = 0; i < y2z2.size(); ++i) {
   newSpring = new Spring(*y2z2[i], *z2y2[i], memYoung);
   springVec.push_back(newSpring);
   }

   printParticle(allParticle);

   }


   void Assembly::TrimPtclBdryByHeight(REAL height,
   const char *iniptclfile,
   const char *ParticleFile)
   {
   readParticle(iniptclfile);
  
   ParticlePArray::iterator itr;
   for (itr = particleVec.begin(); itr != particleVec.end(); ) {
   if ( (*itr)->getType() == 1 ) { // 1-fixed
   Vec center=(*itr)->getCurrPos();
   if(center.getZ() > height)
   {
   delete (*itr); // release memory
   itr = particleVec.erase(itr); 
   }
   else {
   (*itr)->setType(10); // 10-ghost
   ++itr;
   }
   }
   }
  
   printParticle(ParticleFile);
   }


   void Assembly::deGravitation(int   totalSteps,  
   int   snapNum,
   int   interval,
   bool  toRebuild,
   const char *iniptclfile,   
   const char *ParticleFile, 
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   // pre_1: open streams for output.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: deGravitation" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf.precision(OPREC);
   progressInf << std::setw(OWID) << "iteration"
   << std::setw(OWID) << "poss_contact"
   << std::setw(OWID) << "actual_contact"
   << std::setw(OWID) << "penetration"
   << std::setw(OWID) << "avg_normal"
   << std::setw(OWID) << "avg_tangt"
   << std::setw(OWID) << "avg_velocity"
   << std::setw(OWID) << "avg_omga"
   << std::setw(OWID) << "avg_force"
   << std::setw(OWID) << "avg_moment"
   << std::setw(OWID) << "trans_energy"
   << std::setw(OWID) << "rotat_energy"
   << std::setw(OWID) << "kinet_energy"
   << std::endl;
  
   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: deGravitation" << std::endl; exit(-1); }
   debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
   // pre_2. create particles from existing files.
   if (toRebuild) readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
  
   // pre_3. define variables used in iterations
   REAL avgNormal=0;
   REAL avgTangt=0;
   int  stepsnum=0;
   char stepsstr[4];
   char stepsfp[50];

   // iterations starting ...
   iteration=0; 
   do
   {
   // 1. find contacts between particles.
   findContact();
      
   // 2. set particles' forces/moments as zero before each re-calculation,
   clearContactForce();	
      
   // 3. calculate contact forces/moments and apply them to particles.
   internalForce(avgNormal, avgTangt);
      
   // 4. update particles' velocity/omga/position/orientation based on force/moment.
   updateParticle();
      
   // 5. (1) output particles and contacts information as snapNum.
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   time(&timeStamp);
   g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
   ++stepsnum;
   }
      
   // 5. (2) output stress and strain info.
   if (iteration % interval == 0) {
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << getTransEnergy()
   << std::setw(OWID) << getRotatEnergy()
   << std::setw(OWID) << getKinetEnergy()
   << std::endl;
   }
      
   // 7. loop break conditions.
   if (contactVec.size() == 0) break;
      
   } while (++iteration < totalSteps);
  
   // post_1. store the final snapshot of particles & contacts.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);
  
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
   g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
  
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // actual deposit function for the case of fixed particle boundaries
   void Assembly::deposit_p(int   totalSteps,  
   int   snapNum,
   int   interval,
   REAL dimn,
   REAL rsize,
   const char *iniptclfile,   
   const char *ParticleFile, 
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   // pre_1: open streams for output.
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: deposit_p" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "deposit..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: deposit_p" << std::endl; exit(-1); }
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from existing files.
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

   // pre_3: define variables used in iterations.
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
   REAL void_ratio=0;

   // iterations starting ...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles.
   findContact();

   // 2. set particles' forces/moments as zero before each re-calculation,
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles.
   internalForce(avgNormal, avgTangt);

   // 4. update particles' velocity/omga/position/orientation based on force/moment.
   updateParticle();

   // 5. calculate specimen void ratio.
   l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
   l24=dimn*rsize;
   l13=dimn*rsize;
   bulkVolume=l13*l24*l56;
   void_ratio=bulkVolume/getParticleVolume()-1;

   // 6. (1) output particles and contacts information as snapNum.
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	    
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 6. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;
   }

   // 7. loop break conditions.


   } while (++iteration < totalSteps);
    
   // post_1. store the final snapshot of particles & contacts.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);

   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // squeeze paticles inside a container by moving the boundaries
   void Assembly::squeeze(int   totalSteps,  
   int   init_steps,
   int   snapNum,
   int   interval,
   int   flag,
   const char *iniptclfile,   
   const char *inibdryfile,
   const char *ParticleFile, 
   const char *boundaryfile,
   const char *contactfile,
   const char *progressfile, 
   const char *debugfile)
   {
   // pre_1: open streams for output.
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: squeeze" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "deposit..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: squeeze" << std::endl; exit(-1); }
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from existing files.
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
   readBoundary(inibdryfile);   // create boundaries.

   // pre_3: define variables used in iterations.
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];

   int         mid[2]={1,3};    // boundary 1 and 3
   UPDATECTL   midctl[2];
   REAL void_ratio=0;
   REAL bdry_penetr[7];
   int         bdry_cntnum[7];
   for (int i=0;i<7;++i) {
   bdry_penetr[i]=0; bdry_cntnum[i]=0;
   }

   // iterations starting ...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles.
   findContact();
   findBdryContact();

   // 2. set particles' forces/moments as zero before each re-calculation,
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles.
   internalForce(avgNormal, avgTangt);

   // 4. calculate boundary forces/moments and apply them to particles.
   boundaryForce(bdry_penetr, bdry_cntnum);
	
   // 5. update particles' velocity/omga/position/orientation based on force/moment.
   updateParticle();

   // 6. calculate sample void ratio.
   l56=getTopFreeParticlePosition().getZ() -getApt(6).getZ();
   l24=getApt(2).getY()-getApt(4).getY();
   l13=getApt(1).getX()-getApt(3).getX(); bulkVolume=l13*l24*l56;
   void_ratio=bulkVolume/getParticleVolume()-1;

   // displacement control
   if (iteration > init_steps) {
   if (flag==1) // loosen, totally remove the wall
   midctl[0].tran=Vec(timeStep*1.0e+0*flag,0,0);
   else         // squeeze
   midctl[0].tran=Vec(timeStep*5.0e-3*flag,0,0);
   updateBoundary(mid,midctl,2);
   }

   // 7. (1) output particles and contacts information as snapNum.
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	    
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output stress and strain info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*(getActualContactNum()
   +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
   +bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
   << std::endl;
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << bdry_penetr[1]
   << std::setw(OWID) << bdry_penetr[2]
   << std::setw(OWID) << bdry_penetr[3]
   << std::setw(OWID) << bdry_penetr[4]
   << std::setw(OWID) << bdry_penetr[6]
   << std::setw(OWID) << bdry_cntnum[1]
   << std::setw(OWID) << bdry_cntnum[2]
   << std::setw(OWID) << bdry_cntnum[3]
   << std::setw(OWID) << bdry_cntnum[4]
   << std::setw(OWID) << bdry_cntnum[6]
   << std::endl;

   }

   // 8. loop break conditions.

   } while (++iteration < totalSteps);
    
   // post_1. store the final snapshot of particles & contacts.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);

   strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
   printBoundary(stepsfp);

   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   void Assembly::iso_MemBdry(int   totalSteps,  
   int   snapNum, 
   int   interval,
   REAL  sigma3,
   REAL  rRadius,
   bool  toRebuild,
   const char *iniptclfile, 
   const char *ParticleFile,
   const char *contactfile, 
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile);
   if(!progressInf) { debugInf << "stream error: isoMemBdry" << std::endl; exit(-1);}
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf.precision(OPREC);
   progressInf << std::setw(OWID) << "iteration"
   << std::setw(OWID) << "poss_contact"
   << std::setw(OWID) << "actual_contact"
   << std::setw(OWID) << "penetration"
   << std::setw(OWID) << "avg_normal"
   << std::setw(OWID) << "avg_tangt"
   << std::setw(OWID) << "avg_velocity"
   << std::setw(OWID) << "avg_omga"
   << std::setw(OWID) << "avg_force"
   << std::setw(OWID) << "avg_moment"
   << std::setw(OWID) << "trans_energy"
   << std::setw(OWID) << "rotat_energy"
   << std::setw(OWID) << "kinet_energy"
   << std::endl;
  
   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: isoMemBdry" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);
  
   // pre_2. create particles from file and calculate forces caused by hydraulic pressure
   if (toRebuild) readParticle(iniptclfile);
   REAL radius = gradation.getMinPtclRadius();
   if (gradation.getSize().size() == 1 &&
   gradation.getPtclRatioBA() == 1.0 && 
   gradation.getPtclRatioCA() == 1.0)
   radius *= rRadius; // determine how tiny the boundary particles are
   REAL mag  = radius*radius*4*sigma3;
   Vec v1 = allContainer.getMinCorner();
   Vec v2 = allContainer.getMaxCorner();
   REAL x1 = v1.getX();
   REAL y1 = v1.getY();
   REAL z1 = v1.getZ();
   REAL x2 = v2.getX();
   REAL y2 = v2.getY();
   REAL z2 = v2.getZ();
   ParticlePArray::const_iterator  it;
   Vec pos;
   for (it=particleVec.begin();it!=particleVec.end();++it)
   {
   pos = (*it)->getCurrPos();
   if (pos.getX() < x1)
   (*it)->setConstForce( Vec(mag, 0, 0) );
   else if (pos.getX() > x2)
   (*it)->setConstForce( Vec(-mag, 0, 0) );
   else if (pos.getY() < y1)
   (*it)->setConstForce( Vec(0, mag, 0) );
   else if (pos.getY() > y2)
   (*it)->setConstForce( Vec(0, -mag, 0) );
   else if (pos.getZ() < z1)
   (*it)->setConstForce( Vec(0, 0, mag) );
   else if (pos.getZ() > z2)
   (*it)->setConstForce( Vec(0, 0, -mag) );
   }

   // pre_3. define variables used in iterations
   REAL avgNormal=0;
   REAL avgTangt=0;
   int  stepsnum=0;
   char stepsstr[4];
   char stepsfp[50];
  
   // iterations start here...
   iteration=0;
   do 
   {
   // 1. find contacts between particles
   findContact();

   // 2. set particles forces/moments to zero before each re-calculation
   clearContactForce(); // const_force/moment NOT cleared by this call	
      
   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);

   // 4. calculate and apply spring forces to boundary particles
   springForce();
      
   // 5. update particles velocity/omga/position/orientation based on force/moment
   updateParticle();
      
   // 6. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);

   strcpy(stepsfp,"iso_membrane_"); strcat(stepsfp, stepsstr);
   printMemParticle(stepsfp);
   strcpy(stepsfp,"iso_spring_"); strcat(stepsfp, stepsstr); strcat(stepsfp, ".dat");
   plotSpring(stepsfp);
	
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr); 
   printContact(stepsfp);
   time(&timeStamp);
   g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) << std::flush;
   ++stepsnum;
   }
      
   // 6. (2) output stress and strain info
   if (iteration % interval == 0 ) {
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << getTransEnergy()
   << std::setw(OWID) << getRotatEnergy()
   << std::setw(OWID) << getKinetEnergy()
   << std::endl;
   }
	  
   } while (++iteration < totalSteps);
  
   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, "iso_membrane_end");
   printMemParticle(stepsfp);
   strcpy(stepsfp, "iso_spring_end.dat");
   plotSpring(stepsfp);  

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
   g_timeinf << std::setw(4) << "end" << " " << ctime(&timeStamp) << std::flush;
  
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // This function initializes triaxial sample to a certain confining pressure.
   void Assembly::triaxialPtclBdryIni(int   totalSteps,  
   int   snapNum, 
   int   interval,
   REAL  sigma,
   const char *iniptclfile, 
   const char *inibdryfile,
   const char *ParticleFile,
   const char *boundaryfile,
   const char *contactfile, 
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile);
   if(!progressInf) { debugInf << "stream error: triaxialPtclBdryIni" << std::endl; exit(-1);}
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "triaxial..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average        sample            sample     "
   << "     sample          sample          sample          sample          sample          "
   << "sample          sample         sample           sample          sample          sample     "
   << "     sample          sample          sample          void            sample        coordinate"
   << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "          omga            force           moment        density          "
   << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
   << "        ratio          porosity         number"
   << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: triaxialPtclBdryIni" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
   readBoundary(inibdryfile);   // create boundaries

   // pre_3. define variables used in iterations
   REAL H0 = getApt(5).getZ()-getApt(6).getZ();
   REAL l56= 0;
   REAL sigma3_1, sigma3_2;
   REAL epsilon_h;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
    
   int         min[2]={5,6};    // boundary 5 and 6
   UPDATECTL   minctl[2];

   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();
   findBdryContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. calculate boundary forces/moments and apply them to particles
   boundaryForce();

   // 5. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 6. update boundaries' position and orientation
   sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

   // force control
   if (sigma3_1 < sigma)
   minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
   else
   minctl[0].tran=Vec(0,0, timeStep*boundaryRate);

   if (sigma3_2 < sigma)
   minctl[1].tran=Vec(0,0, timeStep*boundaryRate);
   else
   minctl[1].tran=Vec(0,0,-timeStep*boundaryRate);

   updateBoundary(min,minctl,2);
	
   // 7. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	    
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output stress and strain info
   l56=getApt(5).getZ()-getApt(6).getZ();
   epsilon_h = (H0-l56)/H0;
   if (iteration % interval == 0 ) {
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << getDensity()
   << std::setw(OWID) << 0 << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::setw(OWID) << 0
   << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
   << std::setw(OWID) << getAvgPressure()
   << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << l56
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << epsilon_h
   << std::setw(OWID) << epsilon_h
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;

   }

   // 9. loop break condition: through displacement control mechanism
   if (   fabs(sigma3_1-sigma)/sigma < boundaryStressTol && fabs(sigma3_2-sigma)/sigma < boundaryStressTol )
   break;
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);

   strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
   printBoundary(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // This function performs triaxial compression test.
   // Displacement boundaries are used in axial direction.
   void Assembly::triaxialPtclBdry(int   totalSteps,  
   int   snapNum, 
   int   interval,
   const char *iniptclfile, 
   const char *inibdryfile,
   const char *ParticleFile,
   const char *boundaryfile,
   const char *contactfile, 
   const char *progressfile,
   const char *balancedfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile);
   if(!progressInf) { debugInf << "stream error: triaxialPtclBdry" << std::endl; exit(-1);}
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "triaxial..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average        sample            sample     "
   << "     sample          sample          sample          sample          sample          "
   << "sample          sample         sample           sample          sample          sample     "
   << "     sample          sample          sample          void            sample        coordinate"
   << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "          omga            force           moment        density          "
   << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
   << "        ratio          porosity         number"
   << std::endl;

   std::ofstream balancedinf(balancedfile);
   if(!balancedinf) { debugInf << "stream error: triaxialPtclBdry" << std::endl; exit(-1);}
   balancedinf.setf(std::ios::scientific, std::ios::floatfield);
   balancedinf << "triaxial..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average        sample            sample     "
   << "     sample          sample          sample          sample          sample          "
   << "sample          sample         sample           sample          sample          sample     "
   << "     sample          sample          sample          void            sample        coordinate"
   << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "          omga            force           moment        density          "
   << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       epsilon_v"
   << "        ratio          porosity         number"
   << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: triaxialPtclBdry" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
   readBoundary(inibdryfile);   // create boundaries

   // pre_3. define variables used in iterations
   REAL H0 = getApt(5).getZ()-getApt(6).getZ();
   REAL l56= 0;
   REAL sigma3_1, sigma3_2;
   REAL epsilon_h;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
    
   int         min[2]={5,6};    // boundary 5 and 6
   UPDATECTL   minctl[2];

   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();
   findBdryContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. calculate boundary forces/moments and apply them to particles
   boundaryForce();

   // 5. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 6. update boundaries' position and orientation
   sigma3_1=vfabs(getNormalForce(5))/2.5e-3; sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

   // displacement control
   if(iteration < 100001) {
   minctl[0].tran=Vec(0,0,-timeStep*boundaryRate);
   minctl[1].tran=Vec(0,0, timeStep*boundaryRate);

   updateBoundary(min,minctl,2);
   }
   // 7. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
	    
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output stress and strain info
   l56=getApt(5).getZ()-getApt(6).getZ();
   epsilon_h = (H0-l56)/H0;
   if (iteration % interval == 0 ) {
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << getDensity()
   << std::setw(OWID) << 0 << std::setw(OWID) << 0
   << std::setw(OWID) << 0 << std::setw(OWID) << 0
   << std::setw(OWID) << sigma3_1 << std::setw(OWID) << sigma3_2
   << std::setw(OWID) << getAvgPressure()
   << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << l56
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << epsilon_h
   << std::setw(OWID) << epsilon_h
   << std::setw(OWID) << 0
   << std::setw(OWID) << 0
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;

   }

   // 9. loop break condition: through displacement control mechanism
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);

   strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
   printBoundary(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   balancedinf.close();
   debugInf.close();
   }


   // The specimen has been deposited with gravitation within boundaries composed of particles.
   // A rectangular pile is then drived into the particles using displacement control.
   void Assembly::rectPile_Disp(int   totalSteps,  
   int   snapNum, 
   int   interval,
   const char *iniptclfile,  
   const char *inibdryfile,
   const char *ParticleFile, 
   const char *boundaryfile,
   const char *contactfile,  
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error:recPile_Disp" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "pile penetrate..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential        total           sample           sample     "
   << "     sample          sample          sample          sample          sample          "
   << "sample          sample         sample           sample          sample          sample"
   << "          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy          energy          density         "
   << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: recPile_Disp" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);
   debugInf << " iteration    end_bearing     side_friction   total_force" << std::endl;

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 
   readBoundary(inibdryfile);   // create boundaries

   // pre_3. define variables used in iterations
   int    stepsnum=0;
   char   stepsstr[4];
   char   stepsfp[50];
   REAL avgNormal=0;
   REAL avgTangt=0;
    
   int pile[2]={11,12}; // top and down boundaries
   UPDATECTL pilectl[2];

   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();
   findBdryContact();
	
   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. calculate boundary forces/moments and apply them to particles
   boundaryForce();

   // 5. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 6. update boundaries' position and orientation

   // displacement control of the pile
   pilectl[0].tran=Vec(0,0,-timeStep*pileRate);
   pilectl[1].tran=Vec(0,0,-timeStep*pileRate);

   updateBoundary(pile, pilectl, 2); 
   updateRectPile();
   if (iteration % interval == 0) {
   REAL  f7=getShearForce( 7).getZ();
   REAL  f8=getShearForce( 8).getZ();
   REAL  f9=getShearForce( 9).getZ();
   REAL f10=getShearForce(10).getZ();
   REAL  fn=getNormalForce(12).getZ();
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << fn
   << std::setw(OWID) << (f7+f8+f9+f10)
   << std::setw(OWID) << (fn+f7+f8+f9+f10)
   << std::endl;
   }

   // 7. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);
   printRectPile(stepsfp);
	    
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3) << std::endl;
   }

   // 8. loop break condition
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);
   printRectPile(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);

   strcpy(stepsfp, boundaryfile); strcat(stepsfp, "_end");
   printBoundary(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // The specimen has been deposited with gravitation within boundaries composed of particles.
   // An ellipsoidal pile is then drived into the particles using displacement control.
   void Assembly::ellipPile_Disp(int   totalSteps,  
   int   snapNum, 
   int   interval,
   REAL dimn,
   REAL rsize,
   const char *iniptclfile,
   const char *ParticleFile, 
   const char *contactfile,  
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: ellipPile_Disp" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "pile penetrate..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: ellipPile_Disp" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

   // pre_3. define variables used in iterations
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
   REAL void_ratio=0;
    
   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 5. calculate specimen void ratio.
   l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
   l24=dimn*rsize;
   l13=dimn*rsize;
   bulkVolume=l13*l24*l56-ellipPilePeneVol();
   void_ratio=bulkVolume/getParticleVolume()-1;

   // 6. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);

   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 6. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()
   << std::setw(OWID) << ellipPileTipZ()
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
   << std::setw(OWID) << l13*l24*l56
   << std::setw(OWID) << ellipPilePeneVol()
   << std::setw(OWID) << bulkVolume
   << std::endl;
   }

   // 7. loop break condition
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // The specimen has been deposited with gravitation within rigid boundaries.
   // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
   void Assembly::ellipPile_Impact(int   totalSteps,  
   int   snapNum, 
   int   interval,
   REAL dimn,
   const char *iniptclfile,
   const char *inibdryfile,
   const char *ParticleFile, 
   const char *contactfile,  
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: ellipPile_Impact" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "penetrator impact..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: ellipPile_Impact" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles
   readBoundary(inibdryfile);   // create boundaries.

   // pre_3. define variables used in iterations
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
   REAL void_ratio=0;
   REAL bdry_penetr[7];
   int         bdry_cntnum[7];
   for (int i=0;i<7;++i) {
   bdry_penetr[i]=0; bdry_cntnum[i]=0;
   }
    
   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);

   // 4. calculate boundary forces/moments and apply them to particles.
   boundaryForce(bdry_penetr, bdry_cntnum);
	
   // 5. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 6. calculate specimen void ratio.
   l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
   l24=dimn;
   l13=dimn;
   bulkVolume=l13*l24*l56-ellipPilePeneVol();
   void_ratio=bulkVolume/getParticleVolume()-1;

   // 7. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);

   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*(getActualContactNum()
   +bdry_cntnum[1]+bdry_cntnum[2]+bdry_cntnum[3]
   +bdry_cntnum[4]+bdry_cntnum[6])/allParticleVec.size()
   << std::endl;
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << bdry_penetr[1]
   << std::setw(OWID) << bdry_penetr[2]
   << std::setw(OWID) << bdry_penetr[3]
   << std::setw(OWID) << bdry_penetr[4]
   << std::setw(OWID) << bdry_penetr[6]
   << std::setw(OWID) << bdry_cntnum[1]
   << std::setw(OWID) << bdry_cntnum[2]
   << std::setw(OWID) << bdry_cntnum[3]
   << std::setw(OWID) << bdry_cntnum[4]
   << std::setw(OWID) << bdry_cntnum[6]
   << std::endl;

   }

   // 8. loop break condition
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }


   // The specimen has been deposited with gravitation within particle boundaries.
   // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
   void Assembly::ellipPile_Impact_p(int   totalSteps,  
   int   snapNum, 
   int   interval,
   REAL dimn,
   const char *iniptclfile,
   const char *ParticleFile, 
   const char *contactfile,  
   const char *progressfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: ellipPile_Impact_p" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "penetrator impact..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: ellipPile_Impact_p" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles

   // pre_3. define variables used in iterations
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
   REAL void_ratio=0;
    
   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();
	
   // 5. calculate specimen void ratio.
   l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
   l24=dimn;
   l13=dimn;
   bulkVolume=l13*l24*l56-ellipPilePeneVol();
   void_ratio=bulkVolume/getParticleVolume()-1;

   // 6. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);

   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 6. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()
   << std::setw(OWID) << ellipPileTipZ()
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
   << std::setw(OWID) << l13*l24*l56
   << std::setw(OWID) << ellipPilePeneVol()
   << std::setw(OWID) << bulkVolume
   << std::endl;
   }

   // 7. loop break condition
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   debugInf.close();
   }



   // The specimen has been deposited with gravitation within boundaries composed of particles.
   // An ellipsoidal pile is then drived into the particles using force control.
   // Not recommended.
   void Assembly::ellipPile_Force(int   totalSteps,  
   int   snapNum,
   int   interval,
   REAL dimn,
   REAL force,
   int   division,
   const char *iniptclfile,
   const char *ParticleFile, 
   const char *contactfile,  
   const char *progressfile,
   const char *balancedfile,
   const char *debugfile) 
   {
   // pre_1: open streams for output
   // ParticleFile and contactfile are used for snapNum at the end.
   progressInf.open(progressfile); 
   if(!progressInf) { debugInf << "stream error: ellipPile_Force" << std::endl; exit(-1); }
   progressInf.setf(std::ios::scientific, std::ios::floatfield);
   progressInf << "pile penetrate..." << std::endl
   << "     iteration possible  actual      average	    average         average         average"
   << "         average         average         average       translational    rotational       "
   << "kinetic        potential         total           void            sample       coordination"
   << "       sample           sample          sample          sample          sample          sample"
   << "          sample          sample          sample         sample           sample         "
   << " sample          sample          sample          sample          sample" << std::endl
   << "       number  contacts contacts   penetration   contact_normal  contact_tangt     velocity"
   << "         omga            force           moment         energy           energy          "
   << "energy         energy            energy          ratio          porosity         number       "
   << "   density         sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
   << "sigma3_1        sigma3_2           p             width          length           "
   << "height          volume         epsilon_w       epsilon_l       epsilon_h       "
   << "epsilon_v" << std::endl;

   std::ofstream balancedinf(balancedfile);
   if(!balancedinf) { debugInf << "stream error: ellipPile_Force" << std::endl; exit(-1);}
   balancedinf.setf(std::ios::scientific, std::ios::floatfield);
   balancedinf << "pile penetrate..." << std::endl
   << "   iteration   apply_force    pile_tip_pos     pile_force" << std::endl;

   debugInf.open(debugfile);
   if(!debugInf) { debugInf << "stream error: ellipPile_Force" << std::endl; exit(-1);}
   debugInf.setf(std::ios::scientific, std::ios::floatfield);

   // pre_2. create particles and boundaries from files
   readParticle(iniptclfile); // create container and particles, velocity and omga are set zero. 

   // pre_3. define variables used in iterations
   REAL l13, l24, l56;
   REAL avgNormal=0;
   REAL avgTangt=0;
   int         stepsnum=0;
   char        stepsstr[4];
   char        stepsfp[50];
   REAL void_ratio=0;

   REAL zforce_inc=force/division;
   REAL zforce=zforce_inc;

   // iterations start here...
   iteration=0;
   do
   {
   // 1. create possible boundary particles and contacts between particles
   findContact();

   // 2. set particles' forces/moments as zero before each re-calculation
   clearContactForce();	

   // 3. calculate contact forces/moments and apply them to particles
   internalForce(avgNormal, avgTangt);
	
   // 4. update particles' velocity/omga/position/orientation based on force/moment
   updateParticle();

   // 5. calculate specimen void ratio.
   l56=getTopFreeParticlePosition().getZ() - (-dimn/2);
   l24=dimn;
   l13=dimn;
   bulkVolume=l13*l24*l56-ellipPilePeneVol();
   void_ratio=bulkVolume/getParticleVolume()-1;
	
   // 6. update pile external force and position
   if(zforce>ellipPileForce())
   ellipPileUpdate();

   if(fabs(ellipPileForce()-zforce)/zforce < boundaryStressTol ) {
   balancedinf << std::setw(OWID) << iteration
   << std::setw(OWID) << zforce
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
   << std::setw(OWID) << ellipPileForce()
   << std::endl;
   zforce += zforce_inc;
   }

   if( iteration % interval == 0) {
   debugInf << std::setw(OWID) << iteration
   << std::setw(OWID) << zforce
   << std::setw(OWID) << getTopFreeParticlePosition().getZ()-ellipPileTipZ()
   << std::setw(OWID) << ellipPileForce()
   << std::endl;
   }

   // 7. (1) output particles and contacts information
   if (iteration % (totalSteps/snapNum) == 0) {
   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printParticle(stepsfp);

   sprintf(stepsstr, "%03d", stepsnum); 
   strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
   printContact(stepsfp);
   ++stepsnum;
   }

   // 7. (2) output statistics info.
   if (iteration % interval == 0) {
   REAL t1=getTransEnergy();
   REAL t2=getRotatEnergy();
   REAL t3=getPotenEnergy(-0.025);
   progressInf << std::setw(OWID) << iteration
   << std::setw(OWID) << getPossContactNum()
   << std::setw(OWID) << getActualContactNum()
   << std::setw(OWID) << getAvgPenetration()
   << std::setw(OWID) << avgNormal
   << std::setw(OWID) << avgTangt
   << std::setw(OWID) << getAvgVelocity() 
   << std::setw(OWID) << getAvgOmga()
   << std::setw(OWID) << getAvgForce()   
   << std::setw(OWID) << getAvgMoment()
   << std::setw(OWID) << t1
   << std::setw(OWID) << t2
   << std::setw(OWID) << (t1+t2)
   << std::setw(OWID) << t3
   << std::setw(OWID) << (t1+t2+t3)
   << std::setw(OWID) << void_ratio
   << std::setw(OWID) << void_ratio/(1+void_ratio)
   << std::setw(OWID) << 2.0*getActualContactNum()/allParticleVec.size()
   << std::endl;
   }

   // 8. loop break condition
   if (fabs((zforce-force)/force)<0.001)
   break;
	
   } while (++iteration < totalSteps);

   // post_1. store the final snapshot of particles, contacts and boundaries.
   strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_end");
   printParticle(stepsfp);

   strcpy(stepsfp, contactfile); strcat(stepsfp, "_end");
   printContact(stepsfp);
    
   // post_2. close streams
   progressInf.close();
   balancedinf.close();
   debugInf.close();
   }

 */ 
 // rule out

