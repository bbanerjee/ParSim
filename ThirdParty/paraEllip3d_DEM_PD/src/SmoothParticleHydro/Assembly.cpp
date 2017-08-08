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
//    where x1 < x2, y1 < y2, z1 < z2. We also use surface 1, 2, 3, 4, 5, 6
//    accordingly.
//

#include "Assembly.h"
#include "const.h"
#include <cassert>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <utility>

//#define BINNING
//#define TIME_PROFILE

static time_t timeStamp;                // for file timestamping
static struct timeval time_w1, time_w2; // for wall-clock time record
static struct timeval time_p1, time_p2; // for internal wall-clock time
                                        // profiling, can be used on any piece
                                        // of code
static struct timeval time_r1, time_r2; // for internal wall-clock time
                                        // profiling for contact resolution only
                                        // (excluding space search)

namespace dem {

struct timeval
timediff(const struct timeval& time1, const struct timeval& time2)
{
  struct timeval diff;
  diff.tv_sec = time2.tv_sec - time1.tv_sec;
  if ((diff.tv_usec = time2.tv_usec - time1.tv_usec) < 0) {
    diff.tv_usec += 1000000;
    --diff.tv_sec;
  }
  return (diff);
}

long int
timediffmsec(const struct timeval& time1, const struct timeval& time2)
{
  struct timeval diff = timediff(time1, time2);
  return (diff.tv_sec * 1000000 + diff.tv_usec);
}

REAL
timediffsec(const struct timeval& time1, const struct timeval& time2)
{
  return ((REAL)timediffmsec(time1, time2) / 1.0e+6);
}

char*
combineString(char* cstr, const char* str, std::size_t num, std::size_t width)
{
  std::string obj(str);
  std::stringstream ss;
  ss << std::setw(width) << std::setfill('0') << std::right << num;
  obj += ss.str();
  return strcpy(cstr, obj.c_str());
}

// input:   number percentage smaller from data file
// output:  mass percentage smaller to disk file debugInf
// purpose: let mass percentage smaller satisfy particle size distribution curve
// method:  use trial and error method on number percentage until mass
// percentage is satisfied
void
Assembly::tuneMassPercent()
{
  if (mpiRank == 0) {
    REAL minX = dem::Parameter::getSingleton().parameter["minX"];
    REAL minY = dem::Parameter::getSingleton().parameter["minY"];
    REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
    REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
    REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
    REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
    std::size_t particleLayers =
      dem::Parameter::getSingleton().parameter["particleLayers"];

    setContainer(Rectangle(minX, minY, minZ, maxX, maxY, maxZ));

    buildBoundary(5, "deposit_boundary_ini");

    std::size_t sieveNum = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["sieveNum"]);
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL>>& grada =
      dem::Parameter::getSingleton().gradation;
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
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      massPercent[i] = 0;

    for (std::vector<Particle*>::iterator itr = allParticleVec.begin();
         itr != allParticleVec.end(); ++itr)
      for (int i = massPercent.size() - 1; i >= 0;
           --i) { // do not use size_t for descending series
        if ((*itr)->getA() <= massSize[i])
          massPercent[i] += (*itr)->getMass();
      }
    REAL totalMass = massPercent[0];
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      massPercent[i] /= totalMass;
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
  }
}

// input:   particle file with estimated gradation
// output:  particle file with precise gradation
void
Assembly::calcMassPercent()
{
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());

    // statistics of mass distribution
    Gradation massGrad = gradation;
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      massPercent[i] = 0;

    for (std::vector<Particle*>::iterator itr = allParticleVec.begin();
         itr != allParticleVec.end(); ++itr)
      for (int i = massPercent.size() - 1; i >= 0;
           --i) { // do not use size_t for descending series
        if ((*itr)->getA() <= massSize[i])
          massPercent[i] +=
            (*itr)->getMass(); // += 1 for calculating number percentage
      }
    REAL totalMass = massPercent[0];
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      massPercent[i] /= totalMass;
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
  }
}

void
Assembly::depositIntoContainer()
{
  if (mpiRank == 0) {
    REAL minX = dem::Parameter::getSingleton().parameter["minX"];
    REAL minY = dem::Parameter::getSingleton().parameter["minY"];
    REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
    REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
    REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
    REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
    std::size_t particleLayers =
      dem::Parameter::getSingleton().parameter["particleLayers"];

    setContainer(Rectangle(minX, minY, minZ, maxX, maxY, maxZ));

    buildBoundary(5, "deposit_boundary_ini");

    std::size_t sieveNum = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["sieveNum"]);
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL>>& grada =
      dem::Parameter::getSingleton().gradation;
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

  deposit("deposit_boundary_ini", "float_particle_ini");

  if (mpiRank == 0) {
    setContainer(Rectangle(
      allContainer.getMinCorner().getX(), allContainer.getMinCorner().getY(),
      allContainer.getMinCorner().getZ(), allContainer.getMaxCorner().getX(),
      allContainer.getMaxCorner().getY(),
      dem::Parameter::getSingleton().parameter["trimHeight"]));
    buildBoundary(6, "trim_boundary_ini");
    char cstr[50];
    std::size_t endSnap = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["endSnap"]);
    trim(false, combineString(cstr, "deposit_particle_", endSnap, 3),
         "trim_particle_ini");
  }
}

void
Assembly::resumeDepositIntoContainer()
{
  deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());

  if (mpiRank == 0) {
    setContainer(Rectangle(
      allContainer.getMinCorner().getX(), allContainer.getMinCorner().getY(),
      allContainer.getMinCorner().getZ(), allContainer.getMaxCorner().getX(),
      allContainer.getMaxCorner().getY(),
      dem::Parameter::getSingleton().parameter["trimHeight"]));
    buildBoundary(6, "trim_boundary_ini");
    char cstr[50];
    std::size_t endSnap = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["endSnap"]);
    trim(false, combineString(cstr, "deposit_particle_", endSnap, 3),
         "trim_particle_ini");
  }
}

void
Assembly::proceedFromPreset()
{
  deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());
}

void
Assembly::deposit(const char* inputBoundary, const char* inputParticle)
{
  if (mpiRank == 0) {
    readBoundary(inputBoundary);
    readParticle(inputParticle);
    openDepositProg(progressInf, "deposit_progress");
  }
  scatterParticle(); // scatter particles only once; also updates grid for the
                     // first time
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  /**/ REAL timeCount = 0;
  /**/ timeAccrued =
    static_cast<REAL>(dem::Parameter::getSingleton().parameter["timeAccrued"]);
  /**/ REAL timeIncr = timeStep * netStep;
  /**/ REAL timeTotal = timeAccrued + timeStep * netStep;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "deposit_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(combineString(cstr0, "deposit_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    printParticle(combineString(cstr0, "deposit_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "deposit_bdrycntc_", iterSnap - 1, 3));
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "compuT"
             << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%"
             << std::endl;
  /**/ while (timeAccrued < timeTotal) {
    // while (iteration <= endStep) {
    bool toCheckTime = (iteration + 1) % (netStep / netSnap) == 0;

    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    if (toCheckTime)
      time2 = MPI_Wtime();
    commuT = time2 - time0;

    /**/ calcTimeStep(); // use values from last step, must call before
                         // findContact (which clears data)
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    dragForce();

    updateParticle();
    updateGridMaxZ();
    // updateGridMaxZ() is for deposition or explosion. If they go out of side
    // walls, particles are discarded.
    // updateGrid() updates all six directions, thus side walls may "disappear"
    // if particles go far out of side walls
    // and cause some grids to extrude out of side walls.

    /**/ timeCount += timeStep;
    /**/ timeAccrued += timeStep;
    /**/ if (timeCount >= timeIncr / netSnap) {
      // if (iteration % (netStep / netSnap) == 0) {
      if (toCheckTime)
        time1 = MPI_Wtime();
      gatherParticle();
      gatherBdryContact();
      gatherEnergy();
      if (toCheckTime)
        time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "deposit_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(combineString(cstr, "deposit_gridplot_", iterSnap, 3),
                        ".dat"));
        printParticle(combineString(cstr, "deposit_particle_", iterSnap, 3));
        printBdryContact(combineString(cstr, "deposit_bdrycntc_", iterSnap, 3));
        printDepositProg(progressInf);
      }
      printContact(combineString(cstr, "deposit_contact_", iterSnap, 3));

      /**/ timeCount = 0;
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    if (toCheckTime)
      time1 = MPI_Wtime();
    migrateParticle();
    if (toCheckTime)
      time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        toCheckTime) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID)
               << totalT - commuT - migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;
    ++iteration;
  }

  if (mpiRank == 0)
    closeProg(progressInf);
}

void
Assembly::expandCavityParticle()
{
  if (mpiRank == 0) {
    const char* inputParticle =
      dem::Parameter::getSingleton().datafile["particleFile"].c_str();
    REAL percent = dem::Parameter::getSingleton().parameter["expandPercent"];
    readParticle(inputParticle);

    REAL x1 = dem::Parameter::getSingleton().parameter["cavityMinX"];
    REAL y1 = dem::Parameter::getSingleton().parameter["cavityMinY"];
    REAL z1 = dem::Parameter::getSingleton().parameter["cavityMinZ"];
    REAL x2 = dem::Parameter::getSingleton().parameter["cavityMaxX"];
    REAL y2 = dem::Parameter::getSingleton().parameter["cavityMaxY"];
    REAL z2 = dem::Parameter::getSingleton().parameter["cavityMaxZ"];

    std::vector<Particle*> cavityParticleVec;
    std::vector<Particle*>::iterator it;
    Vec center;

    for (it = allParticleVec.begin(); it != allParticleVec.end(); ++it) {
      center = (*it)->getCurrPos();
      if (center.getX() > x1 && center.getX() < x2 && center.getY() > y1 &&
          center.getY() < y2 && center.getZ() > z1 && center.getZ() < z2) {
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

void
Assembly::resumeExpandCavityParticle()
{
  deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());
}

void
Assembly::isotropic()
{
  std::size_t isotropicType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["isotropicType"]);
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openCompressProg(progressInf, "isotropic_progress");
    openCompressProg(balancedInf, "isotropic_balanced");
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL sigmaEnd, sigmaInc, sigmaVar;
  std::size_t sigmaDiv;

  sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
  sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];
  std::vector<REAL>& sigmaPath = dem::Parameter::getSingleton().sigmaPath;
  std::size_t sigma_i = 0;

  if (isotropicType == 1)
    sigmaVar = sigmaEnd;
  else if (isotropicType == 2) {
    REAL sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
    sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
    sigmaVar = sigmaStart;
  } else if (isotropicType == 3) {
    sigmaVar = sigmaPath[sigma_i];
    sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
    sigmaEnd = sigmaPath[sigmaPath.size() - 1];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "isotropic_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(
      combineString(cstr0, "isotropic_gridplot_", iterSnap - 1, 3), ".dat"));
    printParticle(combineString(cstr0, "isotropic_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "isotropic_bdrycntc_", iterSnap - 1, 3));
    printBoundary(combineString(cstr0, "isotropic_boundary_", iterSnap - 1, 3));
    getStartDimension(distX, distY, distZ);
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    calcTimeStep(); // use values from last step, must call before findContact
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    gatherBdryContact(); // must call before updateBoundary
    updateBoundary(sigmaVar, "isotropic");
    updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "isotropic_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(combineString(cstr, "isotropic_gridplot_", iterSnap, 3),
                        ".dat"));
        printParticle(combineString(cstr, "isotropic_particle_", iterSnap, 3));
        printBdryContact(
          combineString(cstr, "isotropic_bdrycntc_", iterSnap, 3));
        printBoundary(combineString(cstr, "isotropic_boundary_", iterSnap, 3));
        printCompressProg(progressInf, distX, distY, distZ);
      }
      printContact(combineString(cstr, "isotropic_contact_", iterSnap, 3));
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

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
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
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
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_i + 1]) {
          sigmaVar = sigmaPath[++sigma_i];
          sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
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

void
Assembly::odometer()
{
  std::size_t odometerType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["odometerType"]);
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openCompressProg(progressInf, "odometer_progress");
    openCompressProg(balancedInf, "odometer_balanced");
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL sigmaEnd, sigmaInc, sigmaVar;
  std::size_t sigmaDiv;

  sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
  sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];
  std::vector<REAL>& sigmaPath = dem::Parameter::getSingleton().sigmaPath;
  std::size_t sigma_i = 0;

  if (odometerType == 1) {
    REAL sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
    sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
    sigmaVar = sigmaStart;
  } else if (odometerType == 2) {
    sigmaVar = sigmaPath[sigma_i];
    sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
    sigmaEnd = sigmaPath[sigmaPath.size() - 1];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "odometer_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(combineString(cstr0, "odometer_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    printParticle(combineString(cstr0, "odometer_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "odometer_bdrycntc_", iterSnap - 1, 3));
    printBoundary(combineString(cstr0, "odometer_boundary_", iterSnap - 1, 3));
    getStartDimension(distX, distY, distZ);
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    calcTimeStep(); // use values from last step, must call before findContact
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    gatherBdryContact(); // must call before updateBoundary
    updateBoundary(sigmaVar, "odometer");
    updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "odometer_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(combineString(cstr, "odometer_gridplot_", iterSnap, 3),
                        ".dat"));
        printParticle(combineString(cstr, "odometer_particle_", iterSnap, 3));
        printBdryContact(
          combineString(cstr, "odometer_bdrycntc_", iterSnap, 3));
        printBoundary(combineString(cstr, "odometer_boundary_", iterSnap, 3));
        printCompressProg(progressInf, distX, distY, distZ);
      }
      printContact(combineString(cstr, "odometer_contact_", iterSnap, 3));
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (odometerType == 1) {
      if (tractionErrorTol(sigmaVar, "odometer")) {
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
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
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_i + 1]) {
          sigmaVar = sigmaPath[++sigma_i];
          sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
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

void
Assembly::triaxial()
{
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openCompressProg(progressInf, "triaxial_progress");
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  REAL sigmaConf = dem::Parameter::getSingleton().parameter["sigmaConf"];
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "triaxial_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(combineString(cstr0, "triaxial_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    printParticle(combineString(cstr0, "triaxial_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "triaxial_bdrycntc_", iterSnap - 1, 3));
    printBoundary(combineString(cstr0, "triaxial_boundary_", iterSnap - 1, 3));
    getStartDimension(distX, distY, distZ);
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    // displacement control relies on constant time step, so do not call
    // calcTimeStep().
    // calcTimeStep(); // use values from last step, must call before
    // findContact
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    gatherBdryContact(); // must call before updateBoundary
    updateBoundary(sigmaConf, "triaxial");
    updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "triaxial_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(combineString(cstr, "triaxial_gridplot_", iterSnap, 3),
                        ".dat"));
        printParticle(combineString(cstr, "triaxial_particle_", iterSnap, 3));
        printBdryContact(
          combineString(cstr, "triaxial_bdrycntc_", iterSnap, 3));
        printBoundary(combineString(cstr, "triaxial_boundary_", iterSnap, 3));
        // printCompressProg(progressInf, distX, distY, distZ); // redundant
      }
      printContact(combineString(cstr, "triaxial_contact_", iterSnap, 3));
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

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

  if (mpiRank == 0)
    closeProg(progressInf);
}

void
Assembly::planeStrain()
{
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openCompressProg(progressInf, "plnstrn_progress");
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  REAL sigmaConf = dem::Parameter::getSingleton().parameter["sigmaConf"];
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "plnstrn_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(combineString(cstr0, "plnstrn_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    printParticle(combineString(cstr0, "plnstrn_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "plnstrn_bdrycntc_", iterSnap - 1, 3));
    printBoundary(combineString(cstr0, "plnstrn_boundary_", iterSnap - 1, 3));
    getStartDimension(distX, distY, distZ);
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    // displacement control relies on constant time step, so do not call
    // calcTimeStep().
    // calcTimeStep(); // use values from last step, must call before
    // findContact
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    gatherBdryContact(); // must call before updateBoundary
    updateBoundary(sigmaConf, "plnstrn");
    updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "plnstrn_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(combineString(cstr, "plnstrn_gridplot_", iterSnap, 3),
                        ".dat"));
        printParticle(combineString(cstr, "plnstrn_particle_", iterSnap, 3));
        printBdryContact(combineString(cstr, "plnstrn_bdrycntc_", iterSnap, 3));
        printBoundary(combineString(cstr, "plnstrn_boundary_", iterSnap, 3));
        // printCompressProg(progressInf, distX, distY, distZ); // redundant
      }
      printContact(combineString(cstr, "plnstrn_contact_", iterSnap, 3));
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

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

  if (mpiRank == 0)
    closeProg(progressInf);
}

void
Assembly::trueTriaxial()
{
  std::size_t trueTriaxialType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["trueTriaxialType"]);
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openCompressProg(progressInf, "trueTriaxial_progress");
    openCompressProg(balancedInf, "trueTriaxial_balanced");
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL sigmaStart, sigmaEndZ, sigmaEndX, sigmaEndY;
  REAL sigmaDiv, sigmaIncZ, sigmaIncX, sigmaIncY, sigmaVarZ, sigmaVarX,
    sigmaVarY;
  REAL sigmaInit[3], sigmaEnd, sigmaInc, sigmaVar;
  std::size_t changeDirc;
  sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];

  if (trueTriaxialType == 1) {
    sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
    sigmaEndZ = dem::Parameter::getSingleton().parameter["sigmaEndZ"];
    sigmaEndX = dem::Parameter::getSingleton().parameter["sigmaEndX"];
    sigmaEndY = dem::Parameter::getSingleton().parameter["sigmaEndY"];
    sigmaIncZ = (sigmaEndZ - sigmaStart) / sigmaDiv;
    sigmaIncX = (sigmaEndX - sigmaStart) / sigmaDiv;
    sigmaIncY = (sigmaEndY - sigmaStart) / sigmaDiv;
    sigmaVarZ = sigmaStart;
    sigmaVarX = sigmaStart;
    sigmaVarY = sigmaStart;
  } else if (trueTriaxialType == 2) {
    sigmaInit[0] = dem::Parameter::getSingleton().parameter["sigmaStartX"];
    sigmaInit[1] = dem::Parameter::getSingleton().parameter["sigmaStartY"];
    sigmaInit[2] = dem::Parameter::getSingleton().parameter["sigmaStartZ"];
    sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
    changeDirc = dem::Parameter::getSingleton().parameter["changeDirc"];
    sigmaInc = (sigmaEnd - sigmaInit[changeDirc]) / sigmaDiv;
    sigmaVar = sigmaInit[changeDirc];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "trueTriaxial_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(
      combineString(cstr0, "trueTriaxial_gridplot_", iterSnap - 1, 3), ".dat"));
    printParticle(
      combineString(cstr0, "trueTriaxial_particle_", iterSnap - 1, 3));
    printBdryContact(
      combineString(cstr0, "trueTriaxial_bdrycntc_", iterSnap - 1, 3));
    printBoundary(
      combineString(cstr0, "trueTriaxial_boundary_", iterSnap - 1, 3));
    getStartDimension(distX, distY, distZ);
  }
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    calcTimeStep(); // use values from last step, must call before findContact
    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();
    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    gatherBdryContact(); // must call before updateBoundary

    if (trueTriaxialType == 1)
      updateBoundary(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY);
    else if (trueTriaxialType == 2) {
      REAL sigmaX, sigmaY, sigmaZ;
      if (changeDirc == 0) {
        sigmaX = sigmaVar;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaVar;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaVar;
      }
      updateBoundary(sigmaZ, "trueTriaxial", sigmaX, sigmaY);
    }

    updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(strcat(
          combineString(cstr, "trueTriaxial_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(strcat(
          combineString(cstr, "trueTriaxial_gridplot_", iterSnap, 3), ".dat"));
        printParticle(
          combineString(cstr, "trueTriaxial_particle_", iterSnap, 3));
        printBdryContact(
          combineString(cstr, "trueTriaxial_bdrycntc_", iterSnap, 3));
        printBoundary(
          combineString(cstr, "trueTriaxial_boundary_", iterSnap, 3));
        printCompressProg(progressInf, distX, distY, distZ);
      }
      printContact(combineString(cstr, "trueTriaxial_contact_", iterSnap, 3));
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (trueTriaxialType == 1) {
      if (tractionErrorTol(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY)) {
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
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
      REAL sigmaX, sigmaY, sigmaZ;
      if (changeDirc == 0) {
        sigmaX = sigmaVar;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaVar;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaVar;
      }
      if (tractionErrorTol(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
        if (mpiRank == 0)
          printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }

      if (changeDirc == 0) {
        sigmaX = sigmaEnd;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaEnd;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaEnd;
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

bool
Assembly::tractionErrorTol(REAL sigma, std::string type, REAL sigmaX,
                           REAL sigmaY)
{
  // sigma implies sigmaZ
  REAL tol = dem::Parameter::getSingleton().parameter["tractionErrorTol"];

  std::map<std::string, REAL> normalForce;
  REAL x1, x2, y1, y2, z1, z2;
  // do not use mergeBoundaryVec because each process calls this function.
  for (std::vector<Boundary*>::const_iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    Vec normal = (*it)->getNormalForce();
    Vec point = (*it)->getPoint();
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
    return (fabs(normalForce["x1"] / areaX - sigma) / sigma <= tol &&
            fabs(normalForce["x2"] / areaX - sigma) / sigma <= tol &&
            fabs(normalForce["y1"] / areaY - sigma) / sigma <= tol &&
            fabs(normalForce["y2"] / areaY - sigma) / sigma <= tol &&
            fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);

  else if (type.compare("odometer") == 0)
    return (fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);

  else if (type.compare("triaxial") == 0)
    return true; // always near equilibrium

  else if (type.compare("trueTriaxial") == 0)
    return (fabs(normalForce["x1"] / areaX - sigmaX) / sigmaX <= tol &&
            fabs(normalForce["x2"] / areaX - sigmaX) / sigmaX <= tol &&
            fabs(normalForce["y1"] / areaY - sigmaY) / sigmaY <= tol &&
            fabs(normalForce["y2"] / areaY - sigmaY) / sigmaY <= tol &&
            fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);
}

void
Assembly::coupleWithGas()
{
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    openDepositProg(progressInf, "couple_progress");
    /*1*/ fluid.initParameter(allContainer, gradation);
    /*2*/ fluid.initialize();
  }
  scatterParticle();
  calcNeighborRanks();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL timeCount = 0;
  timeAccrued =
    static_cast<REAL>(dem::Parameter::getSingleton().parameter["timeAccrued"]);
  REAL timeIncr = timeStep * netStep;
  REAL timeTotal = timeAccrued + timeIncr;
  if (mpiRank == 0) {
    plotBoundary(strcat(
      combineString(cstr0, "couple_bdryplot_", iterSnap - 1, 3), ".dat"));
    plotGrid(strcat(combineString(cstr0, "couple_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    printParticle(combineString(cstr0, "couple_particle_", iterSnap - 1, 3));
    printBdryContact(combineString(cstr0, "couple_bdrycntc_", iterSnap - 1, 3));
    /*3*/ fluid.plot(strcat(
      combineString(cstr0, "couple_fluidplot_", iterSnap - 1, 3), ".dat"));
  }
  /*
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" <<
  std::setw(OWID) << "migraT"
             << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" <<
  std::endl;
  */
  while (timeAccrued < timeTotal) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    calcTimeStep(); // use values from last step, must call before findContact

    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();

    /*4*/ fluid.getPtclInfo(particleVec); // not allParticleVec
    /*5*/ fluid.runOneStep(particleVec);
    /*6*/ fluid.calcPtclForce(particleVec); // not allParticleVec
    /*7*/ fluid.penalize(particleVec);
    // fluid.checkMomentum(particleVec);

    internalForce();
    if (isBdryProcess())
      boundaryForce();

    updateParticle();
    updateGridMaxZ();

    timeCount += timeStep;
    // timeAccrued += timeStep; // note fluid.runOneStep() might change timeStep
    // and print timeAccrued
    if (timeCount >= timeIncr / netSnap) {
      time1 = MPI_Wtime();
      gatherParticle();
      gatherBdryContact();
      gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (mpiRank == 0) {
        plotBoundary(
          strcat(combineString(cstr, "couple_bdryplot_", iterSnap, 3), ".dat"));
        plotGrid(
          strcat(combineString(cstr, "couple_gridplot_", iterSnap, 3), ".dat"));
        printParticle(combineString(cstr, "couple_particle_", iterSnap, 3));
        printBdryContact(combineString(cstr, "couple_bdrycntc_", iterSnap, 3));
        printDepositProg(progressInf);
        /*8*/ fluid.plot(strcat(
          combineString(cstr, "couple_fluidplot_", iterSnap, 3), ".dat"));
      }
      printContact(combineString(cstr, "couple_contact_", iterSnap, 3));

      timeCount = 0;
      ++iterSnap;
    }

    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    time1 = MPI_Wtime();
    migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    /*
    if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore
    gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT <<
    std::setw(OWID) << migraT
               << std::setw(OWID) << totalT << std::setw(OWID) << (commuT +
    migraT)/totalT*100 << std::endl;
    */

    ++iteration;
  }

  if (mpiRank == 0)
    closeProg(progressInf);
}

// particleLayers:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
void
Assembly::generateParticle(std::size_t particleLayers, const char* genParticle)
{
  REAL young = dem::Parameter::getSingleton().parameter["young"];
  REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

  REAL x, y, z;
  Particle* newptcl;
  std::size_t particleNum = 0;
  REAL diameter = gradation.getPtclMaxRadius() * 2.0;

  REAL offset = 0;
  REAL edge = diameter;
  if (gradation.getSize().size() == 1 && gradation.getPtclRatioBA() == 1.0 &&
      gradation.getPtclRatioCA() == 1.0) {
    edge = diameter * 2.0;
    offset = diameter * 0.25;
  }

  REAL x1 = allContainer.getMinCorner().getX() + edge;
  REAL y1 = allContainer.getMinCorner().getY() + edge;
  REAL z1 = allContainer.getMinCorner().getZ() + diameter;
  REAL x2 = allContainer.getMaxCorner().getX() - edge;
  REAL y2 = allContainer.getMaxCorner().getY() - edge;
  // REAL z2 = allContainer.getMaxCorner().getZ() - diameter;
  REAL z2 = dem::Parameter::getSingleton().parameter["floatMaxZ"] - diameter;
  REAL x0 = allContainer.getCenter().getX();
  REAL y0 = allContainer.getCenter().getY();
  REAL z0 = allContainer.getCenter().getZ();

  if (particleLayers == 0) { // just one free particle
    newptcl = new Particle(particleNum + 1, 0, Vec(x0, y0, z0), gradation,
                           young, poisson);
    allParticleVec.push_back(newptcl);
    particleNum++;
  } else if (particleLayers == 1) { // a horizontal layer of free particles
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
        newptcl = new Particle(particleNum + 1, 0, Vec(x, y, z0), gradation,
                               young, poisson);
        allParticleVec.push_back(newptcl);
        particleNum++;
      }
  } else if (particleLayers == 2) { // multiple layers of free particles
    for (z = z1; z - z2 < EPS; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
        //	  for (y = y1 + offset; y - y2 < EPS; y += diameter) {
        y = 0;
      newptcl = new Particle(particleNum + 1, 0, Vec(x, y, z), gradation, young,
                             poisson);
      allParticleVec.push_back(newptcl);
      particleNum++;
      //	  }
      offset *= -1;
    }
  }

  printParticle(genParticle);
}

void
Assembly::trimOnly()
{
  if (mpiRank == 0) {
    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    trim(true, dem::Parameter::getSingleton().datafile["particleFile"].c_str(),
         "trim_particle_end");
  }
}

void
Assembly::trim(bool toRebuild, const char* inputParticle,
               const char* trmParticle)
{
  if (toRebuild)
    readParticle(inputParticle);
  trimHistoryNum = allParticleVec.size();

  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.getX();
  REAL y1 = v1.getY();
  REAL z1 = v1.getZ();
  REAL x2 = v2.getX();
  REAL y2 = v2.getY();
  REAL z2 = v2.getZ();
  REAL maxR = gradation.getPtclMaxRadius();

  std::vector<Particle*>::iterator itr;
  Vec center;

  for (itr = allParticleVec.begin(); itr != allParticleVec.end();) {
    center = (*itr)->getCurrPos();
    if (center.getX() < x1 || center.getX() > x2 || center.getY() < y1 ||
        center.getY() > y2 || center.getZ() < z1 || center.getZ() + maxR > z2) {
      delete (*itr); // release memory
      itr = allParticleVec.erase(itr);
    } else
      ++itr;
  }

  printParticle(trmParticle);
}

void
Assembly::removeBySphere()
{
  if (mpiRank == 0) {

    readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    REAL minR = gradation.getPtclMinRadius();

    REAL x0S = dem::Parameter::getSingleton().parameter["x0S"];
    REAL y0S = dem::Parameter::getSingleton().parameter["y0S"];
    REAL z0S = dem::Parameter::getSingleton().parameter["z0S"];
    REAL r0S = dem::Parameter::getSingleton().parameter["r0S"];

    std::vector<Particle*>::iterator itr;
    Vec center;
    REAL dist;

    for (itr = allParticleVec.begin(); itr != allParticleVec.end();) {
      center = (*itr)->getCurrPos();
      dist = sqrt(pow(center.getX() - x0S, 2) + pow(center.getY() - y0S, 2) +
                  pow(center.getZ() - z0S, 2));
      if (dist <= r0S + minR) {
        delete (*itr); // release memory
        itr = allParticleVec.erase(itr);
      } else
        ++itr;
    }

    printParticle("remove_particle_end");
  }
}

void
Assembly::findParticleInRectangle(const Rectangle& container,
                                  const std::vector<Particle*>& inputParticle,
                                  std::vector<Particle*>& foundParticle)
{
  foundParticle.reserve(inputParticle.size());
  Vec v1 = container.getMinCorner();
  Vec v2 = container.getMaxCorner();
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
  std::vector<Particle*>(foundParticle).swap(foundParticle);
}

REAL
Assembly::getPtclMaxX(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL x0 = (*it)->getCurrPos().getX();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getX() > x0)
      x0 = (*it)->getCurrPos().getX();
  }
  return x0;
}

REAL
Assembly::getPtclMinX(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL x0 = (*it)->getCurrPos().getX();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getX() < x0)
      x0 = (*it)->getCurrPos().getX();
  }
  return x0;
}

REAL
Assembly::getPtclMaxY(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL y0 = (*it)->getCurrPos().getY();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getY() > y0)
      y0 = (*it)->getCurrPos().getY();
  }
  return y0;
}

REAL
Assembly::getPtclMinY(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL y0 = (*it)->getCurrPos().getY();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getY() < y0)
      y0 = (*it)->getCurrPos().getY();
  }
  return y0;
}

REAL
Assembly::getPtclMaxZ(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL z0 = (*it)->getCurrPos().getZ();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getZ() > z0)
      z0 = (*it)->getCurrPos().getZ();
  }
  return z0;
}

REAL
Assembly::getPtclMinZ(const std::vector<Particle*>& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  std::vector<Particle*>::const_iterator it = inputParticle.begin();
  REAL z0 = (*it)->getCurrPos().getZ();
  for (; it != inputParticle.end(); ++it) {
    if ((*it)->getCurrPos().getZ() < z0)
      z0 = (*it)->getCurrPos().getZ();
  }
  return z0;
}

void
Assembly::setCommunicator(boost::mpi::communicator& comm)
{
  boostWorld = comm;
  mpiWorld = MPI_Comm(comm);
  mpiProcX =
    static_cast<int>(dem::Parameter::getSingleton().parameter["mpiProcX"]);
  mpiProcY =
    static_cast<int>(dem::Parameter::getSingleton().parameter["mpiProcY"]);
  mpiProcZ =
    static_cast<int>(dem::Parameter::getSingleton().parameter["mpiProcZ"]);

  // create Cartesian virtual topology (unavailable in boost.mpi)
  int ndim = 3;
  int dims[3] = { mpiProcX, mpiProcY, mpiProcZ };
  int periods[3] = { 0, 0, 0 };
  int reorder = 0; // mpiRank not reordered
  MPI_Cart_create(mpiWorld, ndim, dims, periods, reorder, &cartComm);
  MPI_Comm_rank(cartComm, &mpiRank);
  MPI_Comm_size(cartComm, &mpiSize);
  MPI_Cart_coords(cartComm, mpiRank, ndim, mpiCoords);
  mpiTag = 0;
  assert(mpiRank == boostWorld.rank());
  // debugInf << mpiRank << " " << mpiCoords[0] << " " << mpiCoords[1] << " " <<
  // mpiCoords[2] << std::endl;

  for (int iRank = 0; iRank < mpiSize; ++iRank) {
    int ndim = 3;
    int coords[3];
    MPI_Cart_coords(cartComm, iRank, ndim, coords);
    if (coords[0] == 0 || coords[0] == mpiProcX - 1 || coords[1] == 0 ||
        coords[1] == mpiProcY - 1 || coords[2] == 0 ||
        coords[2] == mpiProcZ - 1)
      bdryProcess.push_back(iRank);
  }
}

void
Assembly::scatterParticle()
{
  // partition particles and send to each process
  if (mpiRank == 0) { // process 0
    setGrid(
      Rectangle(grid.getMinCorner().getX(), grid.getMinCorner().getY(),
                grid.getMinCorner().getZ(), grid.getMaxCorner().getX(),
                grid.getMaxCorner().getY(),
                getPtclMaxZ(allParticleVec) + gradation.getPtclMaxRadius()));

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;

    boost::mpi::request* reqs = new boost::mpi::request[mpiSize - 1];
    std::vector<Particle*> tmpParticleVec;
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
                          v1.getZ() +
                            vspan.getZ() / mpiProcZ * (coords[2] + 1));
      findParticleInRectangle(container, allParticleVec, tmpParticleVec);
      if (iRank != 0)
        reqs[iRank - 1] =
          boostWorld.isend(iRank, mpiTag, tmpParticleVec); // non-blocking send
      if (iRank == 0) {
        particleVec.resize(tmpParticleVec.size());
        for (int i = 0; i < particleVec.size(); ++i)
          particleVec[i] = new Particle(
            *tmpParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec);
  }

  // content of allParticleVec may need to be printed, so do not clear it.
  // if (mpiRank == 0) releaseGatheredParticle();

  // broadcast necessary info
  broadcast(boostWorld, gradation, 0);
  broadcast(boostWorld, boundaryVec, 0);
  broadcast(boostWorld, allContainer, 0);
  broadcast(boostWorld, grid, 0);
}

bool
Assembly::isBdryProcess()
{
  return (mpiCoords[0] == 0 || mpiCoords[0] == mpiProcX - 1 ||
          mpiCoords[1] == 0 || mpiCoords[1] == mpiProcY - 1 ||
          mpiCoords[2] == 0 || mpiCoords[2] == mpiProcZ - 1);
}

void
Assembly::calcNeighborRanks()
{ // the ranks of the neighbor blocks just need to be calculated once before
  // time loop

  // find neighboring blocks
  rankX1 = -1;
  rankX2 = -1;
  rankY1 = -1;
  rankY2 = -1;
  rankZ1 = -1;
  rankZ2 = -1;
  rankX1Y1 = -1;
  rankX1Y2 = -1;
  rankX1Z1 = -1;
  rankX1Z2 = -1;
  rankX2Y1 = -1;
  rankX2Y2 = -1;
  rankX2Z1 = -1;
  rankX2Z2 = -1;
  rankY1Z1 = -1;
  rankY1Z2 = -1;
  rankY2Z1 = -1;
  rankY2Z2 = -1;
  rankX1Y1Z1 = -1;
  rankX1Y1Z2 = -1;
  rankX1Y2Z1 = -1;
  rankX1Y2Z2 = -1;
  rankX2Y1Z1 = -1;
  rankX2Y1Z2 = -1;
  rankX2Y2Z1 = -1;
  rankX2Y2Z2 = -1;
  // x1: -x direction
  int neighborCoords[3] = { mpiCoords[0], mpiCoords[1], mpiCoords[2] };
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
  --neighborCoords[0];
  --neighborCoords[1];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1);
  // x1y2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  ++neighborCoords[1];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2);
  // x1z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z1);
  // x1z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Z2);
  // x2y1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  --neighborCoords[1];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1);
  // x2y2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  ++neighborCoords[1];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2);
  // x2z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z1);
  // x2z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Z2);
  // y1z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z1);
  // y1z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankY1Z2);
  // y2z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z1);
  // y2z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankY2Z2);
  // x1y1z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  --neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z1);
  // x1y1z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  --neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y1Z2);
  // x1y2z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  ++neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z1);
  // x1y2z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  --neighborCoords[0];
  ++neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX1Y2Z2);
  // x2y1z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  --neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z1);
  // x2y1z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  --neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y1Z2);
  // x2y2z1
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  ++neighborCoords[1];
  --neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z1);
  // x2y2z2
  neighborCoords[0] = mpiCoords[0];
  neighborCoords[1] = mpiCoords[1];
  neighborCoords[2] = mpiCoords[2];
  ++neighborCoords[0];
  ++neighborCoords[1];
  ++neighborCoords[2];
  MPI_Cart_rank(cartComm, neighborCoords, &rankX2Y2Z2);
}

void
Assembly::updateGrid()
{
  updateGridMinX();
  updateGridMaxX();
  updateGridMinY();
  updateGridMaxY();
  updateGridMinZ();
  updateGridMaxZ();
}

void
Assembly::updateGridMinX()
{
  REAL pMinX = getPtclMinX(particleVec);
  REAL minX = 0;
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Rectangle(minX - gradation.getPtclMaxRadius(),
                    grid.getMinCorner().getY(), grid.getMinCorner().getZ(),
                    grid.getMaxCorner().getX(), grid.getMaxCorner().getY(),
                    grid.getMaxCorner().getZ()));
}

void
Assembly::updateGridMaxX()
{
  REAL pMaxX = getPtclMaxX(particleVec);
  REAL maxX = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setGrid(Rectangle(grid.getMinCorner().getX(), grid.getMinCorner().getY(),
                    grid.getMinCorner().getZ(),
                    maxX + gradation.getPtclMaxRadius(),
                    grid.getMaxCorner().getY(), grid.getMaxCorner().getZ()));
}

void
Assembly::updateGridMinY()
{
  REAL pMinY = getPtclMinY(particleVec);
  REAL minY = 0;
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Rectangle(grid.getMinCorner().getX(),
                    minY - gradation.getPtclMaxRadius(),
                    grid.getMinCorner().getZ(), grid.getMaxCorner().getX(),
                    grid.getMaxCorner().getY(), grid.getMaxCorner().getZ()));
}

void
Assembly::updateGridMaxY()
{
  REAL pMaxY = getPtclMaxY(particleVec);
  REAL maxY = 0;
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setGrid(Rectangle(grid.getMinCorner().getX(), grid.getMinCorner().getY(),
                    grid.getMinCorner().getZ(), grid.getMaxCorner().getX(),
                    maxY + gradation.getPtclMaxRadius(),
                    grid.getMaxCorner().getZ()));
}

void
Assembly::updateGridMinZ()
{
  REAL pMinZ = getPtclMinZ(particleVec);
  REAL minZ = 0;
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Rectangle(grid.getMinCorner().getX(), grid.getMinCorner().getY(),
                    minZ - gradation.getPtclMaxRadius(),
                    grid.getMaxCorner().getX(), grid.getMaxCorner().getY(),
                    grid.getMaxCorner().getZ()));
}

void
Assembly::updateGridMaxZ()
{
  // update compute grids adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setGrid(Rectangle(grid.getMinCorner().getX(), grid.getMinCorner().getY(),
                    grid.getMinCorner().getZ(), grid.getMaxCorner().getX(),
                    grid.getMaxCorner().getY(),
                    maxZ + gradation.getPtclMaxRadius()));
}

void
Assembly::printBdryContact(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: printBdryContact" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    (*it)->printContactInfo(ofs);
  }

  ofs.close();
}

void
Assembly::gatherEnergy()
{
  calcTransEnergy();
  calcRotatEnergy();
  calcKinetEnergy();
  calcGraviEnergy(allContainer.getMinCorner().getZ());
  calcMechaEnergy();
}

void
Assembly::closeProg(std::ofstream& ofs)
{
  ofs.close();
}

void
Assembly::openDepositProg(std::ofstream& ofs, const char* str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openDepositProg" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "normal_x1"
      << std::setw(OWID) << "normal_x2" << std::setw(OWID) << "normal_y1"
      << std::setw(OWID) << "normal_y2" << std::setw(OWID) << "normal_z1"
      << std::setw(OWID) << "normal_z2"

      << std::setw(OWID) << "contact_x1" << std::setw(OWID) << "contact_x2"
      << std::setw(OWID) << "contact_y1" << std::setw(OWID) << "contact_y2"
      << std::setw(OWID) << "contact_z1" << std::setw(OWID) << "contact_z2"
      << std::setw(OWID) << "contact_inside"

      << std::setw(OWID) << "penetr_x1" << std::setw(OWID) << "penetr_x2"
      << std::setw(OWID) << "penetr_y1" << std::setw(OWID) << "penetr_y2"
      << std::setw(OWID) << "penetr_z1" << std::setw(OWID) << "penetr_z2"

      << std::setw(OWID) << "avgNormal" << std::setw(OWID) << "avgShear"
      << std::setw(OWID) << "avgPenetr"

      << std::setw(OWID) << "transEnergy" << std::setw(OWID) << "rotatEnergy"
      << std::setw(OWID) << "kinetEnergy" << std::setw(OWID) << "graviEnergy"
      << std::setw(OWID) << "mechaEnergy"

      << std::setw(OWID) << "vibra_est_dt" << std::setw(OWID) << "impact_est_dt"
      << std::setw(OWID) << "actual_dt" << std::setw(OWID) << "accruedTime"

      << std::endl;
}

void
Assembly::printDepositProg(std::ofstream& ofs)
{
  REAL var[6];

  // normalForce
  for (std::size_t i = 0; i < 6; ++i)
    var[i] = 0;
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
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
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getContactNum();
  }
  for (std::size_t i = 0; i < 6; ++i)
    ofs << std::setw(OWID) << static_cast<std::size_t>(var[i]);
  ofs << std::setw(OWID) << allContactNum;

  // avgPenetr
  for (std::size_t i = 0; i < 6; ++i)
    var[i] = 0;
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getAvgPenetr();
  }
  for (std::size_t i = 0; i < 6; ++i)
    ofs << std::setw(OWID) << var[i];

  // average data
  ofs << std::setw(OWID) << avgNormal << std::setw(OWID) << avgShear
      << std::setw(OWID) << avgPenetr;

  // energy
  ofs << std::setw(OWID) << transEnergy << std::setw(OWID) << rotatEnergy
      << std::setw(OWID) << kinetEnergy << std::setw(OWID) << graviEnergy
      << std::setw(OWID) << mechaEnergy;

  // time
  ofs << std::setw(OWID) << vibraTimeStep << std::setw(OWID) << impactTimeStep
      << std::setw(OWID) << timeStep << std::setw(OWID) << timeAccrued;

  ofs << std::endl;
}

void
Assembly::getStartDimension(REAL& distX, REAL& distY, REAL& distZ)
{
  REAL x1, x2, y1, y2, z1, z2;
  // use boundaryVec
  for (std::vector<Boundary*>::const_iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it) {
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

void
Assembly::openCompressProg(std::ofstream& ofs, const char* str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openCompressProg" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "traction_x1"
      << std::setw(OWID) << "traction_x2" << std::setw(OWID) << "traction_y1"
      << std::setw(OWID) << "traction_y2" << std::setw(OWID) << "traction_z1"
      << std::setw(OWID) << "traction_z2" << std::setw(OWID) << "mean_stress"

      << std::setw(OWID) << "bulk_volume" << std::setw(OWID) << "density"
      << std::setw(OWID) << "epsilon_x" << std::setw(OWID) << "epsilon_y"
      << std::setw(OWID) << "epsilon_z" << std::setw(OWID) << "epsilon_v"
      << std::setw(OWID) << "void_ratio" << std::setw(OWID) << "porosity"

      << std::setw(OWID) << "velocity_x1" << std::setw(OWID) << "velocity_x2"
      << std::setw(OWID) << "velocity_y1" << std::setw(OWID) << "velocity_y2"
      << std::setw(OWID) << "velocity_z1" << std::setw(OWID) << "velocity_z2"

      << std::setw(OWID) << "contact_x1" << std::setw(OWID) << "contact_x2"
      << std::setw(OWID) << "contact_y1" << std::setw(OWID) << "contact_y2"
      << std::setw(OWID) << "contact_z1" << std::setw(OWID) << "contact_z2"
      << std::setw(OWID) << "contact_inside"

      << std::setw(OWID) << "penetr_x1" << std::setw(OWID) << "penetr_x2"
      << std::setw(OWID) << "penetr_y1" << std::setw(OWID) << "penetr_y2"
      << std::setw(OWID) << "penetr_z1" << std::setw(OWID) << "penetr_z2"

      << std::setw(OWID) << "avgNormal" << std::setw(OWID) << "avgShear"
      << std::setw(OWID) << "avgPenetr"

      << std::setw(OWID) << "transEnergy" << std::setw(OWID) << "rotatEnergy"
      << std::setw(OWID) << "kinetEnergy" << std::setw(OWID) << "graviEnergy"
      << std::setw(OWID) << "mechaEnergy"

      << std::setw(OWID) << "vibra_est_dt" << std::setw(OWID) << "impact_est_dt"
      << std::setw(OWID) << "actual_dt"

      << std::endl;
}

void
Assembly::printCompressProg(std::ofstream& ofs, REAL distX, REAL distY,
                            REAL distZ)
{
  REAL x1, x2, y1, y2, z1, z2;
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
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
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    Vec normal = (*it)->getNormalForce();
    Vec veloc = (*it)->getVeloc();
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
  ofs << std::setw(OWID) << avg / 6;

  // volume
  ofs << std::setw(OWID) << bulkVolume << std::setw(OWID)
      << getMass() / bulkVolume << std::setw(OWID) << 1 - (x2 - x1) / distX
      << std::setw(OWID) << 1 - (y2 - y1) / distY << std::setw(OWID)
      << 1 - (z2 - z1) / distZ << std::setw(OWID)
      << 3 - (x2 - x1) / distX - (y2 - y1) / distY - (z2 - z1) / distZ
      << std::setw(OWID) << voidRatio << std::setw(OWID)
      << voidRatio / (1 + voidRatio);

  // velocity
  for (std::size_t i = 0; i < 6; ++i)
    ofs << std::setw(OWID) << vel[i];

  // contactNum
  for (std::size_t i = 0; i < 6; ++i)
    var[i] = 0;
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getContactNum();
  }
  for (std::size_t i = 0; i < 6; ++i)
    ofs << std::setw(OWID) << static_cast<std::size_t>(var[i]);
  ofs << std::setw(OWID) << allContactNum;

  // avgPenetr
  for (std::size_t i = 0; i < 6; ++i)
    var[i] = 0;
  for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getAvgPenetr();
  }
  for (std::size_t i = 0; i < 6; ++i)
    ofs << std::setw(OWID) << var[i];

  // average data
  ofs << std::setw(OWID) << avgNormal << std::setw(OWID) << avgShear
      << std::setw(OWID) << avgPenetr;

  // energy
  ofs << std::setw(OWID) << transEnergy << std::setw(OWID) << rotatEnergy
      << std::setw(OWID) << kinetEnergy << std::setw(OWID) << graviEnergy
      << std::setw(OWID) << mechaEnergy;

  // time
  ofs << std::setw(OWID) << vibraTimeStep << std::setw(OWID) << impactTimeStep
      << std::setw(OWID) << timeStep;

  ofs << std::endl;
}

void
Assembly::readParticle(const char* inputParticle)
{

  REAL young = dem::Parameter::getSingleton().parameter["young"];
  REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

  std::ifstream ifs(inputParticle);
  if (!ifs) {
    debugInf << "stream error: readParticle" << std::endl;
    exit(-1);
  }
  std::size_t particleNum;
  ifs >> particleNum;
  std::string str;
  ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str;

  std::vector<Particle*>::iterator it;
  for (it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    delete (*it);
  allParticleVec.clear();

  std::size_t id, type;
  REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
  for (std::size_t i = 0; i < particleNum; ++i) {
    ifs >> id >> type >> a >> b >> c >> px >> py >> pz >> dax >> day >> daz >>
      dbx >> dby >> dbz >> dcx >> dcy >> dcz >> vx >> vy >> vz >> omx >> omy >>
      omz >> fx >> fy >> fz >> mx >> my >> mz;
    if (px < 30)
      continue;
    Particle* pt =
      new Particle(id, type, Vec(a, b, c), Vec(px, py, pz), Vec(dax, day, daz),
                   Vec(dbx, dby, dbz), Vec(dcx, dcy, dcz), young, poisson);

    // optional settings for a particle's initial status
    if ((static_cast<std::size_t>(
          dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1) {
      pt->setPrevVeloc(Vec(vx, vy, vz));
      pt->setCurrVeloc(Vec(vx, vy, vz));
      pt->setPrevOmega(Vec(omx, omy, omz));
      pt->setCurrOmega(Vec(omx, omy, omz));
      pt->setForce(Vec(fx, fy, fz));  // initial force
      pt->setMoment(Vec(mx, my, mz)); // initial moment
    }

    // pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial force
    // pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial moment

    allParticleVec.push_back(pt);

    /*
          //////////////////////////////////
          // + 7m
          Particle* pt1= new Particle(id+186, type, Vec(a,b,c),
       Vec(px,py,pz+15), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz),
       young, poisson);

          // optional settings for a particle's initial status
          if ( (static_cast<std::size_t>
       (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
            pt1->setPrevVeloc(Vec(vx,vy,vz));
            pt1->setCurrVeloc(Vec(vx,vy,vz));
            pt1->setPrevOmega(Vec(omx,omy,omz));
            pt1->setCurrOmega(Vec(omx,omy,omz));
            pt1->setForce(Vec(fx,fy,fz));  // initial force
            pt1->setMoment(Vec(mx,my,mz)); // initial moment
          }

          //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial
       force
          //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
       moment

          allParticleVec.push_back(pt1);


          //////////////////////////////////
          // + 7*2m
          Particle* pt2= new Particle(id+31*2, type, Vec(a,b,c),
       Vec(px,py,pz+7*2), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz),
       young, poisson);

          // optional settings for a particle's initial status
          if ( (static_cast<std::size_t>
       (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
            pt2->setPrevVeloc(Vec(vx,vy,vz));
            pt2->setCurrVeloc(Vec(vx,vy,vz));
            pt2->setPrevOmega(Vec(omx,omy,omz));
            pt2->setCurrOmega(Vec(omx,omy,omz));
            pt2->setForce(Vec(fx,fy,fz));  // initial force
            pt2->setMoment(Vec(mx,my,mz)); // initial moment
          }

          //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial
       force
          //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
       moment

          allParticleVec.push_back(pt2);


          //////////////////////////////////
          // + 7*3m
          Particle* pt3= new Particle(id+31*3, type, Vec(a,b,c),
       Vec(px,py,pz+7*3), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz),
       young, poisson);

          // optional settings for a particle's initial status
          if ( (static_cast<std::size_t>
       (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
            pt3->setPrevVeloc(Vec(vx,vy,vz));
            pt3->setCurrVeloc(Vec(vx,vy,vz));
            pt3->setPrevOmega(Vec(omx,omy,omz));
            pt3->setCurrOmega(Vec(omx,omy,omz));
            pt3->setForce(Vec(fx,fy,fz));  // initial force
            pt3->setMoment(Vec(mx,my,mz)); // initial moment
          }

          //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial
       force
          //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
       moment

          allParticleVec.push_back(pt3);


          //////////////////////////////////
          // + 7*4m
          Particle* pt4= new Particle(id+31*4, type, Vec(a,b,c),
       Vec(px,py,pz+7*4), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz),
       young, poisson);

          // optional settings for a particle's initial status
          if ( (static_cast<std::size_t>
       (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
            pt4->setPrevVeloc(Vec(vx,vy,vz));
            pt4->setCurrVeloc(Vec(vx,vy,vz));
            pt4->setPrevOmega(Vec(omx,omy,omz));
            pt4->setCurrOmega(Vec(omx,omy,omz));
            pt4->setForce(Vec(fx,fy,fz));  // initial force
            pt4->setMoment(Vec(mx,my,mz)); // initial moment
          }

          //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial
       force
          //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
       moment

          allParticleVec.push_back(pt4);


          //////////////////////////////////
          // + 7*5m
          Particle* pt5= new Particle(id+31*5, type, Vec(a,b,c),
       Vec(px,py,pz+7*5), Vec(dax,day,daz), Vec(dbx,dby,dbz), Vec(dcx,dcy,dcz),
       young, poisson);

          // optional settings for a particle's initial status
          if ( (static_cast<std::size_t>
       (dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1 ) {
            pt5->setPrevVeloc(Vec(vx,vy,vz));
            pt5->setCurrVeloc(Vec(vx,vy,vz));
            pt5->setPrevOmega(Vec(omx,omy,omz));
            pt5->setCurrOmega(Vec(omx,omy,omz));
            pt5->setForce(Vec(fx,fy,fz));  // initial force
            pt5->setMoment(Vec(mx,my,mz)); // initial moment
          }

          //pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial
       force
          //pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
       moment

          allParticleVec.push_back(pt5);
    */
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

void
Assembly::readParticleMiddleLayers(const char* inputParticle)
{

  REAL young = dem::Parameter::getSingleton().parameter["young"];
  REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];

  std::ifstream ifs(inputParticle);
  if (!ifs) {
    debugInf << "stream error: readParticle" << std::endl;
    exit(-1);
  }
  std::size_t particleNum;
  ifs >> particleNum;
  std::string str;
  ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str;

  std::vector<Particle*>::iterator it;
  for (it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    delete (*it);
  allParticleVec.clear();

  std::size_t id, type;
  REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
  int ptcl_num = 0;
  REAL DEMZmin = dem::Parameter::getSingleton().parameter["DEMZmin"];
  REAL DEMZmax = dem::Parameter::getSingleton().parameter["DEMZmax"];
  for (std::size_t i = 0; i < particleNum; ++i) {
    ifs >> id >> type >> a >> b >> c >> px >> py >> pz >> dax >> day >> daz >>
      dbx >> dby >> dbz >> dcx >> dcy >> dcz >> vx >> vy >> vz >> omx >> omy >>
      omz >> fx >> fy >> fz >> mx >> my >> mz;
    if (pz >= DEMZmin && pz <= DEMZmax) {
      ptcl_num++;
      Particle* pt =
        new Particle(id, 1, Vec(a, b, c), Vec(px, py, pz), Vec(dax, day, daz),
                     Vec(dbx, dby, dbz), Vec(dcx, dcy, dcz), young, poisson);

      // optional settings for a particle's initial status
      if ((static_cast<std::size_t>(
            dem::Parameter::getSingleton().parameter["toInitParticle"])) == 1) {
        pt->setPrevVeloc(Vec(vx, vy, vz));
        pt->setCurrVeloc(Vec(vx, vy, vz));
        pt->setPrevOmega(Vec(omx, omy, omz));
        pt->setCurrOmega(Vec(omx, omy, omz));
        pt->setForce(Vec(fx, fy, fz));  // initial force
        pt->setMoment(Vec(mx, my, mz)); // initial moment
      }

      // pt->setConstForce(Vec(fx,fy,fz));  // constant force, not initial force
      // pt->setConstMoment(Vec(mx,my,mz)); // constant moment, not initial
      // moment

      allParticleVec.push_back(pt);
    }
  }
  particleNum = ptcl_num;

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

void
Assembly::printParticle(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: printParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << allParticleVec.size() << std::endl;
  ofs << std::setw(OWID) << "id" << std::setw(OWID) << "type" << std::setw(OWID)
      << "radius_a" << std::setw(OWID) << "radius_b" << std::setw(OWID)
      << "radius_c" << std::setw(OWID) << "position_x" << std::setw(OWID)
      << "position_y" << std::setw(OWID) << "position_z" << std::setw(OWID)
      << "axle_a_x" << std::setw(OWID) << "axle_a_y" << std::setw(OWID)
      << "axle_a_z" << std::setw(OWID) << "axle_b_x" << std::setw(OWID)
      << "axle_b_y" << std::setw(OWID) << "axle_b_z" << std::setw(OWID)
      << "axle_c_x" << std::setw(OWID) << "axle_c_y" << std::setw(OWID)
      << "axle_c_z" << std::setw(OWID) << "velocity_x" << std::setw(OWID)
      << "velocity_y" << std::setw(OWID) << "velocity_z" << std::setw(OWID)
      << "omga_x" << std::setw(OWID) << "omga_y" << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x" << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z" << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y" << std::setw(OWID) << "moment_z"
      << std::endl;

  Vec vObj;
  std::vector<Particle*>::const_iterator it;
  for (it = allParticleVec.begin(); it != allParticleVec.end(); ++it) {
    ofs << std::setw(OWID) << (*it)->getId() << std::setw(OWID)
        << (*it)->getType() << std::setw(OWID) << (*it)->getA()
        << std::setw(OWID) << (*it)->getB() << std::setw(OWID) << (*it)->getC();

    vObj = (*it)->getCurrPos();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrVeloc();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrOmga();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getForce();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getMoment();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ() << std::endl;
  }

  std::size_t sieveNum = gradation.getSieveNum();
  std::vector<REAL> percent = gradation.getPercent();
  std::vector<REAL> size = gradation.getSize();
  ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
  for (std::size_t i = 0; i < sieveNum; ++i)
    ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i]
        << std::endl;
  ofs << std::endl
      << std::setw(OWID) << gradation.getPtclRatioBA() << std::setw(OWID)
      << gradation.getPtclRatioCA() << std::endl;

  ofs.close();
}

void
Assembly::printParticle(const char* str,
                        std::vector<Particle*>& particleVec) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: printParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << particleVec.size() << std::endl;
  ofs << std::setw(OWID) << "id" << std::setw(OWID) << "type" << std::setw(OWID)
      << "radius_a" << std::setw(OWID) << "radius_b" << std::setw(OWID)
      << "radius_c" << std::setw(OWID) << "position_x" << std::setw(OWID)
      << "position_y" << std::setw(OWID) << "position_z" << std::setw(OWID)
      << "axle_a_x" << std::setw(OWID) << "axle_a_y" << std::setw(OWID)
      << "axle_a_z" << std::setw(OWID) << "axle_b_x" << std::setw(OWID)
      << "axle_b_y" << std::setw(OWID) << "axle_b_z" << std::setw(OWID)
      << "axle_c_x" << std::setw(OWID) << "axle_c_y" << std::setw(OWID)
      << "axle_c_z" << std::setw(OWID) << "velocity_x" << std::setw(OWID)
      << "velocity_y" << std::setw(OWID) << "velocity_z" << std::setw(OWID)
      << "omga_x" << std::setw(OWID) << "omga_y" << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x" << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z" << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y" << std::setw(OWID) << "moment_z"
      << std::endl;

  Vec vObj;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    ofs << std::setw(OWID) << (*it)->getId() << std::setw(OWID)
        << (*it)->getType() << std::setw(OWID) << (*it)->getA()
        << std::setw(OWID) << (*it)->getB() << std::setw(OWID) << (*it)->getC();

    vObj = (*it)->getCurrPos();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrVeloc();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getCurrOmga();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getForce();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = (*it)->getMoment();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ() << std::endl;
  }

  ofs.close();
}

void
Assembly::readBoundary(const char* str)
{
  std::ifstream ifs(str);
  if (!ifs) {
    debugInf << "stream error: readBoundary" << std::endl;
    exit(-1);
  }

  REAL x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  setContainer(Rectangle(x1, y1, z1, x2, y2, z2));
  // compute grid assumed to be the same as container, change in
  // scatterParticle() if necessary.
  setGrid(Rectangle(x1, y1, z1, x2, y2, z2));

  boundaryVec.clear();
  Boundary* bptr;
  std::size_t boundaryNum;
  std::size_t type;
  ifs >> boundaryNum;
  for (std::size_t i = 0; i < boundaryNum; ++i) {
    ifs >> type;
    if (type == 1) // plane boundary
      bptr = new planeBoundary(type, ifs);
    else if (type == 2) // cylindrical boundary
      bptr = new cylinderBoundary(type, ifs);

    boundaryVec.push_back(bptr);
  }

  ifs.close();
}

void
Assembly::printBoundary(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: printBoundary" << std::endl;
    exit(-1);
  }
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
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl
      << std::endl
      << std::setw(OWID) << boundaryVec.size() << std::endl;

  for (std::vector<Boundary*>::const_iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it)
    (*it)->print(ofs);

  ofs.close();
}

void
Assembly::plotBoundary(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: plotBoundary" << std::endl;
    exit(-1);
  }
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
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2
      << std::endl;
  ofs << "1 2 3 4 5 6 7 8" << std::endl;

  ofs.close();
}

void
Assembly::plotGrid(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: plotGrid" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = v2 - v1;

  ofs << "ZONE N=" << (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1)
      << ", E=" << mpiProcX * mpiProcY * mpiProcZ
      << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;

  std::vector<Vec> coords((mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1));
  std::size_t index = 0;
  for (std::size_t i = 0; i < mpiProcX + 1; ++i)
    for (std::size_t j = 0; j < mpiProcY + 1; ++j)
      for (std::size_t k = 0; k < mpiProcZ + 1; ++k)
        coords[index++] = Vec(v1.getX() + vspan.getX() / mpiProcX * i,
                              v1.getY() + vspan.getY() / mpiProcY * j,
                              v1.getZ() + vspan.getZ() / mpiProcZ * k);

  for (std::size_t i = 0; i < (mpiProcX + 1) * (mpiProcY + 1) * (mpiProcZ + 1);
       ++i)
    ofs << std::setw(OWID) << coords[i].getX() << std::setw(OWID)
        << coords[i].getY() << std::setw(OWID) << coords[i].getZ() << std::endl;

  for (int iRank = 0; iRank < mpiSize; ++iRank) {
    int coords[3];
    MPI_Cart_coords(cartComm, iRank, 3, coords);

    int id4 = 1 + coords[0] * (mpiProcZ + 1) * (mpiProcY + 1) +
              coords[1] * (mpiProcZ + 1) + coords[2];
    int id1 = 1 + (coords[0] + 1) * (mpiProcZ + 1) * (mpiProcY + 1) +
              coords[1] * (mpiProcZ + 1) + coords[2];
    int id3 = 1 + coords[0] * (mpiProcZ + 1) * (mpiProcY + 1) +
              (coords[1] + 1) * (mpiProcZ + 1) + coords[2];
    int id2 = 1 + (coords[0] + 1) * (mpiProcZ + 1) * (mpiProcY + 1) +
              (coords[1] + 1) * (mpiProcZ + 1) + coords[2];

    int id8 = 1 + coords[0] * (mpiProcZ + 1) * (mpiProcY + 1) +
              coords[1] * (mpiProcZ + 1) + (coords[2] + 1);
    int id5 = 1 + (coords[0] + 1) * (mpiProcZ + 1) * (mpiProcY + 1) +
              coords[1] * (mpiProcZ + 1) + (coords[2] + 1);
    int id7 = 1 + coords[0] * (mpiProcZ + 1) * (mpiProcY + 1) +
              (coords[1] + 1) * (mpiProcZ + 1) + (coords[2] + 1);
    int id6 = 1 + (coords[0] + 1) * (mpiProcZ + 1) * (mpiProcY + 1) +
              (coords[1] + 1) * (mpiProcZ + 1) + (coords[2] + 1);

    ofs << std::setw(8) << id1 << std::setw(8) << id2 << std::setw(8) << id3
        << std::setw(8) << id4 << std::setw(8) << id5 << std::setw(8) << id6
        << std::setw(8) << id7 << std::setw(8) << id8 << std::endl;
  }

  ofs.close();
}

void
Assembly::findContact()
{ // various implementations
  int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];

  if (ompThreads == 1) { // non-openmp single-thread version, time complexity
                         // bigO(n x n), n is the number of particles.
    contactVec.clear();

#ifdef TIME_PROFILE
    REAL time_r =
      0; // time consumed in contact resolution, i.e., tmpContact.isOverlapped()
    gettimeofday(&time_p1, NULL);
#endif

    std::size_t num1 = particleVec.size();      // particles inside container
    std::size_t num2 = mergeParticleVec.size(); // paticles inside container (at
                                                // front) + particles from
                                                // neighboring blocks (at end)
    for (std::size_t i = 0; i < num1;
         ++i) { // NOT (num1 - 1), in parallel situation where one particle
                // could contact received particles!
      Vec u = particleVec[i]->getCurrPos();
      for (std::size_t j = i + 1; j < num2; ++j) {
        Vec v = mergeParticleVec[j]->getCurrPos();
        if ((vfabs(v - u) <
             particleVec[i]->getA() + mergeParticleVec[j]->getA()) &&
            (particleVec[i]->getType() != 1 ||
             mergeParticleVec[j]->getType() !=
               1) // not both are fixed particles
            && (particleVec[i]->getType() != 5 ||
                mergeParticleVec[j]->getType() !=
                  5) // not both are free boundary particles
            && (particleVec[i]->getType() != 10 ||
                mergeParticleVec[j]->getType() !=
                  10)) { // not both are ghost particles
          Contact tmpContact(
            particleVec[i],
            mergeParticleVec[j]); // a local and temparory object
#ifdef TIME_PROFILE
          gettimeofday(&time_r1, NULL);
#endif
          if (tmpContact.isOverlapped())
            contactVec.push_back(tmpContact); // containers use value semantics,
                                              // so a "copy" is pushed back.
#ifdef TIME_PROFILE
          gettimeofday(&time_r2, NULL);
          time_r += timediffsec(time_r1, time_r2);
#endif
        }
      }
    }

#ifdef TIME_PROFILE
    gettimeofday(&time_p2, NULL);
    debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID)
             << timediffsec(time_p1, time_p2) << std::setw(OWID)
             << "isOverlapped=" << std::setw(OWID) << time_r;
#endif

  }

  else if (ompThreads > 1) { // openmp implementation: various loop scheduling -
                             // (static), (static,1), (dynamic), (dynamic,1)
    contactVec.clear();

#ifdef TIME_PROFILE
    gettimeofday(&time_p1, NULL);
#endif

    std::size_t i, j;
    Vec u, v;
    std::size_t num1 = particleVec.size();
    std::size_t num2 = mergeParticleVec.size();
    int ompThreads = dem::Parameter::getSingleton().parameter["ompThreads"];

#pragma omp parallel for num_threads(ompThreads) private(i, j, u, v) shared(   \
  num1, num2) schedule(dynamic)
    for (i = 0; i < num1; ++i) {
      u = particleVec[i]->getCurrPos();
      for (j = i + 1; j < num2; ++j) {
        v = mergeParticleVec[j]->getCurrPos();
        if ((vfabs(v - u) <
             particleVec[i]->getA() + mergeParticleVec[j]->getA()) &&
            (particleVec[i]->getType() != 1 ||
             mergeParticleVec[j]->getType() !=
               1) // not both are fixed particles
            && (particleVec[i]->getType() != 5 ||
                mergeParticleVec[j]->getType() !=
                  5) // not both are free boundary particles
            && (particleVec[i]->getType() != 10 ||
                mergeParticleVec[j]->getType() !=
                  10)) { // not both are ghost particles
          Contact tmpContact(
            particleVec[i],
            mergeParticleVec[j]); // a local and temparory object
          if (tmpContact.isOverlapped())
#pragma omp critical
            contactVec.push_back(tmpContact); // containers use value semantics,
                                              // so a "copy" is pushed back.
        }
      }
    }

#ifdef TIME_PROFILE
    gettimeofday(&time_p2, NULL);
    debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID)
             << timediffsec(time_p1, time_p2);
#endif

  } // end of openmp implementation
}

void
Assembly::internalForce()
{
  REAL pAvg[3], sum[3];
  for (std::size_t i = 0; i < 3; ++i) {
    pAvg[i] = 0;
    sum[i] = 0;
  }

  if (contactVec.size() > 0) {
    for (std::vector<Contact>::iterator it = contactVec.begin();
         it != contactVec.end(); ++it)
      it->checkinPrevTgt(
        contactTgtVec); // checkin previous tangential force and displacment

#ifdef TIME_PROFILE
    gettimeofday(&time_p1, NULL);
#endif

    contactTgtVec
      .clear(); // contactTgtVec must be cleared before filling in new values.
    for (std::vector<Contact>::iterator it = contactVec.begin();
         it != contactVec.end(); ++it) {
      it->contactForce(); // cannot be parallelized as it may change a
                          // particle's force simultaneously.
      it->checkoutTgt(
        contactTgtVec); // checkout current tangential force and displacment
      pAvg[0] += it->getNormalForce();
      pAvg[1] += it->getTgtForce();
      pAvg[2] += it->getPenetration();
    }
    for (std::size_t i = 0; i < 3; ++i)
      pAvg[i] /= contactVec.size();

#ifdef TIME_PROFILE
    gettimeofday(&time_p2, NULL);
    debugInf << std::setw(OWID) << "internalForce=" << std::setw(OWID)
             << timediffsec(time_p1, time_p2) << std::endl;
#endif
  }

  MPI_Reduce(pAvg, sum, 3, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
  avgNormal = sum[0] / mpiSize;
  avgShear = sum[1] / mpiSize;
  avgPenetr = sum[2] / mpiSize;
}

void
Assembly::dragForce()
{
  for (std::vector<Particle*>::iterator it = particleVec.begin();
       it != particleVec.end(); ++it)
    (*it)->dragForce();
}

void
Assembly::updateParticle()
{
  for (std::vector<Particle*>::iterator it = particleVec.begin();
       it != particleVec.end(); ++it)
    (*it)->update();
}

void
Assembly::updateBoundary(REAL sigma, std::string type, REAL sigmaX, REAL sigmaY)
{
  if (mpiRank == 0) {
    REAL x1, x2, y1, y2, z1, z2;
    for (std::vector<Boundary*>::const_iterator it = mergeBoundaryVec.begin();
         it != mergeBoundaryVec.end(); ++it) {
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
      for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
           it != mergeBoundaryVec.end(); ++it)
        (*it)->updateIsotropic(sigma, areaX, areaY, areaZ);
    } else if (type.compare("odometer") == 0) {
      for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
           it != mergeBoundaryVec.end(); ++it)
        (*it)->updateOdometer(sigma, areaX, areaY, areaZ);
    } else if (type.compare("triaxial") == 0) {
      for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
           it != mergeBoundaryVec.end(); ++it)
        (*it)->updateTriaxial(sigma, areaX, areaY, areaZ);
    } else if (type.compare("plnstrn") == 0) {
      for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
           it != mergeBoundaryVec.end(); ++it)
        (*it)->updatePlaneStrain(sigma, areaX, areaY, areaZ);
    } else if (type.compare("trueTriaxial") == 0) {
      for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
           it != mergeBoundaryVec.end(); ++it)
        (*it)->updateTrueTriaxial(sigma, areaX, areaY, areaZ, sigmaX, sigmaY);
    }

    // update boundaryVec from mergeBoundaryVec and remove contactInfo to reduce
    // MPI transmission
    boundaryVec = mergeBoundaryVec;
    for (std::vector<Boundary*>::iterator it = boundaryVec.begin();
         it != boundaryVec.end(); ++it)
      (*it)->clearContactInfo();

    // update allContainer
    for (std::vector<Boundary*>::const_iterator it = boundaryVec.begin();
         it != boundaryVec.end(); ++it) {
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

void
Assembly::clearContactForce()
{
  for (std::vector<Particle*>::iterator it = particleVec.begin();
       it != particleVec.end(); ++it)
    (*it)->clearContactForce();
}

void
Assembly::findBdryContact()
{
  for (std::vector<Boundary*>::iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it)
    (*it)->findBdryContact(particleVec);
}

void
Assembly::boundaryForce()
{
  for (std::vector<Boundary*>::iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it)
    (*it)->boundaryForce(boundaryTgtMap);
}

void
Assembly::printContact(char* str) const
{
  // There are two implementions of printContact
  // implementation 1: parallel IO, each process prints to a shared data file
  // using a shared pointer.
  //                   and use post-processing tool to remove redundant info.
  MPI_Status status;
  MPI_File contactFile;
  MPI_File_open(mpiWorld, str, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                &contactFile);
  if (boostWorld.rank() == 0 && !contactFile) {
    debugInf << "stream error: printContact" << std::endl;
    exit(-1);
  }

  std::stringstream inf;
  inf.setf(std::ios::scientific, std::ios::floatfield);

  for (std::vector<Contact>::const_iterator it = contactVec.begin();
       it != contactVec.end(); ++it)
    inf << std::setw(OWID) << it->getP1()->getId() << std::setw(OWID)
        << it->getP2()->getId() << std::setw(OWID) << it->getPoint1().getX()
        << std::setw(OWID) << it->getPoint1().getY() << std::setw(OWID)
        << it->getPoint1().getZ() << std::setw(OWID) << it->getPoint2().getX()
        << std::setw(OWID) << it->getPoint2().getY() << std::setw(OWID)
        << it->getPoint2().getZ() << std::setw(OWID) << it->getRadius1()
        << std::setw(OWID) << it->getRadius2() << std::setw(OWID)
        << it->getPenetration() << std::setw(OWID) << it->getTgtDisp()
        << std::setw(OWID) << it->getContactRadius() << std::setw(OWID)
        << it->getR0() << std::setw(OWID) << it->getE0() << std::setw(OWID)
        << it->getNormalForce() << std::setw(OWID) << it->getTgtForce()
        << std::setw(OWID)
        << (it->getPoint1().getX() + it->getPoint2().getX()) / 2
        << std::setw(OWID)
        << (it->getPoint1().getY() + it->getPoint2().getY()) / 2
        << std::setw(OWID)
        << (it->getPoint1().getZ() + it->getPoint2().getZ()) / 2
        << std::setw(OWID) << it->normalForceVec().getX() << std::setw(OWID)
        << it->normalForceVec().getY() << std::setw(OWID)
        << it->normalForceVec().getZ() << std::setw(OWID)
        << it->tgtForceVec().getX() << std::setw(OWID)
        << it->tgtForceVec().getY() << std::setw(OWID)
        << it->tgtForceVec().getZ() << std::setw(OWID) << it->getVibraTimeStep()
        << std::setw(OWID) << it->getImpactTimeStep() << std::endl;

  int length = (OWID * 28 + 1) * contactVec.size();
  // write a file at a location specified by a shared file pointer (blocking,
  // collective)
  // note MPI_File_write_shared is non-collective
  MPI_File_write_ordered(contactFile, const_cast<char*>(inf.str().c_str()),
                         length, MPI_CHAR, &status);
  MPI_File_close(&contactFile);

  // implementation 2: each process prints to an individual file.
  //                   use post-processing tool to merge files and remove
  //                   redundance.
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

  std::vector<Contact>::const_iterator it;
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
        << std::setw(OWID) << ( it->getPoint1().getX() + it->getPoint2().getX()
  )/2
        << std::setw(OWID) << ( it->getPoint1().getY() + it->getPoint2().getY()
  )/2
        << std::setw(OWID) << ( it->getPoint1().getZ() + it->getPoint2().getZ()
  )/2
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

void
Assembly::calcTransEnergy()
{
  REAL pEngy = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getTransEnergy();
  }
  MPI_Reduce(&pEngy, &transEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcRotatEnergy()
{
  REAL pEngy = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getRotatEnergy();
  }
  MPI_Reduce(&pEngy, &rotatEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcKinetEnergy()
{
  REAL pEngy = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getKinetEnergy();
  }
  MPI_Reduce(&pEngy, &kinetEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcGraviEnergy(REAL ref)
{
  REAL pEngy = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getPotenEnergy(ref);
  }
  MPI_Reduce(&pEngy, &graviEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcMechaEnergy()
{
  mechaEnergy = kinetEnergy + graviEnergy;
}

REAL
Assembly::getMass() const
{
  REAL var = 0;
  for (std::vector<Particle*>::const_iterator it = allParticleVec.begin();
       it != allParticleVec.end(); ++it)
    var += (*it)->getMass();
  return var;
}

REAL
Assembly::getParticleVolume() const
{
  REAL var = 0;
  for (std::vector<Particle*>::const_iterator it = allParticleVec.begin();
       it != allParticleVec.end(); ++it)
    if ((*it)->getType() == 0)
      var += (*it)->getVolume();
  return var;
}

void
Assembly::calcTimeStep()
{
  calcVibraTimeStep();
  calcImpactTimeStep();
  calcContactNum();

  std::valarray<REAL> dt(3);
  dt[0] = dem::Parameter::getSingleton().parameter["timeStep"];
  dt[1] = vibraTimeStep;
  dt[2] = impactTimeStep;

  timeStep = dt.min();
}

void
Assembly::calcContactNum()
{
  std::size_t pContactNum = contactVec.size();
  MPI_Reduce(&pContactNum, &allContactNum, 1, MPI_INT, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcVibraTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (contactVec.size() == 0)
    pTimeStep = 1 / EPS;
  else {
    std::vector<Contact>::const_iterator it = contactVec.begin();
    pTimeStep = it->getVibraTimeStep();
    for (++it; it != contactVec.end(); ++it) {
      REAL val = it->getVibraTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &vibraTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
}

void
Assembly::calcImpactTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (contactVec.size() == 0)
    pTimeStep = 1 / EPS;
  else {
    std::vector<Contact>::const_iterator it = contactVec.begin();
    pTimeStep = it->getImpactTimeStep();
    for (++it; it != contactVec.end(); ++it) {
      REAL val = it->getImpactTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &impactTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
}

REAL
Assembly::getAvgTransVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getCurrVeloc());
      ++count;
    }
  return avgv /= count;
}

REAL
Assembly::getAvgRotatVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getCurrOmga());
      ++count;
    }
  return avgv /= count;
}

REAL
Assembly::getAvgForce() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getForce());
      ++count;
    }
  return avgv / count;
}

REAL
Assembly::getAvgMoment() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  std::vector<Particle*>::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getMoment());
      ++count;
    }
  return avgv /= count;
}

void
Assembly::buildBoundary(std::size_t boundaryNum, const char* boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if (!ofs) {
    debugInf << "stream error: buildBoundary" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);

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

  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl
      << std::endl
      << std::setw(OWID) << boundaryNum << std::endl
      << std::endl;

  if (boundaryNum == 1) { // only a bottom boundary, i.e., boundary 5
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl;

  } else if (boundaryNum == 5) { // no top boundary, i.e., no boundary 6
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 1 << std::setw(OWID) << -1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x1 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 2
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 2 << std::setw(OWID) << 1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x2 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 3
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 3 << std::setw(OWID) << 0 << std::setw(OWID) << -1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y1 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 4
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 4 << std::setw(OWID) << 0 << std::setw(OWID) << 1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y2 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 5
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl;
  } else if (boundaryNum == 6) { // all 6 boundaries
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 1 << std::setw(OWID) << -1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x1 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 2
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 2 << std::setw(OWID) << 1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x2 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 3
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 3 << std::setw(OWID) << 0 << std::setw(OWID) << -1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y1 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 4
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 4 << std::setw(OWID) << 0 << std::setw(OWID) << 1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y2 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 5
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl

        // boundary 6
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 6 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl;
  }

  ofs.close();
}

} // namespace dem ends

/*
// create a specimen from discreate particles through floating and then
gravitation,
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
newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99,
young, poisson);
particleVec.push_back(newptcl);
particleNum++;
}

// particle boundary 2
y=dimn/2*(grid+1)/10;
for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99,
young, poisson);
particleVec.push_back(newptcl);
particleNum++;
}

// particle boundary 3
x=-dimn/2*(grid+1)/10;
for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99,
young, poisson);
particleVec.push_back(newptcl);
particleNum++;
}

// particle boundary 4
y=-dimn/2*(grid+1)/10;
for (x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5)
for(z=-dimn/2; z<-dimn/2 + dimn*wall; z+=dimn/2/5) {
newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99,
young, poisson);
particleVec.push_back(newptcl);
particleNum++;
}

// particle boundary 6
z=-dimn/2;
for (y=-dimn/2*grid/10; y<dimn/2*grid/10*est; y+=dimn/2/5)
for( x=-dimn/2*grid/10; x<dimn/2*grid/10*est; x+=dimn/2/5) {
newptcl = new Particle(particleNum+1, 1, Vec(x,y,z), grad.ptclsize[0]*0.99,
young, poisson);
particleVec.push_back(newptcl);
particleNum++;
}

if (particleLayers == 0) {      // just one free particle
newptcl = new Particle(particleNum+1, 0, Vec(dimn/2/40,dimn/2/20,dimn/2), grad,
young, poisson);
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
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1
  << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1
  << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1
  << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
  << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2
  << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
  << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2
  << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2
  << std::endl;
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
  ofs << "ZONE N=" << totalMemParticle << ", E=" << totalSpring << ",
  DATAPACKING=POINT, ZONETYPE=FELINESEG" << std::endl;
  Particle *pt = NULL;
  Vec vt;
  for (std::size_t i = 0; i < memBoundary.size(); ++i)
  for (std::size_t j = 0; j < memBoundary[i].size(); ++j)
  for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k) {
  pt = memBoundary[i][j][k];
  vt = pt->getCurrPos();
  ofs << std::setw(OWID) << vt.getX() << std::setw(OWID) << vt.getY() <<
  std::setw(OWID) << vt.getZ() << std::endl;
  }
  for (std::size_t i = 0; i < springVec.size(); ++i) {
  ofs << std::setw(OWID) << springVec[i]->getParticleId1() - trimHistoryNum  <<
  std::setw(OWID) << springVec[i]->getParticleId2() - trimHistoryNum <<
  std::endl;
  }

  ofs.close();
  }

  void Assembly::printMemParticle(const char *str) const  {
  std::ofstream ofs(str);
  if(!ofs) { debugInf << "stream error: printMemParticle" << std::endl;
  exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  std::size_t totalMemParticle = 0;
  for (std::size_t i = 0; i < memBoundary.size(); ++i)
  for (std::size_t j = 0; j < memBoundary[i].size(); ++j)
  for (std::size_t k = 0; k < memBoundary[i][j].size(); ++k)
  ++totalMemParticle;

  ofs << std::setw(OWID) << totalMemParticle << std::setw(OWID) << 1 <<
  std::endl;
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
  std::vector<Particle*> vec1d;  // 1-dimension
  std::vector< std::vector<Particle*>  > vec2d; // 2-dimension
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
  // 3: implementation 3, no partition, each thread leaps by ts until num/2 and
  handles two particles.
  // 4: implementation 4, no partition, parallel for, various loop scheduling:
  (static), (static,1), (dynamic), (dynamic,1)

  //start of def OPENMP
  #ifdef OPENMP

  #if OPENMP_IMPL == 0
  // OpenMP implementation 0: ts partitions, each thread handles a partition,
  max diff = n*n*(1-1/ts)/ts
  // implementation is based on linked list, also works for vector but not
  efficient.
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
  std::vector<Particle*>::iterator ot, it, pt;

  #pragma omp parallel num_threads(nThreads) private(tid, ts, tnum, it, pt, i,
  j, u, v) shared(num) reduction(+: possContact)
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
  && ( (*it)->getType() !=  1 || (*pt)->getType() != 1  )      // not both are
  fixed particles
  && ( (*it)->getType() !=  5 || (*pt)->getType() != 5  )      // not both are
  free boundary particles
  && ( (*it)->getType() != 10 || (*pt)->getType() != 10 )  ) { // not both are
  ghost particles
  contact<Particle> tmpContact(*it, *pt); // a local and temparory object
  ++possContact;
  if(tmpContact.isOverlapped())
  #pragma omp critical
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  }
  }
  }
  }

  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2);
  #endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
  } // end of OpenMP implementation 0

  #elif OPENMP_IMPL == 1
  // OpenMP implementation 1: ts partitions, each thread handles a partition,
  max diff = n*n*(1-1/ts)/ts
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

  #pragma omp parallel num_threads(nThreads) private(tid, ts, start, end, i, j,
  u, v) shared(num) reduction(+: possContact)
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
  && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )
  // not both are fixed particles
  && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )
  // not both are free boundary particles
  && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) {
  // not both are ghost particles
  contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and
  temparory object
  ++possContact;
  if(tmpContact.isOverlapped())
  #pragma omp critical
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  }
  }
  }
  }

  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf <<  std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2);
  #endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
  } // end of OpenMP implementation 1

  #elif OPENMP_IMPL == 2
  // OpenMP implementation 2: no partitions, each thread leaps by ts until
  completed, max diff = n*(ts-1)/ts
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

  #pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, u, v)
  shared(num) reduction(+: possContact)
  {
  tid = omp_get_thread_num();
  ts  = omp_get_num_threads();

  // explore each partition
  for (i = tid; i < num; i += ts) {
  u = particleVec[i]->getCurrPos();
  for (j = i + 1; j < num; ++j) {
  v = particleVec[j]->getCurrPos();
  if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
  && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )
  // not both are fixed particles
  && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )
  // not both are free boundary particles
  && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) {
  // not both are ghost particles
  contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and
  temparory object
  ++possContact;
  if(tmpContact.isOverlapped())
  #pragma omp critical
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  }
  }
  }
  }

  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2);
  #endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
  } // end of OpenMP implementation 2

  #elif OPENMP_IMPL == 3
  // OpenMP implementation 3: no partitions, each thread leaps by ts until num/2
  and handles two particles, max diff = 0
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

  #pragma omp parallel num_threads(nThreads) private(tid, ts, i, j, k, u, v)
  shared(num) reduction(+: possContact)
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
  && ( particleVec[k]->getType() !=  1 || particleVec[j]->getType() != 1  )
  // not both are fixed particles
  && ( particleVec[k]->getType() !=  5 || particleVec[j]->getType() != 5  )
  // not both are free boundary particles
  && ( particleVec[k]->getType() != 10 || particleVec[j]->getType() != 10 )  ) {
  // not both are ghost particles
  contact<Particle> tmpContact(particleVec[k], particleVec[j]); // a local and
  temparory object
  ++possContact;
  if(tmpContact.isOverlapped())
  #pragma omp critical
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  }
  }
  }
  }
  }

  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2);
  #endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
  } // end of OpenMP implementation 3

  #elif OPENMP_IMPL == 4
  // OpenMP implementation 4: no partitions, parallel for, various loop
  scheduling: (static), (static,1), (dynamic), (dynamic,1)
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

  #pragma omp parallel for num_threads(nThreads) private(i, j, u, v) shared(num)
  reduction(+: possContact) schedule(dynamic)
  for (i = 0; i < num - 1; ++i) {
  u = particleVec[i]->getCurrPos();
  for (j = i + 1; j < num; ++j) {
  v = particleVec[j]->getCurrPos();
  if (   ( vfabs(v-u) < particleVec[i]->getA() + particleVec[j]->getA() )
  && ( particleVec[i]->getType() !=  1 || particleVec[j]->getType() != 1  )
  // not both are fixed particles
  && ( particleVec[i]->getType() !=  5 || particleVec[j]->getType() != 5  )
  // not both are free boundary particles
  && ( particleVec[i]->getType() != 10 || particleVec[j]->getType() != 10 )  ) {
  // not both are ghost particles
  contact<Particle> tmpContact(particleVec[i], particleVec[j]); // a local and
  temparory object
  ++possContact;
  if(tmpContact.isOverlapped())
  #pragma omp critical
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  }
  }
  }


  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2);
  #endif
  possContactNum   = possContact;
  actualContactNum = contactVec.size();
  } // end of OpenMP implementation 4

  #endif

  #else //else of def OPENMP, i.e., serial versions start here:

  //start of ndef BINNING
  #ifndef BINNING
  void Assembly::findContact() { // serial version, O(n x n), n is the number of
  particles.
  contactVec.clear();
  possContactNum = 0;

  #ifdef TIME_PROFILE
  REAL time_r = 0; // time consumed in contact resolution, i.e.,
  tmpContact.isOverlapped()
  gettimeofday(&time_p1,NULL);
  #endif

  int num1 = particleVec.size();  // particles inside container
  int num2 = mergeParticleVec.size(); // particles inside container (at front) +
  particles from neighboring blocks (at end)
  for (int i = 0; i < num1 - 1; ++i) {
  Vec u = particleVec[i]->getCurrPos();
  for (int j = i + 1; j < num2; ++j) {
  Vec v = mergeParticleVec[j]->getCurrPos();
  if (   ( vfabs(v - u) < particleVec[i]->getA() + mergeParticleVec[j]->getA())
  && ( particleVec[i]->getType() !=  1 || mergeParticleVec[j]->getType() != 1  )
  // not both are fixed particles
  && ( particleVec[i]->getType() !=  5 || mergeParticleVec[j]->getType() != 5  )
  // not both are free boundary particles
  && ( particleVec[i]->getType() != 10 || mergeParticleVec[j]->getType() != 10 )
  ) { // not both are ghost particles
  Contact tmpContact(particleVec[i], mergeParticleVec[j]); // a local and
  temparory object
  ++possContactNum;
  #ifdef TIME_PROFILE
  gettimeofday(&time_r1,NULL);
  #endif
  if(tmpContact.isOverlapped())
  contactVec.push_back(tmpContact);    // containers use value semantics, so a
  "copy" is pushed back.
  #ifdef TIME_PROFILE
  gettimeofday(&time_r2,NULL);
  time_r += timediffsec(time_r1, time_r2);
  #endif
  }
  }
  }

  #ifdef TIME_PROFILE
  gettimeofday(&time_p2,NULL);
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" <<
  std::setw(OWID) << time_r;
  #endif

  actualContactNum = contactVec.size();
  }

  //else of ndef BINNING
  #else
  void Assembly::findContact() { // serial version, binning methods, cell
  slightly larger than maximum particle
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
  typedef std::pair<bool, std::vector<Particle*> > cellT;
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
  //debugInf <<  i << " " << j << " " << k << " " << "m n size=" << m << " " <<
  n << " " <<  cellVec[i][j][k].size() << std::endl;
  pt = cellVec[i][j][k].second[n];
  v  = pt->getCurrPos();
  if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
  ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed
  particles
  ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free
  boundary particles
  ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost
  particles
  contact<Particle> tmpContact(it, pt); // a local and temparory object
  ++possContactNum;
  #ifdef TIME_PROFILE
  gettimeofday(&time_r1,NULL);
  #endif
  if(tmpContact.isOverlapped())
  contactVec.push_back(tmpContact);   // containers use value semantics, so a
  "copy" is pushed back.
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
  if (ci > -1 && ci < nx && cj > -1 && cj < ny && ck > -1 && ck < nz &&
  cellVec[ci][cj][ck].first == false ) {
  //debugInf << "i j k m ncell ci cj ck size contacts= " << i << " " << j << " "
  << k << " " << m  << " " << ncell << " " << ci << " " << cj << " " << ck << "
  " << cellVec[ci][cj][ck].second.size() << " "  << contactVec.size() <<
  std::endl;
  std::vector<Particle*> vt = cellVec[ci][cj][ck].second;
  for (int n = 0; n < vt.size(); ++n) {
  pt = vt[n];
  v  = pt->getCurrPos();
  if ( ( vfabs(u-v) < it->getA() + pt->getA() )  &&
  ( it->getType() !=  1 || pt->getType() != 1 ) &&   // not both are fixed
  particles
  ( it->getType() !=  5 || pt->getType() != 5 ) &&   // not both are free
  boundary particles
  ( it->getType() != 10 || pt->getType() != 10)  ) { // not both are ghost
  particles
  contact<Particle> tmpContact(it, pt); // a local and temparory object
  ++possContactNum;
  #ifdef TIME_PROFILE
  gettimeofday(&time_r1,NULL);
  #endif
  if(tmpContact.isOverlapped())
  contactVec.push_back(tmpContact);   // containers use value semantics, so a
  "copy" is pushed back.
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
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID) <<
  timediffsec(time_p1, time_p2) << std::setw(OWID) << "isOverlapped=" <<
  std::setw(OWID) << time_r;
  #endif

  actualContactNum = contactVec.size();
  }

  //end of ndef BINNING
  #endif

  //end of def OPENMP
  #endif


  Vec Assembly::getTopFreeParticlePosition() const {
  std::vector<Particle*>::const_iterator it,jt,kt;
  it=particleVec.begin();
  while (it!=particleVec.end() && (*it)->getType()!=0)   // find the 1st free
  particle
  ++it;

  if (it==particleVec.end())    // no free particles
  return 0;

  jt=it;
  kt=it;

  // two cases:
  // 1: 1st particle is not free
  // 2: 1st particle is free
  if (++kt!=particleVec.end()) { // case1: more than 2 particles; case 2: more
  than 1 particle
  for(++it;it!=particleVec.end();++it) {
  if ((*it)->getType()==0)
  if ((*it)->getCurrPos().getZ() > (*jt)->getCurrPos().getZ())
  jt=it;
  }
  return (*jt)->getCurrPos();
  }
  else {
  if ((*it)->getType()==0)  // case1: only 2 particles, the 2nd one is free;
  case2: only 1 particle
  return (*it)->getCurrPos();
  else
  return 0;
  }

  }



  REAL Assembly::ellipPileForce() {
  REAL val=0;
  for(std::vector<Particle*>::iterator
  it=particleVec.begin();it!=particleVec.end();++it)
  if ((*it)->getType()==3) {
  val = (*it)->getForce().getZ();
  break;
  }
  return val;
  }

  Vec Assembly::ellipPileDimn() {
  Vec val;
  for(std::vector<Particle*>::iterator
  it=particleVec.begin();it!=particleVec.end();++it)
  if ((*it)->getType()==3) {
  val = Vec((*it)->getA(), (*it)->getB(), (*it)->getC());
  break;
  }
  return val;
  }

  REAL Assembly::ellipPileTipZ() {
  REAL val=0;
  for(std::vector<Particle*>::iterator
  it=particleVec.begin();it!=particleVec.end();++it)
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
  REAL low=ellipPileTipZ() + ellipPileDimn().getX() -
  getTopFreeParticlePosition().getZ();
  REAL lowint=low-pow(low,3)/3.0/pow(ellipPileDimn().getX(),2);
  val = Pi * ellipPileDimn().getY() * ellipPileDimn().getZ()
  *(2.0/3*ellipPileDimn().getX()-lowint);
  }
  return val;
  }

  void Assembly::ellipPileUpdate() {
  for(std::vector<Particle*>::iterator
  it=particleVec.begin();it!=particleVec.end();++it) {
  if ((*it)->getType()==3) {
  (*it)->setCurrVeloc(Vec(0, 0, -pileRate));
  (*it)->setCurrPos( (*it)->getPrevPos() + (*it)->getCurrVeloc() * timeStep);
  }
  }
  }





  void Assembly::springForce() {
  for (vector<Spring*>::iterator it = springVec.begin(); it != springVec.end();
  ++it)
  (*it)->applyForce();
  }

  void Assembly::readCavityBoundary(const char *str) {
  std::ifstream ifs(str);
  if(!ifs) { debugInf << "stream error: readCavityBoundary" << std::endl;
  exit(-1); }

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
  if(!ofs) { debugInf << "stream error: printCavityBoundary" << std::endl;
  exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  ofs << std::setw(OWID) << cavityBoundaryVec.size() << std::endl;
  std::vector<Boundary*>::const_iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
  (*rt)->display(ofs);
  ofs << std::endl;

  ofs.close();
  }



  void Assembly::findCavityContact() {
  std::vector<Boundary*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
  (*rt)->findBdryContact(allParticleVec);
  }


  void Assembly::cavityBoundaryForce() {
  std::vector<Boundary*>::iterator rt;
  for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
  (*rt)->boundaryForce(boundaryTgtMap);
  }

  Vec Assembly::getNormalForce(int bdry) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==bdry)
  return (*it)->getNormalForce();
  }
  return 0;
  }

  Vec Assembly::getShearForce(int bdry) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==bdry)
  return (*it)->getShearForce();
  }
  return 0;
  }

  REAL Assembly::getAvgNormal(int bdry) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==bdry)
  return (*it)->getAvgNormal();
  }
  return 0;
  }

  Vec Assembly::getApt(int bdry) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==bdry)
  return (*it)->getApt();
  }
  return 0;
  }


  Vec Assembly::getDirc(int bdry) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==bdry)
  return (*it)->getDirc();
  }
  return 0;
  }

  REAL Assembly::getArea(int n) const {
  std::vector<Boundary*>::const_iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==n)
  return (*it)->area;
  }
  return 0;
  }

  void Assembly::setArea(int n, REAL a) {
  std::vector<Boundary*>::iterator it;
  for(it=boundaryVec.begin();it!=boundaryVec.end();++it) {
  if((*it)->getBdryID()==n)
  (*it)->area=a;
  }
  }

  REAL Assembly::getAvgPressure() const {
  std::vector<Boundary*>::const_iterator rt;
  REAL avgpres=0;
  for(rt=boundaryVec.begin();rt!=boundaryVec.end();++rt)
  avgpres+=vfabs((*rt)->getNormalForce())/(*rt)->getArea();
  return avgpres/=boundaryVec.size();
  }

  // only update CoefOfLimits[0] for specified boundaries
  void Assembly::updateBoundary(int bn[], UPDATECTL rbctl[], int num) {
  for(int i=0;i<num;i++) {
  for(std::vector<Boundary*>::iterator
  rt=boundaryVec.begin();rt!=boundaryVec.end();++rt) {
  if((*rt)->getBdryID()==bn[i]) {
  (*rt)->update(rbctl[i]);
  break;
  }
  }
  }
  }

  // update CoefOfLimits[1,2,3,4] for all 6 boundaries
  void Assembly::updateBoundary6() {
  for(std::vector<Boundary*>::iterator
  rt=boundaryVec.begin();rt!=boundaryVec.end();++rt) {
  if((*rt)->getBdryID()==1 || (*rt)->getBdryID()==3) {
  for(std::vector<Boundary*>::iterator
  lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
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
  for(std::vector<Boundary*>::iterator
  lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
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
  for(std::vector<Boundary*>::iterator
  lt=boundaryVec.begin();lt!=boundaryVec.end();++lt) {
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
  if(!progressInf) { debugInf << "stream error: angleOfRepose" << std::endl;
  exit(-1); }
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
  if(!debugInf) { debugInf << "stream error: angleOfRepose" << std::endl;
  exit(-1); }
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
  std::vector<Particle*> lastPtcls;
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
  newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i),
  gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr +
  maxDiameter * i), gradation, young, poisson);
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

  lastPtcls.clear(); // do not delete those pointers to release memory;
  particleVec will do it.
  zCurr = getPtclMaxZ(allParticleVec) + maxDiameter;

  for ( int i = 0; i != layers; ++i) {
  newPtcl = new Particle(particleNum+1, 0, Vec(0, 0, zCurr + maxDiameter * i),
  gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(maxDiameter, 0, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(-maxDiameter, 0, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(0, maxDiameter, zCurr +
  maxDiameter * i), gradation, young, poisson);
  particleVec.push_back(newPtcl);
  ++particleNum;
  //lastPtcls.push_back(newPtcl);

  newPtcl = new Particle(particleNum+1, 0, Vec(0, -maxDiameter, zCurr +
  maxDiameter * i), gradation, young, poisson);
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

  // 6. update particles' velocity/omga/position/orientation based on
  force/moment.
  updateParticle();

  // 7. (1) output particles and contacts information as snapNum.
  if (toSnapshot) {
  sprintf(stepsstr, "%03d", stepsnum);
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp,
  stepsstr);
  printParticle(stepsfp);

  sprintf(stepsstr, "%03d", stepsnum);
  strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
  printContact(stepsfp);
  time(&timeStamp);
  g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) <<
  std::flush;
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

  } while (particleNum < 2000); //( zCurr < container.getMaxCorner().getZ() );
  //(++iteration < totalSteps);

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
  buildBoundary(1,              // 1-only bottom boundary; 5-no top
  boundary;6-boxed 6 boundaries
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

  std::vector<Particle*>::iterator itr;
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

  std::vector<Particle*>::iterator itr;
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
  if(!ofs) { debugInf << "stream error: printCavityParticle" << std::endl;
  exit(-1); }
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
  std::vector<Particle*>::const_iterator  it;
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
  if(!ofs) { debugInf << "stream error: buildCavityBoundary" << std::endl;
  exit(-1); }

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

  std::vector<Particle*> vec1d;  // 1-dimension
  std::vector< std::vector<Particle*>  > vec2d; // 2-dimension
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
  std::vector<Particle*> x1y1;
  std::vector<Particle*> x1y2;
  std::vector<Particle*> x1z1;
  std::vector<Particle*> x1z2;

  std::vector<Particle*> x2y1;
  std::vector<Particle*> x2y2;
  std::vector<Particle*> x2z1;
  std::vector<Particle*> x2z2;

  std::vector<Particle*> y1x1;
  std::vector<Particle*> y1x2;
  std::vector<Particle*> y1z1;
  std::vector<Particle*> y1z2;

  std::vector<Particle*> y2x1;
  std::vector<Particle*> y2x2;
  std::vector<Particle*> y2z1;
  std::vector<Particle*> y2z2;

  std::vector<Particle*> z1x1;
  std::vector<Particle*> z1x2;
  std::vector<Particle*> z1y1;
  std::vector<Particle*> z1y2;

  std::vector<Particle*> z2x1;
  std::vector<Particle*> z2x2;
  std::vector<Particle*> z2y1;
  std::vector<Particle*> z2y2;

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

  std::vector<Particle*>::iterator itr;
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
  if(!progressInf) { debugInf << "stream error: deGravitation" << std::endl;
  exit(-1); }
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
  if(!debugInf) { debugInf << "stream error: deGravitation" << std::endl;
  exit(-1); }
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles from existing files.
  if (toRebuild) readParticle(iniptclfile); // create container and particles,
  velocity and omga are set zero.

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

  // 4. update particles' velocity/omga/position/orientation based on
  force/moment.
  updateParticle();

  // 5. (1) output particles and contacts information as snapNum.
  if (iteration % (totalSteps/snapNum) == 0) {
  sprintf(stepsstr, "%03d", stepsnum);
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp,
  stepsstr);
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
  if(!progressInf) { debugInf << "stream error: deposit_p" << std::endl;
  exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "deposit..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: deposit_p" << std::endl; exit(-1);
  }
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from existing files.
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.

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

  // 4. update particles' velocity/omga/position/orientation based on
  force/moment.
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
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp,
  stepsstr);
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
  if(!progressInf) { debugInf << "stream error: squeeze" << std::endl; exit(-1);
  }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "deposit..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: squeeze" << std::endl; exit(-1); }
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from existing files.
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.
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

  // 5. update particles' velocity/omga/position/orientation based on
  force/moment.
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
  strcpy(stepsfp, ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp,
  stepsstr);
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
  if(!progressInf) { debugInf << "stream error: isoMemBdry" << std::endl;
  exit(-1);}
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
  if(!debugInf) { debugInf << "stream error: isoMemBdry" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles from file and calculate forces caused by hydraulic
  pressure
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
  std::vector<Particle*>::const_iterator  it;
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

  // 5. update particles velocity/omga/position/orientation based on
  force/moment
  updateParticle();

  // 6. (1) output particles and contacts information
  if (iteration % (totalSteps/snapNum) == 0) {
  sprintf(stepsstr, "%03d", stepsnum);
  strcpy(stepsfp,ParticleFile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
  printParticle(stepsfp);

  strcpy(stepsfp,"iso_membrane_"); strcat(stepsfp, stepsstr);
  printMemParticle(stepsfp);
  strcpy(stepsfp,"iso_spring_"); strcat(stepsfp, stepsstr); strcat(stepsfp,
  ".dat");
  plotSpring(stepsfp);

  sprintf(stepsstr, "%03d", stepsnum);
  strcpy(stepsfp, contactfile); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
  printContact(stepsfp);
  time(&timeStamp);
  g_timeinf << std::setw(4) << stepsnum << " " << ctime(&timeStamp) <<
  std::flush;
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
  if(!progressInf) { debugInf << "stream error: triaxialPtclBdryIni" <<
  std::endl; exit(-1);}
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "triaxial..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average        sample
  sample     "
  << "     sample          sample          sample          sample
  sample          "
  << "sample          sample         sample           sample          sample
  sample     "
  << "     sample          sample          sample          void
  sample        coordinate"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "          omga            force           moment        density          "
  << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h
  epsilon_v"
  << "        ratio          porosity         number"
  << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: triaxialPtclBdryIni" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from files
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.
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

  // 5. update particles' velocity/omga/position/orientation based on
  force/moment
  updateParticle();

  // 6. update boundaries' position and orientation
  sigma3_1=vfabs(getNormalForce(5))/2.5e-3;
  sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

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
  if (   fabs(sigma3_1-sigma)/sigma < boundaryStressTol &&
  fabs(sigma3_2-sigma)/sigma < boundaryStressTol )
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
  if(!progressInf) { debugInf << "stream error: triaxialPtclBdry" << std::endl;
  exit(-1);}
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "triaxial..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average        sample
  sample     "
  << "     sample          sample          sample          sample
  sample          "
  << "sample          sample         sample           sample          sample
  sample     "
  << "     sample          sample          sample          void
  sample        coordinate"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "          omga            force           moment        density          "
  << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h
  epsilon_v"
  << "        ratio          porosity         number"
  << std::endl;

  std::ofstream balancedinf(balancedfile);
  if(!balancedinf) { debugInf << "stream error: triaxialPtclBdry" << std::endl;
  exit(-1);}
  balancedinf.setf(std::ios::scientific, std::ios::floatfield);
  balancedinf << "triaxial..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average        sample
  sample     "
  << "     sample          sample          sample          sample
  sample          "
  << "sample          sample         sample           sample          sample
  sample     "
  << "     sample          sample          sample          void
  sample        coordinate"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "          omga            force           moment        density          "
  << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h
  epsilon_v"
  << "        ratio          porosity         number"
  << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: triaxialPtclBdry" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from files
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.
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

  // 5. update particles' velocity/omga/position/orientation based on
  force/moment
  updateParticle();

  // 6. update boundaries' position and orientation
  sigma3_1=vfabs(getNormalForce(5))/2.5e-3;
  sigma3_2=vfabs(getNormalForce(6))/2.5e-3;

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


  // The specimen has been deposited with gravitation within boundaries composed
  of particles.
  // A rectangular pile is then drived into the particles using displacement
  control.
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
  if(!progressInf) { debugInf << "stream error:recPile_Disp" << std::endl;
  exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "pile penetrate..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential        total           sample           sample "
  << "     sample          sample          sample          sample
  sample          "
  << "sample          sample         sample           sample          sample
  sample"
  << "          sample          sample          sample" << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy          energy          density         "
  << "sigma1_1        sigma1_2        sigma2_1        sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: recPile_Disp" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);
  debugInf << " iteration    end_bearing     side_friction   total_force" <<
  std::endl;

  // pre_2. create particles and boundaries from files
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.
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

  // 5. update particles' velocity/omga/position/orientation based on
  force/moment
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


  // The specimen has been deposited with gravitation within boundaries composed
  of particles.
  // An ellipsoidal pile is then drived into the particles using displacement
  control.
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
  if(!progressInf) { debugInf << "stream error: ellipPile_Disp" << std::endl;
  exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "pile penetrate..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: ellipPile_Disp" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from files
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.

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

  // 4. update particles' velocity/omga/position/orientation based on
  force/moment
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
  // An ellipsoidal penetrator is then impacted into the particles with initial
  velocity.
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
  if(!progressInf) { debugInf << "stream error: ellipPile_Impact" << std::endl;
  exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "penetrator impact..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: ellipPile_Impact" << std::endl;
  exit(-1);}
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

  // 5. update particles' velocity/omga/position/orientation based on
  force/moment
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


  // The specimen has been deposited with gravitation within particle
  boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial
  velocity.
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
  if(!progressInf) { debugInf << "stream error: ellipPile_Impact_p" <<
  std::endl; exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "penetrator impact..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: ellipPile_Impact_p" << std::endl;
  exit(-1);}
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

  // 4. update particles' velocity/omga/position/orientation based on
  force/moment
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



  // The specimen has been deposited with gravitation within boundaries composed
  of particles.
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
  if(!progressInf) { debugInf << "stream error: ellipPile_Force" << std::endl;
  exit(-1); }
  progressInf.setf(std::ios::scientific, std::ios::floatfield);
  progressInf << "pile penetrate..." << std::endl
  << "     iteration possible  actual      average	    average
  average         average"
  << "         average         average         average       translational
  rotational       "
  << "kinetic        potential         total           void            sample
  coordination"
  << "       sample           sample          sample          sample
  sample          sample"
  << "          sample          sample          sample         sample
  sample         "
  << " sample          sample          sample          sample          sample"
  << std::endl
  << "       number  contacts contacts   penetration   contact_normal
  contact_tangt     velocity"
  << "         omga            force           moment         energy
  energy          "
  << "energy         energy            energy          ratio          porosity
  number       "
  << "   density         sigma1_1        sigma1_2        sigma2_1
  sigma2_2        "
  << "sigma3_1        sigma3_2           p             width          length "
  << "height          volume         epsilon_w       epsilon_l       epsilon_h "
  << "epsilon_v" << std::endl;

  std::ofstream balancedinf(balancedfile);
  if(!balancedinf) { debugInf << "stream error: ellipPile_Force" << std::endl;
  exit(-1);}
  balancedinf.setf(std::ios::scientific, std::ios::floatfield);
  balancedinf << "pile penetrate..." << std::endl
  << "   iteration   apply_force    pile_tip_pos     pile_force" << std::endl;

  debugInf.open(debugfile);
  if(!debugInf) { debugInf << "stream error: ellipPile_Force" << std::endl;
  exit(-1);}
  debugInf.setf(std::ios::scientific, std::ios::floatfield);

  // pre_2. create particles and boundaries from files
  readParticle(iniptclfile); // create container and particles, velocity and
  omga are set zero.

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

  // 4. update particles' velocity/omga/position/orientation based on
  force/moment
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
