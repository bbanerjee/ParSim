#include <Peridynamics/Peridynamics.h>

#include <InputOutput/PeriParticleFileReader.h>
#include <Core/Const/const.h>
#include <Core/Math/IntVec.h>
#include <chrono>

using namespace pd;

using Timer = std::chrono::steady_clock;
using Seconds = std::chrono::seconds;
using IntVec = dem::IntVec;
using Vec = dem::Vec;
using Matrix = dem::Matrix;
using Box = dem::Box;
using ParticlePArray = dem::ParticlePArray;

Peridynamics::Peridynamics() 
{
}

Peridynamics::~Peridynamics()
{
  d_allPeriParticlesInitial.clear();
  allPeriParticleVec.clear(); // this is part of d_allPeriParticlesInitial
  interfacePeriParticleVec.clear();
  outerfacePeriParticleVec.clear();

  // delete memories in current cpus
  periParticleVec.clear();
  fixedPeriParticleVec.clear(); // this is part of periParticleVec
  bottomBoundaryVec.clear();
  frontBoundaryVec.clear();
  leftBoundaryVec.clear();
  topBoundaryInnerVec.clear();
  topBoundaryEdgeVec.clear();
  topBoundaryCornerVec.clear();

  topBoundaryBondVec.clear();
  bottomBoundaryBondVec.clear();
  leftBoundaryBondVec.clear();
  rightBoundaryBondVec.clear();
  frontBoundaryBondVec.clear();
  backBoundaryBondVec.clear();

  recvPeriBondVec.clear();
  periBondVec.clear();
} // ~Peridynamics()

// Pure peridynamics only
// not used in DEM-PD coupling code, i.e. rigidInclusion()
void
Peridynamics::initializePeridynamics(const std::string& inputFile)
{
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Problem Initialization " << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  //std::cout << "Read data file ..." << std::endl;
  readPeriDynamicsData(inputFile);

  //std::cout << "Calculate particle volume ..." << std::endl;
  calcParticleVolume();
  // writeMeshCheckVolume("checkv.dat"); exit(1);

  //std::cout << "Calculate horizon size ..." << std::endl;
  calcHorizonSize();

  //std::cout << "Construct neighor list ..." << std::endl;
  constructNeighbor();

  //std::cout << "Calculate Kinv ..." << std::endl;
  calcParticleKinv();

  for (auto& particle : periParticleVec) {
    particle->initial();
  }

  // prescrible the essential boundary condition
  //std::cout << "Prescribe the boundary condition ..." << std::endl;
  prescribeEssentialBoundaryCondition(0);

  // calculate the stress at each particle
  //std::cout << "Calculate particle stress ..." << std::endl;
  calcParticleStress();

  // calculate the acceleration for each particle
  //std::cout << "Calculate particle acceleration ..." << std::endl;
  calcParticleAcceleration();

  // traction boundary
  applyTractionBoundary(0);

  // apply initial velocity boundary condition, in this case give the
  // cubicTopBoundaryVec particles initial velocities
  for (auto& particle : cubicTopBoundaryVec) {
    particle->setInitVelocity(Vec(0.0, 0.0, 1.0));
  }
} // end initializePeridynamics

// readData - reads particle positions and mesh connectivities
// @param std::string&  - reference of the input file name
// Caveats: should only be called by master cpu
void
Peridynamics::readPeriDynamicsData(const std::string& inputFile)
{ 
  PeriParticleFileReader reader;
  reader.read(inputFile, d_allPeriParticlesInitial, d_connectivity);

  // Compute the maximum distance between adjacent particles
  d_maxDistBetweenParticles = -1.0e16;
  for (auto& element : d_connectivity) {
    auto node0 = element[0]; 
    for (auto node : element.nodes) {
      auto maxDist = vfabs(d_allPeriParticlesInitial[node0-1]->getInitPosition() -
                           d_allPeriParticlesInitial[node-1]->getInitPosition());
      d_maxDistBetweenParticles = (d_maxDistBetweenParticles > maxDist) ? d_maxDistBetweenParticles : maxDist;
    }
  }

} // readPeriDynamicsData()

void
Peridynamics::calcParticleVolume()
{
  if (d_connectivity[0][4] == 0) {
    calcParticleVolumeForTet();
  } else {
    calcParticleVolumeForHex();
  }
}

void
Peridynamics::calcParticleVolumeForHex()
{
  std::vector<int> numOfPieces(nPeriParticle, 0);
  std::vector<REAL> particleVolume(nPeriParticle, 0.0);
  REAL xi, eta, zeta;

  int nip = 2;
  Matrix gp_loc3D;
  Matrix gp_weight3D;
  gauss3D(nip, gp_loc3D, gp_weight3D);

  Matrix xl(3, 8);
  Matrix shp;

  for (const auto& element : d_connectivity) {

    // Compute the node position matrix and the number of 
    // Elements attached to a node
    int nodeIndex = 1;
    for (const auto& node : element.nodes) {
      const auto& pos = d_allPeriParticlesInitial[node-1]->getInitPosition();
      xl(1, nodeIndex) = pos.x();
      xl(2, nodeIndex) = pos.y();
      xl(3, nodeIndex) = pos.z();
      numOfPieces[node-1] += 1;
      nodeIndex++;
    }

    // Compute the volume of the element and the contribution
    // of each node
    for (int ik = 0; ik < 8; ik++) {
      xi = gp_loc3D(1, ik + 1);
      eta = gp_loc3D(2, ik + 1);
      zeta = gp_loc3D(3, ik + 1);
      // call fem function to get the volume
      REAL xsj = 0.0;
      shp3d(xi, eta, zeta, xl, shp, xsj);
      for (const auto& node : element.nodes) {
        particleVolume[node-1] += gp_weight3D(1, ik + 1) * xsj / 8.0;
      }
    }
  }

  int particleIndex = 0;
  for (auto& volume : particleVolume) {
    volume *= 8.0 / (REAL(numOfPieces[particleIndex]));
    particleIndex++;
  }

  // store the particle volume into the object PeriParticle
  particleIndex = 0;
  for (auto& particle : d_allPeriParticlesInitial) {
    particle->setParticleVolume(particleVolume[particleIndex]);
    particleIndex++;
  }

} // end calcParticleVolume

void
Peridynamics::calcParticleVolumeForTet()
{
}

void
Peridynamics::calcHorizonSize()
{
  if (d_connectivity[0][4] == 0) {
    calcHorizonSizeForTet();
  } else {
    calcHorizonSizeForHex();
  }
}

void
Peridynamics::calcHorizonSizeForHex()
{
  d_maxHorizonSize = 0;
  for (const auto& element : d_connectivity) {
    // get initial positions of the particles in this d_connectivity
    int n1 = element[0]; 
    int n2 = element[1];
    int n4 = element[3];
    int n5 = element[4];
    Vec coord1 = d_allPeriParticlesInitial[n1 - 1]->getInitPosition();
    Vec coord2 = d_allPeriParticlesInitial[n2 - 1]->getInitPosition();
    Vec coord4 = d_allPeriParticlesInitial[n4 - 1]->getInitPosition();
    Vec coord5 = d_allPeriParticlesInitial[n5 - 1]->getInitPosition();

    REAL tmp1 = vfabs(coord2 - coord1);
    REAL tmp2 = vfabs(coord4 - coord1);
    REAL tmp3 = vfabs(coord5 - coord1);

    REAL tmpmax = 1.5075 * std::max(std::max(tmp1, tmp2), std::max(tmp2, tmp3));
    d_maxHorizonSize = (tmpmax > d_maxHorizonSize) ? tmpmax : d_maxHorizonSize;

    for (const auto& node : element.nodes) {
      d_allPeriParticlesInitial[node - 1]->replaceHorizonSizeIfLarger(tmpmax);
    }

  }
  //std::cout << "d_maxHorizonSize: " << d_maxHorizonSize << std::endl;

} // end calcHorizonSize()

void
Peridynamics::calcHorizonSizeForTet()
{
}

// this function should be called after
// scattering in each cpu to construct
// peri-bonds
// neighbor - searches and constructs the particle neighborlists
// construct neighborlist for all particles ...
// compute the weighting function for all particles ...
void
Peridynamics::constructNeighbor()
{ 
  if (periParticleVec.empty())
    return;

  for (auto& particle : periParticleVec) {
    particle->clearPeriBonds(); // bondVec should be empty at this time
  }
  periBondVec.clear();

  int ompThreads = util::getParam<int>("ompThreads");
  std::size_t i_nt, j_nt;
  std::size_t num = periParticleVec.size();

#pragma omp parallel for \
        num_threads(ompThreads) \
        private(i_nt, j_nt) \
        shared(num) \
        schedule(dynamic)
  for (i_nt = 0; i_nt < num - 1; i_nt++) {

    Vec coord0_i = periParticleVec[i_nt]->getInitPosition();
    REAL horizonSize_i = periParticleVec[i_nt]->getHorizonSize();

    for (j_nt = i_nt + 1; j_nt < num; j_nt++) {

      Vec coord0_j = periParticleVec[j_nt]->getInitPosition();
      REAL horizonSize_j = periParticleVec[j_nt]->getHorizonSize();

      REAL bond_length = vfabs(coord0_i - coord0_j);

      // This will lead to the fact that horizion is not a sphere!!!
      REAL horizonSize_ij = (horizonSize_i + horizonSize_j) * 0.5; 
      REAL ratio = bond_length / horizonSize_ij;

      // establish the neighbor list
      if (ratio <= 2.0) {

        // create bond
        PeriBondP bond = std::make_shared<PeriBond>(
          bond_length, periParticleVec[i_nt], periParticleVec[j_nt]);

#pragma omp critical
        {
          periParticleVec[i_nt]->pushBackBondVec(bond);
          periParticleVec[j_nt]->pushBackBondVec(bond);
          periBondVec.push_back(bond);
        }

        // for the factor of 3d window function
        REAL factor = 1.5 / (dem::Pi * horizonSize_ij * horizonSize_ij * horizonSize_ij); 

        // weighting function (influence function)
        if (ratio < 1.0) {
          bond->setWeight(
            factor * (2.0 / 3.0 - ratio * ratio + 0.5 * ratio * ratio * ratio));
        } else {
          bond->setWeight(factor * (2.0 - ratio) * (2.0 - ratio) *
                             (2.0 - ratio) / 6.0);
        }
      } // if(ratio<2.0)

    } // end j_nt
  }   // end i_nt

} // end constNeighbor()

// this should be called after
// commuDEMPeriParticle(), and find
// peri-bonds between communicated
// periParticles
// Compute the inverse of the shape tensor K
void
Peridynamics::calcParticleKinv()
{ 
  for (auto& particle : periParticleVec) {
    particle->calcParticleKinv();
  }
  for (auto& particle : recvPeriParticleVec) {
    particle->calcParticleKinv();
  }

} // end calcParticleKinv()

void
Peridynamics::prescribeEssentialBoundaryCondition(const int istep)
{
  // fix z displacement in the z direction
  for (auto& particle : bottomBoundaryVec) {
    particle->prescribeBottomDisplacement(0.0);    
  }

  REAL dispz = 1.8*0.05;
  if (istep <= 200) {
    dispz *= double(istep)/200.;
  }

  for (auto& particle : topBoundaryVec) {
    particle->prescribeTopDisplacement(dispz);    
  }

} // end prescribeEssentialBoundaryCondition()

void
Peridynamics::calcParticleStress()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num = periParticleVec.size();
  int i;

#pragma omp parallel for \
        num_threads(ompThreads) \
        private(i) \
        shared(num) \
        schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->calcParticleStress();
  }
  for (auto& particle : recvPeriParticleVec) {
    particle->calcParticleStress();
  }
} // calcParticleStress()

void
Peridynamics::calcParticleAcceleration()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num = periParticleVec.size();
  int i;

#pragma omp parallel for \
        num_threads(ompThreads) \
        private(i) \
        shared(num) \
        schedule(dynamic)

  for (i = 0; i < num; i++) {
    periParticleVec[i]->calcParticleAcceleration();
  }

} // calcParticleAcceleration()

void
Peridynamics::applyTractionBoundary(int g_iteration)
{
  // set traction boundary
  REAL force = util::getParam<REAL>("periForce");
  REAL rampStep = util::getParam<REAL>("rampStep");
  REAL framp = 0;
  if (g_iteration <= rampStep) {
    framp = g_iteration / rampStep;
  } else {
    framp = 1;
  }
  framp = framp * force;
  dem::Vec tmp_vec = dem::Vec(0, 0, framp);
  for (auto& peri_pt : topBoundaryInnerVec) {
    peri_pt->addAccelerationByForce(tmp_vec);
  }
  for (auto& peri_pt : topBoundaryEdgeVec) {
    peri_pt->addAccelerationByForce(tmp_vec * 0.5);
  }
  for (auto& peri_pt : topBoundaryCornerVec) {
    peri_pt->addAccelerationByForce(tmp_vec * 0.25);
  }
  for (auto& peri_pt : bottomBoundaryVec) {
    peri_pt->addAccelerationByForce(-tmp_vec);
  }

} // applyTractionBoundary()

// delete those peri-points that are inside sand particles
// actually what are deleted are the peri-points that are inside
// a enlarged sand particle to avoid peri-points go to inside the
// sand after a few steps. Here the enlarged sand particle is a little
// larger than the sand particle, the delta = 0.5*maxDistBetweenParticles and
// delta is the difference of principle lengths
void
Peridynamics::removeInsidePeriParticles(const ParticlePArray& allDEMParticleVec)
{
  bool is_inside;
  allPeriParticleVec.clear();

  for (auto& peri_pt : d_allPeriParticlesInitial) {

    Vec coord_peri = peri_pt->currentPosition();

    // if this peri-point is inside sand particles
    is_inside = false; 

    // remove the inside peri-points
    for (auto& dem_pt : allDEMParticleVec) {

      // enlarged sand particle
      REAL a = dem_pt->getA() + 0.5 * d_maxHorizonSize;
      REAL b = dem_pt->getB() + 0.5 * d_maxHorizonSize;
      REAL c = dem_pt->getC() + 0.5 * d_maxHorizonSize;

      // this is important to get local coordinate
      Vec coord_peri_local = 
        dem_pt->globalToLocal(coord_peri - dem_pt->currentPosition()); 

      REAL x_pt = coord_peri_local.x() / a;
      REAL y_pt = coord_peri_local.y() / b;
      REAL z_pt = coord_peri_local.z() / c;

      if (x_pt * x_pt + y_pt * y_pt + z_pt * z_pt < 1) { 
        // this peri-points is inside sand particle
        is_inside = true;
        break; // do not need to check other sand particles
      }
    }

    if (is_inside == false) { 
      // this peri-point is not within any sand particles
      allPeriParticleVec.push_back(peri_pt);
    }
  }
} // end removeInsidePeriParticles()

// this is to scatter the peridynamics particles
// two point: 
// (1) Partition the peri-points before any bonds constructed, all
//     peri-points will be treated as free peri-points.
//     This is also what the coupling model shows, i.e. even the
//     peri-points that are bonded to DEM particles still have properties
//     of free peri-points
// (2) Before partition there are no bonds in PeriParticle, after
//     partition, constructNeighbor() is called, then these peri-bonds,
//     boundary-bonds, and peri-DEM-bonds will be generated
void
Peridynamics::scatterPeriParticle(const Box& allContainer)
{
  if (mpiRank == 0) { // process 0
    setGrid(Box(allContainer, d_maxDistBetweenParticles*0.2));

    Vec v1 = d_periGrid.getMinCorner();
    Vec v2 = d_periGrid.getMaxCorner();
    Vec vspan = v2 - v1;
    vspan /= d_mpiProcs;

    boost::mpi::request reqs[mpiSize - 1];

    PeriParticlePArray insidePeriParticleVec;
    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      insidePeriParticleVec.clear(); // do not release memory!

      int ndim = 3;
      IntVec coords;
      MPI_Cart_coords(cartComm, iRank, ndim, coords.data());

      Vec spanLo = vspan*coords;
      Vec spanHi = vspan*(coords + 1);
      Box container(v1 + spanLo, v1 + spanHi);

      findPeriParticleInBox(container, allPeriParticleVec, insidePeriParticleVec);

      if (iRank != 0) {

        // non-blocking send
        reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag, insidePeriParticleVec); 

      } else {

        // Make a deep copy
        periParticleVec.resize(insidePeriParticleVec.size());
        for (auto i = 0u; i < periParticleVec.size(); ++i)
          periParticleVec[i] = std::make_shared<PeriParticle>(
            *insidePeriParticleVec[i]); 

      } // now particleVec do not share memeory with allParticleVec
    }

    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send

  } else { // other processes except 0

    boostWorld.recv(0, mpiTag, periParticleVec);

  }

  // broadcast necessary info
  broadcast(boostWorld, d_maxDistBetweenParticles, 0);
  broadcast(boostWorld, d_maxHorizonSize, 0);

} // scatterDEMPeriParticle

void
Peridynamics::findPeriParticleInBox(const Box& container,
                                    const PeriParticlePArray& periParticles,
                                    PeriParticlePArray& insideParticles)
{
  for (const auto& particle : periParticles) {
    // it is critical to use EPS
    if (container.inside(particle->currentPosition(), dem::EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

// construct the Matrix members in
// periParticleVec,
// since currently the transfer of pointer array in Matrix is not implemented
// well
void
Peridynamics::constructPeriMatrix()
{ 
  for (auto& pt : periParticleVec) {
    pt->constructMatrixMember(); // construct these Matrix members
  }
} // constructPeriMatrix()

void 
Peridynamics::createPatch(int iteration, 
                          const REAL& ghostWidth) 
{
  // determine container of each process
  Vec v1 = d_periGrid.getMinCorner();
  Vec v2 = d_periGrid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP = std::make_unique<PeriPatch>(cartComm, 
    mpiRank, d_mpiCoords, lower, upper, ghostWidth, dem::EPS);
}

void
Peridynamics::updatePatch(int iteration,
                          const REAL& ghostWidth)
{
  // determine container of each process
  Vec v1 = d_periGrid.getMinCorner();
  Vec v2 = d_periGrid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  Box container(lower, upper);
  d_patchP->update(iteration, lower, upper, ghostWidth);
}

// before the transfer to peri-points, their bondVec should be cleared
void
Peridynamics::commuPeriParticle(int iteration,
                                const double& maxDEMParticleRadius)
{
  // constructNeighbor() is based on 2*horizonSize,
  // refer to PeriParitcle.calcAcceleration(), in this function the
  // deformationGradient and sigma of the other peri-points in the bond are
  // needed,
  // this means that the deformationGradient, sigma and Kinv of the other
  // peri-point should also be calculated even if it is in recvParticleVec, thus
  // we need to transfer 2*cellSize peri-points, the peri-points in the outer
  // cell are used to calculate the deformationGradient, sigma and Kinv
  // of the peri-points in inner cell
  REAL cellSize = std::max(4*d_maxHorizonSize,
                           maxDEMParticleRadius + 3*d_maxDistBetweenParticles); 
  updatePatch(iteration, cellSize);
  
  // Initialize merged particle array (particles + ghosts)
  mergePeriParticleVec.clear();
  mergePeriParticleVec = periParticleVec;

  // Initialize an array for just the received particles during each phase
  //PeriParticlePArray recvParticleVec;

  // Plimpton scheme: x-ghost exchange
  d_patchP->sendRecvGhostXMinus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->sendRecvGhostXPlus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineReceivedParticlesX(iteration, mergePeriParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvPeriParticleVec);

  // Plimpton scheme: y-ghost exchange
  d_patchP->sendRecvGhostYMinus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->sendRecvGhostYPlus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineReceivedParticlesY(iteration, mergePeriParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvPeriParticleVec);

  // Plimpton scheme: z-ghost exchange
  d_patchP->sendRecvGhostZMinus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->sendRecvGhostZPlus(boostWorld, iteration, mergePeriParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineReceivedParticlesZ(iteration, mergePeriParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvPeriParticleVec);

  // Update the bonds for the particles that have been
  // received into a patch
  if (!recvPeriParticleVec.empty()) {
    updatePeridynamicBonds(recvPeriParticleVec, recvPeriBondVec, mergePeriParticleVec);
  }
}

// This function should be called after commuPeriParticle() in each cpu to
// construct peri-bonds and neighbors - 
// searches and constructs the particle neighborlists
// construct neighborlist for all particles ...
// compute the weighting function for all particles ...
void
Peridynamics::updatePeridynamicBonds(PeriParticlePArray& recvParticles,
                                     PeriBondPArray& recvPeriBondVec,
                                     PeriParticlePArray& periParticleVec)
{
  // Clear the new bond vector
  recvPeriBondVec.clear();

  for (auto& recvParticle : recvParticles) {

    // Clear the bonds associated with the received particle
    recvParticle->clearPeriBonds();

    // Get the initial position and horizon of the received particle
    Vec recvPos = recvParticle->getInitPosition();
    REAL recvHorizon = recvParticle->getHorizonSize();

    // Find the neighbors in the current patch (the merged data)
    for (auto& particle : periParticleVec) {

      // Make sure the particles are not identical
      if (recvParticle->getId() == particle->getId()) continue;

      // Get the initial position and horizon of the patch particle
      Vec pos = particle->getInitPosition();
      REAL horizon = particle->getHorizonSize();

      // Compute bond length
      REAL bondLength = vfabs(recvPos - pos);

      // Find average horizon size of the two particles
      REAL horizonSize = recvHorizon + horizon;

      // If the bondlength < horizonSize add neighbor
      if (bondLength < horizonSize) {

        // create bond
        PeriBondP bond =
          std::make_shared<pd::PeriBond>(bondLength, recvParticle, particle);
        
        // Add the bond to the two particles
        recvParticle->pushBackBondVec(bond);
        particle->pushBackBondVec(bond);

        // Add the bond to the list of new bonds
        recvPeriBondVec.push_back(bond);

        // Compute the weighting function (influence function)
        // using the factor for a 3d window function
        REAL ratio =  2.0*bondLength/horizonSize;
        REAL factor = 3.0/(2.0*dem::Pi*horizonSize*horizonSize*horizonSize); 
        if (ratio < 1.0) {
          REAL ratioSq = ratio*ratio;
          REAL ratioCu = ratioSq*ratio;
          bond->setWeight(factor*(2.0/3.0 - ratioSq + 0.5*ratioCu));
        } else {
          bond->setWeight(factor*(2.0-ratio)*(2.0 - ratio)*(2.0-ratio)/6.0);
        }

      }

    }
  }
}

void
Peridynamics::findRecvPeriBonds() 
{
  // Update the bonds for the particles that have been
  // received into a patch
  if (!recvPeriParticleVec.empty()) {
    updatePeridynamicBonds(recvPeriParticleVec, recvPeriBondVec, mergePeriParticleVec);
  }
}

bool
Peridynamics::isBdryProcess()
{
  return (d_mpiCoords.x() == 0 || d_mpiCoords.x() == d_mpiProcs.x() - 1 ||
          d_mpiCoords.y() == 0 || d_mpiCoords.y() == d_mpiProcs.y() - 1 ||
          d_mpiCoords.z() == 0 || d_mpiCoords.z() == d_mpiProcs.z() - 1);
}

void
Peridynamics::releaseRecvPeriParticle()
{
  recvPeriBondVec.clear();
  PeriBondPArray().swap(recvPeriBondVec); // actual memory release

  // release memory of received particles
  recvPeriParticleVec.clear();
  PeriParticlePArray().swap(recvPeriParticleVec); // actual memory release
}

void
Peridynamics::releasePeriBondVec()
{
  periBondVec.clear();
  PeriBondPArray().swap(periBondVec); // actual memory release
}

void
Peridynamics::updatePeriGrid(const PeriParticlePArray& particles)
{
  REAL minX = 0.0, minY = 0.0, minZ = 0.0;
  REAL maxX = 0.0, maxY = 0.0, maxZ = 0.0;

  REAL minval = 1.0e16;
  REAL maxval = -1.0e16;

  if (particles.size() > 0) {
    REAL pMaxX = maxval, pMaxY = maxval, pMaxZ = maxval;
    REAL pMinX = minval, pMinY = minval, pMinZ = minval;

    REAL xPos, yPos, zPos;
    dem::Vec position;
    for (const auto& particle : particles) {
      position = particle->currentPosition();
      xPos = position.x();
      yPos = position.y();
      zPos = position.z();
      pMinX = (pMinX < xPos) ? pMinX : xPos;
      pMinY = (pMinY < yPos) ? pMinY : yPos;
      pMinZ = (pMinZ < zPos) ? pMinZ : zPos;
      pMaxX = (pMaxX > xPos) ? pMaxX : xPos;
      pMaxY = (pMaxY > yPos) ? pMaxY : yPos;
      pMaxZ = (pMaxZ > zPos) ? pMaxZ : zPos;
    }
    MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  } else {
    MPI_Allreduce(&maxval, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&minval, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&maxval, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&minval, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
    MPI_Allreduce(&maxval, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
    MPI_Allreduce(&minval, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  }

  // no need to broadcast grid as it is updated in each process
  setGrid(Box(minX - d_maxDistBetweenParticles * 0.2, 
              minY - d_maxDistBetweenParticles * 0.2,
              minZ - d_maxDistBetweenParticles * 0.2, 
              maxX + d_maxDistBetweenParticles * 0.2,
              maxY + d_maxDistBetweenParticles * 0.2, 
              maxZ + d_maxDistBetweenParticles * 0.2));
}

void
Peridynamics::migratePeriParticle(int iteration)
{
  // Compute the (x,y,z) patch dimensions
  Vec domainWidth = d_periGrid.getMaxCorner() - d_periGrid.getMinCorner();
  Vec patchWidth = domainWidth / d_mpiProcs;

  // Migrate particles in the x-direction
  dem::ParticleIDHashMap sentParticles;
  PeriParticlePArray recvParticles;
  d_patchP->sendRecvMigrateXMinus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->sendRecvMigrateXPlus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineSentParticlesX(iteration, sentParticles);
  d_patchP->combineReceivedParticlesX(iteration, recvParticles);
  d_patchP->deleteSentParticles<PeriParticleP>(iteration, sentParticles, periParticleVec);
  for (auto& particle : recvParticles) {
    particle->constructMatrixMember();
  }
  d_patchP->addReceivedParticles(iteration, recvParticles, periParticleVec);

  // Migrate particles in the y-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateYMinus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->sendRecvMigrateYPlus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineSentParticlesY(iteration, sentParticles);
  d_patchP->combineReceivedParticlesY(iteration, recvParticles);
  d_patchP->deleteSentParticles<PeriParticleP>(iteration, sentParticles, periParticleVec);
  for (auto& particle : recvParticles) {
    particle->constructMatrixMember();
  }
  d_patchP->addReceivedParticles(iteration, recvParticles, periParticleVec);

  // Migrate particles in the z-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateZMinus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->sendRecvMigrateZPlus(boostWorld, iteration, patchWidth, periParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineSentParticlesZ(iteration, sentParticles);
  d_patchP->combineReceivedParticlesZ(iteration, recvParticles);
  d_patchP->deleteSentParticles<PeriParticleP>(iteration, sentParticles, periParticleVec);
  for (auto& particle : recvParticles) {
    particle->constructMatrixMember();
  }
  d_patchP->addReceivedParticles(iteration, recvParticles, periParticleVec);

  // delete outgoing particles
  d_patchP->removeParticlesOutsidePatch<PeriParticleP>(periParticleVec);
}

// update allPeriParticleVec: process 0 collects all updated particles from
// each other process
void
Peridynamics::gatherPeriParticle()
{
  // store Matrix sigma to each value
  for (auto& particle : periParticleVec) {
    particle->assignSigma(); 
  }

  if (mpiRank != 0) {

    boostWorld.send(0, mpiTag, periParticleVec);

  } else { 

    // allPeriParticleVec is cleared before filling with new data
    releaseGatheredPeriParticle();

    // duplicate PeriParticleVec so that it is not destroyed by
    // allPeriParticleVec in next iteration,
    // otherwise it causes memory error.
    PeriParticlePArray dupPatchParticles(periParticleVec.size());
    std::size_t index = 0;
    for (const auto& particle : periParticleVec) {
      dupPatchParticles[index] = std::make_shared<PeriParticle>(*particle);
      dupPatchParticles[index]->releaseBondVec();
      index++;
    }

    // fill allParticleVec with dupParticleVec 
    allPeriParticleVec.insert(allPeriParticleVec.end(),
                              dupPatchParticles.begin(),
                              dupPatchParticles.end());

    // and add received particles from all non-zero ranks
    //long gatherRam = 0;
    for (int iRank = 1; iRank < mpiSize; ++iRank) {
      PeriParticlePArray recvParticles;
      boostWorld.recv(iRank, mpiTag, recvParticles);
      allPeriParticleVec.insert(allPeriParticleVec.end(),
                                recvParticles.begin(),
                                recvParticles.end());
      //gatherRam += recvParticles.size();
    }
    // dem::debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    // << gatherRam * sizeof(Particle) << std::endl;
  }
}

void
Peridynamics::releaseGatheredPeriParticle()
{
  allPeriParticleVec.clear();
  PeriParticlePArray().swap(allPeriParticleVec); // actual memory release
}

void
Peridynamics::openPeriProgress(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    dem::debugInf << "stream error: openPeriProgress" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(dem::OPREC);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"KE\" \"P\" \"Mises\" \"Volume\" \"horizonSize\""
      << std::endl;
}

void
Peridynamics::printPeriProgress(std::ofstream& ofs, const int iframe) const
{
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  REAL sigma11, sigma12; // sigma13;
  REAL /*sigma21,*/ sigma22, sigma23;
  REAL sigma31, /*sigma32,*/ sigma33;
  for (const auto& pt : allPeriParticleVec) {
    //        if( (*pt)->getInitPosition().x()>
    //        0.5*(util::getParam<REAL>("Xmin")+util::getParam<REAL>("Xmax"))
    //        )
    //        continue;
    sigma11 = pt->getSigma11();
    sigma12 = pt->getSigma12();
    // sigma13=(*pt)->getSigma13();
    // sigma21=(*pt)->getSigma21();
    sigma22 = pt->getSigma22();
    sigma23 = pt->getSigma23();
    sigma31 = pt->getSigma31();
    // sigma32=(*pt)->getSigma32();
    sigma33 = pt->getSigma33();
    pressure = sigma11 + sigma22 + sigma33;
    vonMisesStress =
      sqrt(0.5 * ((sigma11 - sigma22) * (sigma11 - sigma22) +
                  (sigma22 - sigma33) * (sigma22 - sigma33) +
                  (sigma11 - sigma33) * (sigma11 - sigma33)) +
           3 * (sigma12 * sigma12 + sigma23 * sigma23 + sigma31 * sigma31));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::setw(20) << pt->getParticleVolume() << std::setw(20)
        << pt->getHorizonSize() << std::endl;
    ofs.flush();
  }
}

void
Peridynamics::printPeriProgressHalf(std::ofstream& ofs, const int iframe) const
{
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  REAL sigma11, sigma12; // sigma13;
  REAL /*sigma21,*/ sigma22, sigma23;
  REAL sigma31, /*sigma32,*/ sigma33;
  for (const auto& pt : allPeriParticleVec) {
    if (pt->getInitPosition().x() >
        0.5 * (util::getParam<REAL>("Xmin") + util::getParam<REAL>("Xmax")))
      continue;
    sigma11 = pt->getSigma11();
    sigma12 = pt->getSigma12();
    // sigma13=(*pt)->getSigma13();
    // sigma21=(*pt)->getSigma21();
    sigma22 = pt->getSigma22();
    sigma23 = pt->getSigma23();
    sigma31 = pt->getSigma31();
    // sigma32=(*pt)->getSigma32();
    sigma33 = pt->getSigma33();
    pressure = sigma11 + sigma22 + sigma33;
    vonMisesStress =
      sqrt(0.5 * ((sigma11 - sigma22) * (sigma11 - sigma22) +
                  (sigma22 - sigma33) * (sigma22 - sigma33) +
                  (sigma11 - sigma33) * (sigma11 - sigma33)) +
           3 * (sigma12 * sigma12 + sigma23 * sigma23 + sigma31 * sigma31));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::setw(20) << pt->getParticleVolume() << std::setw(20)
        << pt->getHorizonSize() << std::endl;
    ofs.flush();
  }
}

// Not used in the DEM-PD coupling code
void
Peridynamics::solvePurePeridynamics(const std::string& outputFile,
                                    int printInterval)
{
  // open the tecplot file for output
  std::ofstream ofs(outputFile);
  int iframe = 0;
  writeParticleTecplot(ofs,iframe);
  //std::cout << std::string(72, '-') << std::endl;
  //std::cout << "Start of the time loop " << std::endl;
  //std::cout << std::string(72, '-') << std::endl;
  std::ofstream datafile("uxyz.dat");
  datafile.setf(std::ios::scientific, std::ios::floatfield);
  datafile.precision(10);
  datafile << "VARIABLES = \"Time step\", \"UX\", \"UY\", \"UZ\"" << std::endl;

  int nsteps = 10000;    // is not true
  for (int istep = 1; istep <= nsteps; istep++) {

    runFirstHalfStep();
    prescribeEssentialBoundaryCondition(istep);
    checkBondParticleAlive();
    calcParticleStress();
    calcParticleAcceleration();
    ApplyExternalForce(istep);
    runSecondHalfStep();

    if ( istep % printInterval == 0) {
      //std::cout << "*** current time step is    " << istep << std::endl;
      iframe++;
      writeParticleTecplot(ofs,iframe);
      datafile << istep
               << std::setw(20) << periParticleVec[568]->getDisplacement().x()
               << std::setw(20) << periParticleVec[568]->getDisplacement().y()
               << std::setw(20) << periParticleVec[568]->getDisplacement().z()
               << std::endl;
    }
    if ( istep % 200 == 0) {
      writeDisplacementData("ux.dat","uy.dat","uz.dat");
    }

  } // time loop
  ofs.close();
  datafile.close();
  //std::cout << std::string(72, '-') << std::endl;
  //std::cout << "Simulation Finished !" << std::endl;
  //std::cout << std::string(72, '-') << std::endl;
  writeDisplacementData("ux.dat","uy.dat","uz.dat");
} // end solve()

// Not used in the DEM-PD coupling code
void
Peridynamics::writeDisplacementData(const std::string& outputFilex,
                                    const std::string& outputFiley,
                                    const std::string& outputFilez)
{
  // displacement along the x axis
  std::ofstream ofs(outputFilex);
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "VARIABLES = \"X\", \"UX\"" << std::endl;
  for (int index = 0; index < 5; index++) {
    int node = Uxindex[index];
    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().x()
        << std::setw(20) << periParticleVec[node]->getDisplacement().x()
        << std::endl;
  }
  ofs.flush();
  ofs.close();

  //displacement along the y axis
  ofs.open(outputFiley);
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "VARIABLES = \"Y\", \"UY\"" << std::endl;
  for (int index = 0; index < 5; index++) {
    int node = Uyindex[index];
    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().y()
        << std::setw(20) << periParticleVec[node]->getDisplacement().y()
        << std::endl;
  }
  ofs.flush();
  ofs.close();

  //displacement along the z axis
  ofs.open(outputFilez);
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "VARIABLES = \"Z\", \"UZ\"" << std::endl;
  for (int index = 0; index < 39; index++) {
    int node = Uzindex[index];
    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().z()
        << std::setw(20) << periParticleVec[node]->getDisplacement().z()
        << std::endl;
  }
  ofs.flush();
  ofs.close();
} // end writeDisplacementData

void
Peridynamics::runFirstHalfStep()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();

#pragma omp parallel for \
        num_threads(ompThreads) \
        private(i) \
        shared(num) \
        schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->updateDisplacement();
  }
} // end runFirstHalfStep()

void
Peridynamics::runSecondHalfStep()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();

#pragma omp parallel for \
        num_threads(ompThreads) \
        private(i) \
        shared(num) \
        schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->updateVelocity();
  }
} // end runSecondHalfStep()


// Not used in DEM-PD coupling code, i.e. rigidInclustion()
void
Peridynamics::writeMesh(const std::string& outputFile)
{
  std::ofstream ofs(outputFile);
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Mesh Checking\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\"" << std::endl;
  ofs << "ZONE N = " << nPeriParticle << " E = " << nele 
      << ", F = FEPOINT ET = BRICK" << std::endl;
  for (int node = 0; node < nPeriParticle; node++) {
    ofs << std::setw(20) << periParticleVec[node]->getInitPosition().x()
        << std::setw(20) << periParticleVec[node]->getInitPosition().y()
        << std::setw(20) << periParticleVec[node]->getInitPosition().z() 
        << std::endl;
  }
  for (int iel = 0; iel < nele; iel++) {
    for (int node = 0; node < 8; node++) {
      ofs << std::setw(10) << d_connectivity[iel][node];
    }
    ofs << std::endl;
  }
  ofs.close();
} // end writeMesh

// Not used in DEM-PD coupling code, i.e. rigidInclusion()
void
Peridynamics::writeMeshCheckVolume(const std::string& outputFile)
{
  std::ofstream ofs(outputFile);
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Volume Checking\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"V\"" << std::endl;
  ofs << "ZONE N = " << nPeriParticle << " E = " << nele 
      << ", F = FEPOINT ET = BRICK" << std::endl;
  // Output the coordinates and the array information
  for (const auto& particle : periParticleVecInitial) {
    ofs << std::setw(20) << particle->getInitPosition().x()
        << std::setw(20) << particle->getInitPosition().y()
        << std::setw(20) << particle->getInitPosition().z()
        << std::setw(20) << particle->getParticleVolume()
        << std::endl;
  }
  for (int iel = 0; iel < nele; iel++) {
    for (int node = 0; node < 8; node++) {
      ofs << std::setw(10) << d_connectivity[iel][node];
    }
    ofs << std::endl;
  }
  ofs.close();
} // end writeMeshCheckVolume

void
Peridynamics::writeParticleTecplot(std::ofstream& ofs, const int iframe) const
{
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  if (iframe == 0) {
    ofs << "Title = \"Particle Information\"" << std::endl;
    ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
           "\"Vz\" \"KE\" \"P\" \"Mises\""
        << std::endl;
  }
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  for (const auto& pt : allPeriParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::endl;
    ofs.flush();
  }

} // end writeparticleTecplot

// this function has some problems, d_allPeriParticlesInitial does not exist
// anymore after the first gather
void
Peridynamics::printPeriDomain(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" "
         "\"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\""
      << std::endl;
  //    ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" <<
  //    nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" <<
  //    std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  Matrix Kinv_tmp;
  Matrix deformationG;
  for (const auto& pt : periParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    Kinv_tmp = pt->getParticleKinv();
    deformationG = pt->getDeformationGradient();
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << pt->getAcceleration().x()
        << std::setw(20) << pt->getAcceleration().y() << std::setw(20)
        << pt->getAcceleration().z() << std::setw(20)
        << vfabs(pt->getVelocity()) << std::setw(20) << pressure
        << std::setw(20) << vonMisesStress << std::setw(20)
        << pt->getParticleVolume() << std::setw(20) << pt->getHorizonSize()
        << std::setw(20) << pt->getBondsNumber() << std::setw(20)
        << dem::det(Kinv_tmp) << std::setw(20) << dem::det(deformationG)
        << std::endl;
    ofs.flush();
  }
  //    for(int iel = 0; iel < nele; iel++){
  //        for(int node = 0; node < 8; node++){
  //        ofs << std::setw(10) << d_connectivity[iel][node];
  //        }
  //        ofs << std::endl;
  //    }

  ofs.close();

} // end printPeriDomain

// this function has some problems, d_allPeriParticlesInitial does not exist
// anymore after the first gather
void
Peridynamics::printRecvPeriDomain(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" "
         "\"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\""
      << std::endl;
  //    ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" <<
  //    nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" <<
  //    std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  Matrix Kinv_tmp;
  Matrix deformationG;
  for (const auto& pt : recvPeriParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    Kinv_tmp = pt->getParticleKinv();
    deformationG = pt->getDeformationGradient();
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << pt->getAcceleration().x()
        << std::setw(20) << pt->getAcceleration().y() << std::setw(20)
        << pt->getAcceleration().z() << std::setw(20)
        << vfabs(pt->getVelocity()) << std::setw(20) << pressure
        << std::setw(20) << vonMisesStress << std::setw(20)
        << pt->getParticleVolume() << std::setw(20) << pt->getHorizonSize()
        << std::setw(20) << pt->getBondsNumber() << std::setw(20)
        << dem::det(Kinv_tmp) << std::setw(20) << dem::det(deformationG)
        << std::endl;
    ofs.flush();
  }
  //    for(int iel = 0; iel < nele; iel++){
  //        for(int node = 0; node < 8; node++){
  //        ofs << std::setw(10) << d_connectivity[iel][node];
  //        }
  //        ofs << std::endl;
  //    }

  ofs.close();

} // end printRecvPeriDomain

void
Peridynamics::printPeriParticle(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(dem::OPREC);

  ofs << std::setw(dem::OWID) << "Node ID" << std::setw(dem::OWID) << "U1"
      << std::setw(dem::OWID) << "U2" << std::setw(dem::OWID) << "U3" << std::setw(dem::OWID)
      << "S.Mises" << std::setw(dem::OWID) << "S11" << std::setw(dem::OWID) << "S22"
      << std::setw(dem::OWID) << "S33" << std::setw(dem::OWID) << "S12" << std::setw(dem::OWID)
      << "S13" << std::setw(dem::OWID) << "S23" << std::setw(dem::OWID) << "X"
      << std::setw(dem::OWID) << "Y" << std::setw(dem::OWID) << "Z" << std::endl;

  int id = 0;
  Matrix sigma;
  REAL vonMisesStress;
  for (const auto& pt : allPeriParticleVec) {
    sigma = pt->getSigma();
    id++;
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    ofs << std::setw(dem::OWID) << id << std::setw(dem::OWID) << pt->getDisplacement().x()
        << std::setw(dem::OWID) << pt->getDisplacement().y() << std::setw(dem::OWID)
        << pt->getDisplacement().z() << std::setw(dem::OWID) << vonMisesStress
        << std::setw(dem::OWID) << sigma(1, 1) << std::setw(dem::OWID) << sigma(2, 2)
        << std::setw(dem::OWID) << sigma(3, 3) << std::setw(dem::OWID) << sigma(1, 2)
        << std::setw(dem::OWID) << sigma(1, 3) << std::setw(dem::OWID) << sigma(2, 3)
        << std::setw(dem::OWID) << pt->getInitPosition().x() << std::setw(dem::OWID)
        << pt->getInitPosition().y() << std::setw(dem::OWID)
        << pt->getInitPosition().z() << std::endl;
  }

  ofs.close();
}

void
Peridynamics::printPeriDomainSphere(const std::string& str) const
{
  /* // not used in DEM-PD coupling code, i.e. rigidInclusion()
        std::ofstream ofs(str);
        if(!ofs) {
              //std::cout << "stream error!" << endl; exit(-1);
        }

      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "Title = \"Particle Information\"" << std::endl;
      ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Sigma_rr\" \"Sigma_tt\"
     \"Tao_rt\" " << std::endl;
      ofs << "ZONE T = \"periDomain\" "<< std::endl;
      // Output the coordinates and the array information
      REAL pressure, vonMisesStress;
      Matrix sigma;
      // this is not to print all the peri-points, it just print the peri-points
     in the interface between the sphere and the box
      for(PeriParticlePArray::const_iterator pt =
     interfacePeriParticleVec.begin(); pt!= interfacePeriParticleVec.end();
     pt++) {
          sigma = (*pt)->getSigma();
          Vec currPosition = (*pt)->currentPosition();
          REAL R = vfabs(currPosition);
          REAL theta = acos(currPosition.z()/R);

          // for phi should notice that x can be zero
          REAL phi;
          if( abs(currPosition.x())<dem::EPS ){
          if( currPosition.y()>0 ){
              phi = dem::Pi*0.5;
          }
          else{
              phi = dem::Pi*1.5;
          }
          }
          else {
              phi = atan(currPosition.y()/currPosition.x());
          }

          Matrix Qtrans(3,3);    // the transfor tensor from cartesian to
     spherical coordinate
                  // refer to
     http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
          Qtrans(1,1) = sin(theta)*cos(phi); Qtrans(1,2) = sin(theta)*sin(phi);
     Qtrans(1,3) = cos(theta);
          Qtrans(2,1) = cos(theta)*cos(phi); Qtrans(2,2) = cos(theta)*sin(phi);
     Qtrans(2,3) = -sin(theta);
          Qtrans(3,1) = -sin(phi);           Qtrans(3,2) = cos(phi);
     Qtrans(3,3) = 0;

          Matrix sigma_sphere = Qtrans*sigma*trans(Qtrans);

          ofs << std::setw(20) << currPosition.x()
          << std::setw(20) << currPosition.y()
          << std::setw(20) << currPosition.z()
          << std::setw(20) << sigma_sphere(1,1)
          << std::setw(20) << sigma_sphere(2,2)
          << std::setw(20) << sigma_sphere(1,2)
          << std::endl;
          ofs.flush();
      }

      ofs.close();
  */
} // end printPeriDomainSphere

void
Peridynamics::calcDeformationGradient()
{
  // calcDeformationGradient - calculates the deformation gradient for all
  // peri-particles
}

void
Peridynamics::setInitIsv()
{
  REAL isv_tmp;
  if (util::getParam<int>("typeConstitutive") ==
      1) { // 1---implicit, 2---explicit
    isv_tmp = util::getParam<REAL>("Chi");
  } else {
    isv_tmp = util::getParam<REAL>("c");
  }

  for (auto& pt : d_allPeriParticlesInitial) {
    pt->setInitIsv(isv_tmp);
  }

} // setInitIsv()

/*    void Peridynamics::checkBondParticleAlive(){
    // compute the bond length and check for whether a bond or particle is alive
   or not
    for(PeriParticlePArray::iterator pt=periParticleVec.begin();
   pt!=periParticleVec.end(); pt++){
        if( (*pt)->getIsAlive() ){ // particle not alive, then go to next
   particle
        (*pt)->checkParticleAlive(util::getParam<REAL>("bondStretchLimit"));
        }

    } // end particle

    } // end checkBondParticleAlive()
*/

void
Peridynamics::checkBondParticleAlive()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int pt;
  num = periBondVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num)       \
  schedule(dynamic)
  for (pt = 0; pt < num; pt++) {
    periBondVec[pt]->checkIfAlive();
  }

  // compute the bond length and check for whether a bond or particle is alive
  // or not
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num)       \
  schedule(dynamic)
  for (pt = 0; pt < num; pt++) {
    if (periParticleVec[pt]
          ->getIsAlive()) { // particle not alive, then go to next particle
      periParticleVec[pt]->checkParticleAlive();
    }
  } // end particle

} // end checkBondParticleAlive()

void
Peridynamics::calcRecvParticleKinv()
{ // this should be called after
  // commuDEMPeriParticle(), and find
  // peri-bonds between communicated
  // periParticles
  // Compute the inverse of the shape tensor K
  for (auto& pt : recvPeriParticleVec) {
    pt->calcParticleKinv();
  }

} // end calcParticleKinv()

void
Peridynamics::ApplyExternalForce(int istep)
{
  // deal with the external force, applied at the top of the boundary
  REAL factor = 0.0;
  REAL rampStep = util::getParam<REAL>("rampStep");
  if (istep <= rampStep) {
    factor = REAL(istep) / REAL(rampStep);
  } else {
    factor = REAL(1.0);
  }
  //    for(PeriParticlePArray::iterator pt=topBoundaryVec.begin();
  //    pt!=topBoundaryVec.end(); pt++){
  //        (*pt)->setAcceleration(Vec(0.0,0.0,factor*util::getParam<REAL>("bodyDensity")));
  //    }
} // end ApplyExternalForce

void
Peridynamics::clearPeriDEMBonds()
{
  // for (auto bt = periDEMBondVec.begin(); bt != periDEMBondVec.end(); bt++) {
  //  delete (*bt);
  //}
  periDEMBondVec.clear();
  PeriDEMBondPArray().swap(periDEMBondVec); // actual memory release
  dem::plan_gravity = util::getParam<REAL>("periDensity") * d_maxDistBetweenParticles *
                 d_maxDistBetweenParticles * 9.8; // rho*l^2*g, used in PeriDEMBond.cpp
}

void
Peridynamics::eraseBrokenPeriDEMBonds()
{
  // BB: Feb 3, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  periDEMBondVec.erase(std::remove_if(periDEMBondVec.begin(),
                                      periDEMBondVec.end(),
                                      [](const PeriDEMBondP& bond) {
                                        if (bond->getIsAlive()) {
                                          return false;
                                        }
                                        return true;
                                      }),
                       periDEMBondVec.end());

  /*
  for (PeriDEMBondPArray::iterator bt = periDEMBondVec.begin();
       bt != periDEMBondVec.end();) {
    if ((*bt)->getIsAlive()) {  // still alive
      bt++;
    } else {
      delete (*bt);
      bt = periDEMBondVec.erase(bt);  // release memory
    }
  }
  */
}

// if a peri-point is first time (i.e. previous step peri-point is not bonded to
// dem particle) entering the influence domain of a DEM particle,
// then construct peri-DEM bond between this peri-point and this DEM particle
void
Peridynamics::findPeriDEMBonds(dem::ParticlePArray mergeParticleVec)
{
  eraseBrokenPeriDEMBonds();
  // construct sand peri-bonds
  // here mergePeriParticleVec is used, since if we have a DEM particle that is
  // crossing the boundary between two cpus, then sand-peri bonds between
  // the communicated peri-points can provide force on the DEM particle from the
  // peri-points in neighboring cpus
  // mergeParticleVec is used, since if a DEM particle is crossing the boundary
  // between two cpus, then these sand-peri bonds between the communicated
  // DEM particles can provide force on peri-points from the DEM particles in
  // neighboring cpus
  REAL delta = d_maxDistBetweenParticles * 3.2; // 3 layers of peri-points

  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int peri_pt;
  num = mergePeriParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt)              \
  shared(num, delta) schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    Vec xyz_peri = mergePeriParticleVec[peri_pt]->currentPosition();

    for (auto& dem_pt : mergeParticleVec) {
      // check and construct the periDEMBondVec in this particle
      REAL ra = dem_pt->getA();
      REAL rb = dem_pt->getB();
      REAL rc = dem_pt->getC();
      Vec xyz_peri_tmp = dem_pt->globalToLocal(
        xyz_peri - dem_pt->currentPosition()); // this is very important, since
                                          // all calculations below for
                                          // ellipsoid
      REAL x_peri =
        xyz_peri_tmp.x(); // are based on the local coordinate of the ellipsoid
      REAL y_peri = xyz_peri_tmp.y();
      REAL z_peri = xyz_peri_tmp.z();
      REAL x_peri_de =
        x_peri /
        (ra + delta); // are based on the local coordinate of the ellipsoid
      REAL y_peri_de = y_peri / (rb + delta);
      REAL z_peri_de = z_peri / (rc + delta);
      if (x_peri_de * x_peri_de + y_peri_de * y_peri_de +
            z_peri_de * z_peri_de <
          1) {
        // check if (*peri_pt) is first time entering the influence domain of
        // (*dem_pt)
        if (mergePeriParticleVec[peri_pt]->isBonded(dem_pt->getId()))
          continue; // already bonded

        // calcualte the local coordinates of the projector on the surface of
        // the ellipsoid particle
        Vec projector_local;

        // get the projector of this peri-point, *peri_pt, on the surface of
        // this ellipsoid, *dem_pt

        // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in the
        // same octant
        REAL kd = sqrt(ra * ra * rb * rb * rc * rc /
                       (x_peri * x_peri * rb * rb * rc * rc +
                        y_peri * y_peri * ra * ra * rc * rc +
                        z_peri * z_peri * ra * ra * rb * rb));
        REAL dd =
          sqrt(x_peri * x_peri + y_peri * y_peri + z_peri * z_peri) * (1 - kd);

        REAL n1 = 2 * x_peri / (ra + dd) / (ra + dd);
        REAL n2 = 2 * y_peri / (rb + dd) / (rb + dd);
        REAL n3 = 2 * z_peri / (rc + dd) / (rc + dd);

        REAL aa = n1 * n1 * rb * rb * rc * rc + n2 * n2 * ra * ra * rc * rc +
                  n3 * n3 * ra * ra * rb * rb;
        REAL bb = 2 * (n1 * x_peri * rb * rb * rc * rc +
                       n2 * y_peri * ra * ra * rc * rc +
                       n3 * z_peri * ra * ra * rb * rb);
        REAL cc = rb * rb * rc * rc * x_peri * x_peri +
                  ra * ra * rc * rc * y_peri * y_peri +
                  ra * ra * rb * rb * z_peri * z_peri -
                  ra * ra * rb * rb * rc * rc;

        REAL kn = (-bb + sqrt(bb * bb - 4 * aa * cc)) * 0.5 / aa;

        // //std::cout << "kn: " << kn << std::endl;

        projector_local.setX(x_peri + kn * n1);
        projector_local.setY(y_peri + kn * n2);
        projector_local.setZ(z_peri + kn * n3);

        //            // this is to use the cross point instead of the normal
        //            point as the projector point
        //            projector_local.setX(x_peri*kd);
        //            projector_local.setY(y_peri*kd);
        //            projector_local.setZ(z_peri*kd);

        PeriDEMBondP bond_tmp = std::make_shared<PeriDEMBond>(
          projector_local, dem_pt.get(), mergePeriParticleVec[peri_pt].get());

//            // this is used to test the coupled force model, October 10, 2014
//            // in this test model, the sand-peri-points will move along the
//            dem-particle
//            PeriDEMBondP bond_tmp =
//            std::make_shared<PeriDEMBond>(projector_local,
//               (*dem_pt).get(), (*peri_pt).get());

#pragma omp critical
        periDEMBondVec.push_back(bond_tmp);
      }
    }
  }

} // findPeriDEMBonds()

// construct boundary bonds, sand bonds, July 15, 2014
void
Peridynamics::constructBoundarySandPeriBonds()
{
  /*
      // delete previous bonds vector
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=topBoundaryBondVec.begin(); pt!=topBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=bottomBoundaryBondVec.begin(); pt!=bottomBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=leftBoundaryBondVec.begin(); pt!=leftBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=rightBoundaryBondVec.begin(); pt!=rightBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=frontBoundaryBondVec.begin(); pt!=frontBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=backBoundaryBondVec.begin(); pt!=backBoundaryBondVec.end(); pt++){
          delete (*pt);
          }

          topBoundaryBondVec.clear();
          bottomBoundaryBondVec.clear();
          leftBoundaryBondVec.clear();
          rightBoundaryBondVec.clear();
          frontBoundaryBondVec.clear();
          backBoundaryBondVec.clear();

          for(std::vector<dem::PeriDEMBond*>::iterator
  pt=periDEMBondVec.begin(); pt!=periDEMBondVec.end(); pt++){
          delete (*pt);
          }
          periDEMBondVec.clear();


      REAL delta = d_maxDistBetweenParticles*3;    // 3 layers of peri-points

      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();

  // no boundaries in pile penetration
  //    // construct boundary peri-bonds
  //    PeriBoundaryBond* bond_tmp;
  //    for(PeriParticlePArray::iterator pt=periParticleVec.begin();
  pt!=periParticleVec.end(); pt++){ // overlap over peri-points
  //        Vec xyz_pt = (*pt)->currentPosition();
  //        REAL x_pt = xyz_pt.x();
  //        REAL y_pt = xyz_pt.y();
  //        REAL z_pt = xyz_pt.z();
  //
  //        if( x_pt-x_min<delta ){    // back boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_min,y_pt,z_pt), (*pt));
  //        backBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( x_max-x_pt<delta ){    // front boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_max,y_pt,z_pt), (*pt));
  //        frontBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( y_pt-y_min<delta ){    // left boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_min,z_pt), (*pt));
  //        leftBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( y_max-y_pt<delta ){    // right boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_max,z_pt), (*pt));
  //        rightBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( z_pt-z_min<delta ){    // bottom boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_min), (*pt));
  //        bottomBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( z_max-z_pt<delta ){    // top boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_max), (*pt));
  //        topBoundaryBondVec.push_back(bond_tmp);
  //        }
  //
  //    }
  //

      // construct sand peri-bonds
      for(PeriParticlePArray::iterator peri_pt=interfacePeriParticleVec.begin();
  peri_pt!=interfacePeriParticleVec.end(); peri_pt++){
          Vec xyz_peri = (*peri_pt)->currentPosition();

          for(ParticlePArray::iterator dem_pt=ParticleVec.begin();
  dem_pt!=ParticleVec.end(); dem_pt++){
          // check and construct the periDEMBondVec in this particle
          REAL ra = (*dem_pt)->getA();
          REAL rb = (*dem_pt)->getB();
          REAL rc = (*dem_pt)->getC();
              xyz_peri = (*dem_pt)->localVec(xyz_peri-(*dem_pt)->currentPosition());
  // this is very important, since all calculations below for ellipsoid
              REAL x_peri = xyz_peri.x();                        // are based
  on the local coordinate of the ellipsoid
              REAL y_peri = xyz_peri.y();
              REAL z_peri = xyz_peri.z();
  //        if(
  (x_peri/(ra+delta))*(x_peri/(ra+delta))+(y_peri/(rb+delta))*(y_peri/(rb+delta))+(z_peri/(rc+delta))*(z_peri/(rc+delta))
  < 1 ){
                  // calcualte the local coordinates of the projector on the
  surface of the ellipsoid particle
                  Vec projector_local;

                  // get the projector of this peri-point, *peri_pt, on the
  surface of this ellipsoid, *dem_pt

              // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in
  the same octant
              REAL kd =
  sqrt(ra*ra*rb*rb*rc*rc/(x_peri*x_peri*rb*rb*rc*rc+y_peri*y_peri*ra*ra*rc*rc+z_peri*z_peri*ra*ra*rb*rb));
  //            REAL dd =
  sqrt(x_peri*x_peri+y_peri*y_peri+z_peri*z_peri)*(1-kd);
  //
  //            REAL n1 = 2*x_peri/(ra+dd)/(ra+dd);
  //            REAL n2 = 2*y_peri/(rb+dd)/(rb+dd);
  //            REAL n3 = 2*z_peri/(rc+dd)/(rc+dd);
  //
  //            REAL aa = n1*n1*rb*rb*rc*rc+n2*n2*ra*ra*rc*rc+n3*n3*ra*ra*rb*rb;
  //            REAL bb =
  2*(n1*x_peri*rb*rb*rc*rc+n2*y_peri*ra*ra*rc*rc+n3*z_peri*ra*ra*rb*rb);
  //            REAL cc =
  rb*rb*rc*rc*x_peri*x_peri+ra*ra*rc*rc*y_peri*y_peri+ra*ra*rb*rb*z_peri*z_peri-ra*ra*rb*rb*rc*rc;
  //
  //            REAL kn = (-bb+sqrt(bb*bb-4*aa*cc))*0.5/aa;
  //
  ////std::cout << "kn: " << kn << std::endl;
  //
  //
  //            projector_local.setX(x_peri+kn*n1);
  //            projector_local.setY(y_peri+kn*n2);
  //            projector_local.setZ(z_peri+kn*n3);
  //
  //
              // this is to use the cross point instead of the normal point as
  the projector point
              projector_local.setX(x_peri*kd);
              projector_local.setY(y_peri*kd);
              projector_local.setZ(z_peri*kd);

              PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, *dem_pt,
  *peri_pt);

  //            // this is used to test the coupled force model, October 10,
  2014
  //            // in this test model, the sand-peri-points will move along the
  dem-particle
  //            PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local,
  *dem_pt, *peri_pt);

              periDEMBondVec.push_back(bond_tmp);

  //        }
          }
      }

  */
} // end constructBoundarySandPeriBonds

void
Peridynamics::applyCoupledForces()
{
  /* no boundaries in the pile penetration
      // calculate the coupled forces between boundaries
      // calculate current length of the bonds
      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();
      for(PeriBoundaryBondPArray::iterator
     bond_pt=bottomBoundaryBondVec.begin();
     bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_min, 6);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=topBoundaryBondVec.begin();
     bond_pt!=topBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_max, 5);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=backBoundaryBondVec.begin();
     bond_pt!=backBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_min, 3);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=frontBoundaryBondVec.begin();
     bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_max, 1);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=leftBoundaryBondVec.begin();
     bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_min, 4);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=rightBoundaryBondVec.begin();
     bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_max, 2);
          // pass the boundary coordinate and the number of the boundary
      }
  */
  // calculate coupled forces between sand particles and peri-points
  for (auto& bond_pt : periDEMBondVec) {
    bond_pt->applyBondForce();
  }
} // end applyCoupledForces()

// this is used to test the coupling model. In this test model, the
// sand-peri-points will move along the
// dem-particles. And the peri-points have no influence on the dem-particles.
// October 10, 2014
void
Peridynamics::applyCoupledBoundary()
{
  /* no boundaries in the pile penetration
      // calculate the coupled forces between boundaries
      // calculate current length of the bonds
      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();
      for(PeriBoundaryBondPArray::iterator
     bond_pt=bottomBoundaryBondVec.begin();
     bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_min, 6);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=topBoundaryBondVec.begin();
     bond_pt!=topBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_max, 5);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=backBoundaryBondVec.begin();
     bond_pt!=backBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_min, 3);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=frontBoundaryBondVec.begin();
     bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_max, 1);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=leftBoundaryBondVec.begin();
     bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_min, 4);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=rightBoundaryBondVec.begin();
     bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_max, 2);
          // pass the boundary coordinate and the number of the boundary
      }
  */

  // calculate coupled forces between sand particles and peri-points
  for (auto& bond_pt : periDEMBondVec) {
    bond_pt->applyBondBoundary();
  }

} // end applyCoupledBoundary()

void
Peridynamics::constructRecvPeriMatrix()
{ // construct the Matrix members in
  // periParticleVec,
  // since currently the transfer of pointer array in Matrix is not implemented
  // well
  for (auto& pt : recvPeriParticleVec) {
    pt->constructMatrixMember(); // construct these Matrix members
  }

} // constructRecvPeriMatrix()

void
Peridynamics::findBoundaryPeriParticles()
{ // it does not matter since this
  // will be only called once
  bottomBoundaryVec.clear();
  frontBoundaryVec.clear();
  leftBoundaryVec.clear();
  topBoundaryInnerVec.clear();
  topBoundaryEdgeVec.clear();
  topBoundaryCornerVec.clear();

  REAL x1 = util::getParam<REAL>("Xmin");
  REAL x2 = util::getParam<REAL>("Xmax");
  REAL y1 = util::getParam<REAL>("Ymin");
  REAL y2 = util::getParam<REAL>("Ymax");
  REAL z1 = util::getParam<REAL>("Zmin");
  REAL z2 = util::getParam<REAL>("Zmax");
  for (auto& peri_pt : periParticleVec) {
    REAL x_peri = peri_pt->getInitPosition().x();
    REAL y_peri = peri_pt->getInitPosition().y();
    REAL z_peri = peri_pt->getInitPosition().z();
    if (z_peri == z2) { // top boundary
      //        if( (x_peri==x2&&y_peri==y2) || x_peri==x1&&y_peri==y2 ||
      //        x_peri==x2&&y_peri==y1 || x_peri==x1&&y_peri==y1){    // corner
      //        points
      //        topBoundaryCornerVec.push_back(*peri_pt);
      //        }
      //        else if(x_peri==x1 || y_peri==y1 || x_peri==x2 || y_peri==y2){
      //        // edge points
      //        topBoundaryEdgeVec.push_back(*peri_pt);
      //        }
      //        else {    // inner points
      //            topBoundaryInnerVec.push_back(*peri_pt);
      //        }
    } else if (z_peri < z1 + 4.5 * d_maxDistBetweenParticles) { // bottom boundary
      bottomBoundaryVec.push_back(peri_pt);
    }

    if (x_peri < x1 + 4.5 * d_maxDistBetweenParticles ||
        x_peri > x2 - 4.5 * d_maxDistBetweenParticles) {
      leftBoundaryVec.push_back(peri_pt);
    }
    if (y_peri < y1 + 4.5 * d_maxDistBetweenParticles ||
        y_peri > y2 - 4.5 * d_maxDistBetweenParticles) {
      frontBoundaryVec.push_back(peri_pt);
    }
  }

} // findBoundaryPeriParticles()

void
Peridynamics::findFixedPeriParticles()
{ // it does not matter since this will
  // be only called once
  fixedPeriParticleVec.clear();
  REAL radius = util::getParam<REAL>("fixRadius");
  radius = radius +
           d_maxDistBetweenParticles *
             0.2; // in order to find all points on the spherical surface
  REAL x0 = util::getParam<REAL>("periFixCentroidX");
  REAL y0 = util::getParam<REAL>("periFixCentroidY");
  REAL z0 = util::getParam<REAL>("periFixCentroidZ");
  for (auto& peri_pt : periParticleVec) {
    REAL x_peri = peri_pt->getInitPosition().x();
    REAL y_peri = peri_pt->getInitPosition().y();
    REAL z_peri = peri_pt->getInitPosition().z();
    if ((x_peri - x0) * (x_peri - x0) + (y_peri - y0) * (y_peri - y0) +
          (z_peri - z0) * (z_peri - z0) <
        radius * radius + dem::EPS) {
      fixedPeriParticleVec.push_back(peri_pt);
    }
  }

} // findFixedPeriParticles()

void
Peridynamics::applyPeriBoundaryCondition()
{
  dem::Vec zero_vec = dem::Vec(0, 0, 0);
  //    for(PeriParticlePArray::iterator peri_pt=fixedPeriParticleVec.begin();
  //    peri_pt!=fixedPeriParticleVec.end(); peri_pt++){
  //        (*peri_pt)->prescribeDisplacement(zero_vec);
  //    }

  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int peri_pt;
  num = bottomBoundaryVec.size();
  //#pragma omp parallel for num_threads(ompThreads) private(peri_pt)
  // shared(num) schedule(dynamic)
  //    for(peri_pt=0; peri_pt<num; peri_pt++){
  //        bottomBoundaryVec[peri_pt]->prescribeBottomDisplacement(0);
  //    }
  num = leftBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num)  \
  schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    leftBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
  }
  num = frontBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num)  \
  schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    frontBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
  }

} // applyPeriBoundaryCondition()

