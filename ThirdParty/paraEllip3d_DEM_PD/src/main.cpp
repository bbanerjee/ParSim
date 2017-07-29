///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: paraEllip3d //
//                                 Author: Beichuan Yan //
//                                  Email: beichuan.yan@colorado.edu //
//                              Institute: University of Colorado at Boulder //
///////////////////////////////////////////////////////////////////////////////////////////////////////
#include <Core/Const/const.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DiscreteElements.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/Parameter.h>
#include <Simulations/Command.h>
#include <Simulations/CommandHandler.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <ctime>
#include <iostream>
#include <vector>
//#include <boost/timer/timer.hpp>

// serialization of pointers to objects of derived classes
// http://www.boost.org/doc/libs/1_50_0/libs/serialization/doc/serialization.html#registration
// this is called registration, and this method of registration is referred to
// as "key export".
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(dem::PlaneBoundary, "PlaneBoundary")

// optimization on non-template types which have a fixed amount of data stored
// at fixed field positions
BOOST_IS_MPI_DATATYPE(dem::Vec)
BOOST_IS_MPI_DATATYPE(dem::Particle)
BOOST_IS_MPI_DATATYPE(dem::Contact)

using dem::Parameter;

int
main(int argc, char* argv[])
{

  boost::mpi::environment boostEnv(argc, argv);
  boost::mpi::communicator boostWorld;
  // boost::timer::auto_cpu_timer boostTimer;
  REAL time0 = MPI_Wtime();

  if (boostWorld.rank() == 0) {
    if (argc != 2) {
      std::cerr << "Usage: paraEllip3d input.txt" << std::endl;
      return -1;
    }

    Parameter::get().readIn(argv[1]);
    Parameter::get().writeOut();
    int mpiProcX = util::getParam<int>("mpiProcX");;
    int mpiProcY = util::getParam<int>("mpiProcY");;
    int mpiProcZ = util::getParam<int>("mpiProcZ");;
    if (mpiProcX * mpiProcY * mpiProcZ != boostWorld.size()) {
      std::cerr << "Number of MPI processes does not match grid in data file!"
                << std::endl;
      return -1;
    }
  }
  broadcast(boostWorld, Parameter::get(), 0); // broadcast from root process 0

  dem::DiscreteElements dem;
  pd::Peridynamics pd;
  dem.setCommunicator(boostWorld);
  pd.setCommunicator(boostWorld, 
                     dem::DiscreteElements::getMPIWorldComm(), 
                     dem::DiscreteElements::getMPICartComm(), 
                     dem::DiscreteElements::getMPIRank(),
                     dem::DiscreteElements::getMPISize(), 0, 
                     dem::DiscreteElements::getMPIProcs(), 
                     dem::DiscreteElements::getMPICoords());

  // only root process prints to debugInf
  if (boostWorld.rank() == 0) {
    dem::debugInf.open("debugInf");
    if (!dem::debugInf) {
      std::cerr << "stream error: main.cpp debugInf" << std::endl;
      exit(-1);
    }
    dem::debugInf.setf(std::ios::scientific, std::ios::floatfield);
  }

  // parallel IO
  MPI_File_open(MPI_Comm(boostWorld), "overlapInf",
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                &dem::overlapInf);
  if (boostWorld.rank() == 0 && !dem::overlapInf) {
    std::cerr << "stream error: main.cpp overlapInf" << std::endl;
    exit(-1);
  }

  int simuType = util::getParam<int>("simuType");

  // Run the commands for each simulation type
  dem::CommandHandler handler;
  dem::CommandP command = handler.handleCommand(simuType);
  if (simuType > 3000) {
    command->execute(&dem, &pd);
  } else {
    command->execute(&dem);
  }

  dem::debugInf << std::endl
                << "MPI_Wtime: " << MPI_Wtime() - time0 << " seconds"
                << std::endl;
  dem::debugInf.close();
  MPI_File_close(&dem::overlapInf);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// notes:
//
// particleLayers (settings for free particles):
//  0 - a single free particle
//  1 - a layer of free particles
//  2 - multiple layers of free particles
//
// particle type (settings for individual particles):
//  0 - free particle
//  1 - fixed particle
//  2 - special case 2 (pure moment): translate first, then rotate only,
//  MNT_START needs to be defined
//  3 - special case 3 (displacemental ellipsoidal pile): translate in vertical
//  direction only
//  4 - special case 4 (impacting ellipsoidal penetrator): impact with inital
//  velocity in vertical direction only
//  5 - free boundary particle
// 10 - ghost particle
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
// degravitation, no boundary, quasi-static
assemb.deGravitation(1000,               // total_steps
                     5,                  // number of snapshots
                     1,                  // print interval
                     true,               // recreate from input file
                     "dep_particle_end", // input file, initial particles
                     "dgr_particle",     // output file, resulted particles,
including snapshots
                     "dgr_contact",      // output file, resulted contacts,
including snapshots
                     "dgr_progress",     // output file, statistical info
                     "dgr_debug");       // output file, debug info
// container properties
dem::Box container(0.05, 0.05, 0.05, vec(0, 0, 0)); // dimx, dimy, dimz, center
// particle shape, size and percentage
REAL ptcl_ratio_ba = 1;//0.8;  // ratio of radius b to radius a
REAL ptcl_ratio_ca = 1;//0.6;  // ratio of radius c to radius a
std::vector<REAL> percent; // mass percentage of particles smaller than a
certain size
std::vector<REAL> ptclSize;    // particle size
percent.push_back(1.00); ptclSize.push_back(2.5e-3);
//percent.push_back(0.80); ptclSize.push_back(2.3e-3);
//percent.push_back(0.60); ptclSize.push_back(2.0e-3);
//percent.push_back(0.30); ptclSize.push_back(1.5e-3);
//percent.push_back(0.10); ptclSize.push_back(1.0e-3);
dem::Gradation ptclGradation(percent.size(), percent, ptclSize, ptcl_ratio_ba,
ptcl_ratio_ca);
assemb.setContainer(container);
assemb.setGradation(ptclGradation);
///////////////////////////////////////////////////////////////////////////////////
// trim(), createMemParticle() and iso_MemBdry() must be called together
because:
// 1. trim() and createMemParticle() share variable HistoryNum.
// 2. springs in iso_MemBdry() reference particles created in
createMemParticle().
assemb.trim(false,               // recreate from input file or not
            "dgr_particle_end",  // input
            "trm_particle_end"); // output
assemb.createMemParticle(0.25,   // size relative to minimum radius
                         false,  // recreate from input file or not
                         "trm_particle_end",
                         "mem_particle_end");
assemb.iso_MemBdry(100000,
                   100,
                   10,
                   1.0e+3,
                   0.25,        // size relative to minimum radius
                   false,       // recreate from input file or not, must be
false
                   "mem_particle_end",
                   "iso_particle",
                   "iso_contact",
                   "iso_progress",
                   "iso_debug");
*/

/*
// record run time
time(&time2);
dem::g_timeinf << "simulation start time: " << ctime(&time1);
dem::g_timeinf << "simulation  end  time: " << ctime(&time2);
dem::g_timeinf.close();
*/

/*
// for triaxalPtclBdry
assemb.TrimPtclBdryByHeight(0.061,
                            "dep_particle_end",
                            "dep_particle_trimmed");


// for de-fe coupling
assemb.TrimPtclBdryByHeight(0.061,
                            "pile_particle_ini",
                            "pile_particle_trimmed");


assemb.triaxialPtclBdryIni(100000,             // total_steps
                           100,                // number of snapshots
                           10,                 // print interval
                           2.5e+6,             // confining pressure to achieve
                           "ini_particle_500k",// input file, initial particles
                           "ini_boundary_500k",// input file, initial boundaries
                           "ini_particle",     // output file, resulted
particles, including snapshots
                           "ini_boundary",     // output file, resulted
boundaries
                           "ini_contact",      // output file, resulted
contacts, including snapshots
                           "ini_progress",     // output file, statistical info
                           "ini_debug");       // output file, debug info


assemb.triaxialPtclBdryIni(100000,             // total_steps
                           100,                // number of snapshots
                           10,                 // print interval
                           5.0e+5,             // confining pressure to achieve
                           "ini_particle_ini", // input file, initial particles
                           "ini_boundary_ini", // input file, initial boundaries
                           "ini_particle",     // output file, resulted
particles, including snapshots
                           "ini_boundary",     // output file, resulted
boundaries
                           "ini_contact",      // output file, resulted
contacts, including snapshots
                           "ini_progress",     // output file, statistical info
                           "ini_debug");       // output file, debug info


assemb.triaxialPtclBdry(100000,             // total_steps
                        100,                // number of snapshots
                        10,                 // print interval
                        "tri_particle_ini", // input file, initial particles
                        "tri_boundary_ini", // input file, initial boundaries
                        "tri_particle",     // output file, resulted particles,
including snapshots
                        "tri_boundary",     // output file, resulted boundaries
                        "tri_contact",      // output file, resulted contacts,
including snapshots
                        "tri_progress",     // output file, statistical info
                        "tri_balanced",     // output file, balanced status
                        "tri_debug");       // output file, debug info


// size, shape, and Gradation of particles
int rorc      = 1;     // box = 1 or cylindrical = 0
REAL dimn     = 0.05;  // specimen dimension
REAL ratio_ba = 0.8;   // ratio of radius b to radius a
REAL ratio_ca = 0.6;   // ratio of radius c to radius a
std::vector<REAL> percent;  // mass percentage of particles smaller than a
certain size
std::vector<REAL> ptclSize; // particle size
percent.push_back(1.00); ptclSize.push_back(2.5e-3);
//percent.push_back(0.80); ptclSize.push_back(2.0e-3);
//percent.push_back(0.60); ptclSize.push_back(1.6e-3);
//percent.push_back(0.30); ptclSize.push_back(1.0e-3);
//percent.push_back(0.10); ptclSize.push_back(0.5e-3);
dem::Gradation grad(rorc, dimn, ratio_ba, ratio_ca, percent.size(), percent,
ptclSize);
assemb.deposit_PtclBdry(grad,
                        2,                  // particleLayers, setting of free
particles
                        1.0,                // relative container size,
0.8/1.0/1.2---small/medium/large
                        100000,             // total_steps
                        100,                // number of snapshots
                        10,                 // print interval
                        "flo_particle_end", // output file, initial particles
                        "dep_particle",     // output file, resulted particles,
including snapshots
                        "dep_contact",      // output file, resulted contacts,
including snapshots
                        "dep_progress",     // output file, statistical info
                        "dep_debug");       // output file, debug info


assemb.scale_PtclBdry(20000,             // total_steps
                      100,               // number of snapshots
                      10,                // print interval
                      0.05,              // dimension of
particle-composed-boundary
                      1.0,               // relative container size,
0.8/1.0/1.2---small/medium/large
                      "dep_particle_end",// input file, initial particles
                      "scl_particle",    // output file, resulted particles,
including snapshots
                      "scl_contact",     // output file, resulted contacts,
including snapshots
                      "scl_progress",    // output file, statistical info
                      "scl_debug");      // output file, debug info


assemb.ellipPile_Disp(50000,              // total_steps
                      100,                // number of snapshots
                      10,                 // print interval
                      0.05,               // dimension of
particle-composed-boundary
                      1.0,                // relative container size,
0.8/1.0/1.2---small/medium/large
                      "pile_particle_ini",// input file, initial particles, an
ellipsoidal pile info added
                      "pile_particle",    // output file, resulted particles,
including snapshots
                      "pile_contact",     // output file, resulted contacts,
including snapshots
                      "pile_progress",    // output file, statistical info
                      "pile_debug");      // output file, debug info


assemb.rectPile_Disp(50000,              // total_steps
                     100,                // number of snapshots
                     10,                 // print interval
                     "pile_particle_ini",// input file, initial particles
                     "pile_boundary_ini",// input file, initial boundaries,
rectangle pile boundary info added
                     "pile_particle",    // output file, resulted particles,
including snapshots
                     "pile_boundary",    // output file, resulted boundaries
                     "pile_contact",     // output file, resulted contacts,
including snapshots
                     "pile_progress",    // output file, statistical info
                     "pile_debug");      // output file, debug info


assemb.ellipPile_Impact(50000,              // total_steps
                        100,                // number of snapshots
                        10,                 // print interval
                        0.05,               // size of
particle-composed-boundary
                        "ipt_particle_ini", // input file, initial particles, an
ellipsoidal pile info added
                        "dep_boundary_ini", // input file, initial boundaries
                        "ipt_particle",     // output file, resulted particles,
including snapshots
                        "ipt_contact",      // output file, resulted contacts,
including snapshots
                        "ipt_progress",     // output file, statistical info
                        "ipt_debug");       // output file, debug info


assemb.ellipPile_Impact_p(50000,              // total_steps
                          100,                // number of snapshots
                          10,                 // print interval
                          0.05,               // size of
particle-composed-boundary
                          "ipt_particle_ini", // input file, initial particles,
an ellipsoidal pile info added
                          "ipt_particle",     // output file, resulted
particles, including snapshots
                          "ipt_contact",      // output file, resulted contacts,
including snapshots
                          "ipt_progress",     // output file, statistical info
                          "ipt_debug");       // output file, debug info

assemb.squeeze(300000,             // total_steps
               100000,             // initial_steps to reach equilibrium
               100,                // number of snapshots
               10,                 // print interval
               +1,                 // -1 squeeze; +1 loosen
               "flo_particle_end", // input file, initial particles
               "dep_boundary_ini", // input file, initial boundaries
               "dep_particle",     // output file, resulted particles, including
snapshots
               "dep_boundary",     // output file, resulted boundaries
               "dep_contact",      // output file, resulted contacts, including
snapshots
               "dep_progress",     // output file, statistical info
               "dep_debug");       // output file, debug info


assemb.collapse(rorc,
                100000,
                100,
                10,                // print interval
                "cre_particle",    // input file, initial particles
                "clp_boundary",    // output file, initial boundaries
                "clp_particle",    // output file, resulted particles, including
snapshots
                "clp_contact",     // output file, resulted contacts, including
snapshots
                "clp_progress",    // output file, statistical info
                "clp_debug");      // output file, debug info

*/
