/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 *  partextract.cc: Print out a uintah data archive for particle data
 *
 *  Written by:
 *   Biswajit Banerjee
 *   July 2004
 *
 */

#ifdef _WIN32
#define strtoll _strtoi64
#elif defined(__PGI)
// pgCC version 7.1-2 does not define strtoll (in stdlib.h or
// anywhere)... However, this seems to fake the compiler into
// not complaining.
#define _ISOC99_SOURCE
#endif

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/SymmMatrix3.h>
#include <Core/OS/Dir.h>
// #include <Core/Containers/Array3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <stdlib.h> // for strtoll

using namespace std;
using namespace Uintah;

// declarations
void
usage(const std::string& badarg, const std::string& progname);
void
getParticleStrains(DataArchive* da, int mat, long64 particleID, string flag);
void
getParticleStresses(DataArchive* da, int mat, long64 particleID, string flag);
void
printParticleVariable(DataArchive* da,
                      int mat,
                      string particleVariable,
                      long64 particleID,
                      unsigned long time_step_lower,
                      unsigned long time_step_upper,
                      unsigned long time_step_inc);
void
computeEquivStress(const Matrix3& sig, double& sigeqv);
void
computeEquivStrain(const Matrix3& F, double& epseqv);
void
computeTrueStrain(const Matrix3& F, Vector& strain);
void
computeGreenLagrangeStrain(const Matrix3& F, Matrix3& E);
void
computeGreenAlmansiStrain(const Matrix3& F, Matrix3& e);
void
computeStretchRotation(const Matrix3& F, Matrix3& R, Matrix3& U);
void
printCauchyStress(const Matrix3& stress);

// Main
int
main(int argc, char** argv)
{
  /*
   * Default values
   */
  int mat         = -1;
  bool do_partvar = false;
  //  bool do_partid=false;
  bool do_part_stress           = false;
  bool do_part_strain           = false;
  bool do_av_part_stress        = false;
  bool do_av_part_strain        = false;
  bool do_equiv_part_stress     = false;
  bool do_equiv_part_strain     = false;
  bool do_true_part_strain      = false;
  bool do_lagrange_part_strain  = false;
  bool do_euler_part_strain     = false;
  unsigned long time_step_lower = 0;
  unsigned long time_step_upper = 1;
  unsigned long time_step_inc   = 1;
  bool tslow_set                = false;
  bool tsup_set                 = false;
  //  bool tsinc_set = false;
  string filebase;
  string particleVariable;
  long64 particleID = 0;

  // set defaults for cout
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(8);
  /*
   * Parse arguments
   */
  for (int i = 1; i < argc; i++) {
    string s = argv[i];
    if (s == "-mat") {
      mat = atoi(argv[++i]);
    } else if (s == "-partvar") {
      do_partvar       = true;
      particleVariable = argv[++i];
      if (particleVariable[0] == '-') {
        usage("-partvar <particle variable name>", argv[0]);
      }
    } else if (s == "-partid") {
      //      do_partid=true;
      string id = argv[++i];
      if (id[0] == '-') {
        usage("-partid <particle id>", argv[0]);
      }
      particleID = strtoll(argv[i], (char**)nullptr, 10);
    } else if (s == "-part_stress") {
      do_part_stress = true;
      if (++i < argc) {
        s = argv[i];
        if (s == "avg") {
          do_av_part_stress = true;
        } else if (s == "equiv") {
          do_equiv_part_stress = true;
        } else if (s == "all") {
          do_part_stress = true;
        } else {
          usage("-part_stress [avg or equiv or all]", argv[0]);
        }
      }
    } else if (s == "-part_strain") {
      do_part_strain = true;
      if (++i < argc) {
        s = argv[i];
        if (s == "avg") {
          do_av_part_strain = true;
        } else if (s == "true") {
          do_true_part_strain = true;
        } else if (s == "equiv") {
          do_equiv_part_strain = true;
        } else if (s == "all") {
          do_part_strain = true;
        } else if (s == "lagrangian") {
          do_lagrange_part_strain = true;
        } else if (s == "eulerian") {
          do_euler_part_strain = true;
        } else {
          usage("-part_strain [avg / true / all / lagrangian / eulerian]",
                argv[0]);
        }
      }
    } else if (s == "-timesteplow") {
      time_step_lower = strtoul(argv[++i], (char**)nullptr, 10);
      tslow_set       = true;
    } else if (s == "-timestephigh") {
      time_step_upper = strtoul(argv[++i], (char**)nullptr, 10);
      tsup_set        = true;
    } else if (s == "-timestepinc") {
      time_step_inc = strtoul(argv[++i], (char**)nullptr, 10);
      //      tsinc_set = true;
    }
  }
  filebase = argv[argc - 1];

  if (filebase == "" || filebase == argv[0]) {
    std::cerr << "No archive file specified\n";
    usage("", argv[0]);
  }

  try {
    DataArchive* da = scinew DataArchive(filebase);

    // Get the particle stresses
    if (do_part_stress) {
      if (do_av_part_stress) {
        std::cout << "\t Volume average stress = " << std::endl;
        getParticleStresses(da, mat, particleID, "avg");
      } else if (do_equiv_part_stress) {
        getParticleStresses(da, mat, particleID, "equiv");
      } else {
        getParticleStresses(da, mat, particleID, "all");
      }
    }

    // Get the particle strains
    if (do_part_strain) {
      if (do_av_part_strain) {
        std::cout << "\t Volume average strain = " << std::endl;
        getParticleStrains(da, mat, particleID, "avg");
      } else if (do_true_part_strain) {
        getParticleStrains(da, mat, particleID, "true");
      } else if (do_equiv_part_strain) {
        getParticleStrains(da, mat, particleID, "equiv");
      } else if (do_lagrange_part_strain) {
        getParticleStrains(da, mat, particleID, "lagrange");
      } else if (do_euler_part_strain) {
        getParticleStrains(da, mat, particleID, "euler");
      } else {
        getParticleStrains(da, mat, particleID, "all");
      }
    }

    // Print a particular particle variable
    if (do_partvar) {
      std::vector<int> index;
      std::vector<double> times;
      da->queryTimesteps(index, times);
      ASSERTEQ(index.size(), times.size());
      if (!tslow_set) {
        time_step_lower = 0;
      } else if (time_step_lower >= times.size()) {
        std::cerr << "timesteplow must be between 0 and " << times.size() - 1
                  << std::endl;
        abort();
      }
      if (!tsup_set) {
        time_step_upper = times.size() - 1;
      } else if (time_step_upper >= times.size()) {
        std::cerr << "timestephigh must be between 0 and " << times.size() - 1
                  << std::endl;
        abort();
      }
      printParticleVariable(da,
                            mat,
                            particleVariable,
                            particleID,
                            time_step_lower,
                            time_step_upper,
                            time_step_inc);
    }
  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << std::endl;
    abort();
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    abort();
  }
}

void
usage(const std::string& badarg, const std::string& progname)
{
  if (badarg != "") {
    std::cerr << "Error parsing argument: " << badarg << std::endl;
  }

  std::cerr << "Usage: " << progname << " [options] <archive file>\n\n";
  std::cerr << "Valid options are:\n";
  std::cerr << "  -mat <material id>\n";
  std::cerr << "  -partvar <variable name>\n";
  std::cerr << "  -partid <particleid>\n";
  std::cerr << "  -part_stress [avg or equiv or all]\n";
  std::cerr << "  -part_strain [avg/true/equiv/all/lagrangian/eulerian]\n",
    std::cerr << "  -timesteplow [int] (only outputs timestep from int)\n";
  std::cerr << "  -timestephigh [int] (only outputs timesteps upto int)\n";
  std::cerr << "USAGE IS NOT FINISHED\n\n";
  exit(1);
}

////////////////////////////////////////////////////////////////////////////
//
// Get particle strains (avg, true or all)
//
////////////////////////////////////////////////////////////////////////////
void
getParticleStrains(DataArchive* da, int mat, long64 particleID, string flag)
{

  // Parse the flag and check which option is needed
  bool doAverage  = false;
  bool doTrue     = false;
  bool doLagrange = false;
  bool doEuler    = false;
  bool doEquiv    = false;
  //   bool doAll = false;
  if (flag == "avg") {
    doAverage = true;
  } else if (flag == "equiv") {
    doEquiv = true;
  } else if (flag == "true") {
    doTrue = true;
  } else if (flag == "lagrange") {
    doLagrange = true;
  } else if (flag == "euler") {
    doEuler = true;
  }
  //   else doAll = true;

  // Check if all the required variables are there .. for all cases
  // we need p.deformationGradient and for the volume average we need p.volume
  std::vector<std::string> vars;
  std::vector<int> num_matls;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matls, types);
  ASSERTEQ(vars.size(), types.size());

  bool gotVolume = false;
  bool gotDeform = false;
  for (unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == "p.volume") {
      gotVolume = true;
    }
    if (var == "p.deformationGradient") {
      gotDeform = true;
    }
  }
  if (!gotDeform) {
    std::cerr << "\n **Error** getParticleStrains : DataArchiver does not "
              << "contain p.deformationGradient\n";
    exit(1);
  }
  if (doAverage && !gotVolume) {
    std::cerr << "\n **Error** getParticleStrains : DataArchiver does not "
              << "contain p.volume\n";
    exit(1);
  }

  // Now that the variables have been found, get the data for all available
  // time steps  from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  //  std::cerr <<  "There are " << index.size() << " timesteps:\n";

  unsigned long time_step_lower = 0;
  unsigned long time_step_upper = times.size() - 1;

  // Loop thru all time steps and store the volume and variable (stress/strain)
  for (unsigned long t = time_step_lower; t <= time_step_upper; t++) {
    double time = times[t];
    // std::cerr <<  "Time = " << time << std::endl;
    GridP grid = da->queryGrid(t);

    std::vector<double> volumeVector;
    std::vector<Matrix3> deformVector;
    double totVol = 0.0;

    // Loop thru all the levels
    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);

      // Loop thru all the patches
      Level::const_patch_iterator iter = level->patchesBegin();
      int patchIndex                  = 0;
      for (; iter != level->patchesEnd(); iter++) {
        const Patch* patch = *iter;
        ++patchIndex;

        // Loop thru all the variables
        for (int v = 0; v < (int)vars.size(); v++) {
          std::string var = vars[v];

          if (var != "p.volume" && var != "p.deformationGradient") {
            continue;
          }
          const Uintah::TypeDescription* td = types[v];

          // Check if the variable is a ParticleVariable
          if (td->getType() ==
              Uintah::TypeDescription::Type::ParticleVariable) {

            // loop thru all the materials
            ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
            ConsecutiveRangeSet::iterator matlIter = matls.begin();
            for (; matlIter != matls.end(); matlIter++) {
              int matl = *matlIter;

              if (mat != -1 && matl != mat) {
                continue;
              }

              // Find the name of the variable
              if (doAverage) {
                std::cerr << "Finding volume average strains ... " << std::endl;
                if (var == "p.volume") {
                  switch (td->getSubType()->getType()) {
                    case Uintah::TypeDescription::Type::double_type: {
                      ParticleVariable<double> value;
                      da->query(value, var, matl, patch, t);
                      ParticleSubset* pset = value.getParticleSubset();
                      if (pset->numParticles() > 0) {
                        ParticleSubset::iterator iter = pset->begin();
                        for (; iter != pset->end(); iter++) {
                          volumeVector.push_back(value[*iter]);
                          totVol += value[*iter];
                        }
                      }
                    } break;
                    case Uintah::TypeDescription::Type::float_type: {
                      ParticleVariable<float> value;
                      da->query(value, var, matl, patch, t);
                      ParticleSubset* pset = value.getParticleSubset();
                      if (pset->numParticles() > 0) {
                        ParticleSubset::iterator iter = pset->begin();
                        for (; iter != pset->end(); iter++) {
                          volumeVector.push_back((double)(value[*iter]));
                          totVol += value[*iter];
                        }
                      }
                    } break;
                    default:
                      cerr
                        << "Particle Variable of unknown type: "
                        << static_cast<
                             std::underlying_type<TypeDescription::Type>::type>(
                             td->getSubType()->getType())
                        << std::endl;
                      break;
                  }
                } else if (var == "p.deformationGradient") {
                  // std::cerr <<  "Material: " << matl << std::endl;
                  ParticleVariable<Matrix3> value;
                  da->query(value, var, matl, patch, t);
                  ParticleSubset* pset = value.getParticleSubset();
                  if (pset->numParticles() > 0) {
                    ParticleSubset::iterator iter = pset->begin();
                    for (; iter != pset->end(); iter++) {
                      particleIndex idx = *iter;

                      // Find the right stretch tensor
                      Matrix3 deformGrad = value[idx];
                      Matrix3 stretch(0.0), rotation(0.0);
                      deformGrad.polarDecomposition(
                        stretch, rotation, 1.0e-10, true);

                      cout.setf(ios::scientific, ios::floatfield);
                      cout.precision(6);

                      deformVector.push_back(stretch);
                    }
                  } // end of pset > 0
                }   // end of var if
                continue;
              } // end of doAverage

              // If not an average calculation
              if (var == "p.deformationGradient") {
                // std::cerr <<  "Material: " << matl << std::endl;
                ParticleVariable<Matrix3> value;
                da->query(value, var, matl, patch, t);
                ParticleVariable<long64> pid;
                da->query(pid, "p.particleID", matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if (pset->numParticles() > 0) {
                  ParticleSubset::iterator iter = pset->begin();
                  if (particleID == 0) {
                    for (; iter != pset->end(); iter++) {
                      particleIndex idx = *iter;

                      // Get the deformation gradient
                      Matrix3 F = value[idx];
                      cout.setf(ios::scientific, ios::floatfield);
                      cout.precision(6);

                      std::cout << time << " " << pid[idx] << " " << patchIndex
                                << " " << matl;
                      if (doTrue) {
                        Vector trueStrain(0.0);
                        computeTrueStrain(F, trueStrain);
                      } else if (doEquiv) {
                        double equiv_strain = 0.0;
                        computeEquivStrain(F, equiv_strain);
                      } else if (doLagrange) {
                        Matrix3 E(0.0);
                        computeGreenLagrangeStrain(F, E);
                      } else if (doEuler) {
                        Matrix3 e(0.0);
                        computeGreenAlmansiStrain(F, e);
                      } else {
                        Matrix3 U(0.0), R(0.0);
                        computeStretchRotation(F, R, U);
                      }
                    }
                  } else {
                    for (; iter != pset->end(); iter++) {
                      particleIndex idx = *iter;

                      if (particleID != pid[*iter]) {
                        continue;
                      }

                      // Get the deformation gradient
                      Matrix3 F = value[idx];
                      cout.setf(ios::scientific, ios::floatfield);
                      cout.precision(6);

                      std::cout << time << " " << pid[idx] << " " << patchIndex
                                << " " << matl;
                      if (doTrue) {
                        Vector trueStrain(0.0);
                        computeTrueStrain(F, trueStrain);
                      } else if (doEquiv) {
                        double equiv_strain = 0.0;
                        computeEquivStrain(F, equiv_strain);
                      } else if (doLagrange) {
                        Matrix3 E(0.0);
                        computeGreenLagrangeStrain(F, E);
                      } else if (doEuler) {
                        Matrix3 e(0.0);
                        computeGreenAlmansiStrain(F, e);
                      } else {
                        Matrix3 U(0.0), R(0.0);
                        computeStretchRotation(F, R, U);
                      }
                      break;
                    } // end of pset iter loop
                  }   // end of particleId = 0 if
                }     // end of numparticle > 0 if
              }       // end of var compare if
            }         // end of material loop
          }           // end of ParticleVariable if
        }             // end of variable loop
      }               // end of patch loop
    }                 // end of level loop

    if (doAverage) {
      // Now that the volume vector and variable vector are available just
      // do a weighted average
      ASSERTEQ(volumeVector.size(), deformVector.size());
      Matrix3 avVar;
      for (unsigned int ii = 0; ii < volumeVector.size(); ++ii) {
        avVar += ((deformVector[ii] * volumeVector[ii]) / totVol);
      }
      for (int ii = 0; ii < 3; ++ii) {
        for (int jj = 0; jj < 3; ++jj) {
          std::cout << avVar(ii, jj) << "  ";
        }
        std::cout << std::endl;
      }
    }
  } // end of time step loop
}

////////////////////////////////////////////////////////////////////////////
//
// Get particle stresses (avg, equiv or all)
//
////////////////////////////////////////////////////////////////////////////
void
getParticleStresses(DataArchive* da, int mat, long64 particleID, string flag)
{

  // Parse the flag and check which option is needed
  bool doAverage = false;
  bool doEquiv   = false;
  //   bool doAll = false;
  if (flag == "avg") {
    doAverage = true;
  } else if (flag == "equiv") {
    doEquiv = true;
  }
  //   else doAll = true;

  // Check if all the required variables are there .. for all cases
  // we need p.stress and for the volume average we need p.volume
  std::vector<std::string> vars;
  std::vector<int> num_matls;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matls, types);
  ASSERTEQ(vars.size(), types.size());

  bool gotVolume = false;
  bool gotStress = false;
  for (unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == "p.volume") {
      gotVolume = true;
    }
    if (var == "p.stress") {
      gotStress = true;
    }
  }
  if (!gotStress) {
    std::cerr << "\n **Error** getParticleStresses : DataArchiver does not "
              << "contain p.stress\n";
    exit(1);
  }
  if (doAverage && !gotVolume) {
    std::cerr << "\n **Error** getParticleStresses : DataArchiver does not "
              << "contain p.volume\n";
    exit(1);
  }

  // Now that the variables have been found, get the data for all available
  // time steps // from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  // std::cout << "There are " << index.size() << " timesteps:\n";

  unsigned long time_step_lower = 0;
  unsigned long time_step_upper = times.size() - 1;

  // Loop thru all time steps and store the volume and variable (stress/strain)
  for (unsigned long t = time_step_lower; t <= time_step_upper; t++) {
    double time = times[t];
    // std::cout << "Time = " << time << std::endl;
    GridP grid = da->queryGrid(t);

    std::vector<double> volumeVector;
    std::vector<Matrix3> stressVector;
    double totVol = 0.0;

    // Loop thru all the levels
    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);

      // Loop thru all the patches
      Level::const_patch_iterator iter = level->patchesBegin();
      int patchIndex                  = 0;
      for (; iter != level->patchesEnd(); iter++) {
        const Patch* patch = *iter;
        ++patchIndex;

        // Loop thru all the variables
        for (int v = 0; v < (int)vars.size(); v++) {
          std::string var                   = vars[v];
          const Uintah::TypeDescription* td = types[v];

          // Check if the variable is a ParticleVariable
          if (td->getType() ==
              Uintah::TypeDescription::Type::ParticleVariable) {

            // loop thru all the materials
            ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
            ConsecutiveRangeSet::iterator matlIter = matls.begin();
            for (; matlIter != matls.end(); matlIter++) {
              int matl = *matlIter;

              if (mat != -1 && matl != mat) {
                continue;
              }

              // Find the name of the variable
              if (doAverage) {
                if (var == "p.volume") {
                  switch (td->getSubType()->getType()) {
                    case Uintah::TypeDescription::Type::double_type: {
                      ParticleVariable<double> value;
                      da->query(value, var, matl, patch, t);
                      ParticleSubset* pset = value.getParticleSubset();
                      if (pset->numParticles() > 0) {
                        ParticleSubset::iterator iter = pset->begin();
                        for (; iter != pset->end(); iter++) {
                          volumeVector.push_back(value[*iter]);
                          totVol += value[*iter];
                        }
                      }
                    } break;
                    case Uintah::TypeDescription::Type::float_type: {
                      ParticleVariable<float> value;
                      da->query(value, var, matl, patch, t);
                      ParticleSubset* pset = value.getParticleSubset();
                      if (pset->numParticles() > 0) {
                        ParticleSubset::iterator iter = pset->begin();
                        for (; iter != pset->end(); iter++) {
                          volumeVector.push_back((double)(value[*iter]));
                          totVol += value[*iter];
                        }
                      }
                    } break;
                    default:
                      cerr
                        << "Particle Variable of unknown type: "
                        << static_cast<
                             std::underlying_type<TypeDescription::Type>::type>(
                             td->getSubType()->getType());
                      break;
                  }
                } else if (var == "p.stress") {
                  ParticleVariable<Matrix3> value;
                  da->query(value, var, matl, patch, t);
                  ParticleSubset* pset = value.getParticleSubset();
                  if (pset->numParticles() > 0) {
                    ParticleSubset::iterator iter = pset->begin();
                    for (; iter != pset->end(); iter++) {
                      stressVector.push_back(value[*iter]);
                    }
                  }
                }
                continue;
              }

              // If not a volume averaged quantity
              if (var == "p.stress") {
                ParticleVariable<Matrix3> value;
                da->query(value, var, matl, patch, t);
                ParticleVariable<long64> pid;
                da->query(pid, "p.particleID", matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if (pset->numParticles() > 0) {
                  ParticleSubset::iterator iter = pset->begin();
                  for (; iter != pset->end(); iter++) {
                    if (particleID == 0) {
                      Matrix3 stress = value[*iter];

                      cout.setf(ios::scientific, ios::floatfield);
                      cout.precision(6);

                      std::cout << time << " " << pid[*iter] << " "
                                << patchIndex << " " << matl;

                      if (doEquiv) {
                        double sigeff = 0.0;
                        computeEquivStress(stress, sigeff);
                      } else {
                        printCauchyStress(stress);
                      }
                    } else {
                      Matrix3 stress = value[*iter];

                      if (particleID != pid[*iter]) {
                        continue;
                      }

                      cout.setf(ios::scientific, ios::floatfield);
                      cout.precision(6);

                      std::cout << time << " " << pid[*iter] << " "
                                << patchIndex << " " << matl;

                      if (doEquiv) {
                        double sigeff = 0.0;
                        computeEquivStress(stress, sigeff);
                      } else {
                        printCauchyStress(stress);
                      }
                      break;
                    }
                  }
                }
              } // end of var compare if
            }   // end of material loop
          }     // end of ParticleVariable if
        }       // end of variable loop
      }         // end of patch loop
    }           // end of level loop

    if (doAverage) {
      // Now that the volume vector and variable vector are available just
      // do a weighted average
      ASSERTEQ(volumeVector.size(), stressVector.size());
      Matrix3 avVar;
      for (unsigned int ii = 0; ii < volumeVector.size(); ++ii) {
        avVar += ((stressVector[ii] * volumeVector[ii]) / totVol);
      }
      for (int ii = 0; ii < 3; ++ii) {
        for (int jj = 0; jj < 3; ++jj) {
          std::cout << avVar(ii, jj) << "  ";
        }
        std::cout << std::endl;
      }
    }
  } // end of time step loop
}

////////////////////////////////////////////////////////////////////////////
//
// Print a particle variable
//
////////////////////////////////////////////////////////////////////////////
void
printParticleVariable(DataArchive* da,
                      int mat,
                      string particleVariable,
                      long64 particleID,
                      unsigned long time_step_lower,
                      unsigned long time_step_upper,
                      unsigned long time_step_inc)
{

  // Check if the particle variable is available
  std::vector<std::string> vars;
  std::vector<int> num_matls;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matls, types);
  ASSERTEQ(vars.size(), types.size());
  bool variableFound = false;
  for (unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == particleVariable) {
      variableFound = true;
    }
  }
  if (!variableFound) {
    std::cerr << "Variable " << particleVariable << " not found\n";
    exit(1);
  }

  // Now that the variable has been found, get the data for all
  // available time steps // from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  // std::cout << "There are " << index.size() << " timesteps:\n";

  // Loop thru all time steps and store the volume and variable (stress/strain)
  for (unsigned long t = time_step_lower; t <= time_step_upper;
       t += time_step_inc) {
    double time = times[t];
    // std::cout << "Time = " << time << std::endl;
    GridP grid = da->queryGrid(t);

    // Loop thru all the levels
    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);

      // Loop thru all the patches
      int patchIndex = 0;
      for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
           iter++) {
        const Patch* patch = *iter;
        ++patchIndex;

        // Loop thru all the variables
        std::size_t v = 0;
        for (const auto& var : vars) {
          const Uintah::TypeDescription* td      = types[v++];
          const Uintah::TypeDescription* subtype = td->getSubType();

          // Check if the variable is a ParticleVariable
          if (td->getType() ==
              Uintah::TypeDescription::Type::ParticleVariable) {

            // loop thru all the materials
            auto matls = da->queryMaterials(var, patch, t);
            for (const auto& matl : matls) {

              if (mat != -1 && matl != mat) {
                continue;
              }

              // Find the name of the variable
              if (var == particleVariable) {
                // std::cout << "Material: " << matl << std::endl;
                switch (subtype->getType()) {
                  case Uintah::TypeDescription::Type::double_type: {
                    ParticleVariable<double> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    // std::cout << "numParticles = " << pset->numParticles() <<
                    // std::endl;
                    if (pset->numParticles() > 0) {
                      if (particleID == 0) {
                        for (const auto& pidx : *pset) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[pidx];
                          std::cout << " " << value[pidx] << std::endl;
                        }
                      } else {
                        for (const auto& pidx : *pset) {
                          if (particleID != pid[pidx]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[pidx];
                          std::cout << " " << value[pidx] << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::float_type: {
                    ParticleVariable<float> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      if (particleID == 0) {
                        for (; iter != pset->end(); iter++) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter] << std::endl;
                        }
                      } else {
                        for (; iter != pset->end(); iter++) {
                          if (particleID != pid[*iter]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter] << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::int_type: {
                    ParticleVariable<int> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      if (particleID == 0) {
                        for (; iter != pset->end(); iter++) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter] << std::endl;
                        }
                      } else {
                        for (; iter != pset->end(); iter++) {
                          if (particleID != pid[*iter]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter] << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Point: {
                    ParticleVariable<Point> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      if (particleID == 0) {
                        for (; iter != pset->end(); iter++) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter](0) << " "
                                    << value[*iter](1) << " " << value[*iter](2)
                                    << std::endl;
                        }
                      } else {
                        for (; iter != pset->end(); iter++) {
                          if (particleID != pid[*iter]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter](0) << " "
                                    << value[*iter](1) << " " << value[*iter](2)
                                    << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Vector: {
                    ParticleVariable<Vector> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      if (particleID == 0) {
                        for (; iter != pset->end(); iter++) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter][0] << " "
                                    << value[*iter][1] << " " << value[*iter][2]
                                    << std::endl;
                        }
                      } else {
                        for (; iter != pset->end(); iter++) {
                          if (particleID != pid[*iter]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          std::cout << " " << value[*iter][0] << " "
                                    << value[*iter][1] << " " << value[*iter][2]
                                    << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Matrix3: {
                    ParticleVariable<Matrix3> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    da->query(pid, "p.particleID", matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      if (particleID == 0) {
                        for (; iter != pset->end(); iter++) {
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          for (int ii = 0; ii < 3; ++ii) {
                            for (int jj = 0; jj < 3; ++jj) {
                              std::cout << " " << value[*iter](ii, jj);
                            }
                          }
                          std::cout << std::endl;
                        }
                      } else {
                        for (; iter != pset->end(); iter++) {
                          if (particleID != pid[*iter]) {
                            continue;
                          }
                          std::cout << time << " " << patchIndex << " " << matl;
                          std::cout << " " << pid[*iter];
                          for (int ii = 0; ii < 3; ++ii) {
                            for (int jj = 0; jj < 3; ++jj) {
                              std::cout << " " << value[*iter](ii, jj);
                            }
                          }
                          std::cout << std::endl;
                        }
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::long64_type: {
                    ParticleVariable<long64> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        std::cout << " " << value[*iter] << std::endl;
                      }
                    }
                  } break;
                  default:
                    cerr
                      << "Particle Variable of unknown type: "
                      << static_cast<
                           std::underlying_type<TypeDescription::Type>::type>(
                           subtype->getType())
                      << "\n";
                    break;
                }
              } // end of var compare if
            }   // end of material loop
          }     // end of ParticleVariable if
        }       // end of variable loop
      }         // end of patch loop
    }           // end of level loop
  }             // end of time step loop
}

void
computeEquivStress(const Matrix3& stress, double& sigeff)
{
  Matrix3 I;
  I.Identity();
  Matrix3 s = stress - I * stress.Trace() / 3.0;
  sigeff    = sqrt(1.5 * s.NormSquared());

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  std::cout << " " << sigeff << std::endl;
}

void
computeEquivStrain(const Matrix3& F, double& equiv_strain)
{
  Matrix3 I;
  I.Identity();
  Matrix3 E    = (F.Transpose() * F - I) * 0.5;
  equiv_strain = sqrt(E.NormSquared() / 1.5);

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  std::cout << " " << equiv_strain;
  std::cout << std::endl;
}

void
computeTrueStrain(const Matrix3& F, Vector& strain)
{
  // Compute the left Cauchy-Green tensor
  Matrix3 C = F.Transpose() * F;

  // Find the eigenvalues of C
  Vector Lambda;
  Matrix3 direction;
  C.eigen(Lambda, direction);

  // Find the stretches
  Vector lambda;
  for (int ii = 0; ii < 3; ++ii) {
    lambda[ii] = sqrt(Lambda[ii]);
    strain[ii] = log(lambda[ii]);
  }

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  for (int ii = 0; ii < 3; ++ii) {
    std::cout << " " << strain[ii];
  }
  std::cout << std::endl;
}

void
computeGreenLagrangeStrain(const Matrix3& F, Matrix3& E)
{
  Matrix3 I;
  I.Identity();
  Matrix3 C = F.Transpose() * F;
  E         = (C - I) * 0.5;

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  for (int ii = 0; ii < 3; ++ii) {
    std::cout << " " << E(ii, ii);
  }
  std::cout << " " << E(1, 2);
  std::cout << " " << E(2, 0);
  std::cout << " " << E(0, 1);
  std::cout << std::endl;
}

// Calculate the Almansi-Hamel strain tensor
void
computeGreenAlmansiStrain(const Matrix3& F, Matrix3& e)
{
  Matrix3 I;
  I.Identity();
  Matrix3 b    = F * F.Transpose();
  Matrix3 binv = b.Inverse();
  e            = (I - binv) * 0.5;

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  for (int ii = 0; ii < 3; ++ii) {
    std::cout << " " << e(ii, ii);
  }
  std::cout << " " << e(1, 2);
  std::cout << " " << e(2, 0);
  std::cout << " " << e(0, 1);
  std::cout << std::endl;
}

void
computeStretchRotation(const Matrix3& F, Matrix3& R, Matrix3& U)
{
  F.polarDecomposition(U, R, 1.0e-10, true);

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      std::cout << " " << F(ii, jj);
    }
  }
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      std::cout << " " << U(ii, jj);
    }
  }
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      std::cout << " " << R(ii, jj);
    }
  }
  std::cout << std::endl;
}

void
printCauchyStress(const Matrix3& stress)
{
  double sig11 = stress(0, 0);
  double sig12 = stress(0, 1);
  double sig13 = stress(0, 2);
  double sig22 = stress(1, 1);
  double sig23 = stress(1, 2);
  double sig33 = stress(2, 2);

  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(6);

  std::cout << " " << sig11 << " " << sig22 << " " << sig33 << " " << sig23
            << " " << sig13 << " " << sig12 << std::endl;
}
