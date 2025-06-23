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

/********************************************************************************
    Crack.cc
    PART SEVEN: UPDATE CRACK FRONT AND CALCULATE CRACK-FRONT NORMALS

    Created by Yajun Guo in 2002-2005.
********************************************************************************/

#include "Crack.h"
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <vector>

using namespace Uintah;

using std::string;
using std::vector;

void
Crack::addComputesAndRequiresCrackFrontNodeSubset(
  Task* /*t*/,
  const PatchSet* /*patches*/,
  const MaterialSet* /*matls*/) const
{
  // Currently do nothing
}

void
Crack::CrackFrontNodeSubset(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* /*matls*/,
                            DataWarehouse* /*old_dw*/,
                            DataWarehouse* /*new_dw*/)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int pid, patch_size;
    MPI_Comm_rank(d_mpi_crack_comm, &pid);
    MPI_Comm_size(d_mpi_crack_comm, &patch_size);

    int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      if (d_calFractParametersStep || d_doCrackPropagationStep) {
        // cfnset - subset of crack-front nodes in each patch
        d_cfnset[m][pid].clear();
        for (int j = 0; j < (int)d_cfSegNodes[m].size(); j++) {
          int node = d_cfSegNodes[m][j];
          if (patch->containsPointInExtraCells(d_cx[m][node])) {
            d_cfnset[m][pid].push_back(j);
          }
        }
        MPI_Barrier(d_mpi_crack_comm);

        // Broadcast cfnset to all the ranks
        for (int i = 0; i < patch_size; i++) {
          int num; // number of crack-front nodes in patch i
          if (i == pid) {
            num = d_cfnset[m][i].size();
          }
          MPI_Bcast(&num, 1, MPI_INT, i, d_mpi_crack_comm);
          d_cfnset[m][i].resize(num);
          MPI_Bcast(&d_cfnset[m][i][0], num, MPI_INT, i, d_mpi_crack_comm);
        }
        MPI_Barrier(d_mpi_crack_comm);

        // cfsset - subset of crack-front segment center in each patch
        d_cfsset[m][pid].clear();
        int ncfSegs = (int)d_cfSegNodes[m].size() / 2;
        for (int j = 0; j < ncfSegs; j++) {
          int n1     = d_cfSegNodes[m][2 * j];
          int n2     = d_cfSegNodes[m][2 * j + 1];
          Point cent = d_cx[m][n1] + (d_cx[m][n2] - d_cx[m][n1]) / 2.;
          if (patch->containsPointInExtraCells(cent)) {
            d_cfsset[m][pid].push_back(j);
          }
        }

        MPI_Barrier(d_mpi_crack_comm);

        // Broadcast cfsset to all the ranks
        for (int i = 0; i < patch_size; i++) {
          int num; // number of crack-front segments in patch i
          if (i == pid) {
            num = d_cfsset[m][i].size();
          }
          MPI_Bcast(&num, 1, MPI_INT, i, d_mpi_crack_comm);
          d_cfsset[m][i].resize(num);
          MPI_Bcast(&d_cfsset[m][i][0], num, MPI_INT, i, d_mpi_crack_comm);
        }
      }
    } // End of loop over matls
  }
}

void
Crack::addComputesAndRequiresRecollectCrackFrontSegments(
  Task* t,
  const PatchSet* /*patches*/,
  const MaterialSet* /*matls*/) const
{
  Ghost::GhostType gac = Ghost::AroundCells;
  int NGC              = 2 * d_NGN;
  t->needs(Task::NewDW, lb->gMassLabel, gac, NGC);
  t->needs(Task::NewDW, lb->GMassLabel, gac, NGC);
  t->needs(Task::OldDW, lb->pSizeLabel, Ghost::None);
  t->needs(Task::OldDW, lb->pDefGradLabel, Ghost::None);
}

void
Crack::RecollectCrackFrontSegments(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* /*matls*/,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    auto interpolator = d_flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    Vector dx = patch->dCell();

    int pid, patch_size;
    MPI_Comm_rank(d_mpi_crack_comm, &pid);
    MPI_Comm_size(d_mpi_crack_comm, &patch_size);

    delt_vartype curTimestep;
    old_dw->get(curTimestep, lb->delTLabel, getLevel(patches));

    int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

      // Cell mass of the material
      double d_cell_mass =
        mpm_matl->getInitialDensity() * dx.x() * dx.y() * dx.z();

      // Get nodal mass information
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      Ghost::GhostType gac = Ghost::AroundCells;
      constNCVariable<double> gMass, Gmass;
      int NGC = 2 * d_NGN;
      new_dw->get(gMass, lb->gMassLabel, dwi, patch, gac, NGC);
      new_dw->get(Gmass, lb->GMassLabel, dwi, patch, gac, NGC);

      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      old_dw->get(pSize, lb->pSizeLabel, pset);
      old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

      if (d_doCrackPropagationStep) {

        // Task 1: Detect if crack-front nodes are inside the material

        std::vector<short> cfSegNodesInMat;
        cfSegNodesInMat.resize(d_cfSegNodes[m].size());
        for (int i = 0; i < (int)d_cfSegNodes[m].size(); i++) {
          cfSegNodesInMat[i] = YES;
        }

        for (int i = 0; i < patch_size; i++) {
          int num      = d_cfnset[m][i].size();
          short* inMat = new short[num];

          if (pid == i) { // Rank i does it
            for (int j = 0; j < num; j++) {
              int idx  = d_cfnset[m][i][j];
              int node = d_cfSegNodes[m][idx];
              Point pt = d_cx[m][node];
              inMat[j] = YES;

              interpolator->findCellAndWeights(
                pt, ni, S, pSize[j], pDefGrad[j]);

              for (int k = 0; k < d_n8or27; k++) {
                double totalMass = gMass[ni[k]] + Gmass[ni[k]];
                if (totalMass < d_cell_mass / 32.) {
                  inMat[j] = NO;
                  break;
                }
              }
            }
          }

          MPI_Bcast(&inMat[0], num, MPI_SHORT, i, d_mpi_crack_comm);

          for (int j = 0; j < num; j++) {
            int idx              = d_cfnset[m][i][j];
            cfSegNodesInMat[idx] = inMat[j];
          }
          delete[] inMat;
        } // End of loop over patches

        MPI_Barrier(d_mpi_crack_comm);

        // Task 2: Detect if the centers of crack-front segments
        //         are inside the material

        std::vector<short> cfSegCenterInMat;
        cfSegCenterInMat.resize(d_cfSegNodes[m].size() / 2);
        for (int i = 0; i < (int)d_cfSegNodes[m].size() / 2; i++) {
          cfSegCenterInMat[i] = YES;
        }

        for (int i = 0; i < patch_size; i++) {
          int num      = d_cfsset[m][i].size();
          short* inMat = new short[num];

          if (pid == i) { // Rank i does it
            for (int j = 0; j < num; j++) {
              int idx    = d_cfsset[m][i][j];
              int node1  = d_cfSegNodes[m][2 * idx];
              int node2  = d_cfSegNodes[m][2 * idx + 1];
              Point cent = d_cx[m][node1] + (d_cx[m][node2] - d_cx[m][node1]) / 2.;
              inMat[j]   = YES;

              interpolator->findCellAndWeights(
                cent, ni, S, pSize[j], pDefGrad[j]);

              for (int k = 0; k < d_n8or27; k++) {
                double totalMass = gMass[ni[k]] + Gmass[ni[k]];
                if (totalMass < d_cell_mass / 32.) {
                  inMat[j] = NO;
                  break;
                }
              }
            }
          }

          MPI_Bcast(&inMat[0], num, MPI_SHORT, i, d_mpi_crack_comm);

          for (int j = 0; j < num; j++) {
            int idx               = d_cfsset[m][i][j];
            cfSegCenterInMat[idx] = inMat[j];
          }
          delete[] inMat;
        } // End of loop over patches

        MPI_Barrier(d_mpi_crack_comm);

        // Task 3: Recollect crack-front segments, discarding the
        //         dead ones. A segment is regarded dead if both
        //         ends of it are outside the material

        // Store crack-front parameters in temporary arraies
        int old_size  = (int)d_cfSegNodes[m].size();
        int* copyData = scinew int[old_size];

        for (int i = 0; i < old_size; i++) {
          copyData[i] = d_cfSegNodes[m][i];
        }
        d_cfSegNodes[m].clear();

        // Collect the active crack-front segments
        for (int i = 0; i < old_size / 2; i++) {
          int nd1 = copyData[2 * i];
          int nd2 = copyData[2 * i + 1];

          short thisSegActive = NO;
          if (cfSegNodesInMat[2 * i] || cfSegNodesInMat[2 * i + 1] ||
              cfSegCenterInMat[i]) {
            thisSegActive = YES;
          }

          if (thisSegActive) {
            d_cfSegNodes[m].push_back(nd1);
            d_cfSegNodes[m].push_back(nd2);
          }
        }
        delete[] copyData;

        if (d_cfSegNodes[m].size() > 0) { // New crack front is still in material
          // Seek the start crack point (sIdx), re-arrange crack-front nodes
          int node0 = d_cfSegNodes[m][0];
          int segs[2];
          FindSegsFromNode(m, node0, segs);

          if (segs[R] >= 0) {
            int numNodes = (int)d_cfSegNodes[m].size();
            int sIdx     = 0;
            for (int i = 0; i < numNodes; i++) {
              int segsT[2];
              FindSegsFromNode(m, d_cfSegNodes[m][i], segsT);
              if (segsT[R] < 0) {
                sIdx = i;
                break;
              }
            }

            int* copyData = new int[numNodes];
            int rIdx      = 2 * segs[R] + 1;
            for (int i = 0; i < numNodes; i++) {
              int oldIdx = sIdx + i;
              if (oldIdx > rIdx) {
                oldIdx -= (rIdx + 1);
              }
              copyData[i] = d_cfSegNodes[m][oldIdx];
            }

            for (int i = 0; i < numNodes; i++) {
              d_cfSegNodes[m][i] = copyData[i];
            }

            delete[] copyData;
          }

          // Task 4: Get previous index, and minimum & maximum indexes
          //         for crack-front nodes

          FindCrackFrontNodeIndexes(m);

          // Task 5: Calculate outer normals, tangential normals and binormals
          //         of crack plane at crack-front nodes

          if (d_smoothCrackFront) {
            short smoothSuccessfully = SmoothCrackFrontAndCalculateNormals(m);
            if (!smoothSuccessfully) {
              CalculateCrackFrontNormals(m);
            }
          } else { // Calculate crack-front normals directly
            CalculateCrackFrontNormals(m);
          }
        } // End if(d_cfSegNodes[m].size()>0)

        else { // Crack has penetrated the material
          // If all crack-front segments dead, the material is broken.
          if (d_ce[m].size() > 0) { // for the material with crack(s) initially
            if (pid == 0) {
              std::cout << "!!! Material " << m << " is broken." << std::endl;
            }
          }
        }

      } // End of if(d_doCrackPropagationStep!="false")

      // Save crack elements, crack nodes and crack-front nodes
      // for crack geometry visualization
      if (d_saveCrackGeometry) {
        if (pid == 0) {
          OutputCrackGeometry(m, curTimestep);
        }
      }

    } // End of loop over matls
    // delete interpolator;
  }
}

// Output cracks for visualization
void
Crack::OutputCrackGeometry(const int& m, const int& timestep)
{
  if (d_ce[m].size() > 0) { // for the materials with cracks
    bool timeToDump = d_dataArchiver->isOutputTimestep();
    if (timeToDump) {
      // Create output files in format:
      // ce.matXXX.timestepYYYYY (crack elems)
      // cx.matXXX.timestepYYYYY (crack points)
      // cf.matXXX.timestepYYYYY (crack front nodes)
      // Those files are stored in .uda.XXX/tXXXXX/crackData/

      char timestepbuf[10], matbuf[10];
      sprintf(timestepbuf, "%d", timestep);
      sprintf(matbuf, "%d", m);

      // Create output directories
      char crackDir[200] = "";
      strcat(crackDir, d_udaDir.c_str());
      strcat(crackDir, "/t");
      if (timestep < 10) {
        strcat(crackDir, "0000");
      } else if (timestep < 100) {
        strcat(crackDir, "000");
      } else if (timestep < 1000) {
        strcat(crackDir, "00");
      } else if (timestep < 10000) {
        strcat(crackDir, "0");
      }
      strcat(crackDir, timestepbuf);
      strcat(crackDir, "/crackData");

      MKDIR(crackDir, 0777);

      // Specify output file names
      char ceFileName[200] = "";
      strcat(ceFileName, crackDir);
      strcat(ceFileName, "/ce.mat");
      if (m < 10) {
        strcat(ceFileName, "00");
      } else if (m < 100) {
        strcat(ceFileName, "0");
      }
      strcat(ceFileName, matbuf);

      char cxFileName[200] = "";
      strcat(cxFileName, crackDir);
      strcat(cxFileName, "/cx.mat");
      if (m < 10) {
        strcat(cxFileName, "00");
      } else if (m < 100) {
        strcat(cxFileName, "0");
      }
      strcat(cxFileName, matbuf);

      char cfFileName[200] = "";
      strcat(cfFileName, crackDir);
      strcat(cfFileName, "/cf.mat");
      if (m < 10) {
        strcat(cfFileName, "00");
      } else if (m < 100) {
        strcat(cfFileName, "0");
      }
      strcat(cfFileName, matbuf);

      std::ofstream outputCE(ceFileName, std::ios::out);
      std::ofstream outputCX(cxFileName, std::ios::out);
      std::ofstream outputCF(cfFileName, std::ios::out);

      if (!outputCE || !outputCX || !outputCF) {
        std::cout << "Error: failure to open files for storing crack geometry"
                  << std::endl;
        exit(1);
      }

      // Output crack elems
      for (int i = 0; i < (int)d_ce[m].size(); i++) {
        outputCE << d_ce[m][i].x() << " " << d_ce[m][i].y() << " " << d_ce[m][i].z()
                 << std::endl;
      }

      // Output crack nodes
      for (int i = 0; i < (int)d_cx[m].size(); i++) {
        outputCX << d_cx[m][i].x() << " " << d_cx[m][i].y() << " " << d_cx[m][i].z()
                 << std::endl;
      }

      // Output crack-front nodes
      for (int i = 0; i < (int)d_cfSegNodes[m].size() / 2; i++) {
        outputCF << d_cfSegNodes[m][2 * i] << " " << d_cfSegNodes[m][2 * i + 1]
                 << std::endl;
      }
    }
  } // End if(ce[m].size()>0)
}
