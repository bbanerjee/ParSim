/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConduction.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/ImpMPMSolvers/PetscSolver.h>
#include <CCA/Components/MPM/ImpMPMSolvers/SimpleSolver.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Short27.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/DebugStream.h>
#include <vector>

using namespace Uintah;

static DebugStream cout_doing("ImplicitHeatConduction", false);

ImplicitHeatConduction::ImplicitHeatConduction(
  const MaterialManagerP& mat_manager,
  const MPMLabel* labels,
  const MPMFlags* flags)
  : d_mat_manager{ mat_manager }
  , d_mpm_labels{ labels }
  , d_mpm_flags{ flags }
{
  do_IHC         = d_mpm_flags->d_doImplicitHeatConduction;
  d_HC_transient = d_mpm_flags->d_doTransientImplicitHeatConduction;

  d_one_matl = scinew MaterialSubset();
  d_one_matl->add(0);
  d_one_matl->addReference();

  if (d_mpm_flags->d_8or27 == 8) {
    d_num_ghost_particles = 1;
    d_num_ghost_nodes     = 1;
  } else if (d_mpm_flags->d_8or27 == 27 || flags->d_8or27 == 64) {
    d_num_ghost_particles = 2;
    d_num_ghost_nodes     = 2;
  }
}

ImplicitHeatConduction::~ImplicitHeatConduction()
{
  if (d_one_matl && d_one_matl->removeReference()) {
    delete d_one_matl;
  }

  if (do_IHC) {
    d_perproc_patches->removeReference();
    std::cout << "Freeing patches!!\n";
  }
}

void
ImplicitHeatConduction::problemSetup(string solver_type)
{
  if (solver_type == "petsc") {
    d_HC_solver = std::make_unique<MPMPetscSolver>();
    d_HC_solver->initialize();
  } else {
    d_HC_solver = std::make_unique<SimpleSolver>();
    d_HC_solver->initialize();
  }
}

void
ImplicitHeatConduction::scheduleDestroyHCMatrix(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task("ImpMPM::destroyHCMatrix",
                          this,
                          &ImplicitHeatConduction::destroyHCMatrix);
    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleCreateHCMatrix(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task(
      "ImpMPM::createHCMatrix", this, &ImplicitHeatConduction::createHCMatrix);

    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);

    d_perproc_patches = patches;
    d_perproc_patches->addReference();

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleApplyHCBoundaryConditions(
  SchedulerP& schd,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task("ImpMPM::applyHCBoundaryCondition",
                          this,
                          &ImplicitHeatConduction::applyHCBoundaryConditions);

    t->computes(d_mpm_labels->gTemperatureStarLabel, d_one_matl);
    t->needs(
      Task::NewDW, d_mpm_labels->gExternalHeatFluxLabel, Ghost::None, 0);

    t->setType(Task::OncePerProc);
    schd->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleFindFixedHCDOF(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task(
      "ImpMPM::findFixedHCDOF", this, &ImplicitHeatConduction::findFixedHCDOF);

    t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None, 0);

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleFormHCStiffnessMatrix(SchedulerP& sched,
                                                      const PatchSet* patches,
                                                      const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task("ImpMPM::formHCStiffnessMatrix",
                          this,
                          &ImplicitHeatConduction::formHCStiffnessMatrix);

    t->needs(Task::OldDW, d_mpm_labels->delTLabel);
    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);
    t->needs(Task::OldDW, d_mpm_labels->pVolumeLabel, Ghost::AroundNodes, 1);
    t->needs(
      Task::OldDW, d_mpm_labels->pTemperatureLabel, Ghost::AroundNodes, 1);
    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleFormHCQ(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t =
      scinew Task("ImpMPM::formHCQ", this, &ImplicitHeatConduction::formHCQ);

    t->needs(Task::OldDW, d_mpm_labels->delTLabel);
    t->needs(Task::NewDW,
                d_mpm_labels->gTemperatureLabel,
                d_one_matl,
                Ghost::AroundCells,
                1);
    t->needs(
      Task::NewDW, d_mpm_labels->gExternalHeatRateLabel, Ghost::AroundCells, 1);
    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);
    t->needs(Task::OldDW, d_mpm_labels->pVolumeLabel, Ghost::AroundNodes, 1);

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleAdjustHCQAndHCKForBCs(SchedulerP& sched,
                                                      const PatchSet* patches,
                                                      const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task("ImpMPM::adjustHCQAndHCKForBCs",
                          this,
                          &ImplicitHeatConduction::adjustHCQAndHCKForBCs);

    Ghost::GhostType gnone = Ghost::None;

    t->needs(
      Task::NewDW, d_mpm_labels->gTemperatureStarLabel, d_one_matl, gnone, 0);

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleSolveForTemp(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task(
      "ImpMPM::solveForTemp", this, &ImplicitHeatConduction::solveForTemp);

#if 0
  Ghost::GhostType  gnone = Ghost::None;
  t->needs(Task::NewDW, d_mpm_labels->gTemperatureLabel,d_one_matl,gnone,0);
#endif

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::scheduleGetTemperatureIncrement(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  if (do_IHC) {
    Task* t = scinew Task("ImpMPM::getTemperatureIncrement",
                          this,
                          &ImplicitHeatConduction::getTemperatureIncrement);

    t->needs(Task::OldDW, d_mpm_labels->delTLabel);
    t->needs(
      Task::NewDW, d_mpm_labels->gTemperatureLabel, d_one_matl, Ghost::None, 0);
    t->computes(d_mpm_labels->gTemperatureRateLabel, d_one_matl);

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  } else {
    Task* t = scinew Task("ImpMPM::fillgTemperatureRate",
                          this,
                          &ImplicitHeatConduction::fillgTemperatureRate);

    t->computes(d_mpm_labels->gTemperatureRateLabel, d_one_matl);

    sched->addTask(t, patches, matls);
  }
}

void
ImplicitHeatConduction::destroyHCMatrix(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  if (cout_doing.active()) {
    cout_doing << "Doing destroyHCMatrix "
               << "\t\t\t\t\t IMPM"
               << "\n"
               << "\n";
  }

  d_HC_solver->destroyMatrix(false);
}

void
ImplicitHeatConduction::createHCMatrix(const ProcessorGroup* pg,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  int numMatls = d_mat_manager->getNumMaterials("MPM");

  std::map<int, int> dof_diag;
  d_HC_solver->createLocalToGlobalMapping(
    pg, d_perproc_patches, patches, 1, d_mpm_flags->d_8or27);
  int global_offset = 0;
  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    if (cout_doing.active()) {
      cout_doing << "Doing createHCMatrix on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }

    Array3<int> l2g(lowIndex, highIndex);
    d_HC_solver->copyL2G(l2g, patch);
    // set global offset if this is the first patch
    if (pp == 0) {
      global_offset = l2g[lowIndex];
    }
    CCVariable<int> visited;
    new_dw->allocateTemporary(visited, patch, Ghost::AroundCells, 1);
    visited.initialize(0);

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      constParticleVariable<Point> px;
      ParticleSubset* pset;

      pset = old_dw->getParticleSubset(
        dwi, patch, Ghost::AroundNodes, 1, d_mpm_labels->pXLabel);
      old_dw->get(px, d_mpm_labels->pXLabel, pset);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        IntVector cell, ni[8];
        patch->findCell(px[idx], cell);
        if (visited[cell] == 0) {
          visited[cell] = 1;
          patch->findNodesFromCell(cell, ni);
          std::vector<int> dof(0);
          int l2g_node_num;
          for (int k = 0; k < 8; k++) {
            if (patch->containsNode(ni[k])) {
              l2g_node_num =
                l2g[ni[k]] - global_offset; // subtract global offset to map
                                            // into array correctly
              dof.push_back(l2g_node_num);
            }
          }
          for (int I = 0; I < (int)dof.size(); I++) {
            int dofi = dof[I];
            for (int J = 0; J < (int)dof.size(); J++) {
              dof_diag[dofi] += 1;
            }
          }
        }
      }
    }
  }
  d_HC_solver->createMatrix(pg, dof_diag);
}

void
ImplicitHeatConduction::findNeighbors(IntVector n,
                                      std::vector<int>& neigh,
                                      Array3<int>& l2g)
{
  int dof             = l2g[n];
  IntVector lowIndex  = l2g.getLowIndex();
  IntVector highIndex = l2g.getHighIndex();
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        IntVector neighbor = n + IntVector(i, j, k);
        if (neighbor.x() >= lowIndex.x() && neighbor.x() < highIndex.x() &&
            neighbor.y() >= lowIndex.y() && neighbor.y() < highIndex.y() &&
            neighbor.z() >= lowIndex.z() && neighbor.z() < highIndex.z()) {
          int dof_neighbor = l2g[neighbor];
          if (dof_neighbor > 0 && dof_neighbor != dof) {
#if 0
            std::cout << "neighbor: " << neighbor << " dof: " 
                 << dof_neighbor << std::endl;
#endif
            neigh.push_back(dof_neighbor);
          }
        }
      }
    }
  }
}

void
ImplicitHeatConduction::applyHCBoundaryConditions(const ProcessorGroup*,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset*,
                                                  DataWarehouse* old_dw,
                                                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing applyHCBoundaryConditions "
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }
    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    // Apply grid boundary conditions to the temperature before storing the
    // data -- all data is for material 0 -- single temperature field.

    d_HC_solver->copyL2G(l2g, patch);
    int matl = 0;

    NCVariable<double> gtemp;

    new_dw->allocateAndPut(
      gtemp, d_mpm_labels->gTemperatureStarLabel, matl, patch);
    gtemp.initialize(0.);

    for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
         face                 = Patch::nextFace(face)) {

      if (patch->getBCType(face) == Patch::None) {
        int numChildren = patch->getBCDataArray(face)->getNumberChildren(matl);
        for (int child = 0; child < numChildren; child++) {
          Iterator nbound_ptr;
          Iterator nu; // not used

          BoundCondBaseSP temp_bcs = patch->getArrayBCValues(
            face, matl, "Temperature", nu, nbound_ptr, child);
          BoundCond<double>::BoundCondP bc =
            std::dynamic_pointer_cast<BoundCond<double>>(temp_bcs);
          if (bc != 0) {
            if (bc->getBCType() == "Dirichlet") {
              for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                gtemp[*nbound_ptr] = bc->getValue();
              }
              IntVector l, h;
              patch->getFaceNodes(face, 0, l, h);
              for (NodeIterator it(l, h); !it.done(); it++) {
                IntVector n = *it;
                int dof     = l2g[n];
                std::vector<int> neigh;
                // Find the neighbors for this DOF
                findNeighbors(n, neigh, l2g);
                d_HC_solver->d_DOFNeighbors[dof] = neigh;
                d_HC_solver->d_DOF.insert(dof);
              }
            }
          }
        }
      } else {
        continue;
      }
    } // faces

    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matlindex = mpm_matl->getDWIndex();

      constNCVariable<double> gheatflux;
      new_dw->get(gheatflux,
                  d_mpm_labels->gExternalHeatFluxLabel,
                  matlindex,
                  patch,
                  Ghost::None,
                  0);

      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector n = *iter;
        int dof     = l2g[n];
        if (!compare(gheatflux[n], 0.)) {
          d_HC_solver->d_DOFFlux.insert(dof);
          d_HC_solver->fillFluxVector(dof, gheatflux[n]);
        }
      }
    }
  } // patches
}

void
ImplicitHeatConduction::findFixedHCDOF(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing findFixedHCDOF on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    size_t numMatls = d_mat_manager->getNumMaterials("MPM");
    d_HC_solver->copyL2G(l2g, patch);
    NCVariable<double> GMASS;
    new_dw->allocateTemporary(GMASS, patch, Ghost::None, 0);
    GMASS.initialize(0.);

    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matlindex = mpm_matl->getDWIndex();
      constNCVariable<double> gMass;
      new_dw->get(
        gMass, d_mpm_labels->gMassLabel, matlindex, patch, Ghost::None, 0);

      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector n = *iter;
        GMASS[n] += gMass[n];
      }
    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector n = *iter;
      int dof     = l2g[n];
      if (compare(GMASS[n], 0.)) {
#if 0
        std::vector<int> neigh;
        findNeighbors(n,neigh,l2g);
        d_HC_solver->d_DOFNeighbors[dof] = neigh;
#endif
        d_HC_solver->d_DOFZero.insert(dof);
      }
    } // node iterator
  }   // patches
}

void
ImplicitHeatConduction::formHCStiffnessMatrix(const ProcessorGroup*,
                                              const PatchSubset* patches,
                                              const MaterialSubset*,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing formHCStiffnessMatrix " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    Vector dx = patch->dCell();

    LinearInterpolator* interpolator = scinew LinearInterpolator(patch);

    d_HC_solver->copyL2G(l2g, patch);
    size_t numMatls = d_mat_manager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matlindex = mpm_matl->getDWIndex();

      delt_vartype dt;
      old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());

      ParticleSubset* pset;
      pset = old_dw->getParticleSubset(matlindex, patch);
      constParticleVariable<Point> px;
      constParticleVariable<double> pVolume;
      constParticleVariable<double> ptemperature;

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pVolume, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(ptemperature, d_mpm_labels->pTemperatureLabel, pset);

      double v[64];
      double kHC[8][8];
      int dof[8];
      double K   = mpm_matl->getThermalConductivity();
      double Cp  = mpm_matl->getSpecificHeat();
      double rho = mpm_matl->getInitialDensity();

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        // Get the node indices that surround the cell
        std::vector<IntVector> ni(8);
        std::vector<double> S(8);
        std::vector<Vector> d_S(8);

        interpolator->findCellAndWeightsAndShapeDerivatives(
          px[idx], ni, S, d_S);

        for (int ii = 0; ii < 8; ii++) {
          for (int jj = 0; jj < 8; jj++) {
            kHC[ii][jj] = 0.;
          }
        }

        for (int ii = 0; ii < 8; ii++) {
          int l2g_node_num = l2g[ni[ii]];
          dof[ii]          = l2g_node_num;
          for (int jj = 0; jj < 8; jj++) {
            // Thermal inertia terms
            if (d_HC_transient) {
              kHC[ii][jj] += S[jj] * S[ii] * Cp * rho * pVolume[idx] / dt;
            }
            for (int dir = 0; dir < 3; dir++) {
              kHC[ii][jj] += d_S[jj][dir] * d_S[ii][dir] *
                             (K / (dx[dir] * dx[dir])) * pVolume[idx];
            }
          }
        }

        for (int I = 0; I < 8; I++) {
          for (int J = 0; J < 8; J++) {
            v[8 * I + J] = kHC[I][J];
          }
        }
        d_HC_solver->fillMatrix(8, dof, 8, dof, v);
      }
    } // matls
    // delete interpolator;
  }
}

void
ImplicitHeatConduction::formHCQ(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing formHCQ on patch " << patch->getID()
                 << "\t\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    LinearInterpolator* interpolator = scinew LinearInterpolator(patch);
    d_HC_solver->copyL2G(l2g, patch);

    constNCVariable<double> temperature;
    new_dw->get(temperature,
                d_mpm_labels->gTemperatureLabel,
                0,
                patch,
                Ghost::AroundCells,
                1);

#if 0
    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++){
      IntVector n = *iter;
      std::cout << "temperature[" << n << "]= " << temperature[n] << std::endl;
    }
#endif

    int numMatls = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      delt_vartype dt;

      constNCVariable<double> gextheatrate;

      old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());

      new_dw->get(gextheatrate,
                  d_mpm_labels->gExternalHeatRateLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  1);

      ParticleSubset* pset;
      pset = old_dw->getParticleSubset(dwi, patch);
      constParticleVariable<Point> px;
      constParticleVariable<double> pVolume;

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pVolume, d_mpm_labels->pVolumeLabel, pset);

      int dof[8];
      double Cp  = mpm_matl->getSpecificHeat();
      double rho = mpm_matl->getInitialDensity();

      std::vector<IntVector> ni(8);
      std::vector<double> S(8);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        // Get the node indices that surround the cell

        interpolator->findCellAndWeights(px[idx], ni, S);

        bool add = true;

        for (int ii = 0; ii < 8; ii++) {
          //  std::cout << "ni = " << ni[ii] << std::endl;
          int l2g_node_num = l2g[ni[ii]];
          double v         = 0;
          dof[ii]          = l2g_node_num;
          for (int jj = 0; jj < 8; jj++) {
            // Thermal inertia terms
            if (d_HC_transient) {
              v += S[jj] * S[ii] * Cp * rho * pVolume[idx] *
                     temperature[ni[jj]] / dt +
                   gextheatrate[ni[ii]];
            }
          }
          // std::cout << "v[" << l2g_node_num << "]= " << v << std::endl;
          d_HC_solver->fillVector(dof[ii], v, add);
        }
      }
    } // matls
    // delete interpolator;
  } // patches
}

void
ImplicitHeatConduction::adjustHCQAndHCKForBCs(const ProcessorGroup*,
                                              const PatchSubset* patches,
                                              const MaterialSubset*,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing adjustHCQAndHCKForBCs on patch " << patch->getID()
                 << "\t\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }
    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    d_HC_solver->copyL2G(l2g, patch);
    int dwi = 0;

    Ghost::GhostType gnone = Ghost::None;

    constNCVariable<double> temperature;

    // gTemperatureStar are the boundary condition values for the temperature
    // determined earlier in applyHCBoundaryConditions()

    new_dw->get(
      temperature, d_mpm_labels->gTemperatureStarLabel, dwi, patch, gnone, 0);

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector n = *iter;
      int dof[1];
      dof[0]   = l2g[n];
      double v = -temperature[n];
      d_HC_solver->fillTemporaryVector(dof[0], v);
    }
  } // patches
  d_HC_solver->assembleTemporaryVector();
  d_HC_solver->assembleVector();
  d_HC_solver->finalizeMatrix();
  d_HC_solver->applyBCSToRHS();
  d_HC_solver->removeFixedDOFHeat();
}

void
ImplicitHeatConduction::solveForTemp(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  std::vector<double> guess;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing solveForTemp on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

#if 0
    IntVector lowIndex = patch->getNodeLowIndex();
    IntVector highIndex = patch->getNodeHighIndex()+IntVector(1,1,1);
    Array3<int> l2g(lowIndex,highIndex);

    constNCVariable<double> temperature;
    int dwi = 0;
    new_dw->get(temperature,d_mpm_labels->gTemperatureLabel,dwi,patch,Ghost::None,0);

    for (NodeIterator iter = patch->getNodeIterator();!iter.done();iter++){
      guess.push_back(temperature[*iter]);
    }
#endif
  }

  d_HC_solver->solve(guess);
}

void
ImplicitHeatConduction::getTemperatureIncrement(const ProcessorGroup*,
                                                const PatchSubset* patches,
                                                const MaterialSubset*,
                                                DataWarehouse* old_dw,
                                                DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing getTemperatureIncrement on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    delt_vartype dt;
    old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());

    IntVector lowIndex, highIndex;
    if (d_mpm_flags->d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_mpm_flags->d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }
    Array3<int> l2g(lowIndex, highIndex);

    NCVariable<double> tempRate;
    constNCVariable<double> temp;

    int dwi = 0;

    new_dw->allocateAndPut(
      tempRate, d_mpm_labels->gTemperatureRateLabel, dwi, patch);
    new_dw->get(
      temp, d_mpm_labels->gTemperatureLabel, dwi, patch, Ghost::None, 0);

    tempRate.initialize(0.0);

    d_HC_solver->copyL2G(l2g, patch);

    std::vector<double> x;
    int begin = d_HC_solver->getSolution(x);

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector n = *iter;
      int dof     = l2g[n] - begin;
      tempRate[n] = (x[dof] - temp[n]) / dt;
    }
  }
}

void
ImplicitHeatConduction::fillgTemperatureRate(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset*,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing fillgTemperatureRate on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    int dwi = 0;
    NCVariable<double> tempRate;
    new_dw->allocateAndPut(
      tempRate, d_mpm_labels->gTemperatureRateLabel, dwi, patch);
    tempRate.initialize(0.0);
  }
}
