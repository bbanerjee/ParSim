/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>

#include <CCA/Components/MPM/Core/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Short27.h>
#include <Core/Util/DebugStream.h>

using namespace Uintah;

#define EROSION
#undef EROSION

static DebugStream cout_doing("HeatConduction", false);
static DebugStream cout_heat("MPMHeat", false);

HeatConduction::HeatConduction(const MaterialManagerP& mat_manager,
                               const MPMLabel* labels,
                               const MPMFlags* flags)
  : d_mat_manager{ mat_manager }
  , d_mpm_labels{ labels }
  , d_mpm_flags{ flags }
{
  if (flags->d_8or27 == 8) {
    d_num_ghost_particles = 1;
    d_num_ghost_nodes     = 1;
  } else {
    d_num_ghost_particles = 2;
    d_num_ghost_nodes     = 2;
  }
}

void
HeatConduction::scheduleComputeInternalHeatRate(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  Task* t = scinew Task("MPM::computeInternalHeatRate",
                        this,
                        &HeatConduction::computeInternalHeatRate);

  Ghost::GhostType gan   = Ghost::AroundNodes;
  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::OldDW, d_mpm_labels->pXLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pSizeLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pMassLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pVolumeLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pDefGradLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::NewDW, d_mpm_labels->gTemperatureLabel, gan, 2 * d_num_ghost_nodes);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, gnone);
  t->computes(d_mpm_labels->gdTdtLabel);

  if (d_mpm_flags->d_fracture) { // for FractureMPM
    t->requires(
      Task::NewDW, d_mpm_labels->pgCodeLabel, gan, d_num_ghost_particles);
    t->requires(
      Task::NewDW, d_mpm_labels->GTemperatureLabel, gac, 2 * d_num_ghost_nodes);
    t->requires(Task::NewDW, d_mpm_labels->GMassLabel, gnone);
    t->computes(d_mpm_labels->GdTdtLabel);
  }

  sched->addTask(t, patches, matls);
}
//__________________________________
//
void
HeatConduction::scheduleComputeNodalHeatFlux(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (d_mpm_flags->d_computeNodalHeatFlux == false) {
    return;
  }

  // This task only exists to compute the diagnostic gHeatFluxLabel
  // which is not used in any of the subsequent calculations

  Task* t = scinew Task(
    "MPM::computeNodalHeatFlux", this, &HeatConduction::computeNodalHeatFlux);

  Ghost::GhostType gan   = Ghost::AroundNodes;
  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::OldDW, d_mpm_labels->pXLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pSizeLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pDefGradLabel, gan, d_num_ghost_particles);
  t->requires(
    Task::OldDW, d_mpm_labels->pMassLabel, gan, d_num_ghost_particles);
  t->requires(Task::NewDW,
              d_mpm_labels->gTemperatureLabel,
              gac,
              2 * d_num_ghost_particles);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, gnone);
  t->computes(d_mpm_labels->gHeatFluxLabel);

  sched->addTask(t, patches, matls);
}

void
HeatConduction::scheduleSolveHeatEquations(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  /* solveHeatEquations
   *   in(G.MASS, G.INTERNALHEATRATE, G.EXTERNALHEATRATE)
   *   out(G.TEMPERATURERATE) */

  Task* t = scinew Task(
    "MPM::solveHeatEquations", this, &HeatConduction::solveHeatEquations);

  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, gnone);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, gnone);
  t->requires(Task::NewDW, d_mpm_labels->gExternalHeatRateLabel, gnone);
  t->requires(Task::NewDW, d_mpm_labels->gdTdtLabel, gnone);
  t->requires(
    Task::NewDW, d_mpm_labels->gThermalContactTemperatureRateLabel, gnone);
  t->modifies(d_mpm_labels->gTemperatureRateLabel);

  if (d_mpm_flags->d_fracture) { // for FractureMPM
    t->requires(Task::NewDW, d_mpm_labels->GMassLabel, gnone);
    t->requires(Task::NewDW, d_mpm_labels->GVolumeLabel, gnone);
    t->requires(Task::NewDW, d_mpm_labels->GExternalHeatRateLabel, gnone);
    t->requires(Task::NewDW, d_mpm_labels->GdTdtLabel, gnone);
    t->requires(
      Task::NewDW, d_mpm_labels->GThermalContactTemperatureRateLabel, gnone);
    t->computes(d_mpm_labels->GTemperatureRateLabel);
  }

  sched->addTask(t, patches, matls);
}

void
HeatConduction::scheduleIntegrateTemperatureRate(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  /* integrateTemperatureRate
   *   in(G.TEMPERATURE, G.TEMPERATURERATE)
   *   operation(t* = t + t_rate * dt)
   *   out(G.TEMPERATURE_STAR) */

  Task* t = scinew Task("MPM::integrateTemperatureRate",
                        this,
                        &HeatConduction::integrateTemperatureRate);

  const MaterialSubset* mss = matls->getUnion();

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);

  t->requires(Task::NewDW, d_mpm_labels->gTemperatureLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gTemperatureNoBCLabel, Ghost::None);
  t->modifies(d_mpm_labels->gTemperatureRateLabel, mss);
  t->computes(d_mpm_labels->gTemperatureStarLabel);

  if (d_mpm_flags->d_fracture) { // for FractureMPM
    t->requires(Task::NewDW, d_mpm_labels->GTemperatureLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpm_labels->GTemperatureNoBCLabel, Ghost::None);
    t->modifies(d_mpm_labels->GTemperatureRateLabel, mss);
    t->computes(d_mpm_labels->GTemperatureStarLabel);
  }

  sched->addTask(t, patches, matls);
}

void
HeatConduction::computeInternalHeatRate(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing computeInternalHeatRate on patch " << patch->getID()
                 << "\t\t MPM" << std::endl;
    }
    if (cout_heat.active()) {
      cout_heat << " Patch = " << patch->getID() << std::endl;
    }

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();

    Ghost::GhostType gac   = Ghost::AroundCells;
    Ghost::GhostType gnone = Ghost::None;
    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

      if (cout_heat.active()) {
        cout_heat << "  Material = " << m << std::endl;
      }

      int dwi      = mpm_matl->getDWIndex();
      double kappa = mpm_matl->getThermalConductivity();
      double Cv    = mpm_matl->getSpecificHeat();

      constParticleVariable<Point> px;
      constParticleVariable<double> pvol, pMass;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      ParticleVariable<Vector> pTemperatureGradient;
      constNCVariable<double> gTemperature, gMass;
      NCVariable<double> gdTdt;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_num_ghost_particles,
                                                       d_mpm_labels->pXLabel);

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pvol, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
      new_dw->get(gTemperature,
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  gac,
                  2 * d_num_ghost_nodes);
      new_dw->get(gMass, d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->allocateAndPut(gdTdt, d_mpm_labels->gdTdtLabel, dwi, patch);
      new_dw->allocateTemporary(pTemperatureGradient, pset);

      gdTdt.initialize(0.);

      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      constNCVariable<double> GTemperature;
      constNCVariable<double> GMass;
      NCVariable<double> GdTdt;
      if (d_mpm_flags->d_fracture) {
        new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);
        new_dw->get(GTemperature,
                    d_mpm_labels->GTemperatureLabel,
                    dwi,
                    patch,
                    gac,
                    2 * d_num_ghost_nodes);
        new_dw->get(GMass, d_mpm_labels->GMassLabel, dwi, patch, gnone, 0);
        new_dw->allocateAndPut(GdTdt, d_mpm_labels->GdTdtLabel, dwi, patch);
        GdTdt.initialize(0.);
      }

      // Compute the temperature gradient at each particle and project
      // the particle plastic work temperature rate to the grid
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndShapeDerivatives(
          px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

        pTemperatureGradient[idx] = Vector(0.0, 0.0, 0.0);
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          for (int j = 0; j < 3; j++) {
            pTemperatureGradient[idx][j] +=
              gTemperature[ni[k]] * d_S[k][j] * oodx[j];

            if (cout_heat.active()) {
              cout_heat << "   node = " << ni[k]
                        << " gTemp = " << gTemperature[ni[k]]
                        << " idx = " << idx
                        << " pTempGrad = " << pTemperatureGradient[idx][j]
                        << std::endl;
            }
          }
          // Project the mass weighted particle plastic work temperature
          // rate to the grid
        } // Loop over local nodes
      }   // Loop over particles

      if (d_mpm_flags->d_fracture) { // for FractureMPM
        // Compute the temperature gradient at each particle and project
        // the particle plastic work temperature rate to the grid
        for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
             iter++) {
          particleIndex idx = *iter;

          // Get the node indices that surround the cell
          interpolator->findCellAndShapeDerivatives(
            px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

          pTemperatureGradient[idx] = Vector(0.0, 0.0, 0.0);
          for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
            for (int j = 0; j < 3; j++) {
              if (pgCode[idx][k] == 1) { // above crack
                pTemperatureGradient[idx][j] +=
                  gTemperature[ni[k]] * d_S[k][j] * oodx[j];
              } else if (pgCode[idx][k] == 2) { // below crack
                pTemperatureGradient[idx][j] +=
                  GTemperature[ni[k]] * d_S[k][j] * oodx[j];
              }
            }
          } // Loop over local nodes
        }   // Loop over particles
      }     // if fracture

      // Compute rate of temperature change at the grid due to conduction
      // and plastic work
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndShapeDerivatives(
          px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

        // Calculate k/(rho*Cv)
        double alpha     = kappa * pvol[idx] / Cv;
        Vector dT_dx     = pTemperatureGradient[idx];
        double Tdot_cond = 0.0;
        IntVector node(0, 0, 0);

        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          node = ni[k];
          if (patch->containsNode(node)) {
            Vector div(
              d_S[k].x() * oodx[0], d_S[k].y() * oodx[1], d_S[k].z() * oodx[2]);
            Tdot_cond = Dot(div, dT_dx) * (alpha / gMass[node]);
            gdTdt[node] -= Tdot_cond;

            if (cout_heat.active()) {
              cout_heat << "   node = " << node << " div = " << div
                        << " dT_dx = " << dT_dx << " alpha = " << alpha * Cv
                        << " Tdot_cond = " << Tdot_cond * Cv * gMass[node]
                        << " gdTdt = " << gdTdt[node] << std::endl;
            } // cout_heat
          }   // if patch contains node
        }     // Loop over local nodes

        if (d_mpm_flags->d_fracture) { // for FractureMPM
          for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
            node = ni[k];
            if (patch->containsNode(node)) {
              Vector div(d_S[k].x() * oodx[0],
                         d_S[k].y() * oodx[1],
                         d_S[k].z() * oodx[2]);
              if (pgCode[idx][k] == 1) { // above crack
                Tdot_cond = Dot(div, dT_dx) * (alpha / gMass[node]);
                gdTdt[node] -= Tdot_cond;
              } else if (pgCode[idx][k] == 2) { // below crack
                Tdot_cond = Dot(div, dT_dx) * (alpha / GMass[node]);
                GdTdt[node] -= Tdot_cond;
              }
            } // if patch contains node
          }   // Loop over local nodes
        }

      } // Loop over particles
    }   // End of loop over materials
    // delete interpolator;
  } // End of loop over patches
}
//______________________________________________________________________
//
void
HeatConduction::computeNodalHeatFlux(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    // This task only exists to compute the diagnostic gHeatFluxLabel
    // which is not used in any of the subsequent calculations

    if (cout_doing.active()) {
      cout_doing << "Doing computeNodalHeatFlux on patch " << patch->getID()
                 << "\t\t MPM" << std::endl;
    }
    if (cout_heat.active()) {
      cout_heat << " Patch = " << patch->getID() << std::endl;
    }

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();

    Ghost::GhostType gac   = Ghost::AroundCells;
    Ghost::GhostType gnone = Ghost::None;

    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

      if (cout_heat.active()) {
        cout_heat << "  Material = " << m << std::endl;
      }

      int dwi      = mpm_matl->getDWIndex();
      double kappa = mpm_matl->getThermalConductivity();

      NCVariable<Vector> gHeatFlux;
      constNCVariable<double> gTemperature, gMass;
      constParticleVariable<Point> px;
      constParticleVariable<double> pMass;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_num_ghost_particles,
                                                       d_mpm_labels->pXLabel);

      new_dw->get(gTemperature,
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  gac,
                  2 * d_num_ghost_nodes);
      new_dw->get(gMass, d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

      new_dw->allocateAndPut(
        gHeatFlux, d_mpm_labels->gHeatFluxLabel, dwi, patch);
      gHeatFlux.initialize(Vector(0.0));

      //__________________________________
      // Create a temporary variables for the mass weighted nodal
      // temperature gradient
      NCVariable<Vector> gpdTdx;
      ParticleVariable<Vector> pdTdx;
      new_dw->allocateTemporary(gpdTdx, patch, gnone, 0);
      new_dw->allocateTemporary(pdTdx, pset);

      gpdTdx.initialize(Vector(0., 0., 0.));

      // Compute the temperature gradient at each particle
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        pdTdx[idx]        = Vector(0, 0, 0);

        interpolator->findCellAndShapeDerivatives(
          px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          for (int j = 0; j < 3; j++) {
            pdTdx[idx][j] += gTemperature[ni[k]] * d_S[k][j] * oodx[j];
          }
        }
      } // particles

      // project the mass weighted particle temperature gradient to the grid
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(
          px[idx], ni, S, pSize[idx], pDefGrad[idx]);

        Vector pdTdx_massWt = pdTdx[idx] * pMass[idx];

        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          if (patch->containsNode(ni[k])) {
            gpdTdx[ni[k]] += (pdTdx_massWt * S[k]);
          }
        }
      } // particles

      // compute the nodal temperature gradient by dividing
      // gpdTdx by the grid mass
      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector n  = *iter;
        gHeatFlux[n] = -kappa * gpdTdx[n] / gMass[n];
      }
    } // End of loop over materials
    // delete interpolator;
  } // End of loop over patches
}

void
HeatConduction::solveHeatEquations(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* /*old_dw*/,
                                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing solveHeatEquations on patch " << patch->getID()
                 << "\t\t\t MPM" << std::endl;
    }

    string interp_type = d_mpm_flags->d_interpolatorType;
    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi   = mpm_matl->getDWIndex();
      double Cv = mpm_matl->getSpecificHeat();

      // Get required variables for this patch
      constNCVariable<double> mass, externalHeatRate, gVolume;
      constNCVariable<double> thermalContactTemperatureRate, gdTdt;

      new_dw->get(mass, d_mpm_labels->gMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(
        gVolume, d_mpm_labels->gVolumeLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(externalHeatRate,
                  d_mpm_labels->gExternalHeatRateLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(gdTdt, d_mpm_labels->gdTdtLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(thermalContactTemperatureRate,
                  d_mpm_labels->gThermalContactTemperatureRateLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      // for FractureMPM
      constNCVariable<double> Gmass, GexternalHeatRate, Gvolume;
      constNCVariable<double> GthermalContactTemperatureRate, GdTdt;
      if (d_mpm_flags->d_fracture) {
        new_dw->get(
          Gmass, d_mpm_labels->GMassLabel, dwi, patch, Ghost::None, 0);
        new_dw->get(
          Gvolume, d_mpm_labels->GVolumeLabel, dwi, patch, Ghost::None, 0);
        new_dw->get(GexternalHeatRate,
                    d_mpm_labels->GExternalHeatRateLabel,
                    dwi,
                    patch,
                    Ghost::None,
                    0);
        new_dw->get(
          GdTdt, d_mpm_labels->GdTdtLabel, dwi, patch, Ghost::None, 0);
        new_dw->get(GthermalContactTemperatureRate,
                    d_mpm_labels->GThermalContactTemperatureRateLabel,
                    dwi,
                    patch,
                    Ghost::None,
                    0);
      }

      // Create variables for the results
      NCVariable<double> tempRate, GtempRate;
      new_dw->getModifiable(
        tempRate, d_mpm_labels->gTemperatureRateLabel, dwi, patch);

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        tempRate[c] = gdTdt[c] * ((mass[c] - 1.e-200) / mass[c]) +
                      (externalHeatRate[c]) / (mass[c] * Cv) +
                      thermalContactTemperatureRate[c];
      } // End of loop over iter

      if (d_mpm_flags->d_fracture) { // for FractureMPM
        new_dw->allocateAndPut(
          GtempRate, d_mpm_labels->GTemperatureRateLabel, dwi, patch);
        GtempRate.initialize(0.0);
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c  = *iter;
          GtempRate[c] = GdTdt[c] * ((Gmass[c] - 1.e-200) / Gmass[c]) +
                         (GexternalHeatRate[c]) / (Gmass[c] * Cv) +
                         GthermalContactTemperatureRate[c];
        } // End of loop over iter
      }
    }
  }
}

void
HeatConduction::integrateTemperatureRate(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset*,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing integrateTemperatureRate on patch " << patch->getID()
                 << "\t\t MPM" << std::endl;
    }

    Ghost::GhostType gnone = Ghost::None;
    string interp_type     = d_mpm_flags->d_interpolatorType;
    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      constNCVariable<double> temp_old, temp_oldNoBC;
      NCVariable<double> temp_rate, tempStar;
      delt_vartype delT;
      old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

      new_dw->get(
        temp_old, d_mpm_labels->gTemperatureLabel, dwi, patch, gnone, 0);
      new_dw->get(temp_oldNoBC,
                  d_mpm_labels->gTemperatureNoBCLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(
        temp_rate, d_mpm_labels->gTemperatureRateLabel, dwi, patch);
      new_dw->allocateAndPut(
        tempStar, d_mpm_labels->gTemperatureStarLabel, dwi, patch);
      tempStar.initialize(0.0);

      // for FractureMPM
      constNCVariable<double> Gtemp_old, Gtemp_oldNoBC;
      NCVariable<double> Gtemp_rate, GtempStar;
      if (d_mpm_flags->d_fracture) {
        new_dw->get(
          Gtemp_old, d_mpm_labels->GTemperatureLabel, dwi, patch, gnone, 0);
        new_dw->get(Gtemp_oldNoBC,
                    d_mpm_labels->GTemperatureNoBCLabel,
                    dwi,
                    patch,
                    gnone,
                    0);
        new_dw->getModifiable(
          Gtemp_rate, d_mpm_labels->GTemperatureRateLabel, dwi, patch);
        new_dw->allocateAndPut(
          GtempStar, d_mpm_labels->GTemperatureStarLabel, dwi, patch);
        GtempStar.initialize(0.0);
      }

      MPMBoundCond bc;

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        tempStar[c] = temp_old[c] + temp_rate[c] * delT;
        /*
        if (c == IntVector(24,47,28)) {
          std::cout << "c = " << c << " tempStar = " << tempStar[c]
                    << " temp_old = " << temp_old[c]
                    << " temp_rate = " << temp_rate[c] << std::endl;
        }
        */
      }
      // Apply grid boundary conditions to the temperature
      bc.setBoundaryCondition(patch, dwi, "Temperature", tempStar, interp_type);

      if (d_mpm_flags->d_fracture) { // for FractureMPM
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c  = *iter;
          GtempStar[c] = Gtemp_old[c] + Gtemp_rate[c] * delT;
        }
        // Apply grid boundary conditions to the temperature
        bc.setBoundaryCondition(
          patch, dwi, "Temperature", GtempStar, interp_type);
      }

      // Now recompute temp_rate as the difference between the temperature
      // interpolated to the grid (no bcs applied) and the new tempStar
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c  = *iter;
        temp_rate[c] = (tempStar[c] - temp_oldNoBC[c]) / delT;
        /*
        if (c == IntVector(24,47,28)) {
          std::cout << "c = " << c << " tempStar = " << tempStar[c]
                    << " temp_oldNoBC = " << temp_oldNoBC[c]
                    << " temp_rate = " << temp_rate[c] << std::endl;
        }
        */
      }

      if (d_mpm_flags->d_fracture) { // for FractureMPM
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c   = *iter;
          Gtemp_rate[c] = (GtempStar[c] - Gtemp_oldNoBC[c]) / delT;
        }
      } // fracture
    }   // matls
  }     // patches
}
