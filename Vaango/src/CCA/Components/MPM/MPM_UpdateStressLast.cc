/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Parresia Research Limited
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

#include <CCA/Components/MPM/MPM_UpdateStressLast.h>

#include <CCA/Components/MPM/CohesiveZone/CohesiveZoneTasks.h>
#include <CCA/Components/MPM/HeatConduction/HeatConductionTasks.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>

#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Util/DebugStream.h>

#define DEBUG_WITH_PARTICLE_ID

using namespace Uintah;

static DebugStream cout_doing("MPM_USL", false);
static DebugStream cout_heat("MPM_USL_Heat", false);

MPM_UpdateStressLast::MPM_UpdateStressLast(const ProcessorGroup* myworld,
                                           const MaterialManagerP& mat_manager)
  : SerialMPM(myworld, mat_manager)
{
}

void
MPM_UpdateStressLast::scheduleTimeAdvance(const LevelP& level,
                                          SchedulerP& sched)
{
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  const PatchSet* patches  = level->eachPatch();
  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  scheduleComputeParticleBodyForce(sched, patches, matls);
  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleInterpolateParticlesToGrid(sched, patches, matls);
  scheduleMomentumExchangeInterpolated(sched, patches, matls);
  d_cohesiveZoneTasks->scheduleUpdate(sched, patches, matls);
  if (d_boundaryTractionFaces.size() > 0) {
    scheduleComputeContactArea(sched, patches, matls);
  }
  scheduleComputeInternalForce(sched, patches, matls);
  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  scheduleMomentumExchangeIntegrated(sched, patches, matls);
  scheduleSetGridBoundaryConditions(sched, patches, matls);
  if (d_mpm_flags->d_prescribeDeformation) {
    scheduleSetPrescribedMotion(sched, patches, matls);
  }
  if (d_mpm_flags->d_useXPIC) {
    scheduleComputeXPICVelocities(sched, patches, matls);
  }
  d_heatConductionTasks->scheduleCompute(sched, patches, matls);
  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  scheduleComputeDeformationGradient(sched, patches, matls);
  scheduleComputeStressTensor(sched, patches, matls);
  scheduleComputeBasicDamage(sched, patches, matls);
  scheduleUpdateErosionParameter(sched, patches, matls);
  scheduleFindRogueParticles(sched, patches, matls);
  scheduleInsertParticles(sched, patches, matls);
  if (d_mpm_flags->d_refineParticles) {
    scheduleAddParticles(sched, patches, matls);
  }
  if (d_mpm_flags->d_computeScaleFactor) {
    scheduleComputeParticleScaleFactor(sched, patches, matls);
  }
  scheduleParticleRelocation(sched, level, patches, matls);
}

void
MPM_UpdateStressLast::scheduleInterpolateToParticlesAndUpdate(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM_USL::scheduleInterpolateToParticlesAndUpdate");

  Task* t = scinew Task("MPM_USL::interpolateToParticlesAndUpdate",
                        this,
                        &MPM_UpdateStressLast::interpolateToParticlesAndUpdate);

  t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->needs(Task::NewDW,
              d_mpm_labels->gAccelerationLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gVelocityStarLabel,
              gac,
              d_numGhostNodes);
  if (d_mpm_flags->d_useXPIC) {
    t->needs(Task::NewDW,
                d_mpm_labels->gVelocityXPICLabel,
                gac,
                d_numGhostNodes);
    t->needs(Task::OldDW, d_mpm_labels->pVelocityXPICLabel, gnone);
  }
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureRateLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->frictionalWorkLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pTemperatureLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pDispLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pVolumeLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pDefGradLabel, gnone);
  if (d_mpm_flags->d_useLoadCurves) {
    t->needs(Task::OldDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
  }
  if (d_mpm_flags->d_withICE) {
    t->needs(Task::NewDW, d_mpm_labels->dTdt_NCLabel, gac, d_numGhostNodes);
    t->needs(Task::NewDW,
                d_mpm_labels->massBurnFractionLabel,
                gac,
                d_numGhostNodes);
  }

  t->computes(d_mpm_labels->pDispLabel_preReloc);
  t->computes(d_mpm_labels->pVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pXLabel_preReloc);
  t->computes(d_mpm_labels->pParticleIDLabel_preReloc);
  t->computes(d_mpm_labels->pTemperatureLabel_preReloc);
  t->computes(d_mpm_labels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpm_labels->pMassLabel_preReloc);
  t->computes(d_mpm_labels->pSizeLabel_preReloc);

  //__________________________________
  //  reduction variables
  if (d_mpm_flags->d_reductionVars->momentum) {
    t->computes(d_mpm_labels->TotalMomentumLabel);
  }
  if (d_mpm_flags->d_reductionVars->KE) {
    t->computes(d_mpm_labels->KineticEnergyLabel);
  }
  if (d_mpm_flags->d_reductionVars->thermalEnergy) {
    t->computes(d_mpm_labels->ThermalEnergyLabel);
  }
  if (d_mpm_flags->d_reductionVars->centerOfMass) {
    t->computes(d_mpm_labels->CenterOfMassPositionLabel);
  }
  if (d_mpm_flags->d_reductionVars->mass) {
    t->computes(d_mpm_labels->TotalMassLabel);
  }
  if (d_mpm_flags->d_withColor) {
    t->needs(Task::OldDW, d_mpm_labels->pColorLabel, Ghost::None);
    t->computes(d_mpm_labels->pColorLabel_preReloc);
  }

  // Carry Forward particle refinement flag
  if (d_mpm_flags->d_refineParticles) {
    t->needs(Task::OldDW, d_mpm_labels->pRefinedLabel, Ghost::None);
    t->computes(d_mpm_labels->pRefinedLabel_preReloc);
  }

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->needs(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(d_mpm_labels->NC_CCweightLabel, z_matl);

  sched->addTask(t, patches, matls);

  // The task will have a reference to z_matl
  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}

void
MPM_UpdateStressLast::interpolateToParticlesAndUpdate(
  const ProcessorGroup*,
  const PatchSubset* patches,
  const MaterialSubset*,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac   = Ghost::AroundCells;

  for (auto patch : *patches) {
    printTask(patches,
              patch,
              cout_doing,
              "Doing interpolateToParticlesAndUpdate");

    auto interpolator        = d_mpm_flags->d_interpolator->clone(patch);
    auto num_influence_nodes = interpolator->size();
    std::vector<IntVector> ni(num_influence_nodes);
    std::vector<double> S(num_influence_nodes);

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    double totalmass      = 0;
    double partvoldef     = 0.;
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke = 0;

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    /*
    Material* reactant;
    reactant = d_materialManager->getMaterialByName("reactant");
    bool combustion_problem=false;
    int RMI = -99;
    if(reactant != 0){
      RMI = reactant->getDWIndex();
      combustion_problem=true;
    }
    */

    double move_particles = 1.;
    if (!d_mpm_flags->d_doGridReset) {
      move_particles = 0.;
    }

    // Copy NC_CCweight (only material 0)
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_new;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NC_CCweight_new,
                           d_mpm_labels->NC_CCweightLabel,
                           0,
                           patch);
    NC_CCweight_new.copyData(NC_CCweight);

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Copy particle IDs
      constParticleVariable<long64> pParticleID;
      ParticleVariable<long64> pParticleID_new;
      old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pParticleID_new,
                             d_mpm_labels->pParticleIDLabel_preReloc,
                             pset);
      pParticleID_new.copyData(pParticleID);

      // Copy color
      if (d_mpm_flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_mpm_labels->pColorLabel, pset);
        new_dw->allocateAndPut(pColor_new,
                               d_mpm_labels->pColorLabel_preReloc,
                               pset);
        pColor_new.copyData(pColor);
      }

      // Copy refined
      if (d_mpm_flags->d_refineParticles) {
        constParticleVariable<int> pRefined;
        ParticleVariable<int> pRefined_new;
        old_dw->get(pRefined, d_mpm_labels->pRefinedLabel, pset);
        new_dw->allocateAndPut(pRefined_new,
                               d_mpm_labels->pRefinedLabel_preReloc,
                               pset);
        pRefined_new.copyData(pRefined);
      }

      // Get particle variables
      constParticleVariable<double> pMass, pTemperature, pVolume;
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pVolume, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpm_labels->pXLabel, pset);

      constParticleVariable<Vector> pVelocity, pDisp;
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

      constParticleVariable<Vector> pVelocityXPIC;
      if (d_mpm_flags->d_useXPIC) {
        new_dw->get(pVelocityXPIC, d_mpm_labels->pVelocityXPICLabel, pset);
      }

      constParticleVariable<Matrix3> pDefGrad, pSize;
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);

      // Allocate updated particle variables
      ParticleVariable<double> pMass_new, pTemp_new, pTempPrev_new;
      new_dw->allocateAndPut(pMass_new, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pTemp_new,
                             d_mpm_labels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempPrev_new,
                             d_mpm_labels->pTempPreviousLabel_preReloc,
                             pset);

      ParticleVariable<Point> pX_new;
      new_dw->allocateAndPut(pX_new, d_mpm_labels->pXLabel_preReloc, pset);

      ParticleVariable<Vector> pVelocity_new, pDisp_new;
      new_dw->allocateAndPut(pVelocity_new,
                             d_mpm_labels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pDisp_new, d_mpm_labels->pDispLabel_preReloc, pset);

      ParticleVariable<Matrix3> pSize_new;
      new_dw->allocateAndPut(pSize_new, d_mpm_labels->pSizeLabel_preReloc, pset);

      // Get grid variables
      constNCVariable<double> gTemperatureRate, frictionTempRate;
      new_dw->get(gTemperatureRate,
                  d_mpm_labels->gTemperatureRateLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpm_labels->frictionalWorkLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);

      constNCVariable<double> dTdt, massBurnFrac;
      if (d_mpm_flags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpm_labels->dTdt_NCLabel,
                    dwi,
                    patch,
                    gac,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpm_labels->massBurnFractionLabel,
                    dwi,
                    patch,
                    gac,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create, patch, gac, d_numGhostParticles);
        new_dw->allocateTemporary(massBurnFrac_create,
                                  patch,
                                  gac,
                                  d_numGhostParticles);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);

        dTdt         = dTdt_create;         // reference created data
        massBurnFrac = massBurnFrac_create; // reference created data
      }

      constNCVariable<Vector> gVelocityStar, gAcceleration;
      new_dw->get(gVelocityStar,
                  d_mpm_labels->gVelocityStarLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gAcceleration,
                  d_mpm_labels->gAccelerationLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);

      constNCVariable<Vector> gVelocityXPIC;
      if (d_mpm_flags->d_useXPIC) {
        new_dw->get(gVelocityXPIC,
                    d_mpm_labels->gVelocityXPICLabel,
                    dwi,
                    patch,
                    gac,
                    d_numGhostParticles);
      }

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      double Cp = mpm_matl->getSpecificHeat();
      /*
      double rho_frac_min = 0.;
      if (m == RMI) {
        rho_frac_min = 0.1;
      }
      */

      // Loop over particles
      for (auto idx : *pset) {

        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad[idx]);

        Vector velocity(0.0, 0.0, 0.0);
        Vector acceleration(0.0, 0.0, 0.0);
        Vector velocityXPIC(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0.0;
        double burnFraction = 0.0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < num_influence_nodes; k++) {
          IntVector node = ni[k];
          velocity += gVelocityStar[node] * S[k];
          acceleration += gAcceleration[node] * S[k];
          if (d_mpm_flags->d_useXPIC) {
            velocityXPIC += gVelocityXPIC[node] * S[k];
          }
          fricTempRate = frictionTempRate[node] * d_mpm_flags->d_addFrictionWork;
          tempRate +=
            (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
          burnFraction += massBurnFrac[node] * S[k];
        }

        // Update the particle's position and velocity
        if (d_mpm_flags->d_useXPIC) {
          pX_new[idx] = pX[idx] + velocity * delT -
                        0.5 *
                          (acceleration * delT + pVelocity[idx] -
                           2.0 * pVelocityXPIC[idx] + velocityXPIC) *
                          delT;
          pVelocity_new[idx] =
            2.0 * pVelocityXPIC[idx] - velocityXPIC + acceleration * delT;
          pDisp_new[idx] = pDisp[idx] + (pX_new[idx] - pX[idx]);
        } else {
          pX_new[idx]        = pX[idx] + velocity * delT * move_particles;
          pDisp_new[idx]     = pDisp[idx] + velocity * delT;
          pVelocity_new[idx] = pVelocity[idx] + acceleration * delT;
        }

        pTemp_new[idx]     = pTemperature[idx] + tempRate * delT;
        pTempPrev_new[idx] = pTemperature[idx]; // for thermal stress
        pMass_new[idx]     = std::max(pMass[idx] * (1.0 - burnFraction), 0.0);
        pSize_new[idx]     = (pMass_new[idx] / pMass[idx]) * pSize[idx];

        if (cout_heat.active()) {
          cout_heat << "MPM::Particle = " << pParticleID[idx]
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate << " dT = " << (tempRate * delT)
                    << " T_new = " << pTemp_new[idx] << std::endl;
        }

        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
        ke += .5 * pMass[idx] * pVelocity_new[idx].length2();
        CMX = CMX + (pX_new[idx] * pMass[idx]).asVector();
        totalMom += pVelocity_new[idx] * pMass[idx];
        totalmass += pMass_new[idx];
        partvoldef += pVolume[idx];

      } // End loop over particles

      // If load curves are being used with VelocityBC then apply
      // these BCs to the boundary particles
      if (d_mpm_flags->d_useLoadCurves) {

        std::vector<VelocityBC*> vbcP;
        bool do_VelocityBCs = false;
        for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
             ii++) {
          string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
          if (bcs_type == "Velocity") {
            do_VelocityBCs  = true;
            VelocityBC* vbc = dynamic_cast<VelocityBC*>(
              MPMPhysicalBCFactory::mpmPhysicalBCs[ii].get());
            vbcP.push_back(vbc);
          }
        }

        // std::cout << "do_VelocityBCs = " << do_VelocityBCs << std::endl;
        if (do_VelocityBCs) {

          // Get the current time
          simTime_vartype simTimeVar;
          old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
          double time = simTimeVar;

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          // Iterate over the particles
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (!(loadCurveID < 0)) {
              VelocityBC* vbc = vbcP[loadCurveID];
              pVelocity_new[idx] =
                vbc->getVelocityVector(pX[idx], pDisp[idx], time);
              pDisp_new[idx] = pDisp[idx] + pVelocity_new[idx] * delT;
              pX_new[idx] =
                pX[idx] + pVelocity_new[idx] * delT * move_particles;
              // std::cout << " Load curve ID = " << loadCurveID
              //           << " V = " << pVelocity_new[idx]
              //           << " U = " << pDisp_new[idx]
              //           << " x = " << pX_new[idx]
              //           << " num = " << pset->numParticles() << std::endl;
            }
          }
        }
      }

      // Delete particles that have left the domain
      // This is only needed if extra cells are being used.
      // Also delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (auto idx : *pset) {
        if ((pMass_new[idx] <= d_mpm_flags->d_minPartMass) ||
            (pTemp_new[idx] < 0.0)) {
          delset->addParticle(idx);
#ifdef CHECK_PARTICLE_DELETION
          proc0cout << "In " << __FILE__ << ":" << __LINE__ << std::endl;
          proc0cout << "Material = " << m
                    << " Deleted Particle = " << pParticleID_new[idx]
                    << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
                    << " vold = " << pVelocity[idx]
                    << " vnew = " << pVelocity_new[idx]
                    << " massold = " << pMass[idx]
                    << " massnew = " << pMass_new[idx]
                    << " tempold = " << pTemperature[idx]
                    << " tempnew = " << pTemp_new[idx]
                    << " volnew = " << pVolume[idx] << std::endl;
#endif
        }

        if (pVelocity_new[idx].length() > d_mpm_flags->d_maxVel) {
          if (d_mpm_flags->d_deleteRogueParticles) {
            delset->addParticle(idx);
            proc0cout << "\n Warning: particle " << pParticleID[idx]
                      << " hit speed ceiling #1. Deleting particle."
                      << std::endl;
          } else {
            if (pVelocity_new[idx].length() >= pVelocity[idx].length()) {
              pVelocity_new[idx] =
                (pVelocity_new[idx] / pVelocity_new[idx].length()) *
                (d_mpm_flags->d_maxVel * .9);
              proc0cout << "\n Warning: particle " << pParticleID[idx]
                        << " hit speed ceiling #1. Modifying particle velocity "
                           "accordingly."
                        << std::endl;
              // pVelocity_new[idx]=pVelocity[idx];
            }
          }
          proc0cout << "In " << __FILE__ << ":" << __LINE__ << std::endl;
          proc0cout << "Material = " << m
                    << " Deleted Particle = " << pParticleID_new[idx]
                    << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
                    << " vold = " << pVelocity[idx]
                    << " vnew = " << pVelocity_new[idx]
                    << " massold = " << pMass[idx]
                    << " massnew = " << pMass_new[idx]
                    << " tempold = " << pTemperature[idx]
                    << " tempnew = " << pTemp_new[idx]
                    << " vol = " << pVolume[idx] << "\n";
          proc0cout << " F_old = " << pDefGrad[idx] << "\n";
        }
      }

      new_dw->deleteParticles(delset);
    }

    // DON'T MOVE THESE!!!
    //__________________________________
    //  reduction variables
    if (d_mpm_flags->d_reductionVars->mass) {
      new_dw->put(sum_vartype(totalmass), d_mpm_labels->TotalMassLabel);
    }
    if (d_mpm_flags->d_reductionVars->volDeformed) {
      new_dw->put(sum_vartype(partvoldef),
                  d_mpm_labels->TotalVolumeDeformedLabel);
    }
    if (d_mpm_flags->d_reductionVars->momentum) {
      new_dw->put(sumvec_vartype(totalMom), d_mpm_labels->TotalMomentumLabel);
    }
    if (d_mpm_flags->d_reductionVars->KE) {
      new_dw->put(sum_vartype(ke), d_mpm_labels->KineticEnergyLabel);
    }
    if (d_mpm_flags->d_reductionVars->thermalEnergy) {
      new_dw->put(sum_vartype(thermal_energy), d_mpm_labels->ThermalEnergyLabel);
    }
    if (d_mpm_flags->d_reductionVars->centerOfMass) {
      new_dw->put(sumvec_vartype(CMX), d_mpm_labels->CenterOfMassPositionLabel);
    }

    // std::cout << "Solid mass lost this timestep = " << massLost << std::endl;
    // std::cout << "Solid momentum after advection = " << totalMom <<
    // std::endl;

    // std::cout << "THERMAL ENERGY " << thermal_energy << std::endl;
  }
}
