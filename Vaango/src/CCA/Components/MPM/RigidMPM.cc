/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2018-2023 Parresia Research Limited, NZ
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

#include <CCA/Components/MPM/RigidMPM.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

using namespace Uintah;

#undef INTEGRAL_TRACTION

static DebugStream cout_doing("RIGID_MPM", false);

RigidMPM::RigidMPM(const ProcessorGroup* myworld,
                   const MaterialManagerP& mat_manager)
  : SerialMPM(myworld, mat_manager)
{
}

RigidMPM::~RigidMPM() {}

void
RigidMPM::problemSetup(const ProblemSpecP& prob_spec,
                       const ProblemSpecP& restart_prob_spec,
                       GridP& grid,
                       [[maybe_unused]] const std::string& input_ups_dir)
{

  SerialMPM::problemSetup(prob_spec, restart_prob_spec, grid);

  ProblemSpecP cfd_ps = prob_spec->findBlock("CFD");
  if (cfd_ps && UintahParallelComponent::d_myworld->myRank() == 0) {
    std::cout << "\n__________________________________" << std::endl;
    std::cout << "  W A R N I N G :  " << std::endl;
    std::cout << "  You must use stiff MPM material properties" << std::endl;
    std::cout << "  to get the correct pressure solution in rmpmice"
              << std::endl;
    std::cout << "__________________________________\n" << std::endl;
  }
}

void
RigidMPM::computeStressTensor(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  if (cout_doing.active()) {
    cout_doing << "Doing computeStressTensor "
               << "\t\t\t\t RigidMPM" << std::endl;
  }

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->carryForward(patches, mpm_matl, old_dw, new_dw);
  }

  new_dw->put(delt_vartype(999.0), d_mpm_labels->delTLabel, getLevel(patches));
}

void
RigidMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("MPM::computeInternalForce",
                        this,
                        &RigidMPM::computeInternalForce);
  sched->addTask(t, patches, matls);
}

void
RigidMPM::computeInternalForce(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse*,
                               DataWarehouse*)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing computeInternalForce on patch " << patch->getID()
                 << "\t\t\t RigidMPM" << std::endl;
    }
  }
}

void
RigidMPM::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleComputeAndIntegrateAcceleration");

  Task* t = scinew Task("MPM::computeAndIntegrateAcceleration",
                        this,
                        &RigidMPM::computeAndIntegrateAcceleration);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  t->needs(Task::NewDW, d_mpm_labels->gVelocityLabel, Ghost::None);

  t->computes(d_mpm_labels->gVelocityStarLabel);
  t->computes(d_mpm_labels->gAccelerationLabel);

  sched->addTask(t, patches, matls);
}

void
RigidMPM::computeAndIntegrateAcceleration(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset*,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing computeAndIntegrateAcceleration");

    Ghost::GhostType gnone = Ghost::None;
    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<Vector> velocity;
      new_dw->get(velocity, d_mpm_labels->gVelocityLabel, dwi, patch, gnone, 0);

      // Create variables for the results
      NCVariable<Vector> velocity_star, acceleration;
      new_dw->allocateAndPut(velocity_star,
                             d_mpm_labels->gVelocityStarLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(acceleration,
                             d_mpm_labels->gAccelerationLabel,
                             dwi,
                             patch);

      acceleration.initialize(Vector(0., 0., 0.));

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c      = *iter;
        velocity_star[c] = velocity[c];
      }
    } // matls
  }
}

void
RigidMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("MPM::interpolateToParticlesAndUpdate",
                        this,
                        &RigidMPM::interpolateToParticlesAndUpdate);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  Ghost::GhostType gac = Ghost::AroundCells;
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureRateLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureNoBCLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gAccelerationLabel,
              gac,
              d_numGhostNodes);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pTemperatureLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pDispLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel_preReloc, Ghost::None);

  if (d_mpm_flags->d_withICE) {
    t->needs(Task::NewDW, d_mpm_labels->dTdt_NCLabel, gac, d_numGhostNodes);
  }

  t->computes(d_mpm_labels->pDispLabel_preReloc);
  t->computes(d_mpm_labels->pVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pXLabel_preReloc);
  t->computes(d_mpm_labels->pParticleIDLabel_preReloc);
  t->computes(d_mpm_labels->pTemperatureLabel_preReloc);
  t->computes(d_mpm_labels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpm_labels->pMassLabel_preReloc);
  t->computes(d_mpm_labels->pSizeLabel_preReloc);
  t->computes(d_mpm_labels->pXXLabel);

  t->computes(d_mpm_labels->KineticEnergyLabel);
  t->computes(d_mpm_labels->ThermalEnergyLabel);
  t->computes(d_mpm_labels->CenterOfMassPositionLabel);
  t->computes(d_mpm_labels->TotalMomentumLabel);

  // debugging scalar
  if (d_mpm_flags->d_withColor) {
    t->needs(Task::OldDW, d_mpm_labels->pColorLabel, Ghost::None);
    t->computes(d_mpm_labels->pColorLabel_preReloc);
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
RigidMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset*,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing interpolateToParticlesAndUpdate on patch "
                 << patch->getID() << "\t MPM" << std::endl;
    }

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    Vector CMX(0.0, 0.0, 0.0);
    Vector total_mom(0.0, 0.0, 0.0);
    double ke       = 0;
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    /*
    double move_particles=1.;
    if(!d_mpm_flags->d_doGridReset){
      move_particles=0.;
    }
    */

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> px;
      ParticleVariable<Point> pxnew, pxx;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Vector> pVelocitynew;
      ParticleVariable<Matrix3> pSizeNew;
      constParticleVariable<double> pMass, pTemperature;
      ParticleVariable<double> pMassNew, pTempNew;
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pids_new;
      constParticleVariable<Vector> pDisplacement;
      ParticleVariable<Vector> pDisplacementnew;
      constParticleVariable<Matrix3> pDeformationMeasure;

      // for thermal stress analysis
      ParticleVariable<double> pTempPreNew;

      // Get the arrays of grid data on which the new part. values depend
      constNCVariable<Vector> gVelocity_star, gacceleration;
      constNCVariable<double> gTemperatureRate, gTemperature, gTemperatureNoBC;
      constNCVariable<double> dTdt;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pDisplacement, d_mpm_labels->pDispLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pids, d_mpm_labels->pParticleIDLabel, pset);
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
      new_dw->allocateAndPut(pVelocitynew,
                             d_mpm_labels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pxnew, d_mpm_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pxx, d_mpm_labels->pXXLabel, pset);
      new_dw->allocateAndPut(pDisplacementnew, d_mpm_labels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(pMassNew, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pids_new,
                             d_mpm_labels->pParticleIDLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempNew,
                             d_mpm_labels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->get(pDeformationMeasure,
                  d_mpm_labels->pDefGradLabel_preReloc,
                  pset);

      // for thermal stress analysis
      new_dw->allocateAndPut(pTempPreNew,
                             d_mpm_labels->pTempPreviousLabel_preReloc,
                             pset);

      pids_new.copyData(pids);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(pSizeNew, d_mpm_labels->pSizeLabel_preReloc, pset);
      pSizeNew.copyData(pSize);

      // Carry forward NC_CCweight
      constNCVariable<double> NC_CCweight;
      NCVariable<double> NC_CCweight_new;
      Ghost::GhostType gnone = Ghost::None;
      old_dw
        ->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
      new_dw->allocateAndPut(NC_CCweight_new,
                             d_mpm_labels->NC_CCweightLabel,
                             0,
                             patch);
      NC_CCweight_new.copyData(NC_CCweight);

      Ghost::GhostType gac = Ghost::AroundCells;
      new_dw->get(gTemperatureRate,
                  d_mpm_labels->gTemperatureRateLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gTemperature,
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gTemperatureNoBC,
                  d_mpm_labels->gTemperatureNoBCLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gacceleration,
                  d_mpm_labels->gAccelerationLabel,
                  dwi,
                  patch,
                  gac,
                  d_numGhostParticles);
      if (d_mpm_flags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpm_labels->dTdt_NCLabel,
                    dwi,
                    patch,
                    gac,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create;
        new_dw->allocateTemporary(dTdt_create, patch, gac, d_numGhostParticles);
        dTdt_create.initialize(0.);
        dTdt = dTdt_create; // reference created data
      }

      // Get the constitutive model (needed for plastic temperature
      // update) and get the plastic temperature from the plasticity
      // model
      // For RigidMPM, this isn't done

      double Cp = mpm_matl->getSpecificHeat();

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(
          px[idx],
          ni,
          S,
          d_S,
          pSize[idx],
          pDeformationMeasure[idx]);

        double tempRate = 0.0;
        Vector acc(0.0, 0.0, 0.0);

        // Accumulate the contribution from each surrounding vertex
        // All we care about is the temperature field, everything else
        // should be zero.
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          tempRate += (gTemperatureRate[ni[k]] + dTdt[ni[k]]) * S[k];
          acc += gacceleration[ni[k]] * S[k];
        }
        pTempNew[idx]     = pTemperature[idx] + tempRate * delT;
        pVelocitynew[idx] = pVelocity[idx] + acc * delT;
        // If there is no adiabatic heating, add the plastic temperature
        // to the particle temperature

        // Update the particle's position and velocity
        pxnew[idx]    = px[idx] + pVelocity[idx] * delT;
        pDisplacementnew[idx] = pVelocity[idx] * delT;
        pMassNew[idx] = pMass[idx];
        // pxx is only useful if we're not in normal grid resetting mode.
        pxx[idx]         = px[idx];
        pTempPreNew[idx] = pTemperature[idx]; // for thermal stress

        thermal_energy += pTempNew[idx] * pMass[idx] * Cp;
        ke += .5 * pMass[idx] * pVelocitynew[idx].length2();
        CMX = CMX + (pxnew[idx] * pMass[idx]).asVector();
        total_mom += pVelocitynew[idx] * pMass[idx];
      }

      //__________________________________
      //  particle debugging label-- carry forward
      if (d_mpm_flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_mpm_labels->pColorLabel, pset);
        new_dw->allocateAndPut(pColor_new,
                               d_mpm_labels->pColorLabel_preReloc,
                               pset);
        pColor_new.copyData(pColor);
      }

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);
      new_dw->deleteParticles(delset);
    }

    // DON'T MOVE THESE!!!
    new_dw->put(sum_vartype(ke), d_mpm_labels->KineticEnergyLabel);
    new_dw->put(sumvec_vartype(CMX), d_mpm_labels->CenterOfMassPositionLabel);
    new_dw->put(sumvec_vartype(total_mom), d_mpm_labels->TotalMomentumLabel);
    new_dw->put(sum_vartype(thermal_energy), d_mpm_labels->ThermalEnergyLabel);
  }
}
