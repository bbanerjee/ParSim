/*
 * The MIT License
 *
 * Copyright (c) 2015-2018 Parresia Research Limited
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
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <Core/Grid/Variables/VarTypes.h>

//#define XPIC2_UPDATE
#define DEBUG_WITH_PARTICLE_ID

using namespace Uintah;

static DebugStream cout_doing("MPM_USL", false);
static DebugStream cout_heat("MPM_USL_Heat", false);

MPM_UpdateStressLast::MPM_UpdateStressLast(const ProcessorGroup* myworld)
  : SerialMPM(myworld)
{
}

MPM_UpdateStressLast::~MPM_UpdateStressLast()
{
}

void
MPM_UpdateStressLast::scheduleTimeAdvance(const LevelP& level,
                                          SchedulerP& sched)
{
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()))
    return;

  const PatchSet* patches = level->eachPatch();
  const MaterialSet* matls = d_sharedState->allMPMMaterials();
  const MaterialSet* cz_matls = d_sharedState->allCZMaterials();
  const MaterialSet* all_matls = d_sharedState->allMaterials();

  const MaterialSubset* mpm_matls_sub = (matls ? matls->getUnion() : nullptr);
  const MaterialSubset* cz_matls_sub =
    (cz_matls ? cz_matls->getUnion() : nullptr);

  scheduleComputeParticleBodyForce(sched, patches, matls);
  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleInterpolateParticlesToGrid(sched, patches, matls);
  scheduleExMomInterpolated(sched, patches, matls);
  if (flags->d_useCohesiveZones) {
    scheduleUpdateCohesiveZones(sched, patches, mpm_matls_sub, cz_matls_sub,
                                all_matls);
    scheduleAddCohesiveZoneForces(sched, patches, mpm_matls_sub, cz_matls_sub,
                                  all_matls);
  }
  if (d_boundaryTractionFaces.size() > 0) {
    scheduleComputeContactArea(sched, patches, matls);
  }
  scheduleComputeInternalForce(sched, patches, matls);
  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  scheduleExMomIntegrated(sched, patches, matls);
  scheduleSetGridBoundaryConditions(sched, patches, matls);
  if (flags->d_prescribeDeformation) {
    scheduleSetPrescribedMotion(sched, patches, matls);
  }
  #ifdef XPIC2_UPDATE
    scheduleComputeXPICVelocities(sched, patches, matls);
  #endif
  if (flags->d_doExplicitHeatConduction) {
    scheduleComputeHeatExchange(sched, patches, matls);
    scheduleComputeInternalHeatRate(sched, patches, matls);
    // scheduleComputeNodalHeatFlux(sched, patches, matls);
    scheduleSolveHeatEquations(sched, patches, matls);
    scheduleIntegrateTemperatureRate(sched, patches, matls);
  }
  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  scheduleComputeDeformationGradient(sched, patches, matls);
  scheduleComputeStressTensor(sched, patches, matls);
  scheduleComputeBasicDamage(sched, patches, matls);
  scheduleUpdateErosionParameter(sched, patches, matls);
  scheduleFindRogueParticles(sched, patches, matls);
  scheduleInsertParticles(sched, patches, matls);
  if (flags->d_refineParticles) {
    scheduleAddParticles(sched, patches, matls);
  }
  if (flags->d_computeScaleFactor) {
    scheduleComputeParticleScaleFactor(sched, patches, matls);
  }

  sched->scheduleParticleRelocation(
    level, lb->pXLabel_preReloc, d_sharedState->d_particleState_preReloc, lb->pXLabel,
    d_sharedState->d_particleState, lb->pParticleIDLabel, matls, 1);

  if (flags->d_useCohesiveZones) {
    sched->scheduleParticleRelocation(
      level, lb->pXLabel_preReloc, d_sharedState->d_cohesiveZoneState_preReloc, lb->pXLabel,
      d_sharedState->d_cohesiveZoneState, lb->czIDLabel, cz_matls, 2);
  }
}

void
MPM_UpdateStressLast::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                              const PatchSet* patches,
                                                              const MaterialSet* matls)

{
  if (!flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                           getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM_USL::scheduleInterpolateToParticlesAndUpdate");

  Task* t=scinew Task("MPM_USL::interpolateToParticlesAndUpdate",
                      this, &MPM_UpdateStressLast::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, lb->gAccelerationLabel,      gac,  NGN);
  t->requires(Task::NewDW, lb->gVelocityStarLabel,      gac,  NGN);
  #ifdef XPIC2_UPDATE
    t->requires(Task::NewDW, lb->gVelocityXPICLabel,      gac,  NGN);
  #endif
  t->requires(Task::NewDW, lb->gTemperatureRateLabel,   gac,  NGN);
  t->requires(Task::NewDW, lb->frictionalWorkLabel,     gac,  NGN);
  t->requires(Task::OldDW, lb->pXLabel,                 gnone);
  t->requires(Task::OldDW, lb->pMassLabel,              gnone);
  t->requires(Task::OldDW, lb->pParticleIDLabel,        gnone);
  t->requires(Task::OldDW, lb->pTemperatureLabel,       gnone);
  t->requires(Task::OldDW, lb->pVelocityLabel,          gnone);
  #ifdef XPIC2_UPDATE
    t->requires(Task::OldDW, lb->pVelocityXPICLabel,      gnone);
  #endif
  t->requires(Task::OldDW, lb->pDispLabel,              gnone);
  t->requires(Task::OldDW, lb->pSizeLabel,              gnone);
  t->requires(Task::OldDW, lb->pVolumeLabel,            gnone);
  t->requires(Task::OldDW, lb->pDefGradLabel,           gnone);
  if (flags->d_useLoadCurves) {
    t->requires(Task::OldDW, lb->pLoadCurveIDLabel,     Ghost::None);
  }
  if(flags->d_withICE){
    t->requires(Task::NewDW, lb->dTdt_NCLabel,         gac,NGN);
    t->requires(Task::NewDW, lb->massBurnFractionLabel,gac,NGN);
  }

  t->computes(lb->pDispLabel_preReloc);
  t->computes(lb->pVelocityLabel_preReloc);
  t->computes(lb->pXLabel_preReloc);
  t->computes(lb->pParticleIDLabel_preReloc);
  t->computes(lb->pTemperatureLabel_preReloc);
  t->computes(lb->pTempPreviousLabel_preReloc); // for thermal stress 
  t->computes(lb->pMassLabel_preReloc);
  t->computes(lb->pSizeLabel_preReloc);

  //__________________________________
  //  reduction variables
  if(flags->d_reductionVars->momentum){
    t->computes(lb->TotalMomentumLabel);
  }
  if(flags->d_reductionVars->KE){
    t->computes(lb->KineticEnergyLabel);
  }
  if(flags->d_reductionVars->thermalEnergy){
    t->computes(lb->ThermalEnergyLabel);
  }
  if(flags->d_reductionVars->centerOfMass){
    t->computes(lb->CenterOfMassPositionLabel);
  }
  if(flags->d_reductionVars->mass){
    t->computes(lb->TotalMassLabel);
  }
  if(flags->d_withColor) {
    t->requires(Task::OldDW, lb->pColorLabel,  Ghost::None);
    t->computes(lb->pColorLabel_preReloc);
  }

  // Carry Forward particle refinement flag
  if(flags->d_refineParticles){
    t->requires(Task::OldDW, lb->pRefinedLabel,                Ghost::None);
    t->computes(             lb->pRefinedLabel_preReloc);
  }

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->requires(Task::OldDW, lb->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(             lb->NC_CCweightLabel, z_matl);

  sched->addTask(t, patches, matls);

  // The task will have a reference to z_matl
  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
}

void 
MPM_UpdateStressLast::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                                      const PatchSubset* patches,
                                                      const MaterialSubset* ,
                                                      DataWarehouse* old_dw,
                                                      DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac = Ghost::AroundCells;

  for (auto patch : *patches) {
    printTask(patches, patch, cout_doing, "Doing interpolateToParticlesAndUpdate");

    auto interpolator = flags->d_interpolator->clone(patch);
    auto num_influence_nodes = interpolator->size();
    vector<IntVector> ni(num_influence_nodes);
    vector<double> S(num_influence_nodes);

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    double totalmass = 0;
    double partvoldef = 0.;
    Vector CMX(0.0,0.0,0.0);
    Vector totalMom(0.0,0.0,0.0);
    double ke=0;

    delt_vartype delT;
    old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

    /*
    Material* reactant;
    reactant = d_sharedState->getMaterialByName("reactant");
    bool combustion_problem=false;
    int RMI = -99;
    if(reactant != 0){
      RMI = reactant->getDWIndex();
      combustion_problem=true;
    }
    */

    double move_particles=1.;
    if(!flags->d_doGridReset){
      move_particles=0.;
    }

    // Copy NC_CCweight (only material 0)
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_new;
    old_dw->get(NC_CCweight,                lb->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NC_CCweight_new, lb->NC_CCweightLabel, 0, patch);
    NC_CCweight_new.copyData(NC_CCweight);

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int dwi = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Copy particle IDs
      constParticleVariable<long64>  pParticleID;
      ParticleVariable<long64>  pParticleID_new;
      old_dw->get(pParticleID,                lb->pParticleIDLabel,          pset);
      new_dw->allocateAndPut(pParticleID_new, lb->pParticleIDLabel_preReloc, pset);
      pParticleID_new.copyData(pParticleID);

      // Copy color
      if (flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double>pColor_new;
        old_dw->get(pColor,                lb->pColorLabel,          pset);
        new_dw->allocateAndPut(pColor_new, lb->pColorLabel_preReloc, pset);
        pColor_new.copyData(pColor);
      }    

      // Copy refined
      if(flags->d_refineParticles){
        constParticleVariable<int> pRefined;
        ParticleVariable<int> pRefined_new;
        old_dw->get(pRefined,                lb->pRefinedLabel,          pset);
        new_dw->allocateAndPut(pRefined_new, lb->pRefinedLabel_preReloc, pset);
        pRefined_new.copyData(pRefined);
      }

      // Get particle variables
      constParticleVariable<double>  pMass, pTemperature, pVolume;
      old_dw->get(pMass,             lb->pMassLabel,                pset);
      old_dw->get(pVolume,           lb->pVolumeLabel,              pset);
      old_dw->get(pTemperature,      lb->pTemperatureLabel,         pset);

      constParticleVariable<Point>   pX;
      old_dw->get(pX,                lb->pXLabel,                   pset);

      constParticleVariable<Vector>  pVelocity, pDisp;
      old_dw->get(pVelocity,         lb->pVelocityLabel,            pset);
      old_dw->get(pDisp,             lb->pDispLabel,                pset);

      #ifdef XPIC2_UPDATE
        constParticleVariable<Vector>  pVelocityXPIC;
        new_dw->get(pVelocityXPIC,     lb->pVelocityXPICLabel,        pset);
      #endif

      constParticleVariable<Matrix3> pDefGrad, pSize;
      old_dw->get(pDefGrad,          lb->pDefGradLabel,             pset);
      old_dw->get(pSize,             lb->pSizeLabel,                pset);

      // Allocate updated particle variables
      ParticleVariable<double>  pMass_new, pTemp_new, pTempPrev_new;
      new_dw->allocateAndPut(pMass_new,     lb->pMassLabel_preReloc,         pset);
      new_dw->allocateAndPut(pTemp_new,     lb->pTemperatureLabel_preReloc,  pset);
      new_dw->allocateAndPut(pTempPrev_new, lb->pTempPreviousLabel_preReloc, pset);

      ParticleVariable<Point>   pX_new;
      new_dw->allocateAndPut(pX_new,        lb->pXLabel_preReloc,            pset);

      ParticleVariable<Vector>  pVelocity_new, pDisp_new;
      new_dw->allocateAndPut(pVelocity_new, lb->pVelocityLabel_preReloc,     pset);
      new_dw->allocateAndPut(pDisp_new,     lb->pDispLabel_preReloc,         pset);

      ParticleVariable<Matrix3> pSize_new;
      new_dw->allocateAndPut(pSize_new,     lb->pSizeLabel_preReloc,         pset);

      // Get grid variables
      constNCVariable<double> gTemperatureRate, frictionTempRate;
      new_dw->get(gTemperatureRate, lb->gTemperatureRateLabel, dwi, patch, gac, NGP);
      new_dw->get(frictionTempRate, lb->frictionalWorkLabel,   dwi, patch, gac, NGP);

      constNCVariable<double> dTdt, massBurnFrac;
      if (flags->d_withICE) {
        new_dw->get(dTdt,          lb->dTdt_NCLabel,          dwi, patch, gac, NGP);
        new_dw->get(massBurnFrac,  lb->massBurnFractionLabel, dwi, patch, gac, NGP);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create,         patch, gac, NGP);
        new_dw->allocateTemporary(massBurnFrac_create, patch, gac, NGP);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);

        dTdt = dTdt_create;                         // reference created data
        massBurnFrac = massBurnFrac_create;         // reference created data
      }

      constNCVariable<Vector> gVelocityStar, gAcceleration;
      new_dw->get(gVelocityStar,  lb->gVelocityStarLabel, dwi, patch, gac, NGP);
      new_dw->get(gAcceleration,  lb->gAccelerationLabel, dwi, patch, gac, NGP);

      #ifdef XPIC2_UPDATE
        constNCVariable<Vector> gVelocityXPIC;
        new_dw->get(gVelocityXPIC,  lb->gVelocityXPICLabel, dwi, patch, gac, NGP);
      #endif

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      double Cp = mpm_matl->getSpecificHeat();
      /*
      double rho_frac_min = 0.;
      if (m == RMI) {
        rho_frac_min = 0.1;
      }
      */

      // Loop over particles
      for(auto idx : *pset) {

        interpolator->findCellAndWeights(pX[idx], ni, S, pSize[idx], pDefGrad[idx]);

        Vector velocity(0.0,0.0,0.0);
        Vector acceleration(0.0,0.0,0.0);
        #ifdef XPIC2_UPDATE
          Vector velocityXPIC(0.0, 0.0, 0.0);
        #endif
        double fricTempRate = 0.0;
        double tempRate = 0.0;
        double burnFraction = 0.0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < num_influence_nodes; k++) {
          IntVector node = ni[k];
          velocity      += gVelocityStar[node] * S[k];
          acceleration  += gAcceleration[node] * S[k];
          #ifdef XPIC2_UPDATE
            velocityXPIC  += gVelocityXPIC[node] * S[k];
          #endif
          fricTempRate = frictionTempRate[node] * flags->d_addFrictionWork;
          tempRate += (gTemperatureRate[node] + dTdt[node] +
                       fricTempRate) * S[k];
          burnFraction += massBurnFrac[node] * S[k];
        }

        // Update the particle's position and velocity
        #ifdef XPIC2_UPDATE
          pX_new[idx] = pX[idx] + velocity * delT
                        - 0.5 * (acceleration * delT  
                                 + pVelocity[idx] - 2.0 * pVelocityXPIC[idx]
                                 + velocityXPIC) * delT;
          pVelocity_new[idx]  = 2.0*pVelocityXPIC[idx] - velocityXPIC + 
                                acceleration*delT;
          pDisp_new[idx] = pDisp[idx] + (pX_new[idx] - pX[idx]);
        #else
          pX_new[idx]        = pX[idx]        + velocity * delT * move_particles;
          pDisp_new[idx]     = pDisp[idx]     + velocity * delT;
          pVelocity_new[idx] = pVelocity[idx] + acceleration * delT;
        #endif

        pTemp_new[idx]     = pTemperature[idx] + tempRate * delT;
        pTempPrev_new[idx] = pTemperature[idx]; // for thermal stress
        pMass_new[idx]     = std::max(pMass[idx] * (1.0 - burnFraction), 0.0);
        pSize_new[idx]     = (pMass_new[idx] / pMass[idx]) * pSize[idx];

        if (cout_heat.active()) {
          cout_heat << "MPM::Particle = " << pParticleID[idx]
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate
                    << " dT = " << (tempRate*delT)
                    << " T_new = " << pTemp_new[idx] << endl;
        }

        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
        ke += .5*pMass[idx]*pVelocity_new[idx].length2();
        CMX         = CMX + (pX_new[idx]*pMass[idx]).asVector();
        totalMom   += pVelocity_new[idx]*pMass[idx];
        totalmass  += pMass_new[idx];
        partvoldef += pVolume[idx];

      } // End loop over particles

      // If load curves are being used with VelocityBC then apply 
      // these BCs to the boundary particles
      if (flags->d_useLoadCurves) {

        std::vector<VelocityBC*> vbcP;
        bool do_VelocityBCs = false;
        for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size(); ii++) {
          string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
          if (bcs_type == "Velocity") {
            do_VelocityBCs = true;
            VelocityBC* vbc =
              dynamic_cast<VelocityBC*>(MPMPhysicalBCFactory::mpmPhysicalBCs[ii].get());
            vbcP.push_back(vbc);
          }
        }

        //std::cout << "do_VelocityBCs = " << do_VelocityBCs << std::endl;
        if (do_VelocityBCs) {

          // Get the current time
          double time = d_sharedState->getElapsedTime();

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);

          // Iterate over the particles
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx]-1;
            if (!(loadCurveID < 0)) {
              VelocityBC* vbc = vbcP[loadCurveID];
              pVelocity_new[idx] = vbc->getVelocityVector(pX[idx], pDisp[idx], time);
              pDisp_new[idx] = pDisp[idx] + pVelocity_new[idx]*delT;
              pX_new[idx] = pX[idx] + pVelocity_new[idx]*delT*move_particles;
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
        if ( (pMass_new[idx] <= flags->d_minPartMass) || 
             (pTemp_new[idx] < 0.0) ) {
          delset->addParticle(idx);
          #ifdef CHECK_PARTICLE_DELETION
            proc0cout << "In " << __FILE__ << ":" << __LINE__ << std::endl;
            proc0cout << "Material = " << m << " Deleted Particle = " << pParticleID_new[idx] 
                      << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
                      << " vold = " << pVelocity[idx] << " vnew = "<< pVelocity_new[idx]
                      << " massold = " << pMass[idx] << " massnew = " << pMass_new[idx]
                      << " tempold = " << pTemperature[idx] 
                      << " tempnew = " << pTemp_new[idx]
                      << " volnew = " << pVolume[idx] << endl;
          #endif
        }
        
        if (pVelocity_new[idx].length() > flags->d_maxVel) {
          if (flags->d_deleteRogueParticles) {
            delset->addParticle(idx);
            proc0cout << "\n Warning: particle " << pParticleID[idx] 
                      << " hit speed ceiling #1. Deleting particle." 
                      << std::endl;
          } else {
            if (pVelocity_new[idx].length() >= pVelocity[idx].length()) {
              pVelocity_new[idx] = 
                (pVelocity_new[idx]/pVelocity_new[idx].length())*(flags->d_maxVel*.9);      
              proc0cout << "\n Warning: particle "<< pParticleID[idx] 
                        << " hit speed ceiling #1. Modifying particle velocity accordingly."
                        << std::endl;
              //pVelocity_new[idx]=pVelocity[idx];
            }
          }
          proc0cout << "In " << __FILE__ << ":" << __LINE__ << std::endl;
          proc0cout << "Material = " << m << " Deleted Particle = " << pParticleID_new[idx] 
                    << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
                    << " vold = " << pVelocity[idx] << " vnew = "<< pVelocity_new[idx]
                    << " massold = " << pMass[idx] << " massnew = " << pMass_new[idx]
                    << " tempold = " << pTemperature[idx] << " tempnew = " << pTemp_new[idx]
                    << " vol = " << pVolume[idx] << "\n";
          proc0cout << " F_old = " << pDefGrad[idx] << "\n";
        }
      }

      new_dw->deleteParticles(delset);    
    }

    // DON'T MOVE THESE!!!
    //__________________________________
    //  reduction variables
    if(flags->d_reductionVars->mass){
      new_dw->put(sum_vartype(totalmass),      lb->TotalMassLabel);
    }
    if(flags->d_reductionVars->volDeformed){
      new_dw->put(sum_vartype(partvoldef),     lb->TotalVolumeDeformedLabel);
    }
    if(flags->d_reductionVars->momentum){
      new_dw->put(sumvec_vartype(totalMom),    lb->TotalMomentumLabel);
    }
    if(flags->d_reductionVars->KE){
      new_dw->put(sum_vartype(ke),             lb->KineticEnergyLabel);
    }
    if(flags->d_reductionVars->thermalEnergy){
      new_dw->put(sum_vartype(thermal_energy), lb->ThermalEnergyLabel);
    }
    if(flags->d_reductionVars->centerOfMass){
      new_dw->put(sumvec_vartype(CMX),         lb->CenterOfMassPositionLabel);
    }

    // std::cout << "Solid mass lost this timestep = " << massLost << endl;
    // std::cout << "Solid momentum after advection = " << totalMom << endl;

    // std::cout << "THERMAL ENERGY " << thermal_energy << endl;
  }
  
}


