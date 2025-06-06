/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/CohesiveZone/CohesiveZoneTasks.h>

#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/CZLabel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/Patch.h>

#include <fstream>

using namespace Uintah;

static DebugStream cout_doing("CZ_MPM", false);

//______________________________________________________________________
//  Reference: N. P. Daphalapukar, Hongbing Lu, Demir Coker, Ranga Komanduri,
// " Simulation of dynamic crack growth using the generalized interpolation
// material point (GIMP) method," Int J. Fract, 2007, 143:79-102
//______________________________________________________________________
CohesiveZoneTasks::CohesiveZoneTasks(const ProblemSpecP& ps,
                                     const MaterialManagerP& mat_manager,
                                     const MPMLabel* mpm_labels,
                                     const CZLabel* cz_labels,
                                     const MPMFlags* mpm_flags)
{
  d_mat_manager = mat_manager;
  d_mpm_labels  = mpm_labels;
  d_cz_labels   = cz_labels;
  d_mpm_flags   = mpm_flags;

  if (mpm_flags->d_8or27 == 8) {
    d_num_ghost_particles = 1;
    d_num_ghost_nodes     = 1;
  } else {
    d_num_ghost_particles = 2;
    d_num_ghost_nodes     = 2;
  }

  if (mpm_flags->d_useCohesiveZones) {
    cohesiveZoneProblemSetup(ps);
  }
}

void
CohesiveZoneTasks::cohesiveZoneProblemSetup(const ProblemSpecP& prob_spec)
{
  // Search for the MaterialProperties block and then get the MPM section
  ProblemSpecP mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");
  ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");
  for (ProblemSpecP ps = mpm_mat_ps->findBlock("cohesive_zone"); ps != nullptr;
       ps              = ps->findNextBlock("cohesive_zone")) {

    string index("");
    ps->getAttribute("index", index);
    std::stringstream id(index);
    const int DEFAULT_VALUE = -1;
    int index_val           = DEFAULT_VALUE;

    id >> index_val;

    if (!id) {
      // stringstream parsing failed... on many (most) systems, the
      // original value assigned to index_val would be left
      // intact... but on some systems it inserts garbage,
      // so we have to manually restore the value.
      index_val = DEFAULT_VALUE;
    }
    // cout << "Material attribute = " << index_val << ", " << index << ", " <<
    // id << "\n";

    // Create and register as an MPM material
    std::shared_ptr<CZMaterial> mat =
      std::make_shared<CZMaterial>(ps, d_mat_manager, d_mpm_flags);

    mat->registerParticleState(d_cohesiveZoneState,
                               d_cohesiveZoneState_preReloc);

    // When doing restart, we need to make sure that we load the materials
    // in the same order that they were initially created.  Restarts will
    // ALWAYS have an index number as in <material index = "0">.
    // Index_val = -1 means that we don't register the material by its
    // index number.
    if (index_val > -1) {
      d_mat_manager->registerMaterial("CohesiveZone", mat, index_val);
    } else {
      d_mat_manager->registerMaterial("CohesiveZone", mat);
    }
  }
}

void
CohesiveZoneTasks::outputProblemSpec(ProblemSpecP& ps)
{
  for (size_t i = 0; i < d_mat_manager->getNumMaterials("CohesizeZone"); i++) {
    CZMaterial* mat =
      static_cast<CZMaterial*>(d_mat_manager->getMaterial("CohesiveZone", i));
    ProblemSpecP cm_ps = mat->outputProblemSpec(ps);
  }
}

void
CohesiveZoneTasks::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  int numCZM = d_mat_manager->getNumMaterials("CohesizeZone");
  for (int m = 0; m < numCZM; m++) {
    CZMaterial* cz_matl =
      static_cast<CZMaterial*>(d_mat_manager->getMaterial("CohesiveZone", m));
    CohesiveZone* ch = cz_matl->getCohesiveZone();
    ch->scheduleInitialize(level, sched, cz_matl);
  }
}

void
CohesiveZoneTasks::scheduleUpdate(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls)
{
  if (d_mpm_flags->d_useCohesiveZones) {
    const MaterialSet* cz_matls  = d_mat_manager->allMaterials("CohesiveZone");
    const MaterialSet* all_matls = d_mat_manager->allMaterials();
    const MaterialSubset* mpm_matls_sub = (matls ? matls->getUnion() : nullptr);
    const MaterialSubset* cz_matls_sub =
      (cz_matls ? cz_matls->getUnion() : nullptr);

    scheduleUpdateCohesiveZones(
      sched, patches, mpm_matls_sub, cz_matls_sub, all_matls);
    scheduleAddCohesiveZoneForces(
      sched, patches, mpm_matls_sub, cz_matls_sub, all_matls);
  }
}

void
CohesiveZoneTasks::scheduleUpdateCohesiveZones(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSubset* mpm_matls,
                                               const MaterialSubset* cz_matls,
                                               const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(
    patches, cout_doing, "CohesiveZoneTasks::scheduleUpdateCohesiveZones");

  Task* t = scinew Task("CohesiveZoneTasks::updateCohesiveZones",
                        this,
                        &CohesiveZoneTasks::updateCohesiveZones);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  Ghost::GhostType gnone = Ghost::None;
  t->needs(Task::NewDW,
              d_mpm_labels->gVelocityLabel,
              mpm_matls,
              Ghost::AroundCells,
              d_num_ghost_nodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gMassLabel,
              mpm_matls,
              Ghost::AroundCells,
              d_num_ghost_nodes);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czAreaLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czNormLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czTangLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czDispTopLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czDispBottomLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czSeparationLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czForceLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czTopMatLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czBotMatLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czFailedLabel, cz_matls, gnone);
  t->needs(Task::OldDW, d_cz_labels->czIDLabel, cz_matls, gnone);

  t->computes(d_mpm_labels->pXLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czAreaLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czNormLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czTangLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czDispTopLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czDispBottomLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czSeparationLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czForceLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czTopMatLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czBotMatLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czFailedLabel_preReloc, cz_matls);
  t->computes(d_cz_labels->czIDLabel_preReloc, cz_matls);

  sched->addTask(t, patches, matls);
}

void
CohesiveZoneTasks::updateCohesiveZones(const ProcessorGroup*,
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
              "Doing CohesiveZoneTasks::updateCohesiveZones");

    // The following is adapted from "Simulation of dynamic crack growth
    // using the generalized interpolation material point (GIMP) method"
    // Daphalapurkar, N.P., et al., Int. J. Fracture, 143, 79-102, 2007.

    auto interpolator = std::make_unique<LinearInterpolator>(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    unsigned int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    std::vector<constNCVariable<Vector>> gVelocity(numMPMMatls);
    std::vector<constNCVariable<double>> gMass(numMPMMatls);

    for (unsigned int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      new_dw->get(gVelocity[m],
                  d_mpm_labels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_num_ghost_nodes);
      new_dw->get(gMass[m],
                  d_mpm_labels->gMassLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_num_ghost_nodes);
    }

    size_t numCZMatls = d_mat_manager->getNumMaterials("CohesiveZone");
    for (size_t m = 0; m < numCZMatls; m++) {
      CZMaterial* cz_matl =
        static_cast<CZMaterial*>(d_mat_manager->getMaterial("CohesizeZone", m));
      int dwi = cz_matl->getDWIndex();

      // Not populating the delset, but we need this to satisfy Relocate
      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);
      new_dw->deleteParticles(delset);

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> czx;
      ParticleVariable<Point> czx_new;
      constParticleVariable<double> czarea;
      ParticleVariable<double> czarea_new;
      constParticleVariable<long64> czids;
      ParticleVariable<long64> czids_new;
      constParticleVariable<Vector> cznorm, pCZTangent, czDispTop;
      ParticleVariable<Vector> cznorm_new, pCZTangent_new, czDispTop_new;
      constParticleVariable<Vector> czDispBot, czsep, czforce;
      ParticleVariable<Vector> czDispBot_new, czsep_new, czforce_new;
      constParticleVariable<int> pCZTopMat, pCZBotMat, pCZFailed;
      ParticleVariable<int> pCZTopMat_new, pCZBotMat_new, pCZFailed_new;

      old_dw->get(czx, d_mpm_labels->pXLabel, pset);
      old_dw->get(czarea, d_cz_labels->czAreaLabel, pset);
      old_dw->get(cznorm, d_cz_labels->czNormLabel, pset);
      old_dw->get(pCZTangent, d_cz_labels->czTangLabel, pset);
      old_dw->get(czDispTop, d_cz_labels->czDispTopLabel, pset);
      old_dw->get(czDispBot, d_cz_labels->czDispBottomLabel, pset);
      old_dw->get(czsep, d_cz_labels->czSeparationLabel, pset);
      old_dw->get(czforce, d_cz_labels->czForceLabel, pset);
      old_dw->get(czids, d_cz_labels->czIDLabel, pset);
      old_dw->get(pCZTopMat, d_cz_labels->czTopMatLabel, pset);
      old_dw->get(pCZBotMat, d_cz_labels->czBotMatLabel, pset);
      old_dw->get(pCZFailed, d_cz_labels->czFailedLabel, pset);

      new_dw->allocateAndPut(czx_new, d_mpm_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(
        czarea_new, d_cz_labels->czAreaLabel_preReloc, pset);
      new_dw->allocateAndPut(
        cznorm_new, d_cz_labels->czNormLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pCZTangent_new, d_cz_labels->czTangLabel_preReloc, pset);
      new_dw->allocateAndPut(
        czDispTop_new, d_cz_labels->czDispTopLabel_preReloc, pset);
      new_dw->allocateAndPut(
        czforce_new, d_cz_labels->czForceLabel_preReloc, pset);
      new_dw->allocateAndPut(czids_new, d_cz_labels->czIDLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pCZTopMat_new, d_cz_labels->czTopMatLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pCZBotMat_new, d_cz_labels->czBotMatLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pCZFailed_new, d_cz_labels->czFailedLabel_preReloc, pset);
      new_dw->allocateAndPut(
        czDispBot_new, d_cz_labels->czDispBottomLabel_preReloc, pset);
      new_dw->allocateAndPut(
        czsep_new, d_cz_labels->czSeparationLabel_preReloc, pset);

      czarea_new.copyData(czarea);
      czids_new.copyData(czids);
      pCZTopMat_new.copyData(pCZTopMat);
      pCZBotMat_new.copyData(pCZBotMat);

      double sig_max      = cz_matl->getCohesiveNormalStrength();
      double delta_n      = cz_matl->getCharLengthNormal();
      double delta_t      = cz_matl->getCharLengthTangential();
      double tau_max      = cz_matl->getCohesiveTangentialStrength();
      double delta_s      = delta_t;
      double delta_n_fail = cz_matl->getNormalFailureDisplacement();
      double delta_t_fail = cz_matl->getTangentialFailureDisplacement();
      bool rotate_CZs     = cz_matl->getDoRotation();

      double phi_n = M_E * sig_max * delta_n;
      double phi_t = sqrt(M_E / 2) * tau_max * delta_t;
      double q     = phi_t / phi_n;
      // From the text following Eq. 15 in Nitin's paper it is a little hard
      // to tell what r should be, but zero seems like a reasonable value
      // based on the example problem in that paper
      double r = 0.;

      // Loop over particles
      for (auto idx : *pset) {

        Matrix3 size(0.1, 0., 0., 0., 0.1, 0., 0., 0., 0.1);
        Matrix3 defGrad(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(czx[idx], ni, S, size, defGrad);
        size_t NN = ni.size();

        Vector velTop(0.0, 0.0, 0.0);
        Vector velBot(0.0, 0.0, 0.0);
        double massTop = 0.0;
        double massBot = 0.0;
        int TopMat     = pCZTopMat[idx];
        int BotMat     = pCZBotMat[idx];
        double sumSTop = 0.;
        double sumSBot = 0.;

        // Accumulate the contribution from each surrounding vertex
        for (size_t k = 0; k < NN; k++) {
          IntVector node = ni[k];
          if (gMass[TopMat][node] > 2.e-200) {
            velTop += gVelocity[TopMat][node] * S[k];
            sumSTop += S[k];
          }
          if (gMass[BotMat][node] > 2.e-200) {
            velBot += gVelocity[BotMat][node] * S[k];
            sumSBot += S[k];
          }
          massTop += gMass[TopMat][node] * S[k];
          massBot += gMass[BotMat][node] * S[k];
        }
        velTop /= (sumSTop + 1.e-100);
        velBot /= (sumSBot + 1.e-100);

#if 0
        // I'm not sure what this was here for in the first place,
        // but it is disabled for now
        double mass_ratio = 0.0;
        if (massBot > 0.0) {
          mass_ratio = massTop/massBot;
          mass_ratio = min(mass_ratio,1.0/mass_ratio);
        }
        else {
          mass_ratio = 0.0;
        }

        double mass_correction_factor = mass_ratio;
#endif

        // Update the cohesive zone's position and displacements
        czx_new[idx]       = czx[idx] + .5 * (velTop + velBot) * delT;
        czDispTop_new[idx] = czDispTop[idx] + velTop * delT;
        czDispBot_new[idx] = czDispBot[idx] + velBot * delT;
        czsep_new[idx]     = czDispTop_new[idx] - czDispBot_new[idx];

        double disp_old = czsep[idx].length();
        double disp     = czsep_new[idx].length();
        if (disp > 0.0 && rotate_CZs) {
          Matrix3 Rotation;
          Matrix3 Rotation_tang;
          cz_matl->computeRotationMatrix(
            Rotation, Rotation_tang, cznorm[idx], czsep_new[idx]);

          cznorm_new[idx]     = Rotation * cznorm[idx];
          pCZTangent_new[idx] = Rotation_tang * pCZTangent[idx];
        } else {
          cznorm_new[idx]     = cznorm[idx];
          pCZTangent_new[idx] = pCZTangent[idx];
        }

        Vector pCZTangent2 = Cross(pCZTangent_new[idx], cznorm_new[idx]);

        double D_n  = Dot(czsep_new[idx], cznorm_new[idx]);
        double D_t1 = Dot(czsep_new[idx], pCZTangent_new[idx]);
        double D_t2 = Dot(czsep_new[idx], pCZTangent2);

        // Determine if a CZ has failed.
        double czf = 0.0;
        if (pCZFailed[idx] > 0) {
          if (disp >= disp_old) {
            pCZFailed_new[idx] = std::min(pCZFailed[idx] + 1, 1000);
          } else {
            pCZFailed_new[idx] = pCZFailed[idx];
          }
          czf = .001 * ((double)pCZFailed_new[idx]);
        } else if (fabs(D_n) > delta_n_fail) {
          std::cout << "pCZFailed, D_n =  " << std::endl;
          pCZFailed_new[idx] = 1;
        } else if (fabs(D_t1) > delta_t_fail) {
          pCZFailed_new[idx] = 1;
        } else if (fabs(D_t2) > delta_t_fail) {
          pCZFailed_new[idx] = 1;
        } else {
          pCZFailed_new[idx] = 0;
        }

        double normal_stress =
          (phi_n / delta_n) * exp(-D_n / delta_n) *
          ((D_n / delta_n) * exp((-D_t1 * D_t1) / (delta_t * delta_t)) +
           ((1. - q) / (r - 1.)) *
             (1. - exp(-D_t1 * D_t1 / (delta_t * delta_t))) *
             (r - D_n / delta_n));

        double tang1_stress =
          (phi_n / delta_n) * (2. * delta_n / delta_t) * (D_t1 / delta_t) *
          (q + ((r - q) / (r - 1.)) * (D_n / delta_n)) * exp(-D_n / delta_n) *
          exp(-D_t1 * D_t1 / (delta_t * delta_t));

        double tang2_stress =
          (phi_n / delta_n) * (2. * delta_n / delta_s) * (D_t2 / delta_s) *
          (q + ((r - q) / (r - 1.)) * (D_n / delta_n)) * exp(-D_n / delta_n) *
          exp(-D_t2 * D_t2 / (delta_s * delta_s));

        czforce_new[idx] =
          ((normal_stress * cznorm_new[idx] +
            tang1_stress * pCZTangent_new[idx] + tang2_stress * pCZTangent2) *
           czarea_new[idx]) *
          (1.0 - czf);
        /*
                dest << time << " " << czsep_new[idx].x() << " " <<
           czsep_new[idx].y() << " " << czforce_new[idx].x() << " " <<
           czforce_new[idx].y() << endl; if(fabs(normal_force) >= 0.0){ cout <<
           "czx_new " << czx_new[idx] << endl; cout << "czforce_new " <<
           czforce_new[idx] << endl; cout << "czsep_new " << czsep_new[idx] <<
           endl; cout << "czDispTop_new " << czDispTop_new[idx] << endl; cout <<
           "czDispBot_new " << czDispBot_new[idx] << endl; cout << "velTop " <<
           velTop << endl; cout << "velBot " << velBot << endl; cout << "delT "
           << delT << endl;
                }
        */
      }
    }
  }
}

void
CohesiveZoneTasks::scheduleAddCohesiveZoneForces(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSubset* mpm_matls,
  const MaterialSubset* cz_matls,
  const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(
    patches, cout_doing, "CohesiveZoneTasks::scheduleAddCohesiveZoneForces");

  Task* t = scinew Task("CohesiveZoneTasks::addCohesiveZoneForces",
                        this,
                        &CohesiveZoneTasks::addCohesiveZoneForces);

  t->needs(Task::OldDW,
              d_mpm_labels->pXLabel,
              cz_matls,
              Ghost::AroundNodes,
              d_num_ghost_particles);
  t->needs(Task::NewDW,
              d_cz_labels->czForceLabel_preReloc,
              cz_matls,
              Ghost::AroundNodes,
              d_num_ghost_particles);
  t->needs(Task::NewDW,
              d_cz_labels->czTopMatLabel_preReloc,
              cz_matls,
              Ghost::AroundNodes,
              d_num_ghost_particles);
  t->needs(Task::NewDW,
              d_cz_labels->czBotMatLabel_preReloc,
              cz_matls,
              Ghost::AroundNodes,
              d_num_ghost_particles);
  t->needs(Task::NewDW,
              d_mpm_labels->gMassLabel,
              mpm_matls,
              Ghost::AroundCells,
              d_num_ghost_nodes);

  t->modifies(d_mpm_labels->gExternalForceLabel, mpm_matls);

  sched->addTask(t, patches, matls);
}

void
CohesiveZoneTasks::addCohesiveZoneForces(const ProcessorGroup*,
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
              "Doing CohesiveZoneTasks::addCohesiveZoneForces");

    auto interpolator = std::make_unique<LinearInterpolator>(patch);
    vector<IntVector> ni(interpolator->size());
    vector<double> S(interpolator->size());

    unsigned int numMPMMatls = d_mat_manager->getNumMaterials("MPM");

    std::vector<NCVariable<Vector>> gext_force(numMPMMatls);
    std::vector<constNCVariable<double>> gMass(numMPMMatls);
    for (unsigned int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      new_dw->getModifiable(
        gext_force[m], d_mpm_labels->gExternalForceLabel, dwi, patch);
      new_dw->get(gMass[m],
                  d_mpm_labels->gMassLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_num_ghost_nodes);
    }

    unsigned int numCZMatls = d_mat_manager->getNumMaterials("CohesiveZone");
    for (unsigned int m = 0; m < numCZMatls; m++) {
      CZMaterial* cz_matl =
        static_cast<CZMaterial*>(d_mat_manager->getMaterial("CohesiveZone", m));
      int dwi = cz_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_num_ghost_particles,
                                                       d_mpm_labels->pXLabel);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> czx;
      constParticleVariable<Vector> czforce;
      constParticleVariable<int> pCZTopMat, pCZBotMat;

      old_dw->get(czx, d_mpm_labels->pXLabel, pset);
      new_dw->get(czforce, d_cz_labels->czForceLabel_preReloc, pset);
      new_dw->get(pCZTopMat, d_cz_labels->czTopMatLabel_preReloc, pset);
      new_dw->get(pCZBotMat, d_cz_labels->czBotMatLabel_preReloc, pset);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        Matrix3 size(0.1, 0., 0., 0., 0.1, 0., 0., 0., 0.1);
        Matrix3 defGrad(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(czx[idx], ni, S, size, defGrad);
        size_t NN = ni.size();

        int TopMat = pCZTopMat[idx];
        int BotMat = pCZBotMat[idx];

        double totMassTop = 0.;
        double totMassBot = 0.;
        //        double sumSTop = 0.;
        //        double sumSBot = 0.;

        for (size_t k = 0; k < NN; k++) {
          IntVector node = ni[k];
          totMassTop += S[k] * gMass[TopMat][node];
          totMassBot += S[k] * gMass[BotMat][node];
#if 0
          if(gMass[TopMat][node]>d_SMALL_NUM_MPM){
            sumSTop     += S[k];
          }
          if(gMass[BotMat][node]>d_SMALL_NUM_MPM){
            sumSBot     += S[k];
          }
#endif
        }

        // This currently contains three methods for distributing the CZ force
        // to the nodes.
        // The first of these distributes the force from the CZ
        // to the nodes based on a distance*mass weighting.
        // The second distributes the force to the nodes that have mass,
        // but only uses distance weighting.  So, a node that is near the CZ
        // but relatively far from particles may get a large acceleration
        // compared to other nodes, thereby inducing a velocity gradient.
        // The third simply does a distance weighting from the CZ to the nodes.
        // For this version, it is possible that nodes with no material mass
        // will still acquire force from the CZ, leading to ~infinite
        // acceleration, and thus, badness.

        // Accumulate the contribution from each surrounding vertex
        for (size_t k = 0; k < NN; k++) {
          IntVector node = ni[k];
          if (patch->containsNode(node)) {
            // Distribute force according to material mass on the nodes
            // to get an approximately equal contribution to the acceleration
            gext_force[BotMat][node] +=
              czforce[idx] * S[k] * gMass[BotMat][node] / totMassBot;
            gext_force[TopMat][node] -=
              czforce[idx] * S[k] * gMass[TopMat][node] / totMassTop;

            //            gext_force[BotMat][node] += czforce[idx]*S[k]/sumSBot;
            //            gext_force[TopMat][node] -= czforce[idx]*S[k]/sumSTop;

            //            gext_force[BotMat][node] = gext_force[BotMat][node]
            //                                     + czforce[idx] * S[k];
            //            gext_force[TopMat][node] = gext_force[TopMat][node]
            //                                     - czforce[idx] * S[k];
          }
        }
      }
#if 0
      // This is debugging output which is being left in for now (5/10/18)
      // as it may be helpful in generating figures for reports and papers.
      Vector sumForceTop = Vector(0.);
      Vector sumForceBot = Vector(0.);
      for(NodeIterator iter=patch->getExtraNodeIterator();
                       !iter.done();iter++){
        IntVector c = *iter;
        if(gext_force[1][c].length() > 1.e-100){
           cout << "gEF_BM[" << c << "] = " << gext_force[1][c] 
                << ", " << gext_force[1][c]/gMass[1][c] << endl;
           sumForceBot += gext_force[1][c];
        }
        if(gext_force[2][c].length() > 1.e-100){
           cout << "gEF_BM[" << c << "] = " << gext_force[2][c]
                << ", " << gext_force[2][c]/gMass[2][c] << endl;
           sumForceTop += gext_force[2][c];
        }
        if(gext_force[1][c].length() > 1.e-100 &&
           gext_force[2][c].length() > 1.e-100){
           cout << "ratio = " << (gext_force[1][c].x()/gMass[1][c])/
                                 (gext_force[2][c].x()/gMass[2][c]) << endl;
        }
      }
      cout << "SFB = " << sumForceBot << endl;
      cout << "SFT = " << sumForceTop << endl;
#endif
    }
  }
}

void
CohesiveZoneTasks::scheduleParticleRelocation(MaterialSubset* new_mss,
                                              constVarLabel2DArray& old_labels,
                                              constVarLabel2DArray& new_labels)
{
  if (d_mpm_flags->d_useCohesiveZones) {

    const MaterialSet* cz_matls = d_mat_manager->allMaterials("CohesiveZone");

    // update the mss
    const MaterialSubset* mss = cz_matls->getSubset(0);
    new_mss->addSubset(mss);

    // update the labels
    size_t numLabels = d_cohesiveZoneState_preReloc.size();
    for (size_t i = 0; i < numLabels; i++) {
      old_labels.push_back(d_cohesiveZoneState_preReloc[i]);
      new_labels.push_back(d_cohesiveZoneState[i]);
    }
  }
}