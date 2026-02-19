/*
 * The MIT License
 *
 * Copyright (c) 2013-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/DEM/DEMTasks.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/Util/DebugStream.h>

#include <iostream>
#include <map>
#include <vector>

using namespace Uintah;

static DebugStream cout_dem("DEM", false);

DEMTasks::DEMTasks([[maybe_unused]] ProblemSpecP& ps,
                   const MaterialManagerP& mat_manager,
                   const MPMLabel* mpm_labels,
                   const MPMFlags* mpm_flags,
                   int numGhostParticles,
                   int numGhostNodes)
  : d_mat_manager(mat_manager)
  , d_mpm_labels(mpm_labels)
  , d_mpm_flags(mpm_flags)
  , d_num_ghost_particles(numGhostParticles)
  , d_num_ghost_nodes(numGhostNodes)
{
}

void
DEMTasks::scheduleInitialize([[maybe_unused]] SchedulerP& sched,
                             [[maybe_unused]] const PatchSet* patches,
                             [[maybe_unused]] const MaterialSet* matls)
{
  // Initialization of DEM labels is done within SerialMPM::actuallyInitialize
  // through the hook actuallyInitialize() called per material/patch.
}

void
DEMTasks::actuallyInitialize(const Patch* patch,
                            const MPMMaterial* mpm_matl,
                            DataWarehouse* new_dw)
{
  if (d_mpm_flags->d_enableDEM) {
    int matID            = mpm_matl->getDWIndex();
    ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
    ParticleVariable<Vector> pAngularVelocity, pTorque;
    ParticleVariable<Matrix3> pOrientation, pInertiaTensor;
    ParticleVariable<double> pRadius;

    new_dw->getModifiable(pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, pset);
    new_dw->getModifiable(pTorque, d_mpm_labels->pTorqueLabel, pset);
    new_dw->getModifiable(pOrientation, d_mpm_labels->pOrientationLabel, pset);
    new_dw->getModifiable(pRadius, d_mpm_labels->pRadiusLabel, pset);
    new_dw->getModifiable(pInertiaTensor, d_mpm_labels->pInertiaTensorLabel, pset);

    double radius = mpm_matl->getParticleRadius();
    constParticleVariable<double> pMass;
    new_dw->get(pMass, d_mpm_labels->pMassLabel, pset);

    for (auto idx : *pset) {
      pAngularVelocity[idx] = Vector(0, 0, 0);
      pTorque[idx]          = Vector(0, 0, 0);
      pOrientation[idx]     = Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
      if (pMass[idx] > 0) {
        if (pRadius[idx] == 0) {
          pRadius[idx] = radius;
        }
        if (pInertiaTensor[idx].Determinant() == 0) {
          double I = 0.4 * pMass[idx] * pRadius[idx] * pRadius[idx];
          pInertiaTensor[idx] = Matrix3(I, 0, 0, 0, I, 0, 0, 0, I);
        }
      } else {
        pRadius[idx] = 0.0;
        pInertiaTensor[idx] = Matrix3(0.0);
      }
    }
  }
}

void
DEMTasks::computeStableTimestep(const Patch* patch,
                               const MPMMaterial* mpm_matl,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw,
                               double& dt_dem)
{
  double beta = 0.2; // Safety factor
  int matID = mpm_matl->getDWIndex();

  ParticleSubset* pset = nullptr;
  if (new_dw->haveParticleSubset(matID, patch)) {
    pset = new_dw->getParticleSubset(matID, patch);
  } else if (old_dw && old_dw->haveParticleSubset(matID, patch)) {
    pset = old_dw->getParticleSubset(matID, patch);
  }

  if (pset) {
    constParticleVariable<double> pMass;
    if (new_dw->exists(d_mpm_labels->pMassLabel, matID, patch)) {
      new_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
    } else {
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
    }

    double kn = mpm_matl->getDEMNormalStiffness();
    for (auto idx : *pset) {
      if (pMass[idx] > 0) {
        double dt_particle = beta * std::sqrt(pMass[idx] / kn);
        if (dt_particle < dt_dem) {
          dt_dem = dt_particle;
        }
      }
    }
  }
}

void
DEMTasks::scheduleComputeDEMForces(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls)
{
  Task* t = scinew Task("DEMTasks::computeDEMForces", this, &DEMTasks::computeDEMForces);

  Ghost::GhostType gan = Ghost::AroundNodes;
  int num_ghost = d_num_ghost_particles;

  t->needs(Task::OldDW, d_mpm_labels->pXLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pX0Label, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pRadiusLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pAngularVelocityLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pOrientationLabel, gan, num_ghost);
  t->needs(Task::OldDW, d_mpm_labels->pRigidBodyIDLabel, gan, num_ghost);

  t->modifies(d_mpm_labels->pExtForceLabel_preReloc);
  t->computes(d_mpm_labels->pTorqueLabel_preReloc);
  t->computes(d_mpm_labels->pRigidBodyIDLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
DEMTasks::computeDEMForces(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  Ghost::GhostType gan = Ghost::AroundNodes;
  int num_ghost        = d_num_ghost_particles;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int numMPMMatls = d_mat_manager->getNumMaterials("MPM");

    std::vector<constParticleVariable<Point>> pX_all(numMPMMatls);
    std::vector<constParticleVariable<Point>> pX0_all(numMPMMatls);
    std::vector<constParticleVariable<double>> pMass_all(numMPMMatls);
    std::vector<constParticleVariable<Matrix3>> pSize_all(numMPMMatls);
    std::vector<constParticleVariable<double>> pRadius_all(numMPMMatls);
    std::vector<constParticleVariable<Vector>> pVelocity_all(numMPMMatls);
    std::vector<constParticleVariable<Vector>> pAngVel_all(numMPMMatls);
    std::vector<constParticleVariable<Matrix3>> pOrientation_all(numMPMMatls);
    std::vector<constParticleVariable<long64>> pRigidBodyID_all(numMPMMatls);
    std::vector<ParticleSubset*> psets_all(numMPMMatls);

    std::vector<ParticleVariable<Vector>> pExtForce_new(numMPMMatls);
    std::vector<ParticleVariable<Vector>> pTorque_new(numMPMMatls);
    std::vector<ParticleSubset*> psets_real(numMPMMatls);
    std::vector<ParticleVariable<long64>> pRigidBodyID_new(numMPMMatls);

    std::map<long64, int> master_particles;

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      psets_all[m] = old_dw->getParticleSubset(matID, patch, gan, num_ghost, d_mpm_labels->pXLabel);
      if (psets_all[m]) {
        old_dw->get(pX_all[m], d_mpm_labels->pXLabel, psets_all[m]);
        old_dw->get(pX0_all[m], d_mpm_labels->pX0Label, psets_all[m]);
        old_dw->get(pMass_all[m], d_mpm_labels->pMassLabel, psets_all[m]);
        old_dw->get(pSize_all[m], d_mpm_labels->pSizeLabel, psets_all[m]);
        old_dw->get(pRadius_all[m], d_mpm_labels->pRadiusLabel, psets_all[m]);
        old_dw->get(pVelocity_all[m], d_mpm_labels->pVelocityLabel, psets_all[m]);
        old_dw->get(pAngVel_all[m], d_mpm_labels->pAngularVelocityLabel, psets_all[m]);
        old_dw->get(pOrientation_all[m], d_mpm_labels->pOrientationLabel, psets_all[m]);
        old_dw->get(pRigidBodyID_all[m], d_mpm_labels->pRigidBodyIDLabel, psets_all[m]);
      }

      psets_real[m] = old_dw->getParticleSubset(matID, patch);
      if (psets_real[m]) {
        new_dw->getModifiable(pExtForce_new[m], d_mpm_labels->pExtForceLabel_preReloc, psets_real[m]);
        new_dw->allocateAndPut(pTorque_new[m], d_mpm_labels->pTorqueLabel_preReloc, psets_real[m]);
        new_dw->allocateAndPut(pRigidBodyID_new[m], d_mpm_labels->pRigidBodyIDLabel_preReloc, psets_real[m]);

        pRigidBodyID_new[m].copyData(pRigidBodyID_all[m]);
        for (auto idx : *psets_real[m]) {
          pTorque_new[m][idx] = Vector(0.0, 0.0, 0.0);
          if (mpm_matl->isDiscrete() && pMass_all[m][idx] > 0) {
            master_particles[pRigidBodyID_all[m][idx]] = idx;
          }
        }
      }
    }

    for (int m_i = 0; m_i < numMPMMatls; m_i++) {
      if (!psets_real[m_i]) continue;
      MPMMaterial* matl_i = static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m_i));

      for (int m_j = m_i; m_j < numMPMMatls; m_j++) {
        if (!psets_all[m_j]) continue;
        MPMMaterial* matl_j = static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m_j));

        double kn    = 0.5 * (matl_i->getDEMNormalStiffness() + matl_j->getDEMNormalStiffness());
        double kt    = 0.5 * (matl_i->getDEMTangentialStiffness() + matl_j->getDEMTangentialStiffness());
        double gamma = 0.5 * (matl_i->getDEMDampingCoefficient() + matl_j->getDEMDampingCoefficient());
        double mu    = std::min(matl_i->getDEMFrictionCoefficient(), matl_j->getDEMFrictionCoefficient());

        for (auto idx_i : *psets_real[m_i]) {
          Point pos_i      = pX_all[m_i][idx_i];
          double rad_i     = pRadius_all[m_i][idx_i];
          Vector vel_i     = pVelocity_all[m_i][idx_i];
          Vector omega_i   = pAngVel_all[m_i][idx_i];
          Matrix3 orient_i = pOrientation_all[m_i][idx_i];

          std::map<long64, int> unique_bodies;
          if (matl_j->isDiscrete()) {
            for (auto idx_j : *psets_all[m_j]) {
              if (m_i == m_j && idx_i == idx_j) continue;
              if (m_i == m_j && pRigidBodyID_all[m_i][idx_i] == pRigidBodyID_all[m_j][idx_j]) continue;
              if (m_i == m_j && idx_j <= idx_i) continue;

              long64 rbID = pRigidBodyID_all[m_j][idx_j];
              if (unique_bodies.find(rbID) == unique_bodies.end()) {
                unique_bodies[rbID] = idx_j;
              } else {
                double d2_new = (pos_i - pX_all[m_j][idx_j]).length2();
                double d2_old = (pos_i - pX_all[m_j][unique_bodies[rbID]]).length2();
                if (d2_new < d2_old) unique_bodies[rbID] = idx_j;
              }
            }
          } else {
            for (auto idx_j : *psets_all[m_j]) {
              if (m_i == m_j && idx_i == idx_j) continue;
              if (m_i == m_j && idx_j <= idx_i) continue;
              unique_bodies[(long64)idx_j] = idx_j;
            }
          }

          for (auto const& [rbID_j, idx_j] : unique_bodies) {
            Vector totalForce(0, 0, 0);
            Vector arm_i(0, 0, 0);
            Vector arm_j_center(0, 0, 0);
            bool collision = false;

            bool i_is_rigid = matl_i->isDiscrete() && pMass_all[m_i][idx_i] > 0;
            bool j_is_rigid = matl_j->isDiscrete();

            if (i_is_rigid && !j_is_rigid) {
              int objIdx_i = (int)(pRigidBodyID_all[m_i][idx_i] & 0xFFFFFFFF);
              const GeometryObject* obj_i = matl_i->getGeometryObject(objIdx_i);
              if (obj_i) {
                const GeometryPieceP piece_i = obj_i->getPiece();
                Point pos0_i = pX0_all[m_i][idx_i];
                Point pos_j = pX_all[m_j][idx_j];
                double rad_j = pRadius_all[m_j][idx_j];
                Matrix3 rotT = orient_i.Transpose();
                Point localP_j = pos0_i + rotT * (pos_j - pos_i);
                double phi = piece_i->getSDF(localP_j);
                if (rad_j == 0) {
                  Vector worldGrad = orient_i * piece_i->getSDFGradient(localP_j);
                  Matrix3 size_j = pSize_all[m_j][idx_j];
                  rad_j = 0.5 * (std::abs(Dot(size_j.getColumn(0), worldGrad)) +
                                 std::abs(Dot(size_j.getColumn(1), worldGrad)) +
                                 std::abs(Dot(size_j.getColumn(2), worldGrad)));
                }
                if (phi < rad_j) {
                   double overlap = rad_j - phi;
                   Vector worldGrad = orient_i * piece_i->getSDFGradient(localP_j);
                   Vector contact_p = pos_j.asVector() - worldGrad * rad_j;
                   arm_i = contact_p - pos_i.asVector();
                   Vector arm_j = contact_p - pos_j.asVector();
                   Vector omega_j = pAngVel_all[m_j][idx_j];
                   Vector v_contact_i = vel_i + Cross(omega_i, arm_i);
                   Vector v_contact_j = pVelocity_all[m_j][idx_j] + Cross(omega_j, arm_j);
                   Vector v_rel = v_contact_j - v_contact_i;
                   double v_rel_n = Dot(v_rel, worldGrad);
                   double f_n_mag = kn * overlap - gamma * v_rel_n;
                   if (f_n_mag < 0) f_n_mag = 0;
                   Vector force_n = worldGrad * (-f_n_mag);
                   Vector v_rel_t = v_rel - worldGrad * v_rel_n;
                   Vector force_t = v_rel_t * kt * overlap;
                   double ft_mag  = force_t.length();
                   if (ft_mag > mu * f_n_mag) force_t = force_t * (mu * f_n_mag / ft_mag);
                   totalForce = force_n + force_t;
                   collision = true;
                }
              }
            } else if (!i_is_rigid && j_is_rigid) {
               int master_idx_j = master_particles.count(rbID_j) ? master_particles[rbID_j] : idx_j;
               int objIdx_j = (int)(rbID_j & 0xFFFFFFFF);
               const GeometryObject* obj_j = matl_j->getGeometryObject(objIdx_j);
               if (obj_j) {
                 const GeometryPieceP piece_j = obj_j->getPiece();
                 Point pos0_master_j = pX0_all[m_j][master_idx_j];
                 Point pos_master_j = pX_all[m_j][master_idx_j];
                 Matrix3 orient_master_j = pOrientation_all[m_j][master_idx_j];
                 Vector vel_master_j = pVelocity_all[m_j][master_idx_j];
                 Vector omega_master_j = pAngVel_all[m_j][master_idx_j];
                 Matrix3 rotT = orient_master_j.Transpose();
                 Point localP_i = pos0_master_j + rotT * (pos_i - pos_master_j);
                 double phi = piece_j->getSDF(localP_i);
                 double rad_eff_i = rad_i;
                 if (rad_eff_i == 0) {
                    Vector worldGrad = orient_master_j * piece_j->getSDFGradient(localP_i);
                    Matrix3 size_i = pSize_all[m_i][idx_i];
                    rad_eff_i = 0.5 * (std::abs(Dot(size_i.getColumn(0), worldGrad)) +
                                       std::abs(Dot(size_i.getColumn(1), worldGrad)) +
                                       std::abs(Dot(size_i.getColumn(2), worldGrad)));
                 }
                 if (phi < rad_eff_i) {
                    double overlap = rad_eff_i - phi;
                    Vector worldGrad = orient_master_j * piece_j->getSDFGradient(localP_i);
                    Vector contact_p = pos_i.asVector() - worldGrad * rad_eff_i;
                    arm_j_center = contact_p - pos_master_j.asVector();
                    arm_i = contact_p - pos_i.asVector();
                    Vector v_contact_i = vel_i + Cross(omega_i, arm_i);
                    Vector v_contact_j = vel_master_j + Cross(omega_master_j, arm_j_center);
                    Vector v_rel = v_contact_j - v_contact_i;
                    double v_rel_n = Dot(v_rel, worldGrad);
                    double f_n_mag = kn * overlap + gamma * v_rel_n;
                    if (f_n_mag < 0) f_n_mag = 0;
                    Vector force_n = worldGrad * (f_n_mag);
                    Vector v_rel_t = v_rel - worldGrad * v_rel_n;
                    Vector force_t = v_rel_t * kt * overlap;
                    double ft_mag  = force_t.length();
                    if (ft_mag > mu * f_n_mag) force_t = force_t * (mu * f_n_mag / ft_mag);
                    totalForce = force_n + force_t;
                    collision = true;
                 }
               }
            } else if (i_is_rigid && j_is_rigid) {
               Point pos_j = pX_all[m_j][idx_j];
               double rad_j = pRadius_all[m_j][idx_j];
               Vector dist = pos_j - pos_i;
               double d = dist.length();
               if (d < (rad_i + rad_j)) {
                  double overlap = (rad_i + rad_j) - d;
                  Vector normal = dist / d;
                  Vector v_rel = (pVelocity_all[m_j][idx_j] + Cross(pAngVel_all[m_j][idx_j], Vector(0,0,0))) - 
                                 (vel_i + Cross(omega_i, Vector(0,0,0)));
                  double v_rel_n = Dot(v_rel, normal);
                  double f_n_mag = kn * overlap - gamma * v_rel_n;
                  if (f_n_mag < 0) f_n_mag = 0;
                  totalForce = normal * (-f_n_mag);
                  collision = true;
               }
            }

            if (collision) {
              if (cout_dem.active()) {
                cout_dem << "Collision detected between mat " << m_i << " particle " << idx_i 
                         << " and mat " << m_j << " body " << rbID_j << std::endl;
              }
              pExtForce_new[m_i][idx_i] += totalForce;
              pTorque_new[m_i][idx_i] += Cross(arm_i, totalForce);
              if (matl_j->isDiscrete()) {
                if (master_particles.count(rbID_j)) {
                  int master_idx = master_particles[rbID_j];
                  if (patch->containsPoint(pX_all[m_j][master_idx])) {
                    pExtForce_new[m_j][master_idx] -= totalForce;
                    pTorque_new[m_j][master_idx] -= Cross(arm_j_center, totalForce);
                  }
                }
              } else {
                if (patch->containsPoint(pX_all[m_j][idx_j])) {
                  pExtForce_new[m_j][idx_j] -= totalForce;
                }
              }
            }
          }
        }
      }
    }
  }
}

void
DEMTasks::scheduleIntegrateDEMRotation(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  Task* t = scinew Task("DEMTasks::integrateDEMRotation", this, &DEMTasks::integrateDEMRotation);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);
  t->needs(Task::OldDW, d_mpm_labels->pRadiusLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pAngularVelocityLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pOrientationLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pInertiaTensorLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pRigidBodyIDLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->pTorqueLabel_preReloc, Ghost::None);

  t->computes(d_mpm_labels->pAngularVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pOrientationLabel_preReloc);
  t->computes(d_mpm_labels->pInertiaTensorLabel_preReloc);
  t->computes(d_mpm_labels->pRadiusLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
DEMTasks::integrateDEMRotation(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      if (!pset) continue;

      constParticleVariable<Vector> pAngularVelocity, pTorque;
      constParticleVariable<Matrix3> pOrientation, pInertiaTensor;
      constParticleVariable<double> pRadius, pMass;
      constParticleVariable<long64> pRigidBodyID;
      ParticleVariable<Vector> pAngularVelocity_new;
      ParticleVariable<Matrix3> pOrientation_new, pInertiaTensor_new;
      ParticleVariable<double> pRadius_new;

      old_dw->get(pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, pset);
      old_dw->get(pOrientation, d_mpm_labels->pOrientationLabel, pset);
      old_dw->get(pInertiaTensor, d_mpm_labels->pInertiaTensorLabel, pset);
      old_dw->get(pRadius, d_mpm_labels->pRadiusLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pRigidBodyID, d_mpm_labels->pRigidBodyIDLabel, pset);
      new_dw->get(pTorque, d_mpm_labels->pTorqueLabel_preReloc, pset);

      new_dw->allocateAndPut(pAngularVelocity_new, d_mpm_labels->pAngularVelocityLabel_preReloc, pset);
      new_dw->allocateAndPut(pOrientation_new, d_mpm_labels->pOrientationLabel_preReloc, pset);
      new_dw->allocateAndPut(pInertiaTensor_new, d_mpm_labels->pInertiaTensorLabel_preReloc, pset);
      new_dw->allocateAndPut(pRadius_new, d_mpm_labels->pRadiusLabel_preReloc, pset);

      struct RotationState { Vector omega; Matrix3 orient; Matrix3 inertia; };
      std::map<long64, RotationState> master_rotations;

      if (mpm_matl->isDiscrete() || mpm_matl->isRigid()) {
        for (auto idx : *pset) {
          if (pRadius[idx] > 0.0 && pInertiaTensor[idx].Determinant() > 0.0) {
            Matrix3 invI              = pInertiaTensor[idx].Inverse();
            Vector angAcc             = invI * pTorque[idx];
            pAngularVelocity_new[idx] = pAngularVelocity[idx] + angAcc * delT;
            Matrix3 omegaSkew(0.0, -pAngularVelocity[idx].z(), pAngularVelocity[idx].y(),
                              pAngularVelocity[idx].z(), 0.0, -pAngularVelocity[idx].x(),
                              -pAngularVelocity[idx].y(), pAngularVelocity[idx].x(), 0.0);
            Matrix3 identity; identity.Identity();
            Matrix3 deltaR = identity + omegaSkew * delT;
            pOrientation_new[idx] = deltaR * pOrientation[idx];
            pInertiaTensor_new[idx] = GeometryPiece::rotateInertiaTensor(pInertiaTensor[idx], deltaR);
            master_rotations[pRigidBodyID[idx]] = {pAngularVelocity_new[idx], pOrientation_new[idx], pInertiaTensor_new[idx]};
          } else {
            pAngularVelocity_new[idx] = pAngularVelocity[idx];
            pOrientation_new[idx]     = pOrientation[idx];
            pInertiaTensor_new[idx] = pInertiaTensor[idx];
          }
          pRadius_new[idx] = pRadius[idx];
        }
        for (auto idx : *pset) {
          if (pRadius[idx] <= 0.0) {
            long64 rbID = pRigidBodyID[idx];
            if (master_rotations.count(rbID)) {
              pAngularVelocity_new[idx] = master_rotations[rbID].omega;
              pOrientation_new[idx] = master_rotations[rbID].orient;
              pInertiaTensor_new[idx] = master_rotations[rbID].inertia;
            }
          }
        }
      } else {
        pAngularVelocity_new.copyData(pAngularVelocity);
        pOrientation_new.copyData(pOrientation);
        pInertiaTensor_new.copyData(pInertiaTensor);
        pRadius_new.copyData(pRadius);
      }
    }
  }
}
