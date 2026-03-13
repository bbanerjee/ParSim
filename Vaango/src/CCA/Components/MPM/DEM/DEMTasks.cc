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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/DEM/DEMTasks.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Util/DOUT.hpp>

#include <iostream>
#include <map>
#include <vector>

using namespace Uintah;

// Debug streams:  To turn on use  export SCI_DEBUG="DEM:+,SerialMPM:+" 
Dout dem_doing("DEMdoing", "MPM", "DEM tasks", false);
Dout dem_doing_fn("DEMdoingfn", "MPM", "DEM subtasks", false);
Dout dem_dbg("DEMdbg", "MPM", "Debug DEM contact forces", false);

namespace Vaango {

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
  DOUT(dem_doing, "[DEMTasks::DEMTasks] Create DEMTasks object");
}

void
DEMTasks::scheduleInitialize(Task* task)
{
  if (d_mpm_flags->d_enableDEM) {
    DOUT(dem_doing, "[DEMTasks::scheduleInitialize]");
    task->computes(d_mpm_labels->pX0Label);
    task->computes(d_mpm_labels->pRigidBodyIDLabel);
    task->computes(d_mpm_labels->pAngularVelocityLabel);
    task->computes(d_mpm_labels->pTorqueLabel);
    task->computes(d_mpm_labels->pOrientationLabel);
    task->computes(d_mpm_labels->pRadiusLabel);
    task->computes(d_mpm_labels->pInertiaTensorLabel);
  }
}

void
DEMTasks::actuallyInitialize(const Patch* patch,
                             const MPMMaterial* mpm_matl,
                             DataWarehouse* new_dw)
{
  if (d_mpm_flags->d_enableDEM) {
    DOUT(dem_doing, "[DEMTasks::actuallyInitialize] Start");

    int matID            = mpm_matl->getDWIndex();
    ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
    ParticleVariable<Vector> pAngularVelocity, pTorque;
    ParticleVariable<Matrix3> pOrientation, pInertiaTensor;
    ParticleVariable<double> pRadius;

    new_dw->getModifiable(
      pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, pset);
    new_dw->getModifiable(pTorque, d_mpm_labels->pTorqueLabel, pset);
    new_dw->getModifiable(pOrientation, d_mpm_labels->pOrientationLabel, pset);
    new_dw->getModifiable(pRadius, d_mpm_labels->pRadiusLabel, pset);
    new_dw->getModifiable(
      pInertiaTensor, d_mpm_labels->pInertiaTensorLabel, pset);

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
          double I            = 0.4 * pMass[idx] * pRadius[idx] * pRadius[idx];
          pInertiaTensor[idx] = Matrix3(I, 0, 0, 0, I, 0, 0, 0, I);
        }
      } else {
        pRadius[idx]        = 0.0;
        pInertiaTensor[idx] = Matrix3(0.0);
      }
    }
    DOUT(dem_doing, "[DEMTasks::actuallyInitialize] Done");
  }
}

void
DEMTasks::computeStableTimestep(const Patch* patch,
                                const MPMMaterial* mpm_matl,
                                DataWarehouse* new_dw,
                                double& dt_dem)
{
  DOUT(dem_doing, "[DEMTasks::computeStableTimestep] Start");
  double beta = 0.2; // Safety factor
  int matID   = mpm_matl->getDWIndex();
  DOUT(dem_dbg, "[DEMTasks::computeStableTimestep] matID = " << matID << " patch = " << patch);
  DOUT(dem_dbg, "[DEMTasks::computeStableTimestep] new_dw = " << new_dw);
  DOUT(dem_dbg, "[DEMTasks::computeStableTimestep] lo = " << patch->getExtraCellLowIndex());
  DOUT(dem_dbg, "[DEMTasks::computeStableTimestep] hi = " << patch->getExtraCellHighIndex());

  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pMass;
  new_dw->get(pMass, d_mpm_labels->pMassLabel, pset);

  double kn = mpm_matl->getDEMNormalStiffness();
  for (auto idx : *pset) {
    if (pMass[idx] > 0) {
      double dt_particle = beta * std::sqrt(pMass[idx] / kn);
      if (dt_particle < dt_dem) {
        dt_dem = dt_particle;
      }
    }
  }

  DOUT(dem_doing, "[DEMTasks::computeStableTimestep] Done");
}

void
DEMTasks::scheduleComputeDEMForces(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  DOUT(dem_doing, "[DEMTasks::scheduleComputeDEMForces]");
  Task* t = scinew Task(
    "DEMTasks::computeDEMForces", this, &DEMTasks::computeDEMForces);

  Ghost::GhostType gan = Ghost::AroundNodes;
  int num_ghost        = d_num_ghost_particles;

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
  DOUT(dem_doing, "[DEMTasks::computeDEMForces]");
  int numMatls = d_mat_manager->getNumMaterials("MPM");

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector cell_size = patch->dCell();

    DEMParticleSets psets(numMatls);
    DEMParticleInputData inputs(numMatls);
    DEMParticleOutputData outputs(numMatls);

    getAndAllocateParticleData(patch, old_dw, new_dw, psets, inputs, outputs);

    auto master_particles = buildMasterParticleMap(numMatls, psets, inputs);

    for (int mat_i = 0; mat_i < numMatls; mat_i++) {

      // Go to next material if there are no non-ghost particles for this material
      if (!psets.psets_real[mat_i]) {
        continue;
      }

      MPMMaterial* matl_i =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", mat_i));
      bool i_is_dem_material = matl_i->isDEMMaterial();

      for (int mat_j = mat_i; mat_j < numMatls; mat_j++) {

        // Go to next material if there are no particles for this material
        if (!psets.psets_all[mat_j]) {
          continue;
        }

        MPMMaterial* matl_j =
          static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", mat_j));
        bool j_is_dem_material = matl_j->isDEMMaterial();

        DEMContactProps props = getContactProps(matl_i, matl_j);

        // Once per (mat_j, timestep) — mirrors the numInCell build pattern
        DEMRigidBodySpatialIndex spatial_index(mat_j, psets, inputs, patch,
                                               j_is_dem_material);

        // Loop over real (non-ghost) particles in material i
        for (auto pidx_i : *psets.psets_real[mat_i]) {

          // Identify master DEM particles and flag them
          bool i_is_rigid = i_is_dem_material && inputs.pMass_old[mat_i][pidx_i] > 0;

          auto unique_bodies = buildUniqueBodies(mat_i, pidx_i, mat_j,
                                                 inputs, spatial_index, cell_size);

          for (auto const& [rbID_j, pidx_j] : unique_bodies) {
            bool j_is_rigid = j_is_dem_material;

            DEMContactResult contact;
            if (i_is_rigid && !j_is_rigid) {
              contact = computeRigidVsParticle(
                pidx_i, mat_i, pidx_j, mat_j, props, inputs, matl_i);
            } else if (!i_is_rigid && j_is_rigid) {
              contact = computeParticleVsRigid(pidx_i,
                                               mat_i,
                                               pidx_j,
                                               mat_j,
                                               rbID_j,
                                               props,
                                               inputs,
                                               matl_j,
                                               master_particles);
            } else if (i_is_rigid && j_is_rigid) {
              contact =
                computeRigidVsRigid(pidx_i, mat_i, pidx_j, mat_j, props, inputs);
            }

            if (contact.collision) {
              applyContactForces(contact,
                                 pidx_i,
                                 mat_i,
                                 rbID_j,
                                 mat_j,
                                 inputs,
                                 matl_j,
                                 patch,
                                 master_particles,
                                 outputs);
            }
          }
        }
      }
    }
  }
}

// Get contact properties
DEMContactProps
DEMTasks::getContactProps(const MPMMaterial* matl_i,
                          const MPMMaterial* matl_j) const
{
  return { 0.5 * (matl_i->getDEMNormalStiffness() +
                  matl_j->getDEMNormalStiffness()),
           0.5 * (matl_i->getDEMTangentialStiffness() +
                  matl_j->getDEMTangentialStiffness()),
           0.5 * (matl_i->getDEMDampingCoefficient() +
                  matl_j->getDEMDampingCoefficient()),
           std::min(matl_i->getDEMFrictionCoefficient(),
                    matl_j->getDEMFrictionCoefficient()) };
}

// Load data for contact force calculations
void
DEMTasks::getAndAllocateParticleData(const Patch* patch,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw,
                                     DEMParticleSets& psets,
                                     DEMParticleInputData& inputs,
                                     DEMParticleOutputData& outputs)
{
  DOUT(dem_doing, "[DEMTasks::getAndAllocateParticleData]");
  Ghost::GhostType gan     = Ghost::AroundNodes;
  int              numGhost = d_num_ghost_particles;
  int              numMatls = d_mat_manager->getNumMaterials("MPM");

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* matl =
      static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
    int matID = matl->getDWIndex();

    // Get particle data for all particles (including ghosts) for this material
    psets.psets_all[m] = old_dw->getParticleSubset(
      matID, patch, gan, numGhost, d_mpm_labels->pXLabel);

    if (psets.psets_all[m]) {
      old_dw->get(inputs.pX_old[m], d_mpm_labels->pXLabel, psets.psets_all[m]);
      old_dw->get(inputs.pX0_old[m], d_mpm_labels->pX0Label, psets.psets_all[m]);
      old_dw->get(inputs.pMass_old[m], d_mpm_labels->pMassLabel, psets.psets_all[m]);
      old_dw->get(inputs.pSize_old[m], d_mpm_labels->pSizeLabel, psets.psets_all[m]);
      old_dw->get(
        inputs.pRadius_old[m], d_mpm_labels->pRadiusLabel, psets.psets_all[m]);
      old_dw->get(
        inputs.pVelocity_old[m], d_mpm_labels->pVelocityLabel, psets.psets_all[m]);
      old_dw->get(inputs.pAngVel_old[m],
                  d_mpm_labels->pAngularVelocityLabel,
                  psets.psets_all[m]);
      old_dw->get(inputs.pOrientation_old[m],
                  d_mpm_labels->pOrientationLabel,
                  psets.psets_all[m]);
      old_dw->get(inputs.pRigidBodyID_old[m],
                  d_mpm_labels->pRigidBodyIDLabel,
                  psets.psets_all[m]);
    }

    // Get particle subset for real (non-ghost) particles for this material
    psets.psets_real[m] = old_dw->getParticleSubset(matID, patch);
    if (psets.psets_real[m]) {
      new_dw->getModifiable(outputs.pExtForce_new[m],
                            d_mpm_labels->pExtForceLabel_preReloc,
                            psets.psets_real[m]);
      new_dw->allocateAndPut(outputs.pTorque_new[m],
                             d_mpm_labels->pTorqueLabel_preReloc,
                             psets.psets_real[m]);
      new_dw->allocateAndPut(outputs.pRigidBodyID_new[m],
                             d_mpm_labels->pRigidBodyIDLabel_preReloc,
                             psets.psets_real[m]);

      outputs.pRigidBodyID_new[m].copyData(inputs.pRigidBodyID_old[m]);
      for (auto idx : *psets.psets_real[m]) {
        outputs.pTorque_new[m][idx] = Vector(0, 0, 0);
      }
    }
  }
}

// Master particle map
ParticleIDToCurrentIdxMap
DEMTasks::buildMasterParticleMap(int numMatls,
                                 const DEMParticleSets& psets,
                                 const DEMParticleInputData& inputs) const
{
  DOUT(dem_doing, "[DEMTasks::buildMasterParticleMap]");
  ParticleIDToCurrentIdxMap master_particles;
  for (int m = 0; m < numMatls; m++) {
    // If there are no non-ghost particles for this material, skip it 
    if (!psets.psets_real[m]) {
      continue;
    }

    MPMMaterial* matl =
      static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

    // Ignore non-DEM materials
    if (!matl->isDEMMaterial()) {
      continue;
    }

    for (auto idx : *psets.psets_real[m]) {
      if (inputs.pMass_old[m][idx] > 0) {
        master_particles[inputs.pRigidBodyID_old[m][idx]] = idx;
        DOUT(dem_dbg, "[DEM: buildMasterParticleMap] m=" << m << " idx=" << idx
                     << " rbID=" << inputs.pRigidBodyID_old[m][idx]);
      }
    }
  }
  return master_particles;
}

// Unique body selection
// Returns a map from rigid-body ID -> representative particle index for
// all j-material bodies that particle pidx_i should interact with.
ParticleIDToCurrentIdxMap
DEMTasks::buildUniqueBodies(int mat_i,
                            particleIndex pidx_i,
                            int mat_j,     
                            const DEMParticleInputData& inputs,
                            const DEMRigidBodySpatialIndex& spatial_index,
                            const Uintah::Vector& cell_size) const
{
  DOUT(dem_doing_fn, "[DEMTasks::buildUniqueBodies]" 
    << " mat_i=" << mat_i << " pidx_i=" << pidx_i << " mat_j=" << mat_j);
  double cell_radius = 1.5 * std::sqrt(3) * cell_size.maxComponent(); // Diagonal of a cell

  return spatial_index.query(inputs.pX_old[mat_i][pidx_i], mat_i, pidx_i,
                             cell_radius);
}

// ─── Per-case contact force routines ─────────────────────────────────────────
// Case A: rigid SDF body (i) vs. MPM particle (j)
DEMContactResult
DEMTasks::computeRigidVsParticle(particleIndex pidx_i,
                                 int mat_i,
                                 particleIndex pidx_j,
                                 int mat_j,
                                 const DEMContactProps& props,
                                 const DEMParticleInputData& inputs,
                                 const MPMMaterial* matl_i) const
{
  DOUT(dem_doing_fn, "[DEMTasks::computeRigidVsParticle]");
  DEMContactResult res;

  int objIdx_i                = (int)(inputs.pRigidBodyID_old[mat_i][pidx_i] & 0xFFFFFFFF);
  const GeometryObject* obj_i = matl_i->getGeometryObject(objIdx_i);
  if (!obj_i) {
    return res;
  }

  const GeometryPieceP piece_i = obj_i->getPiece();
  Point pos_i                  = inputs.pX_old[mat_i][pidx_i];
  Point pos0_i                 = inputs.pX0_old[mat_i][pidx_i];
  Vector vel_i                 = inputs.pVelocity_old[mat_i][pidx_i];
  Vector omega_i               = inputs.pAngVel_old[mat_i][pidx_i];
  Matrix3 orient_i             = inputs.pOrientation_old[mat_i][pidx_i];
  Matrix3 rotT                 = orient_i.Transpose();

  Point pos_j  = inputs.pX_old[mat_j][pidx_j];
  double rad_j = inputs.pRadius_old[mat_j][pidx_j];

  Point localP_j = pos0_i + rotT * (pos_j - pos_i);
  double phi     = piece_i->getSDF(localP_j);

  if (rad_j == 0) {
    Vector worldGrad = orient_i * piece_i->getSDFGradient(localP_j);
    Matrix3 size_j   = inputs.pSize_old[mat_j][pidx_j];
    rad_j            = 0.5 * (std::abs(Dot(size_j.getColumn(0), worldGrad)) +
                   std::abs(Dot(size_j.getColumn(1), worldGrad)) +
                   std::abs(Dot(size_j.getColumn(2), worldGrad)));
  }

  if (phi >= rad_j) {
    return res;
  }

  double overlap   = rad_j - phi;
  Vector worldGrad = orient_i * piece_i->getSDFGradient(localP_j);
  Vector contact_p = pos_j.asVector() - worldGrad * rad_j;
  Vector arm_i     = contact_p - pos_i.asVector();
  Vector arm_j     = contact_p - pos_j.asVector();
  Vector omega_j   = inputs.pAngVel_old[mat_j][pidx_j];

  Vector v_contact_i = vel_i + Cross(omega_i, arm_i);
  Vector v_contact_j = inputs.pVelocity_old[mat_j][pidx_j] + Cross(omega_j, arm_j);
  Vector v_rel       = v_contact_j - v_contact_i;
  double v_rel_n     = Dot(v_rel, worldGrad);

  double f_n_mag = std::max(0.0, props.kn * overlap - props.gamma * v_rel_n);
  Vector force_n = worldGrad * (-f_n_mag);
  Vector v_rel_t = v_rel - worldGrad * v_rel_n;
  Vector force_t = v_rel_t * props.kt * overlap;
  if (double ft = force_t.length(); ft > props.mu * f_n_mag) {
    force_t *= (props.mu * f_n_mag / ft);
  }

  res.totalForce = force_n + force_t;
  res.arm_i      = arm_i;
  res.collision  = true;
  return res;
}

// Case B: MPM particle (i) vs. rigid SDF body (j)
DEMContactResult
DEMTasks::computeParticleVsRigid(
  particleIndex pidx_i,
  int mat_i,
  particleIndex pidx_j,
  int mat_j,
  long64 rbID_j,
  const DEMContactProps& props,
  const DEMParticleInputData& inputs,
  const MPMMaterial* matl_j,
  const ParticleIDToCurrentIdxMap& master_particles) const
{
  DOUT(dem_dbg, "[DEMTasks::computeParticleVsRigid]");
  DEMContactResult res;

  int objIdx_j                = (int)(rbID_j & 0xFFFFFFFF);
  const GeometryObject* obj_j = matl_j->getGeometryObject(objIdx_j);
  if (!obj_j) {
    return res;
  }

  int master_pidx_j =
    master_particles.count(rbID_j) ? master_particles.at(rbID_j) : pidx_j;

  const GeometryPieceP piece_j = obj_j->getPiece();
  Point pos0_master_j          = inputs.pX0_old[mat_j][master_pidx_j];
  Point pos_master_j           = inputs.pX_old[mat_j][master_pidx_j];
  Matrix3 orient_master_j      = inputs.pOrientation_old[mat_j][master_pidx_j];
  Vector vel_master_j          = inputs.pVelocity_old[mat_j][master_pidx_j];
  Vector omega_master_j        = inputs.pAngVel_old[mat_j][master_pidx_j];
  Matrix3 rotT                 = orient_master_j.Transpose();

  Point pos_i    = inputs.pX_old[mat_i][pidx_i];
  Vector vel_i   = inputs.pVelocity_old[mat_i][pidx_i];
  Vector omega_i = inputs.pAngVel_old[mat_i][pidx_i];
  double rad_i   = inputs.pRadius_old[mat_i][pidx_i];

  Point localP_i = pos0_master_j + rotT * (pos_i - pos_master_j);
  double phi     = piece_j->getSDF(localP_i);
  double rad_eff = rad_i;

  if (rad_eff == 0) {
    Vector worldGrad = orient_master_j * piece_j->getSDFGradient(localP_i);
    Matrix3 size_i   = inputs.pSize_old[mat_i][pidx_i];
    rad_eff          = 0.5 * (std::abs(Dot(size_i.getColumn(0), worldGrad)) +
                     std::abs(Dot(size_i.getColumn(1), worldGrad)) +
                     std::abs(Dot(size_i.getColumn(2), worldGrad)));
  }

  if (phi >= rad_eff) {
    return res;
  }

  double overlap   = rad_eff - phi;
  Vector worldGrad = orient_master_j * piece_j->getSDFGradient(localP_i);
  Vector contact_p = pos_i.asVector() - worldGrad * rad_eff;
  Vector arm_i     = contact_p - pos_i.asVector();
  Vector arm_j_c   = contact_p - pos_master_j.asVector();

  Vector v_contact_i = vel_i + Cross(omega_i, arm_i);
  Vector v_contact_j = vel_master_j + Cross(omega_master_j, arm_j_c);
  Vector v_rel       = v_contact_j - v_contact_i;
  double v_rel_n     = Dot(v_rel, worldGrad);

  double f_n_mag = std::max(0.0, props.kn * overlap + props.gamma * v_rel_n);
  Vector force_n = worldGrad * f_n_mag;
  Vector v_rel_t = v_rel - worldGrad * v_rel_n;
  Vector force_t = v_rel_t * props.kt * overlap;
  if (double ft = force_t.length(); ft > props.mu * f_n_mag) {
    force_t *= (props.mu * f_n_mag / ft);
  }

  res.totalForce   = force_n + force_t;
  res.arm_i        = arm_i;
  res.arm_j_center = arm_j_c;
  res.collision    = true;
  return res;
}

// Case C: rigid sphere (i) vs. rigid sphere (j)
DEMContactResult
DEMTasks::computeRigidVsRigid(particleIndex pidx_i,
                              int mat_i,
                              particleIndex pidx_j,
                              int mat_j,
                              const DEMContactProps& props,
                              const DEMParticleInputData& inputs) const
{
  DOUT(dem_dbg, "[DEMTasks::computeRigidVsRigid]");
  DEMContactResult res;

  Point pos_i    = inputs.pX_old[mat_i][pidx_i];
  double rad_i   = inputs.pRadius_old[mat_i][pidx_i];
  Vector vel_i   = inputs.pVelocity_old[mat_i][pidx_i];
  Vector omega_i = inputs.pAngVel_old[mat_i][pidx_i];
  Point pos_j    = inputs.pX_old[mat_j][pidx_j];
  double rad_j   = inputs.pRadius_old[mat_j][pidx_j];
  Vector vel_j   = inputs.pVelocity_old[mat_j][pidx_j];
  Vector omega_j = inputs.pAngVel_old[mat_j][pidx_j];

  Vector dist = pos_j - pos_i;
  double d    = dist.length();
  if (d >= rad_i + rad_j) {
    return res;
  }

  double overlap = (rad_i + rad_j) - d;
  Vector normal  = dist / d;

  Vector v_rel = (vel_j + Cross(omega_j, Vector(0, 0, 0))) -
                 (vel_i + Cross(omega_i, Vector(0, 0, 0)));
  double v_rel_n = Dot(v_rel, normal);

  double f_n_mag = std::max(0.0, props.kn * overlap - props.gamma * v_rel_n);
  if (f_n_mag < 0) {
    f_n_mag = 0;
  }
  res.totalForce = normal * (-f_n_mag);
  res.collision  = true;
  return res;
}

// Force/torque application
void
DEMTasks::applyContactForces(const DEMContactResult& contact,
                             particleIndex pidx_i,
                             int mat_i,
                             long64 rbID_j,
                             int mat_j,
                             const DEMParticleInputData& inputs,
                             const MPMMaterial* matl_j,
                             const Patch* patch,
                             const ParticleIDToCurrentIdxMap& master_particles,
                             DEMParticleOutputData& outputs)
{
  DOUT(dem_doing, "[DEMTasks::applyContactForces]");
  DOUT(dem_dbg, "Collision detected between mat " << mat_i
             << " particle " << pidx_i << " and mat " << mat_j
             << " body " << rbID_j);

  outputs.pExtForce_new[mat_i][pidx_i] += contact.totalForce;
  outputs.pTorque_new[mat_i][pidx_i]   += Cross(contact.arm_i, contact.totalForce);

  if (matl_j->isDEMMaterial()) {
    if (master_particles.count(rbID_j)) {
      int master_idx = master_particles.at(rbID_j);
      if (patch->containsPoint(inputs.pX_old[mat_j][master_idx])) {
        outputs.pExtForce_new[mat_j][master_idx] -= contact.totalForce;
        outputs.pTorque_new[mat_j][master_idx]   -= Cross(contact.arm_j_center, contact.totalForce);
      }
    }
  } else {
    // Find the real particle index for pidx_j
    // (caller passes the representative pidx_j from unique_bodies)
    int pidx_j_real = (int)rbID_j; // non-discrete: rbID == pidx_j (see buildUniqueBodies)
    if (patch->containsPoint(inputs.pX_old[mat_j][pidx_j_real])) {
      outputs.pExtForce_new[mat_j][pidx_j_real] -= contact.totalForce;
    }
  }
}

void
DEMTasks::scheduleIntegrateDEMRotation(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  DOUT(dem_doing, "[DEMTasks::scheduleIntegrateDEMRotation]");
  Task* t = scinew Task(
    "DEMTasks::integrateDEMRotation", this, &DEMTasks::integrateDEMRotation);

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
  DOUT(dem_doing, "[DEMTasks::integrateDEMRotation]");
  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int numMPMMatls    = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      if (!pset) {
        continue;
      }

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

      new_dw->allocateAndPut(pAngularVelocity_new,
                             d_mpm_labels->pAngularVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(
        pOrientation_new, d_mpm_labels->pOrientationLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pInertiaTensor_new, d_mpm_labels->pInertiaTensorLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pRadius_new, d_mpm_labels->pRadiusLabel_preReloc, pset);

      struct RotationState
      {
        Vector omega;
        Matrix3 orient;
        Matrix3 inertia;
      };
      std::map<long64, RotationState> master_rotations;

      if (mpm_matl->isDEMMaterial() || mpm_matl->isRigid()) {
        for (auto idx : *pset) {
          if (pRadius[idx] > 0.0 && pInertiaTensor[idx].Determinant() > 0.0) {
            Matrix3 invI              = pInertiaTensor[idx].Inverse();
            Vector angAcc             = invI * pTorque[idx];
            pAngularVelocity_new[idx] = pAngularVelocity[idx] + angAcc * delT;
            Matrix3 omegaSkew(0.0,
                              -pAngularVelocity[idx].z(),
                              pAngularVelocity[idx].y(),
                              pAngularVelocity[idx].z(),
                              0.0,
                              -pAngularVelocity[idx].x(),
                              -pAngularVelocity[idx].y(),
                              pAngularVelocity[idx].x(),
                              0.0);
            Matrix3 identity;
            identity.Identity();
            Matrix3 deltaR        = identity + omegaSkew * delT;
            pOrientation_new[idx] = deltaR * pOrientation[idx];
            pInertiaTensor_new[idx] =
              GeometryPiece::rotateInertiaTensor(pInertiaTensor[idx], deltaR);
            master_rotations[pRigidBodyID[idx]] = { pAngularVelocity_new[idx],
                                                    pOrientation_new[idx],
                                                    pInertiaTensor_new[idx] };
          } else {
            pAngularVelocity_new[idx] = pAngularVelocity[idx];
            pOrientation_new[idx]     = pOrientation[idx];
            pInertiaTensor_new[idx]   = pInertiaTensor[idx];
          }
          pRadius_new[idx] = pRadius[idx];
        }
        for (auto idx : *pset) {
          if (pRadius[idx] <= 0.0) {
            long64 rbID = pRigidBodyID[idx];
            if (master_rotations.count(rbID)) {
              pAngularVelocity_new[idx] = master_rotations[rbID].omega;
              pOrientation_new[idx]     = master_rotations[rbID].orient;
              pInertiaTensor_new[idx]   = master_rotations[rbID].inertia;
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

} // end namespace Vaango