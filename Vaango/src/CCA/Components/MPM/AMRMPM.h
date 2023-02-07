/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef UINTAH_HOMEBREW_AMRMPM_H
#define UINTAH_HOMEBREW_AMRMPM_H

#include <CCA/Components/MPM/SerialMPM.h>

#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class GeometryObject;
class AMRMPMLabel;

class AMRMPM : public SerialMPM
{

public:
  enum class CoarsenFlag
  {
    coarsenData,
    zeroData,
  };

public:
  AMRMPM(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);

  virtual ~AMRMPM();

  AMRMPM(const AMRMPM&) = delete;
  AMRMPM(AMRMPM&&)      = delete;
  AMRMPM&
  operator=(const AMRMPM&) = delete;
  AMRMPM&
  operator=(AMRMPM&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid);

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

  void
  schedulePrintParticleCount(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  scheduleFinalizeTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleAnalysis(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRefine(const PatchSet* patches, SchedulerP& scheduler);

  virtual void
  scheduleRefineInterface(const LevelP& fineLevel,
                          SchedulerP& scheduler,
                          bool needCoarse,
                          bool needFine);

  virtual void
  scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched);

  /// Schedule to mark flags for AMR regridding
  virtual void
  scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched);

  /// Schedule to mark initial flags for AMR regridding
  void
  scheduleInitialErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched);

protected:
  virtual void
  actuallyInitialize(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  printParticleCount(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  printParticleLabels(std::vector<const VarLabel*> label,
                      DataWarehouse* dw,
                      int dwi,
                      const Patch* patch);

  void
  actuallyComputeStableTimestep(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

  void
  partitionOfUnity(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset*,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  virtual void
  computeZoneOfInfluence(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* matls,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  virtual void
  interpolateParticlesToGrid(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  // At Coarse Fine interface
  void
  interpolateParticlesToGrid_CFI(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw);

  void
  interpolateParticlesToGrid_CFI_GIMP(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* matls,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw);

  void
  coarsenNodalData_CFI(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw,
                       const CoarsenFlag flag);

  void
  coarsenNodalData_CFI2(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  void
  normalizeNodalVelTempConc(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* matls,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  virtual void
  computeStressTensor(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  void
  updateErosionParameter(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  virtual void
  computeInternalForce(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  computeInternalForce_CFI(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  virtual void
  computeAndIntegrateAcceleration(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  void
  setGridBoundaryConditions(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  // Compute Vel. Grad and Def Grad
  void
  computeLAndF(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  virtual void
  interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);
#if 0
  // At Coarse Fine interface
  void interpolateToParticlesAndUpdate_CFI(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset* matls,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw);
#endif

  // Used to compute the particles initial physical size
  // for use in deformed particle visualization
  virtual void
  computeParticleScaleFactor(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  // Update particle quantities only
  virtual void
  finalParticleUpdate(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  //////////
  // Add new particles to the simulation based on criteria TBD:
  virtual void
  addParticles(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  ////////
  // Find the extents of MPMRefineCell for use in generating rectangular
  // patches
  virtual void
  reduceFlagsExtents(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  refineGrid(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse*,
             DataWarehouse* new_dw);

  void
  coarsen(const ProcessorGroup*,
          const PatchSubset* patches,
          const MaterialSubset* matls,
          DataWarehouse*,
          DataWarehouse* new_dw);

  void
  errorEstimate(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                DataWarehouse*,
                DataWarehouse* new_dw);

  void
  initialErrorEstimate(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse*,
                       DataWarehouse* new_dw);

  //______________________________________________________________________
  //
  void
  schedulePartitionOfUnity(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeZoneOfInfluence(SchedulerP&,
                                 const PatchSet*,
                                 const MaterialSet*);

  virtual void
  scheduleInterpolateParticlesToGrid(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  scheduleInterpolateParticlesToGrid_CFI(SchedulerP&,
                                         const PatchSet*,
                                         const MaterialSet*);

  void
  scheduleCoarsenNodalData_CFI(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*,
                               const CoarsenFlag flag);

  void
  scheduleCoarsenNodalData_CFI2(SchedulerP&,
                                const PatchSet*,
                                const MaterialSet*);

  void
  scheduleNormalizeNodalVelTempConc(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet*);

  virtual void
  scheduleComputeStressTensor(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleUpdateErosionParameter(SchedulerP& sched,
                                 const PatchSet* patches,
                                 const MaterialSet* matls);

  virtual void
  scheduleComputeInternalForce(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  void
  scheduleComputeInternalForce_CFI(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls);

  virtual void
  scheduleComputeAndIntegrateAcceleration(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

  void
  scheduleSetGridBoundaryConditions(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet* matls);

  void
  scheduleComputeLAndF(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleInterpolateToParticlesAndUpdate(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

#if 0
  void scheduleInterpolateToParticlesAndUpdate_CFI(SchedulerP&, 
                                                   const PatchSet*,
                                                   const MaterialSet*);
#endif

  virtual void
  scheduleComputeParticleScaleFactor(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  virtual void
  scheduleFinalParticleUpdate(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleAddParticles(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleReduceFlagsExtents(SchedulerP&, const PatchSet*, const MaterialSet*);

  //  count the total number of particles in the domain
  void
  scheduleCountParticles(const PatchSet* patches, SchedulerP& sched);

  void
  countParticles(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  scheduleDebug_CFI(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  debug_CFI(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset*,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

  // input coarse patches and return coarse & fine level patches with CFI
  void
  coarseLevelCFI_Patches(const PatchSubset* coarsePatches,
                         Level::selectType& CFI_coarsePatches,
                         Level::selectType& CFI_finePatches);

  // input fine patches and return coarse & fine level patches with CFI
  void
  fineLevelCFI_Patches(const PatchSubset* finePatches,
                       Level::selectType& CFI_coarsePatches,
                       Level::selectType& CFI_finePatches);

  // remove duplicate entries in array  This belongs in Core/Grid/fixedvector.h
  void
  removeDuplicates(Level::selectType& array);

protected:
  // Number of cells on the coarse level that
  // contain particles and surround a fine patch.
  // Coarse level particles are used in the task
  // interpolateToParticlesAndUpdate_CFI.
  int d_nPaddingCells_Coarse{ -9 };

  // debugging code used to check the answers (acceleration)
  Vector d_acc_ans{ 0.0, 0.0, 0.0 };
  double d_acc_tol{ 1.0e-7 };
  // debugging code used to check the answers (velocity)
  Vector d_vel_ans{ -100.0, 0.0, 0.0 };
  double d_vel_tol{ 1.0e-7 };

  const VarLabel* pDbgLabel; // debugging labels
  const VarLabel* gSumSLabel;
  const VarLabel* gZOINETLabel{ nullptr };
  const VarLabel* gZOISWBLabel{ nullptr };
  const VarLabel* RefineFlagXMaxLabel;
  const VarLabel* RefineFlagXMinLabel;
  const VarLabel* RefineFlagYMaxLabel;
  const VarLabel* RefineFlagYMinLabel;
  const VarLabel* RefineFlagZMaxLabel;
  const VarLabel* RefineFlagZMinLabel;

  std::vector<MPMPhysicalBC*> d_physicalBCs;
  std::unique_ptr<ScalarDiffusionTasks> d_diffusionTasks{ nullptr };
  SwitchingCriteria* d_switchCriteria{ nullptr };

private:
  std::string d_CFI_interpolator; // user can override interpolator at CFI

  std::unique_ptr<AMRMPMLabel> d_amrmpmLabels{ nullptr };

  // matlsubset for zone of influence
  std::unique_ptr<MaterialSubset> d_oneMaterial{ nullptr };

  std::vector<std::unique_ptr<GeometryObject>> d_refineGeomObjs;

  // refinement criteria threshold knobs
  struct thresholdVar
  {
    std::string name;
    int matl;
    double value;
  };
  std::vector<thresholdVar> d_thresholdVars;

  //______________________________________________________________________
  // Bulletproofing machinery to keep track of
  // how many nodes on a CFI have been 'touched'
  struct faceMarks
  {
    int marks[Patch::numFaces];
    int&
    operator[](Patch::FaceType face)
    {
      return marks[static_cast<int>(face)];
    }
    faceMarks()
    {
      marks[0] = 0;
      marks[1] = 0;
      marks[2] = 0;
      marks[3] = 0;
      marks[4] = 0;
      marks[5] = 0;
    }
  };
  std::map<const Patch*, faceMarks> faceMarks_map[2];

  inline void
  clearFaceMarks(const int whichMap, const Patch* patch)
  {
    faceMarks_map[whichMap].erase(patch);
  }

  inline int
  getFaceMark(int whichMap, const Patch* patch, Patch::FaceType face)
  {
    ASSERT(whichMap >= 0 && whichMap < 2);
    return faceMarks_map[whichMap][patch][face];
  };

  inline void
  setFaceMark(int whichMap, const Patch* patch, Patch::FaceType face, int value)
  {
    ASSERT(whichMap >= 0 && whichMap < 2);
    faceMarks_map[whichMap][patch][face] = value;
  };

  inline void
  computeVelocityGradient(Matrix3& velGrad,
                          std::vector<IntVector>& ni,
                          std::vector<Vector>& d_S,
                          const double* oodx,
                          constNCVariable<Vector>& gVelocity)
  {
    for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
      const Vector& gvel = gVelocity[ni[k]];
      for (int j = 0; j < 3; j++) {
        double d_SXoodx = d_S[k][j] * oodx[j];
        for (int i = 0; i < 3; i++) {
          velGrad(i, j) += gvel[i] * d_SXoodx;
        }
      }
    }
  };

};

} // end namespace Uintah

#endif
