/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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

#ifndef VAANGO_CCA_COMPONENTS_MPM_FRACTUREMPM_H
#define VAANGO_CCA_COMPONENTS_MPM_FRACTUREMPM_H

#include <CCA/Components/MPM/SerialMPM.h>

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Crack/Crack.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>

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

class Crack;
class ThermalContact;
class HeatConduction;

class FractureMPM : public SerialMPM
{
public:
  std::unique_ptr<Crack> crackModel{nullptr}; // for Fracture

public:
  FractureMPM(const ProcessorGroup* myworld,
              const MaterialManagerP& mat_manager);

  virtual ~FractureMPM() noexcept(false) override = default;

  // No copy or move allowed
  FractureMPM(const FractureMPM&) = delete;
  FractureMPM(FractureMPM&&)      = delete;
  auto
  operator=(const FractureMPM&) -> FractureMPM& = delete;
  auto
  operator=(FractureMPM&&) -> FractureMPM& = delete;

  void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "") override;

  void
  outputProblemSpec(ProblemSpecP& ps) override;

  void
  scheduleInitialize(const LevelP& level, SchedulerP&) override;

  void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&) override;

  void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&) override;

protected:
  void
  scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls) override;

  void
  scheduleMomentumExchangeInterpolated(SchedulerP&,
                                       const PatchSet*,
                                       const MaterialSet*) override;

  void
  scheduleComputeStressTensor(SchedulerP&, const PatchSet*, const MaterialSet*) override;

  // for thermal stress analysis
  void
  scheduleComputeParticleTempFromGrid(SchedulerP&,
                                      const PatchSet*,
                                      const MaterialSet*);

  void
  scheduleComputeAccStrainEnergy(SchedulerP&,
                                 const PatchSet*,
                                 const MaterialSet*);

  void
  scheduleComputeContactArea(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleComputeInternalForce(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*) override;

  void
  scheduleComputeArtificialViscosity(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  scheduleComputeAndIntegrateAcceleration(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*) override;

  void
  scheduleMomentumExchangeIntegrated(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*) override;

  void
  scheduleSetGridBoundaryConditions(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet* matls);

  void
  scheduleApplyExternalLoads(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleInterpolateToParticlesAndUpdate(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*) override;

  // for Fracture ----------------------------------
  void
  scheduleParticleVelocityField(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls);

  void
  scheduleAdjustCrackContactInterpolated(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls);

  void
  scheduleAdjustCrackContactIntegrated(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls);

  void
  scheduleCalculateFractureParameters(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls);

  void
  scheduleDoCrackPropagation(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls);

  void
  scheduleMoveCracks(SchedulerP& sched,
                     const PatchSet* patches,
                     const MaterialSet* matls);

  void
  scheduleUpdateCrackFront(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSet* matls);

protected:
  friend class MPMICE;

  virtual void
  actuallyInitialize(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  actuallyComputeStableTimestep(const ProcessorGroup*,
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

  virtual void
  computeStressTensor(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  virtual void
  computeParticleTempFromGrid(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  void
  computeAccStrainEnergy(const ProcessorGroup*,
                         const PatchSubset*,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  virtual void
  computeContactArea(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  virtual void
  computeInternalForce(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  computeArtificialViscosity(const ProcessorGroup*,
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

  //////////
  // This task is to be used for setting particle external force
  // and external heat rate.  I'm creating a separate task so that
  // user defined schemes for setting these can be implemented without
  // editing the core routines
  void
  applyExternalLoads(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw) override;

private:

  // Used functions for fracture 
};

} // end namespace Uintah

#endif
