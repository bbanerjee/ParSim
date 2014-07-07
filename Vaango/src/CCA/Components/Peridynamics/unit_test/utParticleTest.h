#ifndef Vaango_Peridynamics_ut_utParticleTest_h
#define Vaango_Peridynamics_ut_utParticleTest_h

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/GradientComputer/PeridynamicsDefGradComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/BondInternalForceComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/ParticleInternalForceComputer.h>
#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>

#include <CCA/Components/MPM/Contact/Contact.h>

#include <Core/Grid/SimulationStateP.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Parallel/UintahParallelComponent.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/SimpleMaterial.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/ParticleInterpolator.h>

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>

namespace Vaango {
  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   utParticleTest
    \brief   unit test for peridynamics particles
    \author  Bryan Smith
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class utParticleTest : public Uintah::UintahParallelComponent, public Uintah::SimulationInterface {
  public:
    utParticleTest(const Uintah::ProcessorGroup* myworld);
    virtual ~utParticleTest();

    virtual void problemSetup(const Uintah::ProblemSpecP& params, 
                              const Uintah::ProblemSpecP& restart_prob_spec, 
                              Uintah::GridP& grid, 
                              Uintah::SimulationStateP& state);
    virtual void scheduleInitialize(const Uintah::LevelP& level,
				                    Uintah::SchedulerP& sched);
    virtual void scheduleComputeStableTimestep(const Uintah::LevelP& level,
					                           Uintah::SchedulerP&);
    virtual void scheduleTimeAdvance( const Uintah::LevelP& level, 
				                      Uintah::SchedulerP&);

  private:
    void initialize(const Uintah::ProcessorGroup*,
		            const Uintah::PatchSubset* patches, 
                    const Uintah::MaterialSubset* matls,
		            Uintah::DataWarehouse* old_dw, 
                    Uintah::DataWarehouse* new_dw);

    void computeStableTimestep(const Uintah::ProcessorGroup*,
			                   const Uintah::PatchSubset* patches,
			                   const Uintah::MaterialSubset* matls,
			                   Uintah::DataWarehouse* old_dw, 
                               Uintah::DataWarehouse* new_dw);

    void timeAdvance(const Uintah::ProcessorGroup*,
            	     const Uintah::PatchSubset* patches,
		             const Uintah::MaterialSubset* matls,
		             Uintah::DataWarehouse* old_dw, 
                     Uintah::DataWarehouse* new_dw);

    Uintah::SimulationStateP d_sharedState;
    PeridynamicsLabel* d_labels;
    PeridynamicsFlags* d_flags;
    Uintah::SimpleMaterial* d_mymat;

    int d_doOutput;
    int d_numGhostCells;
    double d_delT;

    utParticleTest(const utParticleTest&);
    utParticleTest& operator=(const utParticleTest&);
  };
}

#endif
