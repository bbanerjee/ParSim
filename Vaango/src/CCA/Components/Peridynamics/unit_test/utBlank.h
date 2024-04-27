/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef Vaango_Peridynamics_ut_utBlank_h
#define Vaango_Peridynamics_ut_utBlank_h

#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/GradientComputer/PeridynamicsDefGradComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/BondInternalForceComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/ParticleInternalForceComputer.h>
#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>

#include <CCA/Components/MPM/Contact/Contact.h>

#include <Core/Grid/MaterialManagerP.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Parallel/UintahParallelComponent.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>

namespace Vaango {
  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   utBlank
    \brief   unit test for peridynamics particles
    \author  Bryan Smith
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class utBlank : public Uintah::UintahParallelComponent, public Uintah::SimulationInterface {
  public:
    utBlank(const Uintah::ProcessorGroup* myworld);
    virtual ~utBlank();

    virtual void problemSetup(const Uintah::ProblemSpecP& params, 
                              const Uintah::ProblemSpecP& restart_prob_spec, 
                              Uintah::GridP& grid, 
                              Uintah::MaterialManagerP& mat_manager);
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

    Uintah::MaterialManagerP 
 d_mat_manager;
    PeridynamicsLabel* d_labels;
    PeridynamicsFlags* d_flags;
    Uintah::EmptyMaterial* d_mymat;

    int d_doOutput;
    int d_numGhostCells;
    double d_delT;

    utBlank(const utBlank&);
    utBlank& operator=(const utBlank&);
  };
}

#endif
