/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Parent/ComponentFactory.h>

#include <CCA/Components/Parent/Switcher.h>

#include <CCA/Components/PostProcessUda/PostProcessUda.h>

#include <CCA/Components/Examples/AMRWave.h>
#include <CCA/Components/Examples/Benchmark.h>
#include <CCA/Components/Examples/Burger.h>
#include <CCA/Components/Examples/ParticleTest1.h>
#include <CCA/Components/Examples/Poisson1.h>
#include <CCA/Components/Examples/Poisson2.h>
#include <CCA/Components/Examples/Poisson3.h>
#include <CCA/Components/Examples/Poisson4.h>
#include <CCA/Components/Examples/RegridderTest.h>
#include <CCA/Components/Examples/SolverTest1.h>
#include <CCA/Components/Examples/Wave.h>
#include <CCA/Components/ICE/AMRICE.h>
#include <CCA/Components/ICE/ICE.h>
#include <CCA/Components/ICE/impAMRICE.h>
#include <CCA/Components/MPM/AMRMPM.h>
#include <CCA/Components/MPM/FractureMPM.h>
#include <CCA/Components/MPM/ImpMPM.h>
#include <CCA/Components/MPM/MPM_UpdateStressLast.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/ShellMPM.h>
#include <CCA/Components/MPM/UofU_MPM.h>
#include <CCA/Components/MPMICE/MPMICE.h>
#include <CCA/Components/Peridynamics/Peridynamics.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <sci_defs/cuda_defs.h>
#include <sci_defs/uintah_defs.h>

#ifdef HAVE_CUDA
#include <CCA/Components/Examples/GPUSchedulerTest.h>
#include <CCA/Components/Examples/PoissonGPU1.h>
#include <CCA/Components/Examples/UnifiedSchedulerTest.h>
#endif

#include <iosfwd>
#include <string>

#include <locale>

using namespace Uintah;

std::unique_ptr<UintahParallelComponent>
ComponentFactory::create(ProblemSpecP& ps,
                         const ProcessorGroup* world,
                         const MaterialManagerP& mat_manager,
                         string uda)
{
  // Check for an AMR attribute with the grid.
  bool doAMR           = false;
  ProblemSpecP grid_ps = ps->findBlock("Grid");
  if (grid_ps) {
    grid_ps->getAttribute("doAMR", doAMR);
  }
  // If the AMR block is defined default to turning AMR on.
  ProblemSpecP amr_ps = ps->findBlock("AMR");
  if (amr_ps) {
    doAMR = true;
  }

  std::string sim_comp;

  ProblemSpecP sim_ps = ps->findBlock("SimulationComponent");
  if (sim_ps) {
    sim_ps->getAttribute("type", sim_comp);
  } else {
    // This is probably a <subcomponent>, so the name of the type of
    // the component is in a different place:
    ps->getAttribute("type", sim_comp);
  }
  if (sim_comp == "") {
    throw ProblemSetupException(
      "Could not determine the type of SimulationComponent...",
      __FILE__,
      __LINE__);
  }
  std::transform(sim_comp.begin(), sim_comp.end(), sim_comp.begin(), ::tolower);

  proc0cout << "Simulation Component: \t'" << sim_comp << "'\n";

  std::string turned_off_options;

  proc0cout << "Simulation Component: \t'" << sim_comp << "'\n";
  if (sim_comp == "peridynamics") {
    return std::make_unique<Vaango::Peridynamics>(world, mat_manager);
  }

#ifndef NO_MPM
  if (sim_comp == "mpm" || sim_comp == "MPM") {
    return std::make_unique<SerialMPM>(world, mat_manager);
  }
  if (sim_comp == "mpm_usl" || sim_comp == "MPM_USL") {
    return std::make_unique<MPM_UpdateStressLast>(world, mat_manager);
  }
  if (sim_comp == "uofu_mpm" || sim_comp == "UofU_MPM") {
    return std::make_unique<UofU_MPM>(world, mat_manager);
  }
  if (sim_comp == "mpmf" || sim_comp == "fracturempm" ||
      sim_comp == "FRACTUREMPM") {
    return std::make_unique<FractureMPM>(world, mat_manager);
  }
  if (sim_comp == "rmpm" || sim_comp == "rigidmpm" || sim_comp == "RIGIDMPM") {
    return std::make_unique<RigidMPM>(world, mat_manager);
  }
  if (sim_comp == "amrmpm" || sim_comp == "AMRmpm" || sim_comp == "AMRMPM") {
    return std::make_unique<AMRMPM>(world, mat_manager);
  }
  if (sim_comp == "smpm" || sim_comp == "shellmpm" || sim_comp == "SHELLMPM") {
    return std::make_unique<ShellMPM>(world, mat_manager);
  }
  if (sim_comp == "impm" || sim_comp == "IMPM") {
    return std::make_unique<ImpMPM>(world, mat_manager);
  }
#else
  turned_off_options += "MPM ";
#endif
#ifndef NO_ICE
  if (sim_comp == "ice" || sim_comp == "ICE") {
    ProblemSpecP cfd_ps   = ps->findBlock("CFD");
    ProblemSpecP ice_ps   = cfd_ps->findBlock("ICE");
    ProblemSpecP imp_ps   = ice_ps->findBlock("ImplicitSolver");
    bool doImplicitSolver = (imp_ps);

    if (doAMR) {
      if (doImplicitSolver) {
        return std::make_unique<impAMRICE>(world, mat_manager);
      } else {
        return std::make_unique<AMRICE>(world, mat_manager);
      }
    } else {
      return std::make_unique<ICE>(world, mat_manager);
    }
  }
#else
  turned_off_options += "ICE ";
#endif
#if !defined(NO_MPM) && !defined(NO_ICE)
  if (sim_comp == "mpmice" || sim_comp == "MPMICE") {
    return std::make_unique<MPMICE>(world, mat_manager, MPMType::STAND_MPMICE, doAMR);
  }
  if (sim_comp == "smpmice" || sim_comp == "shellmpmice" ||
      sim_comp == "SHELLMPMICE") {
    return std::make_unique<MPMICE>(world, mat_manager, MPMType::SHELL_MPMICE, doAMR);
  }
  if (sim_comp == "rmpmice" || sim_comp == "rigidmpmice" ||
      sim_comp == "RIGIDMPMICE") {
    return std::make_unique<MPMICE>(world, mat_manager, MPMType::RIGID_MPMICE, doAMR);
  }
#else
  turned_off_options += "MPMICE ";
#endif
  if (sim_comp == "burger" || sim_comp == "BURGER") {
    return std::make_unique<Burger>(world, mat_manager);
  }
  if (sim_comp == "wave" || sim_comp == "WAVE") {
    if (doAMR) {
      return std::make_unique<AMRWave>(world, mat_manager);
    } else {
      return std::make_unique<Wave>(world, mat_manager);
    }
  }
  if (sim_comp == "poisson1" || sim_comp == "POISSON1") {
    return std::make_unique<Poisson1>(world, mat_manager);
  }

#ifdef HAVE_CUDA
  if (sim_comp == "poissongpu1" || sim_comp == "POISSONGPU1") {
    return std::make_unique<PoissonGPU1>(world, mat_manager);
  }
  if (sim_comp == "gpuschedulertest" || sim_comp == "GPUSCHEDULERTEST") {
    return std::make_unique<GPUSchedulerTest>(world, mat_manager);
  }
  if (sim_comp == "unifiedschedulertest" ||
      sim_comp == "UNIFIEDSCHEDULERTEST") {
    return std::make_unique<UnifiedSchedulerTest>(world, mat_manager);
  }
#endif

  if (sim_comp == "regriddertest" || sim_comp == "REGRIDDERTEST") {
    return std::make_unique<RegridderTest>(world, mat_manager);
  }
  if (sim_comp == "poisson2" || sim_comp == "POISSON2") {
    return std::make_unique<Poisson2>(world, mat_manager);
  }
  if (sim_comp == "poisson3" || sim_comp == "POISSON3") {
    return std::make_unique<Poisson3>(world, mat_manager);
  }
  if (sim_comp == "poisson4" || sim_comp == "POISSON4") {
    return std::make_unique<Poisson4>(world, mat_manager);
  }
  if (sim_comp == "benchmark" || sim_comp == "BENCHMARK") {
    return std::make_unique<Benchmark>(world, mat_manager);
  }
  if (sim_comp == "particletest" || sim_comp == "PARTICLETEST") {
    return std::make_unique<ParticleTest1>(world, mat_manager);
  }
  if (sim_comp == "solvertest" || sim_comp == "SOLVERTEST") {
    return std::make_unique<SolverTest1>(world, mat_manager);
  }
#ifdef HAVE_HYPRE
  if (sim_comp == "solvertest2" || sim_comp == "SOLVERTEST2") {
    return std::make_unique<SolverTest2>(world, mat_manager);
  }
#endif
  if (sim_comp == "switcher" || sim_comp == "SWITCHER") {
    return std::make_unique<Switcher>(world, mat_manager, ps, uda);
  }
  if (sim_comp == "postprocessUda") {
    return std::make_unique<PostProcessUda>(world, mat_manager, uda);
  }
  throw ProblemSetupException(
    "Unknown simulationComponent ('" + sim_comp +
      "'). Must specify -ice, -mpm, -mpm_usl"
      "-impm, -mpmice, -burger, -wave, -poisson1, -poisson2, -poisson3 or "
      "-benchmark.\n"
      "Note: the following components were turned off at configure time: " +
      turned_off_options +
      "\n"
      "Make sure that the requested component is supported in this build.",
    __FILE__,
    __LINE__);
}
