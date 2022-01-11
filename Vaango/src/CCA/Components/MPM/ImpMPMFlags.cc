/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ImpMPMFlags.h>
#include <Core/Util/DebugStream.h>

using namespace Uintah;

static DebugStream dbg("ImpMPMFlags", false);

ImpMPMFlags::ImpMPMFlags(const ProcessorGroup* myworld)
  : MPMFlags(myworld)
{

  d_convCritDisp         = 1.e-10;
  d_convCritEnergy       = 4.e-10;
  d_forceIncrementFactor = 1.0;
  d_integrator           = Implicit;
  d_dynamic              = true;

  d_doMechanics       = true;
  d_projectHeatSource = false;

  d_tempSolve                              = false;
  d_interpolateParticleTempToGridEveryStep = true;
}

void
ImpMPMFlags::readMPMFlags(ProblemSpecP& ps, Output* dataArchive)
{
  MPMFlags::readMPMFlags(ps, dataArchive);

  ProblemSpecP root        = ps->getRootNode();
  ProblemSpecP mpm_flag_ps = root->findBlock("MPM");

  if (!mpm_flag_ps)
    return;

  mpm_flag_ps->get("project_heat_source", d_projectHeatSource);
  mpm_flag_ps->get("do_mechanics", d_doMechanics);
  mpm_flag_ps->get("convergence_criteria_disp", d_convCritDisp);
  mpm_flag_ps->get("convergence_criteria_energy", d_convCritEnergy);
  mpm_flag_ps->get("dynamic", d_dynamic);
  mpm_flag_ps->getWithDefault(
    "iters_before_timestep_restart", d_maxNumIterations, 25);
  mpm_flag_ps->getWithDefault(
    "num_iters_to_decrease_delT", d_numItersToDecreaseDelT, 12);
  mpm_flag_ps->getWithDefault(
    "num_iters_to_increase_delT", d_numItersToIncreaseDelT, 4);
  mpm_flag_ps->getWithDefault("delT_decrease_factor", d_delTDecreaseFactor, .6);
  mpm_flag_ps->getWithDefault("delT_increase_factor", d_delTIncreaseFactor, 2);

  mpm_flag_ps->get("solver", d_solverType);
  mpm_flag_ps->get("temperature_solve", d_tempSolve);
  mpm_flag_ps->get("interpolate_particle_temp_to_grid_every_step",
                   d_interpolateParticleTempToGridEveryStep);
}

void
ImpMPMFlags::outputProblemSpec(ProblemSpecP& ps)
{

  MPMFlags::outputProblemSpec(ps);

  ps->appendElement("project_heat_source", d_projectHeatSource);
  ps->appendElement("do_mechanics", d_doMechanics);
  ps->appendElement("convergence_criteria_disp", d_convCritDisp);
  ps->appendElement("convergence_criteria_energy", d_convCritEnergy);
  ps->appendElement("dynamic", d_dynamic);
  ps->appendElement("iters_before_timestep_restart", d_maxNumIterations);
  ps->appendElement("num_iters_to_decrease_delT", d_numItersToDecreaseDelT);
  ps->appendElement("num_iters_to_increase_delT", d_numItersToIncreaseDelT);
  ps->appendElement("delT_decrease_factor", d_delTDecreaseFactor);
  ps->appendElement("delT_increase_factor", d_delTIncreaseFactor);

  ps->appendElement("solver", d_solverType);
  ps->appendElement("temperature_solve", d_tempSolve);
  ps->appendElement("interpolate_particle_temp_to_grid_every_step",
                    d_interpolateParticleTempToGridEveryStep);
}
