/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Ports/ModelInterface.h>

#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SimulationInterface.h>

namespace Uintah {

ModelInterface::ModelInterface(const ProcessorGroup* my_world,
                               MaterialManagerP mat_manager)
  : UintahParallelComponent(my_world)
  , d_materialManager(mat_manager)
{
}

void
ModelInterface::setComponents(SimulationInterface* comp)
{
  SimulationInterface* parent = dynamic_cast<SimulationInterface*>(comp);

  setAMR(parent->isAMR());
  setDynamicRegridding(parent->isDynamicRegridding());

  attachPort("simulator", parent);
  attachPort("scheduler", parent->getScheduler());
  attachPort("regridder", parent->getRegridder());
  attachPort("output", parent->getOutput());

  getComponents();
}

void
ModelInterface::getComponents()
{
  d_simulator = dynamic_cast<SimulationInterface*>(getPort("simulator"));

  if (!d_simulator) {
    throw InternalError(
      "dynamic_cast of 'd_simulator' failed!", __FILE__, __LINE__);
  }

  d_scheduler = dynamic_cast<Scheduler*>(getPort("scheduler"));

  if (!d_scheduler) {
    throw InternalError(
      "dynamic_cast of 'd_scheduler' failed!", __FILE__, __LINE__);
  }

  d_regridder = dynamic_cast<Regridder*>(getPort("regridder"));

  if (isDynamicRegridding() && !d_regridder) {
    throw InternalError(
      "dynamic_cast of 'd_regridder' failed!", __FILE__, __LINE__);
  }

  d_output = dynamic_cast<Output*>(getPort("output"));

  if (!d_output) {
    throw InternalError(
      "dynamic_cast of 'd_output' failed!", __FILE__, __LINE__);
  }
}

void
ModelInterface::releaseComponents()
{
  releasePort("simulator");
  releasePort("scheduler");
  releasePort("regridder");
  releasePort("output");

  d_simulator = nullptr;
  d_scheduler = nullptr;
  d_regridder = nullptr;
  d_output    = nullptr;
}

} // namespace Uintah