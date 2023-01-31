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

#include <CCA/Components/MPM/ThermalContact/NullThermalContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Malloc/Allocator.h>

#include <vector>

using namespace Uintah;

NullThermalContact::NullThermalContact(ProblemSpecP&,
                                       const MaterialManagerP& mat_manager,
                                       const MPMLabel* labels,
                                       const MPMFlags* flags)
  : ThermalContact(mat_manager, labels, flags)
{
}

void
NullThermalContact::outputProblemSpec(ProblemSpecP& ps)
{
}

void
NullThermalContact::computeHeatExchange(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse*,
                                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int numMatls = d_mat_manager->getNumMaterials("MPM");

    for (int m = 0; m < numMatls; m++) {
      int dwindex = d_mat_manager->getMaterial("MPM", m)->getDWIndex();

      NCVariable<double> thermalContactTemperatureRate;
      new_dw->allocateAndPut(thermalContactTemperatureRate,
                             d_mpm_labels->gThermalContactTemperatureRateLabel,
                             dwindex,
                             patch);

      thermalContactTemperatureRate.initialize(0);
      NCVariable<double> GthermalContactTemperatureRate;
      if (d_mpm_flags->d_fracture) {
        new_dw->allocateAndPut(
          GthermalContactTemperatureRate,
          d_mpm_labels->GThermalContactTemperatureRateLabel,
          dwindex,
          patch);
        GthermalContactTemperatureRate.initialize(0);
      }
    }
  }
}

void
NullThermalContact::initializeThermalContact(const Patch* /*patch*/,
                                             int /*vfindex*/,
                                             DataWarehouse* /*new_dw*/)
{
}

void
NullThermalContact::addComputesAndRequires(Task* t,
                                           const PatchSet*,
                                           const MaterialSet*) const
{
  t->computes(d_mpm_labels->gThermalContactTemperatureRateLabel);
  if (d_mpm_flags->d_fracture) {
    t->computes(d_mpm_labels->GThermalContactTemperatureRateLabel);
  }
}
