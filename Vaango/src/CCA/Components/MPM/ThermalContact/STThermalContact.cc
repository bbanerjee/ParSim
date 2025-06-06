/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <CCA/Components/MPM/ThermalContact/STThermalContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <vector>

using namespace Uintah;

STThermalContact::STThermalContact(ProblemSpecP&,
                                   const MaterialManagerP& mat_manager,
                                   const MPMLabel* labels,
                                   const MPMFlags* flags)
  : ThermalContact(mat_manager, labels, flags)
{
}

void
STThermalContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP thermal_ps = ps->appendChild("thermal_contact");
}

void
STThermalContact::computeHeatExchange(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int numMatls = d_mat_manager->getNumMaterials("MPM");

    constNCdoubleArray gMass(numMatls);
    constNCdoubleArray gTemp(numMatls);
    NCdoubleArray thermalContactTemperatureRate(numMatls);
    std::vector<double> Cp(numMatls);
    // for Fracture (additional field)-----------------------------------------
    constNCdoubleArray Gmass(numMatls);
    constNCdoubleArray GTemp(numMatls);
    NCdoubleArray GthermalContactTemperatureRate(numMatls);

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      new_dw->get(
        gMass[dwi], d_mpm_labels->gMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(gTemp[dwi],
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->allocateAndPut(thermalContactTemperatureRate[dwi],
                             d_mpm_labels->gThermalContactTemperatureRateLabel,
                             dwi,
                             patch);
      thermalContactTemperatureRate[dwi].initialize(0.);
      Cp[m] = mpm_matl->getSpecificHeat();
      if (d_mpm_flags->d_fracture) {
        // for Fracture (for additional field)----------------------------------
        new_dw->get(
          Gmass[dwi], d_mpm_labels->GMassLabel, dwi, patch, Ghost::None, 0);
        new_dw->get(GTemp[dwi],
                    d_mpm_labels->GTemperatureLabel,
                    dwi,
                    patch,
                    Ghost::None,
                    0);
        new_dw->allocateAndPut(
          GthermalContactTemperatureRate[dwi],
          d_mpm_labels->GThermalContactTemperatureRateLabel,
          dwi,
          patch);
        GthermalContactTemperatureRate[dwi].initialize(0);
      }
      // -------------------------------------------------------------------
    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      double numerator   = 0.0;
      double denominator = 0.0;
      IntVector c        = *iter;
      for (int m = 0; m < numMatls; m++) {
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
        int n = mpm_matl->getDWIndex();
        numerator += (gTemp[n][c] * gMass[n][c] * Cp[m]);
        denominator += (gMass[n][c] * Cp[m]);
        if (d_mpm_flags->d_fracture) {
          numerator += GTemp[n][c] * Gmass[n][c] * Cp[m];
          denominator += Gmass[n][c] * Cp[m]; // add in second field;
        }
      }

      double contactTemperature = numerator / denominator;

      for (int m = 0; m < numMatls; m++) {
        thermalContactTemperatureRate[m][c] =
          (contactTemperature - gTemp[m][c]) / delT;
        if (d_mpm_flags->d_fracture) {
          GthermalContactTemperatureRate[m][c] =
            (contactTemperature - GTemp[m][c]) / delT;
        }
      }
    }
  }
}

void
STThermalContact::initializeThermalContact(const Patch* /*patch*/,
                                           int /*vfindex*/,
                                           DataWarehouse* /*new_dw*/)
{
}

void
STThermalContact::addComputesAndRequires(Task* t,
                                         const PatchSet*,
                                         const MaterialSet*) const
{
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);
  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gTemperatureLabel, Ghost::None);
  t->computes(d_mpm_labels->gThermalContactTemperatureRateLabel);
  if (d_mpm_flags->d_fracture) {
    // for second field, for Fracture ---------------------------------
    t->needs(Task::NewDW, d_mpm_labels->GMassLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->GTemperatureLabel, Ghost::None);
    t->computes(d_mpm_labels->GThermalContactTemperatureRateLabel);
  }
}
