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

#include <CCA/Components/SwitchingCriteria/SimpleBurn.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

namespace Uintah {

extern DebugStream switching_dbg;

SimpleBurnCriteria::SimpleBurnCriteria(ProblemSpecP& ps)
{
  ps->require("reactant_material", d_material);
  ps->require("ThresholdTemperature", d_temperature);

  proc0cout << "Switching criteria:  \tSimpleBurn, reactant matl: "
            << d_material << " Threshold tempterature " << d_temperature
            << std::endl;

  d_mpm_labels    = std::make_unique<MPMLabel>();
  d_mpmice_labels = std::make_unique<MPMICELabel>();
}

void
SimpleBurnCriteria::problemSetup([[maybe_unused]] const ProblemSpecP& ps,
                                 [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                                 MaterialManagerP& mat_manager)
{
  d_mat_manager = mat_manager;
}
//__________________________________
//
void
SimpleBurnCriteria::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level,
                switching_dbg,
                "Switching Criteria:SimpleBurnCriteria::scheduleSwitchTest");

  Task* t = scinew Task("switchTest", this, &SimpleBurnCriteria::switchTest);

  auto* one_matl = scinew MaterialSubset();
  one_matl->addReference();
  one_matl->add(0);

  Ghost::GhostType gac = Ghost::AroundCells;

  if (level->hasFinerLevel() == false) { // only on finest level
    t->needs(Task::NewDW, d_mpm_labels->gMassLabel, gac, 1);
    t->needs(Task::NewDW,
                d_mpm_labels->gTemperatureLabel,
                one_matl,
                gac,
                1);
    t->needs(Task::OldDW,
                d_mpm_labels->NC_CCweightLabel,
                one_matl,
                gac,
                1);
  }

  t->computes(d_switch_label);

  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());

  if (one_matl && one_matl->removeReference()) {
    delete one_matl;
  }
}

//  This task uses similar logic in the HEChem/simpleBurn.cc
//  to determine if the burning criteria has been reached.
void
SimpleBurnCriteria::switchTest([[maybe_unused]] const ProcessorGroup* group,
                               const PatchSubset* patches,
                               [[maybe_unused]] const MaterialSubset* matls,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  double timeToSwitch = 0;
  const Level* level  = getLevel(patches);

  if (level->hasFinerLevel() == false) { // only on finest level

    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      printTask(patches,
                patch,
                switching_dbg,
                "Doing SimpleBurnCriteria::switchTest");

      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        d_mat_manager->getMaterial("MPM", d_material));
      int indx = mpm_matl->getDWIndex();

      constNCVariable<double> gMass, gTempAllMatls;
      constNCVariable<double> NC_CCweight;
      Ghost::GhostType gac = Ghost::AroundCells;

      new_dw->get(gMass, d_mpm_labels->gMassLabel, indx, patch, gac, 1);
      new_dw
        ->get(gTempAllMatls, d_mpm_labels->gTemperatureLabel, 0, patch, gac, 1);
      old_dw
        ->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);

      IntVector nodeIdx[8];
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        patch->findNodesFromCell(*iter, nodeIdx);

        double Temp_CC_mpm = 0.0;
        double cmass       = 1.e-100;

        double MaxMass = d_SMALL_NUM;
        double MinMass = 1.0 / d_SMALL_NUM;

        for (int in = 0; in < 8; in++) {
          double NC_CCw_mass = NC_CCweight[nodeIdx[in]] * gMass[nodeIdx[in]];
          MaxMass            = std::max(MaxMass, NC_CCw_mass);
          MinMass            = std::min(MinMass, NC_CCw_mass);
          cmass += NC_CCw_mass;
          Temp_CC_mpm += gTempAllMatls[nodeIdx[in]] * NC_CCw_mass;
        }
        Temp_CC_mpm /= cmass;

        double ratio = (MaxMass - MinMass) / MaxMass;

        if (ratio > 0.4 && ratio < 1.0 && MaxMass > d_TINY_RHO) {
          if (Temp_CC_mpm >= d_temperature) {
            timeToSwitch = 1;
            std::cout
              << " \n"
              << " The simpleBurn switching criteria is satisfied in cell " << c
              << " (MaxMass-MinMass)/MaxMass = " << ratio
              << ", temp_CC_mpm= " << Temp_CC_mpm << "\n"
              << std::endl;
            break;
          }
        }
      } // iter
    }   // patches
  }     // finest Level

  max_vartype switch_condition(timeToSwitch);

  const Level* allLevels = 0;
  new_dw->put(switch_condition, d_switch_label, allLevels);
}

} // namespace Uintah