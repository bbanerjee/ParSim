/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#include <CCA/Components/SwitchingCriteria/DDT1Criterion.h>

#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

namespace Uintah {

extern DebugStream switching_dbg;

DDT1Criterion::DDT1Criterion(ProblemSpecP& ps)
{
  ps->require("reactant_material", d_material);
  ps->require("ThresholdTemperature", d_temperature);
  ps->require("BoundaryParticles", d_BP);

  proc0cout << "Switching criteria:  \tDDT1, reactant matl: " << d_material
            << " Threshold tempterature " << d_temperature
            << ", Boundary Particles " << d_BP << std::endl;

  d_mpm_labels     = std::make_unique<MPMLabel>();
  d_ice_labels     = std::make_unique<ICELabel>();
  d_mpm_ice_labels = std::make_unique<MPMICELabel>();
}

void
DDT1Criterion::problemSetup([[maybe_unused]] const ProblemSpecP& ps,
                            [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                            MaterialManagerP& materialManager)
{
  d_materialManager = materialManager;
}

void
DDT1Criterion::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level,
                switching_dbg,
                "Switching Criteria:DDT1Criterion::scheduleSwitchTest");
  Task* t = scinew Task("switchTest", this, &DDT1Criterion::switchTest);

  auto* one_matl = scinew MaterialSubset();
  one_matl->addReference();
  one_matl->add(0);

  const MaterialSubset* mpm_matls =
    d_materialManager->allMaterials("MPM")->getUnion();

  Ghost::GhostType gac = Ghost::AroundCells;
  Ghost::GhostType gan = Ghost::AroundNodes;

  if (level->hasFinerLevel() == false) { // only on the finest level
    t->requires(Task::OldDW, d_ice_labels->vol_frac_CCLabel, mpm_matls, gac, 1);
    t->requires(Task::NewDW, d_mpm_labels->gMassLabel, mpm_matls, gan, 2);
    t->requires(Task::OldDW, d_mpm_labels->pXLabel, mpm_matls, gan, 1);
    t->requires(Task::NewDW,
                d_mpm_labels->gTemperatureLabel,
                one_matl,
                gan,
                2);
    t->requires(Task::OldDW,
                d_mpm_labels->NC_CCweightLabel,
                one_matl,
                gan,
                2);
  }

  t->computes(d_switch_label);

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials());

  if (one_matl && one_matl->removeReference()) {
    delete one_matl;
  }
}

//  This task uses similar logic in the HEChem/DDT1.cc
//  to determine if the burning criteria has been reached.
void
DDT1Criterion::switchTest([[maybe_unused]] const ProcessorGroup* group,
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
                "Doing Switching Criteria:DDT1Criterion::switchTest");

      MPMMaterial* mpm_matl =
        (MPMMaterial*)d_materialManager->getMaterial("MPM", d_material);
      int d_indx = mpm_matl->getDWIndex();

      int numAllMatls = d_materialManager->getNumMaterials();
      int numMPMMatls = d_materialManager->getNumMaterials("MPM");

      // mpm matls
      constNCVariable<double> NC_CCweight;
      constNCVariable<double> gTempAllMatls;
      std::vector<constNCVariable<double>> gMass(numMPMMatls);
      std::vector<CCVariable<double>> temp_CC_mpm(numAllMatls);
      std::vector<constCCVariable<double>> vol_frac_mpm(numAllMatls);

      Ghost::GhostType gac = Ghost::AroundCells;
      Ghost::GhostType gan = Ghost::AroundNodes;
      for (int m = 0; m < numMPMMatls; m++) {
        Material* matl = d_materialManager->getMaterial(m);
        int indx       = matl->getDWIndex();
        new_dw->get(gMass[m], d_mpm_labels->gMassLabel, indx, patch, gan, 2);
        old_dw->get(vol_frac_mpm[m],
                    d_ice_labels->vol_frac_CCLabel,
                    indx,
                    patch,
                    gac,
                    1);
        new_dw->allocateTemporary(temp_CC_mpm[m], patch, gac, 1);
        temp_CC_mpm[m].initialize(0.0);
      }
      new_dw
        ->get(gTempAllMatls, d_mpm_labels->gTemperatureLabel, 0, patch, gan, 2);
      old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gan, 2);

      constParticleVariable<Point> px;
      ParticleSubset* pset =
        old_dw->getParticleSubset(d_indx, patch, gan, 1, d_mpm_labels->pXLabel);
      old_dw->get(px, d_mpm_labels->pXLabel, pset);

      // Which cells contain particles
      CCVariable<double> pFlag;
      new_dw->allocateTemporary(pFlag, patch, gac, 1);
      pFlag.initialize(0.0);
      IntVector nodeIdx[8];

      // count how many reactant particles are in each cell
      for (const auto& idx : *pset) {
        IntVector c;
        patch->findCell(px[idx], c);
        pFlag[c] += 1.0;
      }

      // compute temp_CC_mpm in cells that contain some mass
      // The computational domain needs to hit the ghost cells
      int NGC = 1; // number of ghostCells
      for (CellIterator iter = patch->getCellIterator(NGC); !iter.done();
           iter++) {
        IntVector c = *iter;
        patch->findNodesFromCell(*iter, nodeIdx);

        for (int m = 0; m < numMPMMatls; m++) {
          double Temp_CC = 0.0;
          double cmass   = 1.e-100;
          for (int in = 0; in < 8; in++) {
            double NC_CCw_mass =
              NC_CCweight[nodeIdx[in]] * gMass[m][nodeIdx[in]];
            cmass += NC_CCw_mass;
            Temp_CC += gTempAllMatls[nodeIdx[in]] * NC_CCw_mass;
          }
          if (cmass > 1e-100) {
            Temp_CC /= cmass;
            temp_CC_mpm[m][c] = Temp_CC;
          }
        }
      } // cell iterator

      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        patch->findNodesFromCell(*iter, nodeIdx);

        double MaxMass = d_SMALL_NUM;
        double MinMass = 1.0 / d_SMALL_NUM;
        for (int in = 0; in < 8; in++) {
          double NC_CCw_mass =
            NC_CCweight[nodeIdx[in]] * gMass[d_material][nodeIdx[in]];
          MaxMass = std::max(MaxMass, NC_CCw_mass);
          MinMass = std::min(MinMass, NC_CCw_mass);
        }

        if (MinMass / MaxMass < 0.7 && pFlag[c] > 0) {
          // near interface and containing particles
          for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
              for (int k = -1; k <= 1; k++) {
                IntVector cell = c + IntVector(i, j, k);
                if (pFlag[cell] <= d_BP) {
                  for (int m = 0; m < numMPMMatls; m++) {

                    if (vol_frac_mpm[m][cell] > 0.2 &&
                        temp_CC_mpm[m][cell] > d_temperature) {
                      // cout << " The switching criteria satisfied in cell
                      // "<<cell
                      //      << " vol_frac_mpm " << vol_frac_mpm[m][cell]
                      //      << " temp_CC_mpm " << temp_CC_mpm[m][cell]
                      //      << " matl " << m
                      //      << " main cell " << c << std::endl;
                      timeToSwitch = 1;
                      break;
                    }
                  }
                } // endif
              }   // k
            }     // j
          }       // i
        }         // near a surface
      }           // cell iterator
    }             // patches
  }               // on finest level

  // compute on every level
  max_vartype switch_condition(timeToSwitch);

  const Level* allLevels = 0;
  new_dw->put(switch_condition, d_switch_label, allLevels);
}

} // namespace Uintah