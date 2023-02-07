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

#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/CohesiveZone/CohesiveZone.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/CZLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Patch.h>

#include <fstream>
#include <iostream>

using namespace Uintah;

//______________________________________________________________________
//  Reference: N. P. Daphalapukar, Hongbing Lu, Demir Coker, Ranga Komanduri,
// " Simulation of dynamic crack growth using the generalized interpolation
// material point (GIMP) method," Int J. Fract, 2007, 143:79-102
//______________________________________________________________________
CohesiveZone::CohesiveZone(CZMaterial* czmat,
                           const MaterialManagerP& mat_manager,
                           const MPMLabel* labels,
                           const MPMFlags* flags)
  : d_mat_manager{ mat_manager }
  , d_mpm_labels{ labels }
  , d_mpm_flags{ flags }
{
  d_cz_labels = std::make_unique<CZLabel>();
  registerPermanentCohesiveZoneState(czmat);
}

//______________________________________________________________________
ParticleSubset*
CohesiveZone::createCohesiveZones(CZMaterial* matl,
                                  particleIndex numCohesiveZones,
                                  CCVariable<short int>& cellNAPID,
                                  const Patch* patch,
                                  DataWarehouse* new_dw,
                                  const string filename)
{
  int dwi = matl->getDWIndex();
  ParticleSubset* subset =
    allocateVariables(numCohesiveZones, dwi, patch, new_dw);

  particleIndex start = 0;

  if (filename != "") {
    std::ifstream is(filename.c_str());
    if (!is) {
      throw ProblemSetupException("ERROR Opening cohesive zone file " +
                                    filename + " in createCohesiveZones \n",
                                  __FILE__,
                                  __LINE__);
    }

    // needed for bulletproofing
    std::vector<int> mpmMatlIndex;
    int numMPM = d_mat_manager->getNumMaterials("MPM");

    for (int m = 0; m < numMPM; m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      mpmMatlIndex.push_back(dwi);
    }

    // Field for position, normal, tangential and length.
    // Everything else is assumed to be zero.

    double p1, p2, p3, l4, n5, n6, n7, t8, t9, t10;
    int mt, mb;
    while (is >> p1 >> p2 >> p3 >> l4 >> n5 >> n6 >> n7 >> t8 >> t9 >> t10 >>
           mb >> mt) {

      //__________________________________
      // bulletproofing
      //  the top
      int test1 = count(mpmMatlIndex.begin(), mpmMatlIndex.end(), mb);
      int test2 = count(mpmMatlIndex.begin(), mpmMatlIndex.end(), mt);

      if (test1 == 0 || test2 == 0) {
        std::ostringstream warn;
        warn << "ERROR:MPM:createCohesiveZones\n In the cohesive zone file (" +
                  filename + ") either the top/bottom material";
        warn << "(top: " << mt << " bottom: " << mb
             << ") is not a MPM material ";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }

      Point pos = Point(p1, p2, p3);
      IntVector cell_idx;
      if (patch->containsPoint(pos)) {
        particleIndex pidx  = start;
        pCZPosition[pidx]   = pos;
        pCZArea[pidx]       = l4;
        pCZNormal[pidx]     = Vector(n5, n6, n7);
        pCZTangent[pidx]    = Vector(t8, t9, t10);
        pCZDispTop[pidx]    = Vector(0.0, 0.0, 0.0);
        pCZDispBottom[pidx] = Vector(0.0, 0.0, 0.0);
        pCZSeparation[pidx] = Vector(0.0, 0.0, 0.0);
        pCZForce[pidx]      = Vector(0.0, 0.0, 0.0);
        pCZBotMat[pidx]     = mb;
        pCZTopMat[pidx]     = mt;
        pCZFailed[pidx]     = 0;

        // Figure out unique ID for the CZ
        patch->findCell(pos, cell_idx);
        ASSERT(cell_idx.x() <= 0xffff && cell_idx.y() <= 0xffff &&
               cell_idx.z() <= 0xffff);

        long64 cellID = ((long64)cell_idx.x() << 16) |
                        ((long64)cell_idx.y() << 32) |
                        ((long64)cell_idx.z() << 48);

        short int& myCellNAPID = cellNAPID[cell_idx];
        pCZID[pidx]            = (cellID | (long64)myCellNAPID);
        ASSERT(myCellNAPID < 0x7fff);
        myCellNAPID++;
        start++;
      }
    } // while
    is.close();
  }

  return subset;
}

//__________________________________
//
ParticleSubset*
CohesiveZone::allocateVariables(particleIndex numCZs,
                                int dwi,
                                const Patch* patch,
                                DataWarehouse* new_dw)
{

  ParticleSubset* subset = new_dw->createParticleSubset(numCZs, dwi, patch);

  new_dw->allocateAndPut(pCZPosition, d_mpm_labels->pXLabel, subset);
  new_dw->allocateAndPut(pCZArea, d_cz_labels->czAreaLabel, subset);
  new_dw->allocateAndPut(pCZNormal, d_cz_labels->czNormLabel, subset);
  new_dw->allocateAndPut(pCZTangent, d_cz_labels->czTangLabel, subset);
  new_dw->allocateAndPut(pCZDispTop, d_cz_labels->czDispTopLabel, subset);
  new_dw->allocateAndPut(pCZDispBottom, d_cz_labels->czDispBottomLabel, subset);
  new_dw->allocateAndPut(pCZID, d_cz_labels->czIDLabel, subset);
  new_dw->allocateAndPut(pCZSeparation, d_cz_labels->czSeparationLabel, subset);
  new_dw->allocateAndPut(pCZForce, d_cz_labels->czForceLabel, subset);
  new_dw->allocateAndPut(pCZTopMat, d_cz_labels->czTopMatLabel, subset);
  new_dw->allocateAndPut(pCZBotMat, d_cz_labels->czBotMatLabel, subset);
  new_dw->allocateAndPut(pCZFailed, d_cz_labels->czFailedLabel, subset);

  return subset;
}

//__________________________________
//
particleIndex
CohesiveZone::countCohesiveZones(const Patch* patch, const string filename)
{
  particleIndex sum = 0;

  if (filename != "") {
    std::ifstream is(filename.c_str());
    if (!is) {
      throw ProblemSetupException("ERROR Opening cohesive zone file " +
                                    filename + " in countCohesiveZones\n",
                                  __FILE__,
                                  __LINE__);
    }

    // Field for position, normal, tangential and length.
    // Everything else is assumed to be zero.
    double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, mt, mb;
    while (is >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10 >>
           mt >> mb) {
      // cout << f1 << " " << f2 << " " << f3 << std::endl;
      if (patch->containsPoint(Point(f1, f2, f3))) {
        sum++;
      }
    }
    is.close();
  }

  return sum;
}
//__________________________________
//
vector<const VarLabel*>
CohesiveZone::returnCohesiveZoneState()
{
  return d_cz_state;
}
//__________________________________
//
vector<const VarLabel*>
CohesiveZone::returnCohesiveZoneStatePreReloc()
{
  return d_cz_state_preReloc;
}
//__________________________________
//
void
CohesiveZone::registerPermanentCohesiveZoneState(CZMaterial* czmat)
{
  d_cz_state.push_back(d_cz_labels->czAreaLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czAreaLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czNormLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czNormLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czTangLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czTangLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czDispTopLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czDispTopLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czDispBottomLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czDispBottomLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czSeparationLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czSeparationLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czForceLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czForceLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czTopMatLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czTopMatLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czBotMatLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czBotMatLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czFailedLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czFailedLabel_preReloc);

  d_cz_state.push_back(d_cz_labels->czIDLabel);
  d_cz_state_preReloc.push_back(d_cz_labels->czIDLabel_preReloc);
}
//__________________________________
//
void
CohesiveZone::scheduleInitialize(const LevelP& level,
                                 SchedulerP& sched,
                                 CZMaterial* czmat)
{
  Task* t =
    scinew Task("CohesiveZone::initialize", this, &CohesiveZone::initialize);

  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_cz_labels->czAreaLabel);
  t->computes(d_cz_labels->czNormLabel);
  t->computes(d_cz_labels->czTangLabel);
  t->computes(d_cz_labels->czDispTopLabel);
  t->computes(d_cz_labels->czDispBottomLabel);
  t->computes(d_cz_labels->czSeparationLabel);
  t->computes(d_cz_labels->czForceLabel);
  t->computes(d_cz_labels->czTopMatLabel);
  t->computes(d_cz_labels->czBotMatLabel);
  t->computes(d_cz_labels->czFailedLabel);
  t->computes(d_cz_labels->czIDLabel);
  t->computes(d_cz_labels->pCellNACZIDLabel, zeroth_matl);

  std::vector<int> m(1);
  m[0]                     = czmat->getDWIndex();
  MaterialSet* cz_matl_set = scinew MaterialSet();
  cz_matl_set->addAll(m);
  cz_matl_set->addReference();

  sched->addTask(t, level->eachPatch(), cz_matl_set);

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) {
    delete zeroth_matl; // shouln't happen, but...
  }
}

//__________________________________
//
void
CohesiveZone::initialize(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* cz_matls,
                         DataWarehouse*,
                         DataWarehouse* new_dw)
{
  particleIndex totalCZs = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    //  printTask(patches, patch,cout_doing,"Doing initialize for
    //  CohesiveZones\t");

    CCVariable<short int> cellNACZID;
    new_dw->allocateAndPut(cellNACZID, d_cz_labels->pCellNACZIDLabel, 0, patch);
    cellNACZID.initialize(0);

    for (int m = 0; m < cz_matls->size(); m++) {
      CZMaterial* cz_matl =
        static_cast<CZMaterial*>(d_mat_manager->getMaterial("CohesiveZone", m));
      std::string filename = cz_matl->getCohesiveFilename();
      particleIndex numCZs = countCohesiveZones(patch, filename);
      totalCZs += numCZs;

      std::cout << "Total CZs " << totalCZs << std::endl;

      createCohesiveZones(cz_matl, numCZs, cellNACZID, patch, new_dw, filename);
    }
  }
}
