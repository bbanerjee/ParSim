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

/*
 *  compare_uda.cc: compare results of 2 uintah data archive
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   February 2000
 *
 */

#include <StandAlone/Utils/FieldComparator.h>
#include <StandAlone/Utils/MaterialParticleData.h>
#include <StandAlone/Utils/MaterialParticleVarData.h>
#include <StandAlone/Utils/compare_uda_options.h>
#include <StandAlone/Utils/compare_uda_utils.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/OS/Dir.h>
#include <Core/Util/ProgressiveWarning.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

using namespace Uintah;

using FieldComparator         = Vaango::Utils::CompareUda::FieldComparator;
using MaterialParticleData    = Vaango::Utils::CompareUda::MaterialParticleData;
using MaterialParticleDataMap = std::map<int, MaterialParticleData>;

void
addParticleData(MaterialParticleDataMap& matlParticleDataMap,
                DataArchive* da,
                std::vector<std::string> vars,
                std::vector<const Uintah::TypeDescription*> types,
                LevelP level,
                int timestep)
{
  Level::const_patch_iterator iter;
  for (iter = level->patchesBegin(); iter != level->patchesEnd(); iter++) {
    const Patch* patch = *iter;

    for (int v = 0; v < (int)vars.size(); v++) {
      std::string var = vars[v];

      const Uintah::TypeDescription* td      = types[v];
      const Uintah::TypeDescription* subtype = td->getSubType();

      if (td->getType() == Uintah::TypeDescription::Type::ParticleVariable) {
        ConsecutiveRangeSet matls = da->queryMaterials(var, patch, timestep);

        for (ConsecutiveRangeSet::iterator matlIter = matls.begin();
             matlIter != matls.end();
             matlIter++) {

          int matl = *matlIter;
          // Add a new MaterialPatchData for each matl for this next patch.
          MaterialParticleData& data = matlParticleDataMap[matl];
          data.setMatl(matl);
          ParticleVariableBase* pvb = nullptr;
          switch (subtype->getType()) {
            case Uintah::TypeDescription::Type::double_type:
              pvb = scinew ParticleVariable<double>();
              break;
            case Uintah::TypeDescription::Type::float_type:
              pvb = scinew ParticleVariable<float>();
              break;
            case Uintah::TypeDescription::Type::long64_type:
              pvb = scinew ParticleVariable<long64>();
              break;
            case Uintah::TypeDescription::Type::int_type:
              pvb = scinew ParticleVariable<int>();
              break;
            case Uintah::TypeDescription::Type::Point:
              pvb = scinew ParticleVariable<Point>();
              break;
            case Uintah::TypeDescription::Type::Vector:
              pvb = scinew ParticleVariable<Vector>();
              break;
            case Uintah::TypeDescription::Type::Matrix3:
              pvb = scinew ParticleVariable<Matrix3>();
              break;
            default:
              std::cerr
                << "addParticleData: ParticleVariable of unsupported type: "
                << subtype->getName() << '\n';
              Parallel::exitAll(-1);
          }
          da->query(*pvb, var, matl, patch, timestep);
          data[var].add(pvb, patch); // will add one for each patch
        }
      }
    }
  }
}

int
main(int argc, char** argv)
{
  Uintah::Parallel::initializeManager(argc, argv);

  Vaango::Utils::Options::compare_uda_parse(argc, argv);

  // default to 16 digits of precision when using exact comparison (i.e.
  // Vaango::Utils::Options::rel_tolerance() = 0)
  int digits_precision =
    (Vaango::Utils::Options::rel_tolerance() > 0)
      ? (int)std::ceil(-std::log10(Vaango::Utils::Options::rel_tolerance())) + 1
      : 16;
  std::cerr << std::setprecision(digits_precision);
  std::cout << std::setprecision(digits_precision);

  try {
    DataArchive* da1 = scinew DataArchive(Vaango::Utils::Options::filebase_1());
    DataArchive* da2 = scinew DataArchive(Vaango::Utils::Options::filebase_2());

    using VarTypeVec =
      std::vector<std::pair<std::string, const Uintah::TypeDescription*>>;

    std::vector<std::string> vars, vars2;
    std::vector<int> num_matls, num_matls2;
    std::vector<const Uintah::TypeDescription*> types, types2;
    VarTypeVec vartypes1, vartypes2;

    da1->queryVariables(vars, num_matls, types);
    ASSERTEQ(vars.size(), types.size());

    da2->queryVariables(vars2, num_matls2, types2);
    ASSERTEQ(vars2.size(), types2.size());

    vartypes1.resize(vars.size());
    vartypes2.resize(vars2.size());

    // Create a list of variables minus the ignored variables
    auto ignore_vars = Vaango::Utils::Options::ignore_vars();
    for (auto& var : ignore_vars) {
      std::cout << "Ignoring variable: " << var << std::endl;
    }
    // uda 1
    int count = 0;
    for (size_t i = 0; i < vars.size(); i++) {
      auto me = find(ignore_vars.begin(), ignore_vars.end(), vars[i]);
      if (me == ignore_vars.end()) {
        vartypes1[count] = std::make_pair(vars[i], types[i]);
        count++;
      }
    }
    vars.resize(count);
    vartypes1.resize(vars.size());
    // uda 2
    count = 0;
    for (size_t i = 0; i < vars2.size(); i++) {
      auto me = find(ignore_vars.begin(), ignore_vars.end(), vars2[i]);
      if (me == ignore_vars.end()) {
        vartypes2[count] = std::make_pair(vars2[i], types2[i]);
        count++;
      }
    }
    vars2.resize(count);
    vartypes2.resize(vars2.size());

    // Create a list of variables to compare if the user wants
    // to compare a few variables.  Default is to compare all
    auto compare_vars = Vaango::Utils::Options::compare_vars();
    for (auto& var : compare_vars) {
      std::cout << "Variable: " << var << std::endl;
    }
    if (compare_vars.size() > 0) {
      // uda 1
      count = 0;
      for (size_t i = 0; i < vars.size(); i++) {
        auto me = find(compare_vars.begin(), compare_vars.end(), vars[i]);
        if (me != compare_vars.end()) {
          vartypes1[count] = std::make_pair(vars[i], types[i]);
          count++;
        }
      }
      vars.resize(count);
      vartypes1.resize(vars.size());
      // uda 2
      count = 0;
      for (size_t i = 0; i < vars2.size(); i++) {
        auto me = find(compare_vars.begin(), compare_vars.end(), vars2[i]);
        if (me != compare_vars.end()) {
          vartypes2[count] = std::make_pair(vars2[i], types2[i]);
          count++;
        }
      }
      vars2.resize(count);
      vartypes2.resize(vars2.size());
    }

    size_t vars1_size = vars.size(); // needed for bullet proofing
    size_t vars2_size = vars2.size();

    // sort vars so uda's can be compared if their index files have
    // different orders of variables.
    // Assuming that there are no duplicates in the var names, these will
    // sort alphabetically by varname.
    if (Vaango::Utils::Options::sort_variables()) {
      std::sort(vartypes1.begin(), vartypes1.end());
      std::sort(vartypes2.begin(), vartypes2.end());
    }

    //  create vector of vars to compare
    bool do_udas_have_same_nVars = true;
    if (vartypes1.size() == vartypes2.size()) {
      for (size_t i = 0; i < vars.size(); i++) {
        vars[i]   = vartypes1[i].first;
        types[i]  = vartypes1[i].second;
        vars2[i]  = vartypes2[i].first;
        types2[i] = vartypes2[i].second;
      }
    } else {
      // If the number of variables in each uda differs then find a common set
      // of variables
      do_udas_have_same_nVars = false;
      std::cerr << "\nWARNING: The udas contain a different number of "
                   "variables.  Now comparing the common set of variables.\n";
      VarTypeVec commonVars; // common variables
      std::set_intersection(vartypes1.begin(),
                            vartypes1.end(),
                            vartypes2.begin(),
                            vartypes2.end(),
                            std::back_inserter(commonVars));
      size_t size = commonVars.size();
      vars.resize(size);
      vars2.resize(size);
      vartypes1.resize(size);
      vartypes2.resize(size);
      for (size_t i = 0; i < size; i++) {
        vars[i]   = commonVars[i].first;
        types[i]  = commonVars[i].second;
        vars2[i]  = commonVars[i].first;
        types2[i] = commonVars[i].second;
      }
    }

    for (size_t i = 0; i < vars.size(); i++) {
      if (vars[i] != vars2[i]) {
        std::ostringstream warn;
        warn << "    Variable " << vars[i] << " in "
             << Vaango::Utils::Options::filebase_1() << " does not match\n";
        warn << "    variable " << vars2[i] << " in "
             << Vaango::Utils::Options::filebase_2() << "\n";
        Vaango::Utils::Options::abort_uncomparable(warn);
      }

      if (types[i] != types2[i]) {
        std::ostringstream warn;
        warn << "    Variable " << vars[i]
             << " does not have the same type in both uda directories.\n";
        warn << "    In " << Vaango::Utils::Options::filebase_1()
             << " its type is " << types[i]->getName() << endl;
        warn << "    In " << Vaango::Utils::Options::filebase_2()
             << " its type is " << types2[i]->getName() << endl;
        Vaango::Utils::Options::abort_uncomparable(warn);
      }
    }

    std::vector<int> index;
    std::vector<double> times;
    std::vector<int> index2;
    std::vector<double> times2;

    da1->queryTimesteps(index, times);
    ASSERTEQ(index.size(), times.size());

    da2->queryTimesteps(index2, times2);
    ASSERTEQ(index2.size(), times2.size());

    for (unsigned long tstep = 0; tstep < times.size() && tstep < times2.size();
         tstep++) {
      if (!Vaango::Utils::CompareUda::compare(
            times[tstep],
            times2[tstep],
            Vaango::Utils::Options::abs_tolerance(),
            Vaango::Utils::Options::rel_tolerance())) {
        std::ostringstream warn;
        warn << "Timestep at time " << times[tstep] << " in "
             << Vaango::Utils::Options::filebase_1() << " does not match\n";
        warn << "timestep at time " << times2[tstep] << " in "
             << Vaango::Utils::Options::filebase_2()
             << " within the allowable tolerance.\n";
        Vaango::Utils::Options::abort_uncomparable(warn);
      }

      double time1 = times[tstep];
      double time2 = times2[tstep];
      std::cerr << "time = " << time1 << "\n";

      GridP grid  = da1->queryGrid(tstep);
      GridP grid2 = da2->queryGrid(tstep);

      int maxLevels[2];
      maxLevels[0] = grid->numLevels();
      maxLevels[1] = grid2->numLevels();

      int minLevel[2];
      minLevel[0] = 0;
      minLevel[1] = 0;

      // override if user has specified the levels to compare
      if (Vaango::Utils::Options::uda_levels(0) != -9) {
        minLevel[0]  = Vaango::Utils::Options::uda_levels(0);
        minLevel[1]  = Vaango::Utils::Options::uda_levels(1);
        maxLevels[0] = minLevel[0] + 1;
        maxLevels[1] = minLevel[1] + 1;

        if (maxLevels[0] > grid->numLevels() ||
            maxLevels[1] > grid2->numLevels()) {
          std::ostringstream warn;
          warn << "    The level(s) specified (uda:" << maxLevels[0]
               << " , uda2: " << maxLevels[1] << ") are invalid.\n";
          warn << "    The maximum level index that are valid are (uda:"
               << grid->numLevels() << " , uda2: " << grid2->numLevels()
               << ").\n";
          Vaango::Utils::Options::abort_uncomparable(warn);
        }

      } else if (maxLevels[0] != maxLevels[1]) {
        std::ostringstream warn;
        warn << "    Grid at time " << time1 << " in "
             << Vaango::Utils::Options::filebase_1() << " has "
             << grid->numLevels() << " levels.\n";
        warn << "    Grid at time " << time2 << " in "
             << Vaango::Utils::Options::filebase_2() << " has "
             << grid2->numLevels() << " levels.\n";
        Vaango::Utils::Options::abort_uncomparable(warn);
      }

      // do some consistency checking first
      bool hasParticleIDs  = false;
      bool hasParticleData = false;
      bool hasPerPatchData = false;

      for (size_t v = 0; v < vars.size(); v++) {
        std::string var = vars[v];

        if (var == "p.particleID") {
          hasParticleIDs = true;
        }
        if (types[v]->getType() ==
            Uintah::TypeDescription::Type::ParticleVariable) {
          hasParticleData = true;
        }
        if (types[v]->getType() == Uintah::TypeDescription::Type::PerPatch) {
          hasPerPatchData = true;
        }

        for (int l1 = minLevel[0], l2 = minLevel[1];
             (l1 < maxLevels[0] && l2 < maxLevels[1]);
             l1++, l2++) {
          LevelP level  = grid->getLevel(l1);
          LevelP level2 = grid2->getLevel(l2);

          ConsecutiveRangeSet matls;

          bool first = true;
          Level::const_patch_iterator iter;

          //  bulletproofing does the variable exist in both DAs on this
          //  timestep? This problem mainly occurs if <outputInitTimestep> has
          //  been specified.
          bool existsDA1 = true;
          bool existsDA2 = true;
          for (iter = level->patchesBegin(); iter != level->patchesEnd();
               iter++) {
            const Patch* patch = *iter;
            if (!da1->exists(var, patch, tstep)) {
              existsDA1 = false;
            }
          }
          for (iter = level2->patchesBegin(); iter != level2->patchesEnd();
               iter++) {
            const Patch* patch = *iter;
            if (!da2->exists(var, patch, tstep)) {
              existsDA2 = false;
            }
          }
          if (existsDA1 != existsDA2) {
            std::ostringstream warn;
            warn << "    The variable (" << var
                 << ") was not found on timestep (" << index[tstep]
                 << "), Level-" << level->getIndex() << ", in both udas ("
                 << existsDA1 << ", " << existsDA2 << ").\n"
                 << "    If this occurs on timestep 0 then (" << var
                 << ") was not computed in the initialization task.\n";
            Vaango::Utils::Options::abort_uncomparable(warn);
          }

          for (iter = level->patchesBegin(); iter != level->patchesEnd();
               iter++) {
            const Patch* patch = *iter;

            if (first) {
              matls = da1->queryMaterials(var, patch, tstep);
            } else if (matls != da1->queryMaterials(var, patch, tstep)) {
              std::ostringstream warn;
              warn << "    The material set is not consistent for variable "
                   << var << " across patches at time " << time1 << ".\n";
              warn << "    Previously was: " << matls << endl;
              warn << "    But on patch " << patch->getID() << ": "
                   << da1->queryMaterials(var, patch, tstep) << ".\n";
              Vaango::Utils::Options::abort_uncomparable(warn);
            }
            first = false;
          }

          ASSERT(!first); /* More serious problems would show up if this
                             assertion would fail */
          for (auto iter = level2->patchesBegin(); iter != level2->patchesEnd();
               iter++) {
            const Patch* patch = *iter;

            ConsecutiveRangeSet matls2 = da2->queryMaterials(var, patch, tstep);

            if (matls != matls2) {
              std::ostringstream warn;
              warn << "Inconsistent material sets for variable " << var
                   << " on patch = " << patch->getID() << ", time " << time1
                   << "\n";
              warn << "    " << Vaango::Utils::Options::filebase_1()
                   << " (1) has material set: " << matls << ".\n";
              warn << "    " << Vaango::Utils::Options::filebase_2()
                   << " (2) has material set: "
                   << da2->queryMaterials(var, patch, tstep) << ".\n";
              // If this is timestep 0 and one or both of the DAs are using PIDX
              // format:
              if (index[tstep] == 0 &&
                  (da1->isPIDXFormat() || da2->isPIDXFormat())) {
                std::cout
                  << "Ignoring the following warning because this is "
                     "timestep 0 with a PIDX UDA and thus probably not a "
                     "real issue...\n";
                std::cout << warn.str() << "\n";
              } else {
                Vaango::Utils::Options::abort_uncomparable(warn);
              }
            }
          }
        }
      }

      // Compare Particle and PerPatch Variables
      if (hasPerPatchData || (hasParticleData && !hasParticleIDs)) {

        // Compare particle variables without p.particleID -- patches
        // must be consistent.
        if (hasParticleData && !hasParticleIDs) {
          std::cerr << "Particle data exists without p.particleID output.\n";
          std::cerr << "There must be patch consistency in order to do this "
                       "comparison.\n";
          std::cerr
            << "In order to make a comparison between udas with different\n"
            << "number or distribution of patches, you must either output\n"
            << "p.particleID or don't output any particle variables at all.\n";
          std::cerr << std::endl;
        }

        for (size_t v = 0; v < vars.size(); v++) {
          std::string var = vars[v];

          const Uintah::TypeDescription* td      = types[v];
          const Uintah::TypeDescription* subtype = td->getSubType();

          if (td->getType() !=
                Uintah::TypeDescription::Type::ParticleVariable &&
              td->getType() != Uintah::TypeDescription::Type::PerPatch) {
            continue;
          }

          std::cerr << "\tVariable: " << std::left << std::setw(20) << var
                    << ", type " << td->getName() << "\n";

          if (td->getName() == std::string("-- unknown type --")) {
            std::cerr << "\t\tParticleVariable or PerPatch of unknown type";

            if (Vaango::Utils::Options::strict_types()) {
              std::cerr << ".\nQuitting.\n";
              Parallel::exitAll(-1);
            }
            std::cerr << "; skipping comparison...\n";
            continue;
          }

          for (int l1 = minLevel[0], l2 = minLevel[1];
               (l1 < maxLevels[0] && l2 < maxLevels[1]);
               l1++, l2++) {

            LevelP level  = grid->getLevel(l1);
            LevelP level2 = grid2->getLevel(l2);

            if (level->numPatches() != level2->numPatches()) {
              std::ostringstream warn;
              warn << "    Inconsistent number of patches on level " << l1
                   << " at time " << time1 << ":" << endl;
              warn << "    " << Vaango::Utils::Options::filebase_1() << " has "
                   << level->numPatches() << " patches.\n";
              warn << "    " << Vaango::Utils::Options::filebase_2() << " has "
                   << level2->numPatches() << " patches.\n";
              Vaango::Utils::Options::abort_uncomparable(warn);
            }

            auto iter2 = level2->patchesBegin();
            for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
                 iter++, iter2++) {

              const Patch* patch  = *iter;
              const Patch* patch2 = *iter2;

              if (patch->getID() != patch2->getID()) {
                std::ostringstream warn;
                warn << "    Inconsistent patch ids on level " << l1
                     << " at time " << time1 << endl;
                warn << "    " << Vaango::Utils::Options::filebase_1()
                     << " has patch id " << patch->getID() << " where\n";
                warn << "    " << Vaango::Utils::Options::filebase_2()
                     << " has patch id " << patch2->getID() << ".\n";
                Vaango::Utils::Options::abort_uncomparable(warn);
              }

              // cerr << "\t\tPatch: " << patch->getID() << "\n";

              if (!Vaango::Utils::CompareUda::compare(
                    patch->getExtraBox().lower(),
                    patch2->getExtraBox().lower(),
                    Vaango::Utils::Options::abs_tolerance(),
                    Vaango::Utils::Options::rel_tolerance()) ||
                  !Vaango::Utils::CompareUda::compare(
                    patch->getExtraBox().upper(),
                    patch2->getExtraBox().upper(),
                    Vaango::Utils::Options::abs_tolerance(),
                    Vaango::Utils::Options::rel_tolerance())) {

                std::ostringstream warn;
                warn << "    Inconsistent patch bounds on patch "
                     << patch->getID() << " at time " << time1 << ".\n";

                warn << "    " << Vaango::Utils::Options::filebase_1()
                     << " has bounds " << patch->getExtraBox().lower() << " - "
                     << patch->getExtraBox().upper() << ".\n";

                warn << "    " << Vaango::Utils::Options::filebase_2()
                     << " has bounds " << patch2->getExtraBox().lower() << " - "
                     << patch2->getExtraBox().upper() << ".\n";

                warn << "    "
                     << "Difference is: "
                     << patch->getExtraBox().lower() -
                          patch2->getExtraBox().lower()
                     << " - "
                     << patch->getExtraBox().upper() -
                          patch2->getExtraBox().upper()
                     << ".\n";
                Vaango::Utils::Options::abort_uncomparable(warn);
              }

              ConsecutiveRangeSet matls =
                da1->queryMaterials(var, patch, tstep);
              ConsecutiveRangeSet matls2 =
                da2->queryMaterials(var, patch2, tstep);
              ASSERT(matls == matls2); // should have already been checked

              // loop over materials
              for (ConsecutiveRangeSet::iterator matlIter = matls.begin();
                   matlIter != matls.end();
                   matlIter++) {
                int matl = *matlIter;

                // cerr << "\t\t\tMaterial: " << matl << "\n";

                if (td->getType() ==
                    Uintah::TypeDescription::Type::ParticleVariable) {
                  switch (subtype->getType()) {
                    case Uintah::TypeDescription::Type::double_type:
                      Vaango::Utils::CompareUda::compareParticles<double>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::float_type:
                      Vaango::Utils::CompareUda::compareParticles<float>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::int_type:
                      Vaango::Utils::CompareUda::compareParticles<int>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Point:
                      Vaango::Utils::CompareUda::compareParticles<Point>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Vector:
                      Vaango::Utils::CompareUda::compareParticles<Vector>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::IntVector:
                      Vaango::Utils::CompareUda::compareParticles<IntVector>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Matrix3:
                      Vaango::Utils::CompareUda::compareParticles<Matrix3>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    default:
                      std::cerr
                        << "main: ParticleVariable of unsupported type: "
                        << subtype->getName() << '\n';
                      Parallel::exitAll(-1);
                  }
                } else if (td->getType() ==
                           Uintah::TypeDescription::Type::PerPatch) {

                  switch (subtype->getType()) {
                    case Uintah::TypeDescription::Type::double_type:
                      Vaango::Utils::CompareUda::comparePerPatch<double>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::float_type:
                      Vaango::Utils::CompareUda::comparePerPatch<float>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::int_type:
                      Vaango::Utils::CompareUda::comparePerPatch<int>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Point:
                      Vaango::Utils::CompareUda::comparePerPatch<Point>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Vector:
                      Vaango::Utils::CompareUda::comparePerPatch<Vector>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::IntVector:
                      Vaango::Utils::CompareUda::comparePerPatch<IntVector>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    case Uintah::TypeDescription::Type::Matrix3:
                      Vaango::Utils::CompareUda::comparePerPatch<Matrix3>(
                        da1,
                        da2,
                        var,
                        matl,
                        patch,
                        patch2,
                        time1,
                        tstep,
                        Vaango::Utils::Options::abs_tolerance(),
                        Vaango::Utils::Options::rel_tolerance());
                      break;
                    default:
                      std::cerr << "main: PerPatch of unsupported type: "
                                << subtype->getName() << '\n';
                      Parallel::exitAll(-1);
                  }
                }
              }
            }
          }
        }
      } else if (hasParticleIDs) {
        // Compare Particle variables with p.particleID -- patches don't
        // need to be cosistent.  It will gather and sort the particles
        // so they can be compared in particleID order.

        for (int l1 = minLevel[0], l2 = minLevel[1];
             (l1 < maxLevels[0] && l2 < maxLevels[1]);
             l1++, l2++) {

          LevelP level  = grid->getLevel(l1);
          LevelP level2 = grid2->getLevel(l2);

          MaterialParticleDataMap matlParticleDataMap1;
          MaterialParticleDataMap matlParticleDataMap2;

          addParticleData(matlParticleDataMap1, da1, vars, types, level, tstep);
          addParticleData(
            matlParticleDataMap2, da2, vars2, types2, level2, tstep);

          MaterialParticleDataMap::iterator matlIter;
          MaterialParticleDataMap::iterator matlIter2;

          matlIter  = matlParticleDataMap1.begin();
          matlIter2 = matlParticleDataMap2.begin();

          for (; (matlIter != matlParticleDataMap1.end()) &&
                 (matlIter2 != matlParticleDataMap2.end());
               matlIter++, matlIter2++) {

            // This assert should already have been checked above whan comparing
            // material sets.
            ASSERT((*matlIter).first == (*matlIter).first);

            (*matlIter).second.compare((*matlIter2).second,
                                       time1,
                                       time2,
                                       Vaango::Utils::Options::abs_tolerance(),
                                       Vaango::Utils::Options::rel_tolerance());
          }
          // This assert should already have been check above whan comparing
          // material sets.
          ASSERT(matlIter == matlParticleDataMap1.end() &&
                 matlIter2 == matlParticleDataMap2.end());
        }
      }

      for (size_t v = 0; v < vars.size(); v++) {
        std::string var = vars[v];

        const Uintah::TypeDescription* td      = types[v];
        const Uintah::TypeDescription* subtype = td->getSubType();

        if (td->getType() == Uintah::TypeDescription::Type::ParticleVariable ||
            td->getType() == Uintah::TypeDescription::Type::PerPatch) {
          continue;
        }

        std::cerr << "\tVariable: " << var << ", type " << td->getName()
                  << "\n";

        if (td->getName() == string("-- unknown type --")) {
          std::cerr << "\t\tParticleVariable of unknown type";
          if (Vaango::Utils::Options::strict_types()) {
            std::cerr << ".\nQuitting.\n";
            Parallel::exitAll(-1);
          }
          std::cerr << "; skipping comparison...\n";
          continue;
        }

        Patch::VariableBasis basis =
          Patch::translateTypeToBasis(td->getType(), false);

        for (int l1 = minLevel[0], l2 = minLevel[1];
             (l1 < maxLevels[0] && l2 < maxLevels[1]);
             l1++, l2++) {

          LevelP level  = grid->getLevel(l1);
          LevelP level2 = grid2->getLevel(l2);

          // check patch coverage
          std::vector<Region> region1, region2, difference1, difference2;

          for (int i = 0; i < level->numPatches(); i++) {
            const Patch* patch = level->getPatch(i);
            region1.push_back(Region(patch->getExtraCellLowIndex(),
                                     patch->getExtraCellHighIndex()));
          }
          for (int i = 0; i < level2->numPatches(); i++) {
            const Patch* patch = level2->getPatch(i);
            region2.push_back(Region(patch->getExtraCellLowIndex(),
                                     patch->getExtraCellHighIndex()));
          }

          difference1 = Region::difference(region1, region2);
          difference2 = Region::difference(region1, region2);

          if (!difference1.empty() || !difference2.empty()) {
            std::ostringstream warn;
            warn << "    Patches on level:" << l1
                 << " do not cover the same area.\n";
            Vaango::Utils::Options::abort_uncomparable(warn);
          }

          // map nodes to patches in level and level2 respectively
          Array3<const Patch*> patchMap;
          Array3<const Patch*> patch2Map;

          Vaango::Utils::CompareUda::buildPatchMap(
            level,
            Vaango::Utils::Options::filebase_1(),
            patchMap,
            time1,
            basis);
          Vaango::Utils::CompareUda::buildPatchMap(
            level2,
            Vaango::Utils::Options::filebase_2(),
            patch2Map,
            time2,
            basis);

          serial_for(patchMap.range(), [&](int i, int j, int k) {
            // bulletproofing
            if ((patchMap(i, j, k) == nullptr &&
                 patch2Map(i, j, k) != nullptr) ||
                (patch2Map(i, j, k) == nullptr &&
                 patchMap(i, j, k) != nullptr)) {
              std::ostringstream warn;
              warn << "    Inconsistent patch coverage on level " << l1
                   << " at time " << time1 << ".\n";

              if (patchMap(i, j, k) != nullptr) {
                warn << "    " << IntVector(i, j, k) << " is covered by "
                     << Vaango::Utils::Options::filebase_1() << std::endl
                     << " and not " << Vaango::Utils::Options::filebase_2()
                     << ".\n";
              } else {
                warn << "    " << IntVector(i, j, k) << " is covered by "
                     << Vaango::Utils::Options::filebase_2() << std::endl
                     << " and not " << Vaango::Utils::Options::filebase_1()
                     << ".\n";
              }
              Vaango::Utils::Options::abort_uncomparable(warn);
            }
          });

          Level::const_patch_iterator iter;
          for (iter = level->patchesBegin(); iter != level->patchesEnd();
               iter++) {
            const Patch* patch = *iter;

            ConsecutiveRangeSet matls = da1->queryMaterials(var, patch, tstep);

            FieldComparator* comparator = FieldComparator::makeFieldComparator(
              td,
              subtype,
              patch,
              Vaango::Utils::Options::include_extra_cells());

            if (comparator != 0) {
              comparator->compareFields(
                da1,
                da2,
                var,
                matls,
                patch,
                patch2Map,
                time1,
                tstep,
                Vaango::Utils::Options::abs_tolerance(),
                Vaango::Utils::Options::rel_tolerance());
              delete comparator;
            }
          }
        } // end for (l)
      }   // end for (v)
    }     // end for(tstep)

    if (times.size() != times2.size()) {
      std::ostringstream warn;
      warn << "\n";
      warn << "    " << Vaango::Utils::Options::filebase_1() << " has "
           << times.size() << " timesteps\n";
      warn << "    " << Vaango::Utils::Options::filebase_2() << " has "
           << times2.size() << " timesteps\n";
      Vaango::Utils::Options::abort_uncomparable(warn);
    }

    if (!do_udas_have_same_nVars) {
      std::ostringstream warn;
      warn << "    " << Vaango::Utils::Options::filebase_1() << " has "
           << vars1_size << " variables\n";
      warn << "    " << Vaango::Utils::Options::filebase_2() << " has "
           << vars2_size << " variables\n";
      Vaango::Utils::Options::abort_uncomparable(warn);
    }

    delete da1;
    delete da2;
  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << '\n';
    abort();
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    abort();
  }

  if (Vaango::Utils::Options::tolerance_error()) {
    std::cerr << "\nComparison did NOT fully pass.\n";
    Parallel::exitAll(2);
  } else {
    std::cerr << "\nComparison fully passed!\n";
  }

  return 0;
}
