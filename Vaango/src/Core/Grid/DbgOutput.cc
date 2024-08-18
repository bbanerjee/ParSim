/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DOUT.hpp>

namespace Uintah {

void
printSchedule(const PatchSet* patches,
              Uintah::DebugStream& dbg,
              const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);
    dbg << where << "L-" << getLevel(patches)->getIndex() << std::endl;
  }
}

void
printSchedule(const LevelP& level,
              Uintah::DebugStream& dbg,
              const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);
    dbg << where << "L-" << level->getIndex() << std::endl;
  }
}

void
printTask(const PatchSubset* patches,
          const Patch* patch,
          Uintah::DebugStream& dbg,
          const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);
    if (patches->empty()) {
      dbg << "  \tEmpty patch subset\n";
    } else {
      dbg << where << "  \tL-" << getLevel(patches)->getIndex() << " patch "
          << patch->getGridIndex() << std::endl;
    }
  }
}

void
printTask(const PatchSubset* patches,
          Uintah::DebugStream& dbg,
          const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);

    if (patches->empty()) {
      dbg << "  \tEmpty patch subset\n";
    } else {
      dbg << where << "  \tL-" << getLevel(patches)->getIndex() << " patches "
          << *patches << std::endl;
    }
  }
}

void
printTask(const Patch* patch, Uintah::DebugStream& dbg, const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);
    dbg << where << " \tL-" << patch->getLevel()->getIndex() << " patch "
        << patch->getGridIndex() << std::endl;
  }
}

void
printTask(Uintah::DebugStream& dbg, const string& where) {
  if (dbg.active()) {
    dbg << Uintah::Parallel::getMPIRank() << " ";
    dbg << std::left;
    dbg.width(50);
    dbg << where << " \tAll Patches" << std::endl;
  }
}

//______________________________________________________________________
// Output the task name and the level it's executing on and each of the patches
void
printTask(Dout& out, DetailedTask* dtask) {
  if (out) {
    std::ostringstream msg;
    msg << std::left;
    msg.width(70);
    msg << dtask->getTask()->getName();
    if (dtask->getPatches()) {
      msg << " \t on patches ";
      const PatchSubset* patches = dtask->getPatches();
      for (int p = 0; p < patches->size(); p++) {
        if (p != 0) {
          msg << ", ";
        }
        msg << patches->get(p)->getID();
      }

      if (dtask->getTask()->getType() != Task::OncePerProc &&
          dtask->getTask()->getType() != Task::Hypre) {
        const Level* level = getLevel(patches);
        msg << "\t  L-" << level->getIndex();
      }
    }
    DOUT(out, msg.str());
  }
}

void
printTask(const Patch* patch, Dout& out, const std::string& where) {
  if (out) {
    std::ostringstream msg;
    msg << "Rank-" << Uintah::Parallel::getMPIRank() << " ";
    msg << std::left;
    msg.width(50);
    msg << where << " \tL-" << patch->getLevel()->getIndex() << " patch "
        << patch->getGridIndex();
    DOUT(out, msg.str());
  }
}

void
printTask(const PatchSubset* patches,
          const Patch* patch,
          Dout& out,
          const std::string& where) {
  if (out) {
    std::ostringstream msg;
    msg << "Rank-" << Uintah::Parallel::getMPIRank() << " ";
    msg << std::left;
    msg.width(50);
    msg << where;

    if (patches == nullptr || patches->empty()) {
      msg << "  \tEmpty patch subset";
    } else {
      msg << "  \tL-" << getLevel(patches)->getIndex() << " patch "
          << patch->getGridIndex();
    }

    DOUT(out, msg.str());
  }
}
void
printTask(const PatchSubset* patches, Dout& out, const std::string& where) {
  if (out) {
    std::ostringstream msg;
    msg << "Rank-" << Uintah::Parallel::getMPIRank() << " ";
    msg << std::left;
    msg.width(50);
    msg << where;

    if (patches == nullptr || patches->empty()) {
      msg << "  \tEmpty patch subset";
    } else {
      msg << "  \tL-" << getLevel(patches)->getIndex() << " patches "
          << *patches;
    }

    DOUT(out, msg.str());
  }
}

/*  Output the rank, task name, first patch on the level and thelevel index
  Rank-0   ICE::computeThermoTransportProperties                 Patch-0 L-0
  Rank-0   ICE::computeEquilibrationPressure                     Patch-0 L-0
  Rank-0   ICE::computeVel_FC                                    Patch-0 L-0
  */
void
printTaskLevels(const ProcessorGroup* d_myworld,
                Dout& out,
                DetailedTask* dtask) {
  if (!out) {
    return;
  }

  const PatchSubset* taskPatches = dtask->getPatches();
  if (taskPatches) {
    // Are all patches on the same level?
    bool uniqueLevel = true;
    int L_indx_p0    = taskPatches->get(0)->getLevelP()->getIndex();

    for (int i = 1; i < taskPatches->size(); ++i) {
      int L_indx = taskPatches->get(i)->getLevelP()->getIndex();

      if (L_indx != L_indx_p0) {
        uniqueLevel = false;
        break;
      }
    }

    //  Adjust message if patches are on multiple levels
    if (uniqueLevel) {
      const Level* level      = getLevel(taskPatches);
      const Patch* firstPatch = level->getPatch(0);

      if (taskPatches->contains(firstPatch)) {
        std::ostringstream msg;
        msg << "Rank-" << d_myworld->myRank() << "   ";
        msg << std::left;
        msg.width(50);
        msg << dtask->getTask()->getName();
        msg << "\t Patch-" << firstPatch->getGridIndex();
        msg << "\t L-" << level->getIndex();
        DOUT(out, msg.str());
      }
    } else {
      std::ostringstream msg;
      msg << "Rank-" << d_myworld->myRank() << "   ";
      msg << std::left;
      msg.width(50);
      msg << dtask->getTask()->getName();
      DOUT(out, msg.str());
    }
  }
}

void
printSchedule(const PatchSet* patches, Dout& out, const std::string& where) {
  if (out) {
    std::ostringstream msg;
    msg << "Rank-" << Uintah::Parallel::getMPIRank() << " ";
    msg << std::left;
    msg.width(50);
    msg << where;
    if (patches) {
      msg << " L-" << getLevel(patches)->getIndex();
    }
    DOUT(out, msg.str());
  }
}
void
printSchedule(const LevelP& level, Dout& out, const std::string& where) {
  if (out) {
    std::ostringstream msg;
    msg << "Rank-" << Uintah::Parallel::getMPIRank() << " ";
    msg << std::left;
    msg.width(50);
    msg << where;
    if (level) {
      msg << " L-" << level->getIndex();
    }
    DOUT(out, msg.str());
  }
}

void
printTask(const Patch* patch,
          Dout& out,
          const std::string& where,
          const int timestep,
          const int material,
          const std::string varName) {
  std::ostringstream msg;
  msg << std::left;
  msg.width(50);
  msg << where << std::right << " timestep: " << timestep
      << " matl: " << material << " Var: " << varName;

  printTask(patch, out, msg.str());
}

}  // End namespace Uintah
