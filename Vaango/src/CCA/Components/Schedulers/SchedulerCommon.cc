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

#include <CCA/Components/Schedulers/SchedulerCommon.h>

#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/TaskGraph.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>
#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/LocallyComputedPatchVarMap.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Malloc/Allocator.h>
#include <Core/OS/ProcessInfo.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/Timers/Timers.hpp>
#include <time.h>

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
Uintah::Dout g_schedulercommon_dbg("SchedulerCommon_DBG",
                                   "SchedulerCommon",
                                   "general debug information",
                                   false);
Uintah::Dout g_task_graph_compile("TaskGraphCompile",
                                  "SchedulerCommon",
                                  "task graph compilation info",
                                  false);
} // namespace

namespace Uintah {

SchedulerCommon::SchedulerCommon(const ProcessorGroup* myworld)
  : UintahParallelComponent(myworld)
{
  for (int i = 0; i < Task::TotalDWs; i++) {
    d_dw_map[i] = Task::InvalidDW;
  }

  // Default mapping...
  d_dw_map[Task::OldDW] = 0;
  d_dw_map[Task::NewDW] = 1;

  d_local_patch_var_map = std::make_unique<LocallyComputedPatchVarMap>();
}

SchedulerCommon::~SchedulerCommon()
{
  // list of vars used for AMR regridding
  for (auto& label_material : d_label_matls) {
    for (auto& [label, material] : label_material) {
      material->removeReference();
    }
  }
  d_label_matls.clear();

  if (d_monitoring) {
    if (d_dummy_matl) {
      d_dummy_matl->removeReference();
    }

    // Loop through the global (0) and local (1) tasks
    for (unsigned int i = 0u; i < 2; ++i) {
      for (const auto& it : d_monitoring_tasks[i]) {
        VarLabel::destroy(it.second);
      }
      d_monitoring_values[i].clear();
    }
  }
}

void
SchedulerCommon::setComponents(UintahParallelComponent* comp)
{
  SchedulerCommon* parent = dynamic_cast<SchedulerCommon*>(comp);

  attachPort("load balancer", parent->d_load_balancer);
  attachPort("output", parent->d_output);
  attachPort("application", parent->d_simulator);

  getComponents();
}

void
SchedulerCommon::getComponents()
{
  d_load_balancer = dynamic_cast<LoadBalancer*>(getPort("load balancer"));
  if (!d_load_balancer) {
    throw InternalError(
      "dynamic_cast of 'd_load_balancer' failed!", __FILE__, __LINE__);
  }

  d_output = dynamic_cast<Output*>(getPort("output"));
  if (!d_output) {
    throw InternalError(
      "dynamic_cast of 'd_output' failed!", __FILE__, __LINE__);
  }

  d_simulator = dynamic_cast<SimulationInterface*>(getPort("simulator"));
  if (!d_simulator) {
    throw InternalError(
      "dynamic_cast of 'd_simulator' failed!", __FILE__, __LINE__);
  }
}

void
SchedulerCommon::releaseComponents()
{
  releasePort("load balancer");
  releasePort("output");
  releasePort("simulator");

  d_load_balancer   = nullptr;
  d_output          = nullptr;
  d_simulator       = nullptr;
  d_materialManager = nullptr;
}

void
SchedulerCommon::checkMemoryUse(unsigned long& mem_used,
                                unsigned long& high_water,
                                unsigned long& max_mem_used)
{
  high_water = 0;
  mem_used   = 0;

  if (ProcessInfo::isSupported(ProcessInfo::MEM_SIZE)) {
    mem_used = ProcessInfo::getMemoryResident();
    // printf("1) memuse is %ld (on proc %d)\n", mem_used,
    // Uintah::Parallel::getMPIRank() );
  } else {
    mem_used = (char*)sbrk(0) - s_start_addr;
    // printf("2) memuse is %ld (on proc %d)\n", mem_used,
    // Uintah::Parallel::getMPIRank() );
  }

  if (mem_used > d_max_mem_used) {
    // printf("Max memuse increased\n");
    d_max_mem_used = mem_used;
  }
  max_mem_used = d_max_mem_used;
}

void
SchedulerCommon::resetMaxMemValue()
{
  d_max_mem_used = 0;
}

void
SchedulerCommon::makeTaskGraphDoc([[maybe_unused]] const DetailedTasks* dtasks,
                                  int rank)
{
  // This only happens if "-emit_taskgraphs" is passed to vaango
  if (!d_emit_task_graph) {
    return;
  }

  // ARS NOTE: Outputing and Checkpointing may be done out of snyc
  // now, i.e., turned on just before it happens rather than turned on
  // before the task graph execution.  As such, one should also be
  // checking:

  // d_simulator->activeReductionVariable( "outputInterval" );

  // However, if active, the code below would be called regardless if
  // an output or checkpoint time step or not. That is probably not
  // desired. However, given this code is for debuging it probably
  // fine that it does not happen if doing an output of sync.
  if (!d_output->isOutputTimeStep()) {
    return;
  }

  // make sure to release this DOMDocument after finishing emitting the nodes
  d_graph_doc = ProblemSpec::createDocument("Uintah_TaskGraph");

  ProblemSpecP meta = d_graph_doc->appendChild("Meta");
  meta->appendElement("username", getenv("LOGNAME"));
  time_t t = time(nullptr);
  meta->appendElement("date", ctime(&t));

  d_graph_nodes = d_graph_doc->appendChild("Nodes");

  ProblemSpecP edges_element = d_graph_nodes->appendChild("Edges");

  for (unsigned i = 0; i < d_task_graphs.size(); i++) {
    DetailedTasks* dts = d_task_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->emitEdges(edges_element, rank);
    }
  }
}

bool
SchedulerCommon::useInternalDeps()
{
  // keep track of internal dependencies only if it will emit
  // the taskd_graphs (by default).
  return d_emit_task_graph;
}

void
SchedulerCommon::emitNode(const DetailedTask* dtask,
                          double start,
                          double duration,
                          double execution_duration)
{
  // This only happens if "-emit_taskgraphs" is passed to vaango
  // See makeTaskGraphDoc
  if (d_graph_nodes == nullptr) {
    return;
  }

  ProblemSpecP node = d_graph_nodes->appendChild("node");

  node->appendElement("name", dtask->getName());
  node->appendElement("start", start);
  node->appendElement("duration", duration);
  if (execution_duration > 0) {
    node->appendElement("execution_duration", execution_duration);
  }
}

void
SchedulerCommon::finalizeNodes(int process /* = 0*/)
{
  // This only happens if "-emit_taskgraphs" is passed to vaango
  // See makeTaskGraphDoc
  if (d_graph_doc == nullptr) {
    return;
  }

  std::string timestep_dir(d_output->getLastTimeStepOutputLocation());

  std::ostringstream fname;
  fname << "/taskgraph_" << std::setw(5) << std::setfill('0') << process
        << ".xml";
  std::string file_name(timestep_dir + fname.str());
  d_graph_doc->output(file_name.c_str());

  // m_graphDoc->releaseDocument();
  d_graph_doc   = nullptr;
  d_graph_nodes = nullptr;
}

//______________________________________________________________________
//

void
SchedulerCommon::problemSetup(const ProblemSpecP& prob_spec,
                              const MaterialManagerP& mat_manager)
{
  d_materialManager = mat_manager;

  d_tracking_vars_print_location = PRINT_AFTER_EXEC;

  ProblemSpecP params = prob_spec->findBlock("Scheduler");
  if (params) {
    params->getWithDefault("small_messages", d_use_small_messages, true);
    if (d_use_small_messages) {
      proc0cout
        << "Using small, individual MPI messages (no message combining)\n";
    } else {
      proc0cout << "Using large, combined MPI messages\n";
    }

    // set up the variable tracking system
    setupVarTracker(params);

    // set up task monitoring
    setupTaskMonitoring(params);
  }

  d_no_scrub_vars.insert("refineFlag");
  d_no_scrub_vars.insert("refinePatchFlag");
}

void
SchedulerCommon::setupVarTracker(const ProblemSpecP& params)
{
  ProblemSpecP track = params->findBlock("VarTracker");
  if (track) {
    track->require("start_time", d_tracking_start_time);
    track->require("end_time", d_tracking_end_time);
    track->getWithDefault("level", d_tracking_level, -1);
    track->getWithDefault(
      "start_index", d_tracking_start_index, IntVector(-9, -9, -9));
    track->getWithDefault(
      "end_index", d_tracking_end_index, IntVector(-9, -9, -9));
    track->getWithDefault("patchid", d_tracking_patch_id, -1);

    proc0cout << "\n"
              << "-----------------------------------------------------------\n"
              << "-- Initializing VarTracker...\n"
              << "--  Running from time " << d_tracking_start_time << " to "
              << d_tracking_end_time << "\n"
              << "--  for indices: " << d_tracking_start_index << " to "
              << d_tracking_end_index << "\n";

    ProblemSpecP location = track->findBlock("locations");
    if (location) {
      d_tracking_vars_print_location = 0;
      std::map<std::string, std::string> attributes;
      location->getAttributes(attributes);
      if (attributes["before_comm"] == "true") {
        d_tracking_vars_print_location |= PRINT_BEFORE_COMM;
        proc0cout
          << "--  Printing variable information before communication.\n";
      }
      if (attributes["before_exec"] == "true") {
        d_tracking_vars_print_location |= PRINT_BEFORE_EXEC;
        proc0cout
          << "--  Printing variable information before task execution.\n";
      }
      if (attributes["after_exec"] == "true") {
        d_tracking_vars_print_location |= PRINT_AFTER_EXEC;
        proc0cout
          << "--  Printing variable information after task execution.\n";
      }
    } else {
      // "locations" not specified
      proc0cout << "--  Defaulting to printing variable information after "
                   "task execution.\n";
    }

    for (ProblemSpecP var = track->findBlock("var"); var != nullptr;
         var              = var->findNextBlock("var")) {
      std::map<std::string, std::string> attributes;
      var->getAttributes(attributes);
      std::string name = attributes["label"];
      d_tracking_vars.push_back(name);
      std::string dw = attributes["dw"];

      if (dw == "OldDW") {
        d_tracking_dws.push_back(Task::OldDW);
      } else if (dw == "NewDW") {
        d_tracking_dws.push_back(Task::NewDW);
      } else if (dw == "CoarseNewDW") {
        d_tracking_dws.push_back(Task::CoarseNewDW);
      } else if (dw == "CoarseOldDW") {
        d_tracking_dws.push_back(Task::CoarseOldDW);
      } else if (dw == "ParentOldDW") {
        d_tracking_dws.push_back(Task::ParentOldDW);
      } else if (dw == "ParentOldDW") {
        d_tracking_dws.push_back(Task::ParentNewDW);
      } else {
        // This error message most likely can go away once the .ups validation
        // is put into place:
        printf("WARNING: Hit switch statement default... using NewDW... (This "
               "could possibly be"
               "an error in input file specification.)\n");
        d_tracking_dws.push_back(Task::NewDW);
      }
      proc0cout << "--  Tracking variable '" << name << "' in DataWarehouse '"
                << dw << "'\n";
    }

    for (ProblemSpecP task = track->findBlock("task"); task != nullptr;
         task              = task->findNextBlock("task")) {
      std::map<std::string, std::string> attributes;
      task->getAttributes(attributes);
      std::string name = attributes["name"];
      d_tracking_tasks.push_back(name);
      proc0cout << "--  Tracking variables for specific task: " << name << "\n";
    }
    proc0cout
      << "-----------------------------------------------------------\n\n";
  } else { // Tracking not specified
    // This 'else' won't be necessary once the .ups files are validated... but
    // for now.
    proc0cout << "<VarTracker> not specified in .ups file... no variable "
                 "tracking will take place.\n";
  }
}

void
SchedulerCommon::setupTaskMonitoring(const ProblemSpecP& params)
{
  // Task monitoring variables.
  ProblemSpecP task_monitoring = params->findBlock("TaskMonitoring");
  if (task_monitoring) {
    // Record the task runtime attributes on a per cell basis rather
    // than a per patch basis. Default is per patch.
    task_monitoring->getWithDefault("per_cell", d_monitoring_per_cell, false);

    // Maps for the global tasks to be monitored.
    for (ProblemSpecP attr = task_monitoring->findBlock("attribute");
         attr != nullptr;
         attr = attr->findNextBlock("attribute")) {
      std::string attribute = attr->getNodeValue();

      // Set the variable name AllTasks/ plus the attribute name and
      // store in a map for easy lookup by the attribute name.

      // Note: this modifided name will be needed for saving.
      d_monitoring_tasks[0][attribute] = VarLabel::create(
        "AllTasks/" + attribute, PerPatch<double>::getTypeDescription());

      proc0cout << "--  Monitoring attribute " << attribute << " "
                << "for all tasks. "
                << "VarLabel name = 'AllTasks/" << attribute << "'"
                << std::endl;
    }

    // Maps for the specific tasks to be monitored.
    for (ProblemSpecP task = task_monitoring->findBlock("task");
         task != nullptr;
         task = task->findNextBlock("task")) {
      // Get the task and attribute to be monitored.
      std::map<std::string, std::string> attributes;
      task->getAttributes(attributes);
      std::string taskName  = attributes["name"];
      std::string attribute = attributes["attribute"];

      // Strip off the colons and replace with a forward slash so
      // the tasks are divided by component.
      std::string varName = taskName;
      std::size_t found   = varName.find("::");
      if (found != std::string::npos) {
        varName.replace(found, 2, "/");
      }

      // Set the variable name to the task name plus the attribute
      // name and store in a map for easy lookup by the task and
      // attribute name.

      // Note: this modifided name will be needed for saving.
      d_monitoring_tasks[1][taskName + "::" + attribute] = VarLabel::create(
        varName + "/" + attribute, PerPatch<double>::getTypeDescription());

      proc0cout << "--  Monitoring attribute " << attribute << " "
                << "for task: " << taskName << ".  "
                << "VarLabel name = '" << varName << "/" << attribute << "'"
                << std::endl;
    }
  }

  d_monitoring = (d_monitoring_tasks[0].size() || d_monitoring_tasks[1].size());

  d_monitoring = (d_monitoring_tasks[0].size() || d_monitoring_tasks[1].size());

  if (d_monitoring) {
    d_dummy_matl = std::make_unique<MaterialSubset>();
    d_dummy_matl->add(0);
    d_dummy_matl->addReference();
  }
}

//
// handleError()
//
// The following routine is designed to only print out a given error
// once per error type per variable.  handleError is used by
// printTrackedVars() with each type of error ('errorPosition')
// condition specifically enumerated (by an integer running from 0 to 7).
//
// Returns true if the error message is displayed.
//
bool
handleError(int error_position,
            const std::string& error_message,
            const std::string& variable_name)
{
  static std::vector<std::map<std::string, bool>*> errors_reported(8);

  std::map<std::string, bool>* var_to_reported_map =
    errors_reported[error_position];

  bool reported = (*var_to_reported_map)[variable_name];
  if (!reported) {
    (*var_to_reported_map)[variable_name] = true;
    std::cout << error_message << "\n";
    return true;
  }
  return false;
}

template<class T>
void
SchedulerCommon::printTrackedValues(GridVariable<T>* var,
                                    const IntVector& start,
                                    const IntVector& end)
{
  std::ostringstream message;
  for (int z = start.z(); z < end.z() + 1;
       z++) { // add 1 to high to include x+,y+,z+ extraCells
    for (int y = start.y(); y < end.y() + 1; y++) {
      message << d_myworld->myRank() << "  ";

      for (int x = start.x(); x < end.x() + 1; x++) {
        IntVector c(x, y, z);
        message << " " << c << ": " << (*var)[c];
      }
      message << std::endl;
    }
    message << std::endl;
  }
  DOUT(true, message.str());
}

void
SchedulerCommon::printTrackedVars(DetailedTask* dtask, int when)
{
  bool printed_header = false;
  unsigned task_num;
  for (task_num = 0; task_num < d_tracking_tasks.size(); task_num++) {
    if (d_tracking_tasks[task_num] == dtask->getTask()->getName()) {
      break;
    }
  }

  // Print for all tasks unless one is specified (but disclude DataArchiver
  // tasks)
  if ((task_num == d_tracking_tasks.size() && d_tracking_tasks.size() != 0) ||
      ((std::string(dtask->getTask()->getName())).substr(0, 12) ==
       "DataArchiver")) {
    return;
  }

  if (d_tracking_start_time > d_simulator->getSimTime() ||
      d_tracking_end_time < d_simulator->getSimTime()) {
    return;
  }

  for (std::uint32_t i = 0; i < d_tracking_vars.size(); i++) {
    bool printed_var_name = false;

    // that DW may not have been mapped....
    if (dtask->getTask()->mapDataWarehouse(d_tracking_dws[i]) < 0 ||
        dtask->getTask()->mapDataWarehouse(d_tracking_dws[i]) >=
          static_cast<int>(d_dws.size())) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable ("
           << d_tracking_vars[i] << ") DW is out of range.\n";
      handleError(0, mesg.str(), d_tracking_vars[i]);
      continue;
    }

    OnDemandDataWarehouse* dw =
      d_dws[dtask->getTask()->mapDataWarehouse(d_tracking_dws[i])].get();

    if (dw == nullptr) { // old on initialization timestep
      continue;
    }

    // get the level here, as the grid can be different between the old and new
    // DW
    const Grid* grid = dw->getGrid();
    int level_num;

    if (d_tracking_level == -1) {
      level_num = grid->numLevels() - 1;
    } else {
      level_num = d_tracking_level;
      if (level_num >= grid->numLevels()) {
        continue;
      }
    }

    const LevelP level    = grid->getLevel(level_num);
    const VarLabel* label = VarLabel::find(d_tracking_vars[i]);

    std::cout.precision(16);

    if (!label) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable ("
           << d_tracking_vars[i] << ") because label is nullptr.\n";
      handleError(1, mesg.str(), d_tracking_vars[i]);
      continue;
    }

    const PatchSubset* patches = dtask->getPatches();

    // a once-per-proc task is liable to have multiple levels, and thus calls to
    // getLevel(patches) will fail
    // The task could also run on a different level (coarse or fine).
    const Task::TaskType task_type = dtask->getTask()->getType();
    const bool not_once_per_proc   = (task_type != Task::OncePerProc);
    const bool not_hypre           = (task_type != Task::Hypre);
    const int level_index          = getLevel(patches)->getIndex();
    const bool not_right_level     = (!patches || level_index != level_num);

    if (not_once_per_proc && not_hypre && not_right_level) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable ("
           << d_tracking_vars[i] << ") because patch is non-standard.\n";
      handleError(2, mesg.str(), d_tracking_vars[i]);
      continue;
    }

    if (!patches) {
      return;
    }

    for (auto patch : *patches) {
      if (d_tracking_patch_id != -1 && d_tracking_patch_id != patch->getID()) {
        continue;
      }

      // don't print ghost patches (dw->get will yell at you)
      if ((d_tracking_dws[i] == Task::OldDW &&
           d_load_balancer->getOldProcessorAssignment(patch) !=
             d_myworld->myRank()) ||
          (d_tracking_dws[i] == Task::NewDW &&
           d_load_balancer->getPatchwiseProcessorAssignment(patch) !=
             d_myworld->myRank())) {
        continue;
      }

      const TypeDescription* td = label->typeDescription();
      Patch::VariableBasis basis =
        patch->translateTypeToBasis(td->getType(), false);
      IntVector start = Max(patch->getExtraLowIndex(basis, IntVector(0, 0, 0)),
                            d_tracking_start_index);
      IntVector end   = Min(patch->getExtraHighIndex(basis, IntVector(0, 0, 0)),
                          d_tracking_end_index);

      // loop over matls too
      for (size_t m = 0; m < d_materialManager->getNumMaterials(); m++) {
        if (!dw->exists(label, m, patch)) {
          std::ostringstream mesg;
          mesg << "WARNING: VarTracker: Not printing requested variable ("
               << d_tracking_vars[i] << ") because it does not exist in DW.\n"
               << "            Patch is: " << *patch << "\n";
          if (handleError(3, mesg.str(), d_tracking_vars[i])) {
            // std::cout << "         DW contains (material: " << m << ")\n";
            // dw->print();
          }
          continue;
        }
        if (!(start.x() < end.x() && start.y() < end.y() &&
              start.z() < end.z())) {
          continue;
        }
        const TypeDescription::Type sub_type = td->getSubType()->getType();
        if (sub_type != TypeDescription::Type::double_type &&
            sub_type != TypeDescription::Type::float_type &&
            sub_type != TypeDescription::Type::int_type &&
            sub_type != TypeDescription::Type::Vector) {
          std::ostringstream mesg;
          mesg << "WARNING: VarTracker: Not printing requested variable ("
               << d_tracking_vars[i] << ") because its type is not supported:\n"
               << "             " << td->getName() << "\n";
          handleError(4, mesg.str(), d_tracking_vars[i]);
          continue;
        }

        // pending the task that allocates the var, we may not have allocated it
        // yet
        GridVariableBase* v;
        switch (td->getType()) {
          case TypeDescription::Type::CCVariable:
          case TypeDescription::Type::NCVariable:
          case TypeDescription::Type::SFCXVariable:
          case TypeDescription::Type::SFCYVariable:
          case TypeDescription::Type::SFCZVariable:
            v = dynamic_cast<GridVariableBase*>(
              dw->d_var_DB.get(label, m, patch));
            break;
          default:
            throw InternalError(
              "Cannot track var type of non-grid-type", __FILE__, __LINE__);
            break;
        }

        start = Max(start, v->getLow());
        end   = Min(end, v->getHigh());
        if (!(start.x() < end.x() && start.y() < end.y() &&
              start.z() < end.z())) {
          continue;
        }

        if (!printed_header) {
          std::string location;
          switch (when) {
            case PRINT_BEFORE_COMM:
              location = " before communication of ";
              break;
            case PRINT_BEFORE_EXEC:
              location = " before execution of ";
              break;
            case PRINT_AFTER_EXEC:
              location = " after execution of ";
              break;
          }
          std::cout << d_myworld->myRank() << location << *dtask << std::endl;
          printed_header = true;
        }

        if (!printed_var_name) {
          std::cout << d_myworld->myRank()
                    << "  Variable: " << d_tracking_vars[i] << ", DW "
                    << dw->getID() << ", Patch " << patch->getID() << ", Matl "
                    << m << std::endl;
        }

        switch (sub_type) {
          case TypeDescription::Type::double_type: {
            GridVariable<double>* var = dynamic_cast<GridVariable<double>*>(v);
            printTrackedValues<double>(var, start, end);
          } break;
          case TypeDescription::Type::float_type: {
            GridVariable<float>* var = dynamic_cast<GridVariable<float>*>(v);
            printTrackedValues<float>(var, start, end);
          } break;
          case TypeDescription::Type::int_type: {
            GridVariable<int>* var = dynamic_cast<GridVariable<int>*>(v);
            printTrackedValues<int>(var, start, end);
          } break;
          case TypeDescription::Type::Vector: {
            GridVariable<Vector>* var = dynamic_cast<GridVariable<Vector>*>(v);
            printTrackedValues<Vector>(var, start, end);
          } break;
          default:
            break;
        } // end case variable type
      }   // end for m(aterials)
    }     // end for p(atches)
  }       // end for i : trackingVars.size()
} // end printTrackedVars()

void
SchedulerCommon::addTaskGraph(Scheduler::tgType type, int index)
{
  std::unique_ptr<TaskGraph> tg =
    std::make_unique<TaskGraph>(this, d_myworld, type, index);
  tg->initialize();
  d_task_graphs.push_back(std::move(tg));
}

void
SchedulerCommon::addTask(Task* task,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         int tg_num /* = -1 */)
{
  // Save the DW map
  task->setMapping(d_dw_map);

  bool is_init = d_is_init_timestep || d_is_restart_init_timestep;

  DOUTR(g_schedulercommon_dbg,
        " adding Task: " << task->getName()
                         << ",  # patches: " << (patches ? patches->size() : 0)
                         << ",    # matls: " << (matls ? matls->size() : 0)
                         << ", task-graph: "
                         << ((tg_num < 0) ? (is_init ? "init-tg" : "all")
                                          : std::to_string(tg_num)));

  if (!is_init && tg_num >= (int)d_task_graphs.size()) {
    std::ostringstream msg;
    msg << task->getName() << "::addTask(),  taskgraph index (" << tg_num
        << ") >= num_taskgraphs (" << d_task_graphs.size() << ")";
    throw InternalError(msg.str(), __FILE__, __LINE__);
  }

  std::shared_ptr<Task> task_sp(task);
  addTask(task_sp, patches, matls, tg_num);

  // separate out the standard from the distal ghost cell requirements - for
  // loadbalancer This isn't anything fancy, and could be expanded/modified down
  // the road. It just gets a max ghost cell extent for anything less than
  // MAX_HALO_DEPTH, and another max ghost cell extent for anything >=
  // MAX_HALO_DEPTH.  The idea is that later we will create two neighborhoods
  // with max extents for each as determined here.
  for (auto dep = task->getRequires(); dep != nullptr; dep = dep->next) {
    if (dep->num_ghost_cells >= MAX_HALO_DEPTH) {
      if (dep->num_ghost_cells > this->d_max_distal_ghost_cells) {
        this->d_max_distal_ghost_cells = dep->num_ghost_cells;
      }
    } else {
      if (dep->num_ghost_cells > this->d_max_ghost_cells) {
        this->d_max_ghost_cells = dep->num_ghost_cells;
      }
    }
  }

  if (task->m_max_level_offset > this->d_max_level_offset) {
    this->d_max_level_offset = task->m_max_level_offset;
  }

  // add to init-requires.  These are the vars which require from the OldDW that
  // we'll need for checkpointing, switching, and the like. In the case of
  // treatAsOld Vars, we handle them because something external to the taskgraph
  // needs it that way (i.e., Regridding on a restart requires checkpointed
  // refineFlags).
  for (auto dep = task->getRequires(); dep != nullptr; dep = dep->next) {
    if (isOldDW(dep->mapDataWarehouse()) ||
        d_treat_as_old_vars.find(dep->var->getName()) !=
          d_treat_as_old_vars.end()) {
      d_init_requires.push_back(dep);
      d_init_required_vars.insert(dep->var);
    }
  }

  // for the treat-as-old vars, go through the computes and add them.
  // we can (probably) safely assume that we'll avoid duplicates, since if they
  // were inserted in the above, they wouldn't need to be marked as such
  for (auto dep = task->getComputes(); dep != nullptr; dep = dep->next) {
    d_computed_vars.insert(dep->var);
    if (d_treat_as_old_vars.find(dep->var->getName()) !=
        d_treat_as_old_vars.end()) {
      d_init_requires.push_back(dep);
      d_init_required_vars.insert(dep->var);
    }
  }

  // create reduction task if computes included one or more reduction vars
  for (auto dep = task->getComputes(); dep != 0; dep = dep->next) {
    if (dep->var->typeDescription()->isReductionVariable()) {
      int level_idx =
        dep->reduction_level ? dep->reduction_level->getIndex() : -1;
      int dw = dep->mapDataWarehouse();

      if (!dep->var->isReductionTask()) {
        DOUTR(g_schedulercommon_dbg,
              " Skipping Reduction task for multi compute variable: "
                << dep->var->getName() << " on level " << level_idx << ", DW "
                << dw);
        continue;
      }
      DOUTR(g_schedulercommon_dbg,
            " Creating Reduction task for variable: "
              << dep->var->getName() << " on level " << level_idx << ", DW "
              << dw);

      std::ostringstream taskname;
      taskname << "SchedulerCommon::Reduction: " << dep->var->getName()
               << ", level " << level_idx << ", dw " << dw;

      // Pointer is owned by both d_reduction_tasks and the task graph
      std::shared_ptr<Task> reduction_task =
        std::make_shared<Task>(taskname.str(), Task::Reduction);

      int dw_map[Task::TotalDWs];

      for (int i = 0; i < Task::TotalDWs; i++) {
        dw_map[i] = Task::InvalidDW;
      }

      dw_map[Task::OldDW] = Task::NoDW;
      dw_map[Task::NewDW] = dw;
      reduction_task->setMapping(dw_map);

      if (dep->matls != nullptr) {
        reduction_task->modifies(
          dep->var, dep->reduction_level, dep->matls, Task::OutOfDomain);
        const DataWarehouse* const_dw = get_dw(dw);
        for (int i = 0; i < dep->matls->size(); i++) {
          int mat_idx = dep->matls->get(i);
          VarLabelMatl<Level, DataWarehouse> key(
            dep->var, mat_idx, dep->reduction_level, const_dw);
          d_reduction_tasks[key] = reduction_task;
        }
      } else {
        for (int m = 0; m < task->getMaterialSet()->size(); m++) {
          auto mat_subset = task->getMaterialSet()->getSubset(m);
          reduction_task->modifies(
            dep->var, dep->reduction_level, mat_subset, Task::OutOfDomain);
          const DataWarehouse* const_dw = get_dw(dw);
          for (int i = 0; i < mat_subset->size(); i++) {
            int mat_idx = mat_subset->get(i);
            VarLabelMatl<Level, DataWarehouse> key(
              dep->var, mat_idx, dep->reduction_level, const_dw);

            // For reduction variables there may be multiple computes
            // each of which will create reduction task. The last
            // reduction task should be kept. This is because the
            // tasks do not get sorted.
            if (d_reduction_tasks.find(key) == d_reduction_tasks.end()) {
              DOUTR(g_schedulercommon_dbg,
                    " 2) Excluding previous reduction task for variable: "
                      << dep->var->getName() << " on level " << level_idx
                      << ", DW " << dw << " dep->m_reduction_level "
                      << dep->reduction_level << " material index " << mat_idx);
            }
            d_reduction_tasks[key] = reduction_task;
          }
        }
      }

      // add reduction task to the task graphs
      addTask(reduction_task, nullptr, task->getMaterialSet(), tg_num);
    }
  }
}

void
SchedulerCommon::addTask(std::shared_ptr<Task> task,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const int tg_num)
{
  // During initialization or restart, there is only one task graph.
  if (d_is_init_timestep || d_is_restart_init_timestep) {
    d_task_graphs[d_task_graphs.size() - 1]->addTask(task, patches, matls);
    d_num_tasks++;
  } else {
    // Add it to all "Normal" task graphs (default value == -1, from public
    // addTask() method).
    if (tg_num < 0) {
      for (unsigned int i = 0; i < d_task_graphs.size(); i++) {
        d_task_graphs[i]->addTask(task, patches, matls);
        d_num_tasks++;
      }
    }
    // Otherwise, add this task to a specific task graph.
    else {
      d_task_graphs[tg_num]->addTask(task, patches, matls);
      d_num_tasks++;
    }
  }
}

//
void
SchedulerCommon::initialize(int num_old_dw /* =1 */, int num_new_dw /* =1 */)
{
  // doesn't really do anything except initialize/clear the taskgraph
  //   if the default parameter values are used
  int num_dw  = num_old_dw + num_new_dw;
  int old_num = static_cast<int>(d_dws.size());

  // in AMR cases we will often need to move from many new DWs to one.  In those
  // cases, move the last NewDW to be the next new one.
  if (old_num - d_num_old_dws > 1) {
    d_dws[num_dw - 1] = std::move(d_dws[old_num - 1]);
  }

  // Clear out the data warehouse so that memory will be freed
  for (int i = num_dw; i < old_num; i++) {
    d_dws[i] = nullptr;
  }
  d_dws.resize(num_dw);
  for (; old_num < num_dw; old_num++) {
    d_dws[old_num] = nullptr;
  }
  d_num_old_dws = num_old_dw;

  // clear the taskgraphs, and set the first one
  d_task_graphs.clear();

  d_init_requires.clear();
  d_init_required_vars.clear();
  d_computed_vars.clear();

  d_num_tasks              = 0;
  d_max_ghost_cells        = 0;
  d_max_distal_ghost_cells = 0;
  d_max_level_offset       = 0;

  d_reduction_tasks.clear();

  // During initialization or restart, use only one task graph
  bool is_init           = d_is_init_timestep || d_is_restart_init_timestep;
  size_t num_task_graphs = (is_init) ? 1 : d_num_task_graphs;

  for (size_t i = 0; i < num_task_graphs; ++i) {
    addTaskGraph(NormalTaskGraph, i);
  }
}

void
SchedulerCommon::setParentDWs(DataWarehouse* parent_old_dw,
                              DataWarehouse* parent_new_dw)
{
  OnDemandDataWarehouse* p_old =
    dynamic_cast<OnDemandDataWarehouse*>(parent_old_dw);
  OnDemandDataWarehouse* p_new =
    dynamic_cast<OnDemandDataWarehouse*>(parent_new_dw);
  if (parent_old_dw && parent_new_dw) {
    ASSERT(p_old != 0);
    ASSERT(p_new != 0);
    ASSERT(d_num_old_dws > 2);

    // Transfer ownership to d_dws
    d_dws[0].reset(p_old);
    d_dws[1].reset(p_new);
  }
}

void
SchedulerCommon::clearMappings()
{
  for (int i = 0; i < Task::TotalDWs; i++) {
    d_dw_map[i] = -1;
  }
}

void
SchedulerCommon::mapDataWarehouse(Task::WhichDW which, int dw_tag)
{
  ASSERTRANGE(which, 0, Task::TotalDWs);
  ASSERTRANGE(dw_tag, 0, static_cast<int>(d_dws.size()));
  d_dw_map[which] = dw_tag;
}

DataWarehouse*
SchedulerCommon::get_dw(int idx)
{
  if (0 <= idx && idx < static_cast<int>(d_dws.size())) {
    return d_dws[idx].get();
  }
  return nullptr;
}

//______________________________________________________________________
//
DataWarehouse*
SchedulerCommon::getLastDW(void)
{
  return get_dw(static_cast<int>(d_dws.size()) - 1);
}

void
SchedulerCommon::advanceDataWarehouse(const GridP& grid,
                                      bool initialization /*=false*/)
{
  DOUTR(g_schedulercommon_dbg,
        " advanceDataWarehouse, numDWs = " << d_dws.size());
  ASSERT(d_dws.size() >= 2);

  // TODO: This can cost roughly 1 millisecond of time.  Find a way to reuse
  // data warehouses if possible?  Brad March 6 2018

  // The last becomes last old, and the rest are new
  d_dws[d_num_old_dws - 1] = std::move(d_dws[d_dws.size() - 1]);
  if (d_dws.size() == 2 && d_dws[0] == nullptr) {
    // first datawarehouse -- indicate that it is the "initialization" dw.
    int generation = d_generation++;

    d_dws[1] = std::make_unique<OnDemandDataWarehouse>(
      d_myworld, this, generation, grid, true /* initialization dw */);
  } else {
    for (int i = d_num_old_dws; i < static_cast<int>(d_dws.size()); i++) {
      // in AMR initial cases, you can still be in initialization when you
      // advance again
      replaceDataWarehouse(i, grid, initialization);
    }
  }
}

void
SchedulerCommon::fillDataWarehouses(const GridP& grid)
{
  for (int i = d_num_old_dws; i < static_cast<int>(d_dws.size()); i++) {
    if (!d_dws[i]) {
      replaceDataWarehouse(i, grid);
    }
  }
}

void
SchedulerCommon::replaceDataWarehouse(int index,
                                      const GridP& grid,
                                      bool initialization /*=false*/)
{
  d_dws[index] = std::make_unique<OnDemandDataWarehouse>(
    d_myworld, this, d_generation++, grid, initialization);
  if (initialization) {
    return;
  }
  for (unsigned i = 0; i < d_task_graphs.size(); i++) {
    DetailedTasks* dts = d_task_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->copyoutDWKeyDatabase(d_dws[index].get());
    }
  }
  d_dws[index]->doReserve();
}

const std::vector<const Patch*>*
SchedulerCommon::getSuperPatchExtents(const VarLabel* label,
                                      [[maybe_unused]] int mat_index,
                                      const Patch* patch,
                                      Ghost::GhostType requested_ghost_type,
                                      int requested_num_ghost_cells,
                                      IntVector& required_low,
                                      IntVector& required_high,
                                      IntVector& requested_low,
                                      IntVector& requested_high) const
{
  const SuperPatch* connected_patch_group =
    d_local_patch_var_map->getConnectedPatchGroup(patch);
  if (connected_patch_group == 0) {
    return 0;
  }

  SuperPatch::Region requested_extents = connected_patch_group->getRegion();
  SuperPatch::Region required_extents  = connected_patch_group->getRegion();

  // expand to cover the entire connected patch group
  for (unsigned int i = 0; i < connected_patch_group->getBoxes().size(); i++) {
    // get the minimum extents containing both the expected ghost cells
    // to be needed and the given ghost cells.
    const Patch* member_patch = connected_patch_group->getBoxes()[i];

    Patch::VariableBasis basis =
      Patch::translateTypeToBasis(label->typeDescription()->getType(), true);

    IntVector low_offset = IntVector(0, 0, 0), high_offset = IntVector(0, 0, 0);

    // set requiredLow and requiredHigh as extents without ghost cells
    member_patch->computeExtents(basis,
                                 label->getBoundaryLayer(),
                                 low_offset,
                                 high_offset,
                                 required_low,
                                 required_high);

    // compute ghost cell offsets
    Patch::getGhostOffsets(basis,
                           requested_ghost_type,
                           requested_num_ghost_cells,
                           low_offset,
                           high_offset);

    // set requestedLow and requestedHigh as extents with ghost cells
    member_patch->computeExtents(basis,
                                 label->getBoundaryLayer(),
                                 low_offset,
                                 high_offset,
                                 requested_low,
                                 requested_high);

    SuperPatch::Region required_region =
      SuperPatch::Region(required_low, required_high);
    required_extents = required_extents.enclosingRegion(required_region);
    SuperPatch::Region requested_region =
      SuperPatch::Region(requested_low, requested_high);
    requested_extents = requested_extents.enclosingRegion(requested_region);
  }

  required_low   = required_extents.low_;
  required_high  = required_extents.high_;
  requested_low  = requested_extents.low_;
  requested_high = requested_extents.high_;

  // requested extents must enclose the required extents at lesst.
  ASSERTEQ(Min(required_low, requested_low), requested_low);
  ASSERTEQ(Max(required_high, requested_high), requested_high);

  return &connected_patch_group->getBoxes();
}

void
SchedulerCommon::logMemoryUse()
{
  if (!d_mem_log_file) {
    std::ostringstream fname;
    fname << "uintah_memuse.log.p" << std::setw(5) << std::setfill('0')
          << d_myworld->myRank() << "." << d_myworld->nRanks();
    d_mem_log_file = std::make_unique<std::ofstream>(fname.str().c_str());
    if (!d_mem_log_file) {
      std::cerr << "Error opening file: " << fname.str() << '\n';
    }
  }
  *d_mem_log_file << '\n';
  unsigned long total = 0;
  for (int i = 0; i < static_cast<int>(d_dws.size()); i++) {
    char* name;
    if (i == 0) {
      name = const_cast<char*>("OldDW");
    } else if (i == static_cast<int>(d_dws.size()) - 1) {
      name = const_cast<char*>("NewDW");
    } else {
      name = const_cast<char*>("IntermediateDW");
    }
    if (d_dws[i]) {
      d_dws[i]->logMemoryUse(*d_mem_log_file, total, name);
    }
  }

  for (unsigned i = 0; i < d_task_graphs.size(); i++) {
    DetailedTasks* dts = d_task_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->logMemoryUse(*d_mem_log_file, total, "Taskgraph");
    }
  }
  *d_mem_log_file << "Total: " << total << '\n';
  d_mem_log_file->flush();
}

// Make and return a map that maps strings to VarLabels of
// that name and a list of material indices for which that
// variable is valid (according to d_allcomps in graph).
std::unique_ptr<Scheduler::VarLabelMaterialMap>
SchedulerCommon::makeVarLabelMaterialMap()
{
  std::unique_ptr<VarLabelMaterialMap> result =
    std::make_unique<VarLabelMaterialMap>();
  for (unsigned i = 0; i < d_task_graphs.size(); i++) {
    d_task_graphs[i]->makeVarLabelMaterialMap(result.get());
  }
  return result;
}

void
SchedulerCommon::doEmitTaskGraphDocs()
{
  d_emit_task_graph = true;
}

void
SchedulerCommon::compile()
{
  GridP grid = const_cast<Grid*>(getLastDW()->getGrid());
  GridP old_grid;

  if (d_dws[0]) {
    old_grid = const_cast<Grid*>(get_dw(0)->getGrid());
  }

  if (d_num_tasks > 0) {
    DOUTR(g_schedulercommon_dbg,
          d_myworld->myRank() << " SchedulerCommon starting compile");

    int task_graph_id = 0;
    for (auto& task_graph : d_task_graphs) {
      DOUTR(g_schedulercommon_dbg,
            d_myworld->myRank()
              << "  Compiling task graph#" << task_graph_id << " of "
              << d_task_graphs.size() << " with " << d_num_tasks << "tasks");

      Timers::Simple tg_compile_timer;
      tg_compile_timer.start();

      // check if this TG has any tasks with halo requirements > MAX_HALO_DEPTH
      // (determined in public SchedulerCommon::addTask())
      const bool has_distal_reqs = task_graph->getDistalRequires();

      // NOTE: this single call is where all the TG compilation complexity
      // arises (dependency analysis for auto MPI mesgs)
      [[maybe_unused]] DetailedTasks* dts = task_graph->createDetailedTasks(
        useInternalDeps(), grid, old_grid, has_distal_reqs);

      double compile_time = tg_compile_timer().seconds();

      bool is_init = d_is_init_timestep || d_is_restart_init_timestep;

      DOUT(g_task_graph_compile,
           "Rank-" << std::left << std::setw(5) << d_myworld->myRank()
                   << " time to compile TG-" << std::setw(4)
                   << (is_init ? "init-tg"
                               : std::to_string(task_graph->getIndex()))
                   << ": " << compile_time << " (sec)");
      ++task_graph_id;
    }

    // check scheduler at runtime, that all ranks are executing the same size TG
    // (excluding spatial tasks)
    verifyChecksum();
    DOUTR(g_schedulercommon_dbg, " SchedulerCommon finished compile");
  } else {
    return; // no tasks and nothing to do
  }

  d_local_patch_var_map->reset();

  for (int i = 0; i < grid->numLevels(); i++) {
    const PatchSubset* patches =
      d_load_balancer->getPerProcessorPatchSet(grid->getLevel(i))
        ->getSubset(d_myworld->myRank());
    if (patches->size() > 0) {
      d_local_patch_var_map->addComputedPatchSet(patches);
    }
  }
  for (unsigned int dw = 0; dw < d_dws.size(); dw++) {
    if (d_dws[dw].get()) {
      for (unsigned i = 0; i < d_task_graphs.size(); i++) {
        DetailedTasks* dts = d_task_graphs[i]->getDetailedTasks();
        dts->copyoutDWKeyDatabase(d_dws[dw].get());
      }
      d_dws[dw]->doReserve();
    }
  }

  // create SuperPatch groups - only necessary if
  // OnDemandDataWarehouse::s_combine_memory == true, by default it is false
  d_local_patch_var_map->makeGroups();
}

bool
SchedulerCommon::isOldDW(int idx) const
{
  ASSERTRANGE(idx, 0, static_cast<int>(d_dws.size()));
  return idx < d_num_old_dws;
}

bool
SchedulerCommon::isNewDW(int idx) const
{
  ASSERTRANGE(idx, 0, static_cast<int>(d_dws.size()));
  return idx >= d_num_old_dws;
}

void
SchedulerCommon::finalizeTimestep()
{
  finalizeNodes(d_myworld->myRank());
  for (unsigned int i = d_num_old_dws; i < d_dws.size(); i++) {
    d_dws[i]->finalize();
  }
}

void
SchedulerCommon::scheduleAndDoDataCopy(const GridP& grid)
{
  Timers::Simple timer;
  timer.start();

  // TODO - use the current initReqs and push them back, instead of doing
  // this... clear the old list of vars and matls
  for (auto& label_mat_map : d_label_matls) {
    for (auto& [label, material_subset] : label_mat_map) {
      material_subset->removeReference();
    }
  }

  d_label_matls.clear();
  d_label_matls.resize(grid->numLevels());

  // produce a map from all tasks' requires from the Old DW.  Store the varlabel
  // and matls
  // TODO - only do this ONCE.
  for (auto& tg : d_task_graphs) {
    for (int i = 0; i < tg->getNumTasks(); i++) {
      Task* task = tg->getTask(i);
      if (task->getType() == Task::Output) {
        continue;
      }
      for (auto dep = task->getRequires(); dep != 0; dep = dep->next) {
        bool copy_this_var = (dep->whichdw == Task::OldDW);

        // override to manually copy a var
        if (!copy_this_var) {
          if (d_copy_data_vars.find(dep->var->getName()) !=
              d_copy_data_vars.end()) {
            copy_this_var = true;
          }
        }

        // Overide the logic above.  There are PerPatch variables that
        // cannot/shouldn't be copied to the new grid, for example
        // PerPatch<FileInfo>.
        if (d_not_copy_data_vars.count(dep->var->getName()) > 0) {
          copy_this_var = false;
        }

        if (copy_this_var) {
          // Take care of reduction/sole variables in a different section
          TypeDescription::Type dep_type =
            dep->var->typeDescription()->getType();
          if (dep_type == TypeDescription::Type::ReductionVariable ||
              dep_type == TypeDescription::Type::SoleVariable) {
            continue;
          }

          // check the level on the case where variables are only computed on
          // certain levels
          const PatchSet* ps = task->getPatchSet();
          int level          = -1;
          if (dep->patches) { // just in case the task is over multiple
                              // levels...
            level = getLevel(dep->patches)->getIndex();
          } else if (ps) {
            level = getLevel(ps)->getIndex();
          }

          // we don't want data with an invalid level, or requiring from a
          // different level (remember, we are using an old task graph).  That
          // willbe copied later (and chances are, it's to modify anyway).
          if (level == -1 || level > grid->numLevels() - 1 ||
              dep->patches_dom == Task::CoarseLevel ||
              dep->patches_dom == Task::FineLevel) {
            continue;
          }

          const MaterialSubset* matSubset =
            (dep->matls != 0) ? dep->matls
                              : dep->task->getMaterialSet()->getUnion();

          // if var was already found, make a union of the materials
          std::unique_ptr<MaterialSubset> matls =
            std::make_unique<MaterialSubset>(matSubset->getVector());
          matls->addReference();

          MaterialSubset* union_matls;
          union_matls = (d_label_matls[level][dep->var]).get();

          if (union_matls) {
            for (int i = 0; i < union_matls->size(); i++) {
              if (!matls->contains(union_matls->get(i))) {
                matls->add(union_matls->get(i));
              }
            }
          }
          matls->sort();
          d_label_matls[level][dep->var] = std::move(matls);
        }
      }
    }
  }

  this->initialize(1, 1);
  this->advanceDataWarehouse(grid, true);
  this->clearMappings();
  this->mapDataWarehouse(Task::OldDW, 0);
  this->mapDataWarehouse(Task::NewDW, 1);
  this->mapDataWarehouse(Task::CoarseOldDW, 0);
  this->mapDataWarehouse(Task::CoarseNewDW, 1);

  DataWarehouse* old_dw = this->get_dw(0);
  DataWarehouse* new_dw = this->getLastDW();

  old_dw->setScrubbing(DataWarehouse::ScrubNone);
  new_dw->setScrubbing(DataWarehouse::ScrubNone);

  const Grid* old_grid = old_dw->getGrid();

  std::vector<Task*> data_tasks;
  std::vector<Handle<PatchSet>> refine_patch_sets(grid->numLevels(),
                                                  (PatchSet*)nullptr);
  std::vector<Handle<PatchSet>> copy_patch_sets(grid->numLevels(),
                                                (PatchSet*)nullptr);
  SchedulerP sched(dynamic_cast<Scheduler*>(this));

  d_is_copy_data_timestep = true;

  for (int level = 0; level < grid->numLevels(); level++) {
    LevelP new_level = grid->getLevel(level);
    if (level > 0) {
      if (level >= old_grid->numLevels()) {
        // new level - refine everywhere
        refine_patch_sets[level] =
          const_cast<PatchSet*>(new_level->eachPatch());
        copy_patch_sets[level] = scinew PatchSet;
      } else if (level < old_grid->numLevels()) {
        // find patches with new space - but temporarily, refine everywhere...
        refine_patch_sets[level] = scinew PatchSet;
        copy_patch_sets[level]   = scinew PatchSet;

        std::vector<int> my_patch_ids;
        LevelP old_level = old_dw->getGrid()->getLevel(level);

        // go through the patches, and find if there are patches that weren't
        // entirely covered by patches on the old grid, and interpolate them.
        // then after, copy the data, and if necessary, overwrite interpolated
        // data
        const PatchSubset* ps =
          d_load_balancer->getPerProcessorPatchSet(new_level)->getSubset(
            d_myworld->myRank());

        // for each patch I own
        for (int p = 0; p < ps->size(); p++) {
          const Patch* new_patch = ps->get(p);

          // get the low/high for what we'll need to get
          IntVector low_index, high_index;
          low_index  = new_patch->getCellLowIndex();
          high_index = new_patch->getCellHighIndex();

          // find if area on the new patch was not covered by the old patches
          IntVector dist  = high_index - low_index;
          int total_cells = dist.x() * dist.y() * dist.z();
          int sum         = 0;
          Patch::selectType old_patches;
          old_level->selectPatches(low_index, high_index, old_patches);

          // compute volume of overlapping regions
          for (size_t old = 0; old < old_patches.size(); old++) {
            const Patch* old_patch = old_patches[old];
            IntVector old_low      = old_patch->getCellLowIndex();
            IntVector old_high     = old_patch->getCellHighIndex();

            IntVector low  = Max(old_low, low_index);
            IntVector high = Min(old_high, high_index);
            IntVector dist = high - low;
            sum += dist.x() * dist.y() * dist.z();
          } // for oldPatches

          if (sum != total_cells) {
            my_patch_ids.push_back(new_patch->getID());
          }
        } // for patch

        // Gather size from all processors
        int mycount = my_patch_ids.size();
        std::vector<int> counts(d_myworld->nRanks());
        Uintah::MPI::Allgather(
          &mycount, 1, MPI_INT, &counts[0], 1, MPI_INT, d_myworld->getComm());

        // compute recieve array offset and size
        std::vector<int> displs(d_myworld->nRanks());
        int pos = 0;

        for (int p = 0; p < d_myworld->nRanks(); p++) {
          displs[p] = pos;
          pos += counts[p];
        }

        std::vector<int> all_patch_ids(pos); // receive array;
        Uintah::MPI::Allgatherv(&my_patch_ids[0],
                                counts[d_myworld->myRank()],
                                MPI_INT,
                                &all_patch_ids[0],
                                &counts[0],
                                &displs[0],
                                MPI_INT,
                                d_myworld->getComm());
        // make refine_patch_sets from patch ids
        std::set<int> allPatchIDset(all_patch_ids.begin(), all_patch_ids.end());

        for (auto iter = new_level->patchesBegin();
             iter != new_level->patchesEnd();
             ++iter) {
          Patch* new_patch = *iter;
          if (allPatchIDset.find(new_patch->getID()) != allPatchIDset.end()) {
            refine_patch_sets[level]->add(new_patch);
          } else {
            copy_patch_sets[level]->add(new_patch);
          }
        }
      }

      if (refine_patch_sets[level]->size() > 0) {
        DOUTR(g_schedulercommon_dbg,
              "  Calling scheduleRefine for patches "
                << *refine_patch_sets[level].get_rep());
        d_simulator->scheduleRefine(refine_patch_sets[level].get_rep(), sched);
      }
    } else {
      refine_patch_sets[level] = scinew PatchSet;
      copy_patch_sets[level]   = const_cast<PatchSet*>(new_level->eachPatch());
    }

    //  Scheduling for copyDataToNewGrid
    if (copy_patch_sets[level]->size() > 0) {
      data_tasks.push_back(scinew Task("SchedulerCommon::copyDataToNewGrid",
                                       this,
                                       &SchedulerCommon::copyDataToNewGrid));

      for (auto& [label, mat_subset] : d_label_matls[level]) {
        const VarLabel* var   = label;
        MaterialSubset* matls = mat_subset.get();

        data_tasks.back()->requires(Task::OldDW,
                                    var,
                                    nullptr,
                                    Task::OtherGridDomain,
                                    matls,
                                    Task::NormalDomain,
                                    Ghost::None,
                                    0);
        DOUTR(g_schedulercommon_dbg,
              "Scheduling copy for var "
                << *var << " matl " << *matls
                << " Copies: " << *copy_patch_sets[level].get_rep());
        data_tasks.back()->computes(var, matls);
      }
      addTask(data_tasks.back(),
              copy_patch_sets[level].get_rep(),
              d_materialManager->allMaterials());

      // Monitoring tasks must be scheduled last!!
      scheduleTaskMonitoring(copy_patch_sets[level].get_rep());
    }

    //__________________________________
    //  Scheduling for modifyDataOnNewGrid
    if (refine_patch_sets[level]->size() > 0) {
      data_tasks.push_back(scinew Task("SchedulerCommon::modifyDataOnNewGrid",
                                       this,
                                       &SchedulerCommon::copyDataToNewGrid));

      for (auto& [label, mat_subset] : d_label_matls[level]) {
        const VarLabel* var   = label;
        MaterialSubset* matls = mat_subset.get();

        data_tasks.back()->requires(Task::OldDW,
                                    var,
                                    nullptr,
                                    Task::OtherGridDomain,
                                    matls,
                                    Task::NormalDomain,
                                    Ghost::None,
                                    0);
        DOUTR(g_schedulercommon_dbg,
              "  Scheduling modify for var "
                << *var << " matl " << *matls
                << " Modifies: " << *refine_patch_sets[level].get_rep());
        data_tasks.back()->modifies(var, matls);
      }
      addTask(data_tasks.back(),
              refine_patch_sets[level].get_rep(),
              d_materialManager->allMaterials());

      // Monitoring tasks must be scheduled last!!
      scheduleTaskMonitoring(refine_patch_sets[level].get_rep());
    }

    if (level > 0) {
      d_simulator->scheduleRefineInterface(new_level, sched, 0, 1);
    }
  }

  // set so the load balancer will make an adequate neighborhood, as the default
  // neighborhood isn't good enough for the copy data timestep
  d_is_copy_data_timestep = true;

  this->compile();

  (*d_runtime_stats)[RegriddingCompilationTime] += timer().seconds();

  // save these and restore them, since the next execute will append the
  // scheduler's, and we don't want to.
  double exec_time   = (*d_runtime_stats)[TaskExecTime];
  double local_time  = (*d_runtime_stats)[TaskLocalCommTime];
  double wait_time   = (*d_runtime_stats)[TaskWaitCommTime];
  double reduce_time = (*d_runtime_stats)[TaskReduceCommTime];
  double thread_time = (*d_runtime_stats)[TaskWaitThreadTime];

  timer.reset(true);
  this->execute();

  // copy reduction and sole variables to the new dw
  std::vector<VarLabelMatl<Level>> level_variable_info;
  old_dw->getVarLabelMatlLevelTriples(level_variable_info);

  new_dw->unfinalize();
  for (unsigned int i = 0; i < level_variable_info.size(); i++) {
    VarLabelMatl<Level> current_global_var = level_variable_info[i];

    if (current_global_var.m_label->typeDescription()->isReductionVariable()) {
      const Level* old_level = current_global_var.m_domain;
      const Level* new_level = nullptr;
      if (old_level && old_level->getIndex() < grid->numLevels()) {
        if (old_level->getIndex() >= grid->numLevels()) {
          // the new grid no longer has this level
          continue;
        }
        new_level =
          (new_dw->getGrid()->getLevel(old_level->getIndex())).get_rep();
      }

      // Either both levels need to be null or both need to exist (null levels
      // mean global data)
      if (!old_level || new_level) {
        ReductionVariableBase* v = dynamic_cast<ReductionVariableBase*>(
          current_global_var.m_label->typeDescription()->createInstance());
        old_dw->get(*v,
                    current_global_var.m_label,
                    current_global_var.m_domain,
                    current_global_var.m_matl_index);
        new_dw->put(*v,
                    current_global_var.m_label,
                    new_level,
                    current_global_var.m_matl_index);
        delete v; // copied on the put command
      }
    }
  }

  new_dw->refinalize();

  (*d_runtime_stats)[RegriddingCopyDataTime] += timer().seconds();

  // restore values from before the regrid and data copy
  (*d_runtime_stats)[TaskExecTime]       = exec_time;
  (*d_runtime_stats)[TaskLocalCommTime]  = local_time;
  (*d_runtime_stats)[TaskWaitCommTime]   = wait_time;
  (*d_runtime_stats)[TaskReduceCommTime] = reduce_time;
  (*d_runtime_stats)[TaskWaitThreadTime] = thread_time;

  d_is_copy_data_timestep = false;
}

void
SchedulerCommon::copyDataToNewGrid(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  DOUTR(g_schedulercommon_dbg,
        "SchedulerCommon::copyDataToNewGrid() BGN on patches " << *patches);

  OnDemandDataWarehouse* oldDataWarehouse =
    dynamic_cast<OnDemandDataWarehouse*>(old_dw);
  OnDemandDataWarehouse* newDataWarehouse =
    dynamic_cast<OnDemandDataWarehouse*>(new_dw);

  // For each patch in the patch subset which contains patches in the new grid
  for (int p = 0; p < patches->size(); p++) {
    const Patch* newPatch  = patches->get(p);
    const Level* new_level = newPatch->getLevel();

    // to create once per matl instead of once per matl-var
    std::vector<ParticleSubset*> oldsubsets(
      d_materialManager->getNumMaterials()),
      newsubsets(d_materialManager->getNumMaterials());

    // If there is a level that didn't exist, we don't need to copy it
    if (new_level->getIndex() >= oldDataWarehouse->getGrid()->numLevels()) {
      continue;
    }

    // find old patches associated with this patch
    LevelP old_level =
      oldDataWarehouse->getGrid()->getLevel(new_level->getIndex());

    for (auto& [var_label, material_subset] :
         d_label_matls[old_level->getIndex()]) {
      const VarLabel* label     = var_label;
      MaterialSubset* var_matls = material_subset.get();

      // get the low/high for what we'll need to get
      Patch::VariableBasis basis =
        Patch::translateTypeToBasis(label->typeDescription()->getType(), true);
      IntVector newLowIndex, newHighIndex;
      newPatch->computeVariableExtents(
        basis, IntVector(0, 0, 0), Ghost::None, 0, newLowIndex, newHighIndex);

      //  Loop over materials
      for (int m = 0; m < var_matls->size(); m++) {
        int matl = var_matls->get(m);

        if (!matls->contains(matl)) {
          // std::cout << "We are skipping material " << currentVar.matlIndex_
          // << std::endl;
          continue;
        }

        //  Grid Variables
        switch (label->typeDescription()->getType()) {
          case TypeDescription::Type::PerPatch:
          case TypeDescription::Type::NCVariable:
          case TypeDescription::Type::CCVariable:
          case TypeDescription::Type::SFCXVariable:
          case TypeDescription::Type::SFCYVariable:
          case TypeDescription::Type::SFCZVariable: {
            Patch::selectType oldPatches;
            old_level->selectPatches(newLowIndex, newHighIndex, oldPatches);

            for (size_t oldIdx = 0; oldIdx < oldPatches.size(); oldIdx++) {
              const Patch* oldPatch = oldPatches[oldIdx];

              if (!oldDataWarehouse->exists(label, matl, oldPatch)) {
                continue; // see comment about oldPatchToTest in
                          // ScheduleAndDoDataCopy
              }

              IntVector oldLowIndex;
              IntVector oldHighIndex;

              if (new_level->getIndex() > 0) {
                oldLowIndex  = oldPatch->getLowIndexWithDomainLayer(basis);
                oldHighIndex = oldPatch->getHighIndexWithDomainLayer(basis);
              } else {
                oldLowIndex =
                  oldPatch->getExtraLowIndex(basis, label->getBoundaryLayer());
                oldHighIndex =
                  oldPatch->getExtraHighIndex(basis, label->getBoundaryLayer());
              }

              IntVector copyLowIndex  = Max(newLowIndex, oldLowIndex);
              IntVector copyHighIndex = Min(newHighIndex, oldHighIndex);

              // based on the selectPatches above, we might have patches we
              // don't want to use, so prune them here.
              if (copyLowIndex.x() >= copyHighIndex.x() ||
                  copyLowIndex.y() >= copyHighIndex.y() ||
                  copyLowIndex.z() >= copyHighIndex.z()) {
                continue;
              }

              if (!oldDataWarehouse->exists(label, matl, oldPatch)) {
                SCI_THROW(UnknownVariable(label->getName(),
                                          oldDataWarehouse->getID(),
                                          oldPatch,
                                          matl,
                                          "in copyDataTo GridVariableBase",
                                          __FILE__,
                                          __LINE__));
              }

              if (label->typeDescription()->getType() ==
                  TypeDescription::Type::PerPatch) {
                std::vector<Variable*> varlist;
                oldDataWarehouse->d_var_DB.getlist(
                  label, matl, oldPatch, varlist);
                PerPatchBase* v = nullptr;

                for (auto& var : varlist) {
                  v = dynamic_cast<PerPatchBase*>(var);
                  ASSERT(v->getBasePointer() != nullptr);

                  if (!newDataWarehouse->exists(label, matl, newPatch)) {
                    PerPatchBase* newVariable = v->clone();
                    newDataWarehouse->d_var_DB.put(label,
                                                   matl,
                                                   newPatch,
                                                   newVariable,
                                                   copyTimestep(),
                                                   false);

                  } else {
                    PerPatchBase* newVariable = dynamic_cast<PerPatchBase*>(
                      newDataWarehouse->d_var_DB.get(label, matl, newPatch));

                    if (oldPatch->isVirtual()) {
                      // it can happen where the old patch was virtual and this
                      // is not
                      PerPatchBase* tmpVar = newVariable->clone();
                      tmpVar->copyPointer(*v);
                      newVariable = tmpVar;
                      delete tmpVar;
                    } else {
                      newVariable = v;
                    }
                  }
                }
              } else {
                std::vector<Variable*> varlist;
                oldDataWarehouse->d_var_DB.getlist(
                  label, matl, oldPatch, varlist);
                GridVariableBase* v = nullptr;

                IntVector srclow  = copyLowIndex;
                IntVector srchigh = copyHighIndex;

                for (auto& var : varlist) {
                  v = dynamic_cast<GridVariableBase*>(var);

                  ASSERT(v->getBasePointer() != nullptr);

                  // restrict copy to data range
                  srclow  = Max(copyLowIndex, v->getLow());
                  srchigh = Min(copyHighIndex, v->getHigh());
                  if (srclow.x() >= srchigh.x() || srclow.y() >= srchigh.y() ||
                      srclow.z() >= srchigh.z()) {
                    continue;
                  }

                  if (!newDataWarehouse->exists(label, matl, newPatch)) {
                    GridVariableBase* newVariable = v->cloneType();
                    newVariable->rewindow(newLowIndex, newHighIndex);
                    newVariable->copyPatch(v, srclow, srchigh);
                    newDataWarehouse->d_var_DB.put(label,
                                                   matl,
                                                   newPatch,
                                                   newVariable,
                                                   isCopyDataTimestep(),
                                                   false);
                  } else {
                    GridVariableBase* newVariable =
                      dynamic_cast<GridVariableBase*>(
                        newDataWarehouse->d_var_DB.get(label, matl, newPatch));

                    // make sure it exists in the right region (it might be
                    // ghost data)
                    newVariable->rewindow(newLowIndex, newHighIndex);

                    if (oldPatch->isVirtual()) {
                      // it can happen where the old patch was virtual and
                      // this is not
                      GridVariableBase* tmpVar = newVariable->cloneType();
                      tmpVar->copyPointer(*v);
                      tmpVar->offset(oldPatch->getVirtualOffset());
                      newVariable->copyPatch(tmpVar, srclow, srchigh);
                      delete tmpVar;
                    } else {
                      newVariable->copyPatch(v, srclow, srchigh);
                    }
                  }
                }
              }
            } // end oldPatches
          } break;

          //  Particle Variables
          case TypeDescription::Type::ParticleVariable: {
            ParticleSubset* oldsub = oldsubsets[matl];
            if (!oldsub) {
              // collect the particles from the range encompassing this patch.
              // Use interior cells since extracells aren't collected across
              // processors in the data copy, and they don't matter for
              // particles anyhow (but we will have to reset the bounds to
              // copy the data)
              oldsub = oldDataWarehouse->getParticleSubset(
                matl,
                newPatch->getLowIndexWithDomainLayer(Patch::CellBased),
                newPatch->getHighIndexWithDomainLayer(Patch::CellBased),
                newPatch,
                d_reloc_new_pos_label,
                old_level.get_rep());
              oldsubsets[matl] = oldsub;
              oldsub->addReference();
            }

            ParticleSubset* newsub = newsubsets[matl];
            // it might have been created in Refine
            if (!newsub) {
              if (!newDataWarehouse->haveParticleSubset(matl, newPatch)) {
                newsub = newDataWarehouse->createParticleSubset(
                  oldsub->numParticles(), matl, newPatch);
              } else {
                newsub = newDataWarehouse->getParticleSubset(matl, newPatch);
                ASSERT(newsub->numParticles() == 0);
                newsub->addParticles(oldsub->numParticles());
              }
              newsubsets[matl] = newsub;
            }

            ParticleVariableBase* newv = dynamic_cast<ParticleVariableBase*>(
              label->typeDescription()->createInstance());
            newv->allocate(newsub);
            // don't get and copy if there were no old patches
            if (oldsub->getNeighbors().size() > 0) {
              constParticleVariableBase* var = newv->cloneConstType();
              oldDataWarehouse->get(*var, label, oldsub);

              // reset the bounds of the old var's data so copyData doesn't
              // complain
              ParticleSubset* tempset =
                scinew ParticleSubset(oldsub->numParticles(),
                                      matl,
                                      newPatch,
                                      newPatch->getExtraCellLowIndex(),
                                      newPatch->getExtraCellHighIndex());
              const_cast<ParticleVariableBase*>(&var->getBaseRep())
                ->setParticleSubset(tempset);
              newv->copyData(&var->getBaseRep());
              delete var; // pset and tempset are deleted with it.
            }
            newDataWarehouse->put(*newv, label, true);
            delete newv; // the container is copied
          } break;

          default: {
            SCI_THROW(InternalError("Unknown variable type in copyData: " +
                                      label->getName(),
                                    __FILE__,
                                    __LINE__));
          }
        } // end switch
      }   // end matls
    }     // end label_matls

    for (unsigned i = 0; i < oldsubsets.size(); i++) {
      if (oldsubsets[i] && oldsubsets[i]->removeReference()) {
        delete oldsubsets[i];
      }
    }
  } // end patches

  DOUTR(g_schedulercommon_dbg, "SchedulerCommon::copyDataToNewGrid() END");
}

//______________________________________________________________________
//
void
SchedulerCommon::scheduleParticleRelocation(
  const LevelP& coarsest_level_with_particles,
  const VarLabel* old_pos_label,
  const std::vector<std::vector<const VarLabel*>>& old_other_labels,
  const VarLabel* new_pos_label,
  const std::vector<std::vector<const VarLabel*>>& new_other_labels,
  const VarLabel* particle_id_label,
  const MaterialSet* matls)
{
  if (d_reloc_new_pos_label) {
    ASSERTEQ(d_reloc_new_pos_label, new_pos_label);
  }
  d_reloc_new_pos_label = new_pos_label;

  d_relocate_1.scheduleParticleRelocation(this,
                                          d_myworld,
                                          d_load_balancer,
                                          coarsest_level_with_particles,
                                          old_pos_label,
                                          old_other_labels,
                                          new_pos_label,
                                          new_other_labels,
                                          particle_id_label,
                                          matls);
}

void
SchedulerCommon::scheduleParticleRelocation(
  const LevelP& coarsestLevelwithParticles,
  const VarLabel* pos_label,
  const std::vector<std::vector<const VarLabel*>>& other_labels,
  const MaterialSet* matls)
{
  d_reloc_new_pos_label = pos_label;
  d_relocate_1.scheduleParticleRelocation(this,
                                          d_myworld,
                                          d_load_balancer,
                                          coarsestLevelwithParticles,
                                          pos_label,
                                          other_labels,
                                          matls);
}

void
SchedulerCommon::overrideVariableBehavior(const std::string& var,
                                          bool treatAsOld,
                                          bool copyData,
                                          bool noScrub,
                                          bool notCopyData,
                                          bool noCheckpoint)
{
  // treat variable as an "old" var - will be checkpointed, copied, and only
  // scrubbed from an OldDW
  if (treatAsOld) {
    d_treat_as_old_vars.insert(var);
  }

  // manually copy this variable to the new_dw if regridding occurs
  if (copyData) {
    d_copy_data_vars.insert(var);
    d_no_scrub_vars.insert(var);
  }

  // set variable not to scrub (normally when needed between a normal
  // taskgraph and the regridding phase)
  if (noScrub) {
    d_no_scrub_vars.insert(var);
  }

  // ignore copying this variable between AMR levels
  if (notCopyData) {
    d_not_copy_data_vars.insert(var);
  }

  // do not checkpoint this variable.
  if (noCheckpoint) {
    d_not_checkpoint_vars.insert(var);
  }
}

void
SchedulerCommon::clearTaskMonitoring()
{
  // Loop through the global (0) and local (1) tasks
  for (unsigned int i = 0; i < 2; ++i) {
    d_monitoring_values[i].clear();
  }
}

// Schedule the recording of the task monitoring attribute
// values. This task should be the last task so that the writing is
// done after all task have been executed.
void
SchedulerCommon::scheduleTaskMonitoring(const LevelP& level)
{
  if (!d_monitoring) {
    return;
  }

  // Create and schedule a task that will record each of the
  // tasking monitoring attributes.
  Task* t = scinew Task("SchedulerCommon::recordTaskMonitoring",
                        this,
                        &SchedulerCommon::recordTaskMonitoring);

  // Ghost::GhostType gn = Ghost::None;

  for (unsigned int i = 0; i < 2; ++i) {
    for (const auto& [task_name, var_label] : d_monitoring_tasks[i]) {
      t->computes(var_label, d_dummy_matl.get(), Task::OutOfDomain);

      // treatAsOld copyData noScrub notCopyData noCheckpoint
      overrideVariableBehavior(
        var_label->getName(), false, false, true, true, true);
    }
  }

  addTask(t, level->eachPatch(), d_materialManager->allMaterials());
}

// Schedule the recording of the task monitoring attribute
// values. This task should be the last task so that the writing is
// done after all task have been executed.
void
SchedulerCommon::scheduleTaskMonitoring(const PatchSet* patches)
{
  if (!d_monitoring) {
    return;
  }

  // Create and schedule a task that will record each of the
  // tasking monitoring attributes.
  Task* t = scinew Task("SchedulerCommon::recordTaskMonitoring",
                        this,
                        &SchedulerCommon::recordTaskMonitoring);

  // Ghost::GhostType gn = Ghost::None;

  for (unsigned int i = 0; i < 2; ++i) {
    for (const auto& it : d_monitoring_tasks[i]) {
      t->computes(it.second, d_dummy_matl.get(), Task::OutOfDomain);

      overrideVariableBehavior(
        it.second->getName(), false, false, true, true, true);
      // treatAsOld copyData noScrub notCopyData noCheckpoint
    }
  }

  addTask(t, patches, d_materialManager->allMaterials());
}

// Record the global task monitoring attribute values into the data
// warehouse.
void
SchedulerCommon::recordTaskMonitoring(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* /*matls*/,
                                      [[maybe_unused]] DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  int matlIndex = 0;

  // For all of the patches record the tasking monitoring attribute value.
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    // Loop through the global (0) and local (1) tasks
    for (unsigned int i = 0; i < 2; ++i) {
      for (const auto& it : d_monitoring_tasks[i]) {
        PerPatch<double> value =
          d_monitoring_values[i][it.first][patch->getID()];

        new_dw->put(value, it.second, matlIndex, patch);
      }
    }
  }
}

// Sum the task monitoring attribute values
void
SchedulerCommon::sumTaskMonitoringValues(DetailedTask* dtask)
{
  if (!d_monitoring) {
    return;
  }

  const PatchSubset* patches = dtask->getPatches();

  if (patches && patches->size()) {
    // Compute the cost on a per cell basis so the measured value can
    // be distributed proportionally by cells
    double num_cells = 0;

    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      num_cells += patch->getNumExtraCells();
    }

    double weight;

    // Compute the value on a per cell basis.
    if (d_monitoring_per_cell) {
      weight = num_cells;
    }
    // Compute the value on a per patch basis.
    else {
      weight = 1;
    }
    // Loop through the global (0) and local (1) tasks
    for (auto i = 0; i < 2; ++i) {
      for (const auto& it : d_monitoring_tasks[i]) {
        // Strip off the attribute name from the task name.
        std::string taskName = it.first;

        // For a local task strip off the attribute name.
        if (i == 1) {
          size_t found = taskName.find_last_of("::");
          // std::string attribute = taskName.substr(found + 1);
          taskName = taskName.substr(0, found - 1);
        }

        // Is this task being monitored ?
        if ((i == 0) || // Global monitoring yes, otherwise check.
            (i == 1 && taskName == dtask->getTask()->getName())) {
          bool loadBalancerCost = false;
          double value;

          // Currently the monitoring is limited to the LoadBalancer cost, task
          // exec time, and task wait time.
          if (it.first.find("LoadBalancerCost") != std::string::npos) {
            // The same code is in runTask of the specific scheduler
            // (MPIScheduler and UnifiedScheduler) to use the task
            // execution time which is then weighted by the number of
            // cells in CostModelForecaster::addContribution
            if (!dtask->getTask()->getHasSubScheduler() &&
                !d_is_copy_data_timestep &&
                dtask->getTask()->getType() != Task::Output) {
              value            = dtask->task_exec_time() / num_cells;
              loadBalancerCost = true;
            } else {
              value = 0.0;
            }
          } else if (it.first.find("ExecTime") != std::string::npos) {
            value = dtask->task_exec_time() / weight;
          } else if (it.first.find("WaitTime") != std::string::npos) {
            value = dtask->task_wait_time() / weight;
          } else {
            value = 0.0;
          }

          if (value != 0.0) {
            // Loop through patches and add the contribution.
            for (int p = 0; p < patches->size(); ++p) {
              const Patch* patch = patches->get(p);

              if (d_monitoring_per_cell || loadBalancerCost) {
                d_monitoring_values[i][it.first][patch->getID()] +=
                  patch->getNumExtraCells() * value;
              } else {
                d_monitoring_values[i][it.first][patch->getID()] += value;
              }

              // A cheat ... the only time this task will come here is
              // after the value has been written (and the task is
              // completed) so the value can be overwritten. This
              // allows the monitoring to be monitored.
              if (dtask->getTask()->getName() ==
                  "SchedulerCommon::recordTaskMonitoring") {
                PerPatch<double> value =
                  d_monitoring_values[i][it.first][patch->getID()];
                this->getLastDW()->put(value, it.second, 0, patch);
              }
            }
          }
        }
      }
    }
  }
}

} // namespace Uintah