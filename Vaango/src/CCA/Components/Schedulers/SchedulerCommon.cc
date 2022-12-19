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

#include <CCA/Components/Schedulers/SchedulerCommon.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Components/Schedulers/TaskGraph.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/LocallyComputedPatchVarMap.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Malloc/Allocator.h>
#include <Core/OS/ProcessInfo.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Thread/Time.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cerrno>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>

#include <time.h>

using namespace Uintah;

// Debug: Used to sync cerr so it is readable (when output by
// multiple threads at the same time)  From sus.cc:
extern Uintah::Mutex cerrLock;
extern DebugStream mixedDebug;

static DebugStream dbg("SchedulerCommon", false);

// for calculating memory usage when sci-malloc is disabled.
char * SchedulerCommon::d_start_addr = NULL;


SchedulerCommon::SchedulerCommon(const ProcessorGroup* myworld, const Output* oport)
  : UintahParallelComponent(myworld), m_outPort(oport),
    d_trackingVarsPrintLocation(0), d_maxMemUse(0), m_graphDoc(NULL), m_nodes(NULL)
{
  d_generation = 0;
  d_numOldDWs = 0;

  d_emit_taskgraph = false;
  d_useSmallMessages = true;
  d_memlogfile = 0;
  d_restartable = false;
  for(int i=0;i<Task::TotalDWs;i++)
    d_dwmap[i]=Task::InvalidDW;
  // Default mapping...
  d_dwmap[Task::OldDW]=0;
  d_dwmap[Task::NewDW]=1;

  d_isInitTimestep = false;
  d_isRestartInitTimestep = false;

  m_locallyComputedPatchVarMap = scinew LocallyComputedPatchVarMap;
  d_reloc_new_posLabel = 0;
  d_maxGhost=0;
  d_maxLevelOffset=0;
}

//______________________________________________________________________
//

SchedulerCommon::~SchedulerCommon()
{
  if(d_memlogfile)
    delete d_memlogfile;

  // list of vars used for AMR regridding
  for (auto label_matl_maps : d_label_matls) {
    for (auto& map : label_matl_maps) {
      if (map.second->removeReference()) {
        delete map.second;
      }
    }
  }

  for (unsigned i = 0; i < d_graphs.size(); i++) {
    delete d_graphs[i];
  }

  d_label_matls.clear();

  delete m_locallyComputedPatchVarMap;
}

//______________________________________________________________________
//

void
SchedulerCommon::checkMemoryUse( unsigned long & memuse, 
                                 unsigned long & highwater,
                                 unsigned long & maxMemUse )
{
  highwater = 0; 
  memuse = 0;

#if !defined(DISABLE_SCI_MALLOC)
  size_t nalloc,  sizealloc, nfree,  sizefree, nfillbin,
    nmmap, sizemmap, nmunmap, sizemunmap, highwater_alloc,  
    highwater_mmap, bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse, bytes_inhunks;
  
  GetGlobalStats( DefaultAllocator(),
                  nalloc, sizealloc, nfree, sizefree,
                  nfillbin, nmmap, sizemmap, nmunmap,
                  sizemunmap, highwater_alloc, highwater_mmap,
                  bytes_overhead, bytes_free, bytes_fragmented, bytes_inuse, bytes_inhunks );
  memuse = sizealloc - sizefree;
  highwater = highwater_mmap;

#elif !defined(_WIN32)

  if ( ProcessInfo::IsSupported( ProcessInfo::MEM_SIZE ) ) {
    memuse = ProcessInfo::GetMemoryResident();
    // printf("1) memuse is %d\n", (int)memuse);
  } else {
    memuse = (char*)sbrk(0)-d_start_addr;
    // printf("2) memuse is %d\n", (int)memuse);
  }
#endif

  if( memuse > d_maxMemUse ) {
    // printf("Max memuse increased\n");
    d_maxMemUse = memuse;
  }
  maxMemUse = d_maxMemUse;
}

void
SchedulerCommon::resetMaxMemValue()
{
  d_maxMemUse = 0;
}

//______________________________________________________________________
//

void
SchedulerCommon::makeTaskGraphDoc(const DetailedTasks*/* dt*/, int rank)
{
  if (!d_emit_taskgraph)
    return;

  if (!m_outPort->isOutputTimestep())
    return;
  
  // make sure to release this DOMDocument after finishing emitting the nodes
  m_graphDoc = ProblemSpec::createDocument("Uintah_TaskGraph");
  
  ProblemSpecP meta = m_graphDoc->appendChild("Meta");
  meta->appendElement("username", getenv("LOGNAME"));
  time_t t = time(NULL);
  meta->appendElement("date", ctime(&t));
  
  m_nodes = m_graphDoc->appendChild("Nodes");
  //m_graphDoc->appendChild(m_nodes);
  
  ProblemSpecP edgesElement = m_graphDoc->appendChild("Edges");
  
  for (unsigned i = 0; i < d_graphs.size(); i++) {
    DetailedTasks* dts = d_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->emitEdges(edgesElement, rank);
    }
  }
}

//______________________________________________________________________
//

bool
SchedulerCommon::useInternalDeps()
{
  // keep track of internal dependencies only if it will emit
  // the taskd_graphs (by default).
  return d_emit_taskgraph;
}

//______________________________________________________________________
//

void
SchedulerCommon::emitNode( const DetailedTask* task, 
                           double        start,
                           double        duration,
                           double        execution_duration)
{  
  if (m_nodes == 0)
    return;
   
  ProblemSpecP node = m_nodes->appendChild("node");
  //m_nodes->appendChild(node);

  node->appendElement("name", task->getName());
  node->appendElement("start", start);
  node->appendElement("duration", duration);
  if (execution_duration > 0)
    node->appendElement("execution_duration", execution_duration);
}

//______________________________________________________________________
//

void
SchedulerCommon::finalizeNodes(int process /* = 0*/)
{
  if (m_graphDoc == 0)
    return;

  if (m_outPort->isOutputTimestep()) {
    std::string timestep_dir(m_outPort->getLastTimestepOutputLocation());
      
    std::ostringstream fname;
    fname << "/taskgraph_" << std::setw(5) << std::setfill('0') << process << ".xml";
    std::string file_name(timestep_dir + fname.str());
    m_graphDoc->output(file_name.c_str());
  }
    
  //m_graphDoc->releaseDocument();
  //m_graphDoc = NULL;
  //m_nodes = NULL;
}

//______________________________________________________________________
//

void
SchedulerCommon::problemSetup(const ProblemSpecP& prob_spec,
                              SimulationStateP& state)
{
  d_sharedState = state;

  // Initializing d_trackingStartTime and d_trackingEndTime to default values
  // so that we do not crash when running MALLOC_STRICT.
  d_trackingStartTime = 1;
  d_trackingEndTime = 0;
  d_trackingVarsPrintLocation = PRINT_AFTER_EXEC;
  ProblemSpecP params = prob_spec->findBlock("Scheduler");
  if(params){
    params->getWithDefault("small_messages", d_useSmallMessages, true);
    if (d_useSmallMessages)
      proc0cout << "   Using theoretical scheduler\n";
    ProblemSpecP track = params->findBlock("VarTracker");
    if (track) {
      track->require("start_time", d_trackingStartTime);
      track->require("end_time", d_trackingEndTime);
      track->getWithDefault("level", d_trackingLevel, -1);
      track->getWithDefault("start_index", d_trackingStartIndex, IntVector(-9,-9,-9));
      track->getWithDefault("end_index", d_trackingEndIndex, IntVector(-9,-9,-9));
      track->getWithDefault("patchid", d_trackingPatchID, -1);

      if( d_myworld->myrank() == 0 ) {
        std::cout << "\n";
        std::cout << "-----------------------------------------------------------\n";
        std::cout << "-- Initializing VarTracker...\n";
        std::cout << "--  Running from time " << d_trackingStartTime  << " to " << d_trackingEndTime << "\n";
        std::cout << "--  for indices: " << d_trackingStartIndex << " to " << d_trackingEndIndex << "\n";
      }

      ProblemSpecP location = track->findBlock("locations");
      if (location) {
        d_trackingVarsPrintLocation = 0;
        std::map<std::string,std::string> attributes;
        location->getAttributes(attributes);
        if (attributes["before_comm"] == "true") {
          d_trackingVarsPrintLocation |= PRINT_BEFORE_COMM;
          proc0cout << "--  Printing variable information before communication.\n"; 
        }
        if (attributes["before_exec"] == "true") {
          d_trackingVarsPrintLocation |= PRINT_BEFORE_EXEC;
          proc0cout << "--  Printing variable information before task execution.\n"; 
        }
        if (attributes["after_exec"] == "true") {
          d_trackingVarsPrintLocation |= PRINT_AFTER_EXEC;
          proc0cout << "--  Printing variable information after task execution.\n"; 
        }
      }
      else {
        // "locations" not specified
        proc0cout << "--  Defaulting to printing variable information after task execution.\n"; 
      }

      for (ProblemSpecP var=track->findBlock("var"); var != 0; var = var->findNextBlock("var")) {
        std::map<std::string,std::string> attributes;
        var->getAttributes(attributes);
        std::string name = attributes["label"];
        d_trackingVars.push_back(name);
        std::string dw = attributes["dw"];

        if (dw == "OldDW") {
          d_trackingDWs.push_back(Task::OldDW);
        }
        else if (dw == "NewDW") {
          d_trackingDWs.push_back(Task::NewDW);
        }
        else if (dw == "CoarseNewDW") {
          d_trackingDWs.push_back(Task::CoarseNewDW);
        }
        else if (dw == "CoarseOldDW") {
          d_trackingDWs.push_back(Task::CoarseOldDW);
        }
        else if (dw == "ParentOldDW") {
          d_trackingDWs.push_back(Task::ParentOldDW);
        }
        else if (dw == "ParentOldDW") {
          d_trackingDWs.push_back(Task::ParentNewDW);
        }
        else {
          // This error message most likely can go away once the .ups validation is put into place:
          printf( "WARNING: Hit switch statement default... using NewDW... (This could possibly be"
                  "an error in input file specification.)\n" );
          d_trackingDWs.push_back(Task::NewDW);
        }
        proc0cout << "--  Tracking variable '" << name << "' in DataWarehouse '" << dw << "'\n";
      }

      for (ProblemSpecP task=track->findBlock("task"); task != 0; task = task->findNextBlock("task")) {
        std::map<std::string,std::string> attributes;
        task->getAttributes(attributes);
        std::string name = attributes["name"];
        d_trackingTasks.push_back(name);
        proc0cout << "--  Tracking variables for specific task: " << name << "\n"; 
      }      
      proc0cout << "-----------------------------------------------------------\n\n";
    }
    else { // Tracking not specified
      // This 'else' won't be necessary once the .ups files are validated... but for now.
      proc0cout << "<VarTracker> not specified in .ups file... no variable tracking will take place.\n";
    }
  }

  // If small_messages not specified in UP Scheduler block, still report what's used
  if( d_useSmallMessages ) {
    proc0cout << "   Using small, individual MPI messages (no message combining)\n";
  }
  else {
    proc0cout << "   Using large, combined MPI messages\n";
  }

  d_noScrubVars.insert("refineFlag");
  d_noScrubVars.insert("refinePatchFlag");
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
handleError( int errorPosition, const std::string & errorMessage, const std::string & variableName ) 
{
  static std::vector< std::map<std::string,bool> * > errorsReported( 8 );

  std::map<std::string, bool> * varToReportedMap = errorsReported[ errorPosition ];

  if( varToReportedMap == NULL ) {
    varToReportedMap = new std::map<std::string, bool>;
    errorsReported[ errorPosition ] = varToReportedMap;
  }

  bool reported = (*varToReportedMap)[ variableName ];
  if( !reported ) {
    (*varToReportedMap)[ variableName ] = true;
    std::cout << errorMessage << "\n";
    return true;
  }
  return false;
}

//______________________________________________________________________
//
void
SchedulerCommon::printTrackedVars( DetailedTask* dt, int when )
{
  bool printedHeader = false;

  LoadBalancer* lb = getLoadBalancer();
 
  unsigned taskNum;
  for (taskNum = 0; taskNum < d_trackingTasks.size(); taskNum++) {
    if (d_trackingTasks[taskNum] == dt->getTask()->getName())
      break;
  }

  // Print for all tasks unless one is specified (but disclude DataArchiver tasks)
  if ((taskNum == d_trackingTasks.size() && d_trackingTasks.size() != 0) || 
      ((std::string(dt->getTask()->getName())).substr(0,12) == "DataArchiver"))
    return;

  if( d_sharedState && ( d_trackingStartTime > d_sharedState->getElapsedTime() ||
                         d_trackingEndTime < d_sharedState->getElapsedTime() ) ) {
    return;
  }

  for (int i = 0; i < static_cast<int>(d_trackingVars.size()); i++) {
    bool printedVarName = false;

    // that DW may not have been mapped....
    if (dt->getTask()->mapDataWarehouse(d_trackingDWs[i]) < 0 || 
        dt->getTask()->mapDataWarehouse(d_trackingDWs[i]) >= static_cast<int>(d_dws.size())) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i]
           << ") DW is out of range.\n";
      handleError( 0, mesg.str(), d_trackingVars[i] );

      continue;
    }

    OnDemandDataWarehouseP dw = d_dws[dt->getTask()->mapDataWarehouse(d_trackingDWs[i])];

    if (dw == 0) { // old on initialization timestep
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i] 
           << ") because DW is NULL.  Requested DW was: " 
           << dt->getTask()->mapDataWarehouse(d_trackingDWs[i]) << "\n";
      handleError( 1, mesg.str(), d_trackingVars[i] );
      continue;
    }

    // get the level here, as the grid can be different between the old and new DW
    
    const Grid* grid = dw->getGrid();

    int levelnum;
    
    if (d_trackingLevel == -1) {
      levelnum = grid->numLevels() - 1;
    }
    else {
      levelnum = d_trackingLevel;
      if (levelnum >= grid->numLevels())
        continue;
    }

    const LevelP level = grid->getLevel(levelnum);
    const VarLabel* label = VarLabel::find(d_trackingVars[i]);

    std::cout.precision(16);

    if (!label) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i]
           << ") because label is NULL.\n";
      handleError( 2, mesg.str(), d_trackingVars[i] );
      continue;
    }

    const PatchSubset* patches = dt->getPatches();
    
    // a once-per-proc task is liable to have multiple levels, and thus calls to getLevel(patches) will fail
    if( dt->getTask()->getType() != Task::OncePerProc && (!patches || getLevel(patches)->getIndex() != levelnum) ) {
      std::ostringstream mesg;
      mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i]
           << ") because patch is non-standard.\n";
      handleError( 3, mesg.str(), d_trackingVars[i] );
      continue;
    }

    for (int p = 0; patches && p < patches->size(); p++) {

      const Patch* patch = patches->get(p);
      if (d_trackingPatchID != -1 && d_trackingPatchID != patch->getID()) {
        std::ostringstream mesg;
        mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i]
             << ") because patch ID does not match.\n" 
             << "            (Error only printed once.)\n"
             << "         Tracking Patch ID: " << d_trackingPatchID << ", patch id: " << patch->getID() << "\n";
        handleError( 4, mesg.str(), d_trackingVars[i] );
        continue;
      }

      // don't print ghost patches (dw->get will yell at you)
      if ((d_trackingDWs[i] == Task::OldDW && lb->getOldProcessorAssignment(patch) != d_myworld->myrank()) ||
          (d_trackingDWs[i] == Task::NewDW && lb->getPatchwiseProcessorAssignment(patch) != d_myworld->myrank())) {
        continue;
      }

      const TypeDescription* td = label->typeDescription();
      Patch::VariableBasis basis = patch->translateTypeToBasis(td->getType(), false);
      IntVector start = 
        Max(patch->getExtraLowIndex(basis, IntVector(0,0,0)), d_trackingStartIndex);
      IntVector end = 
        Min(patch->getExtraHighIndex(basis, IntVector(0,0,0)), d_trackingEndIndex);

      // loop over matls too
      for (int m = 0; m < d_sharedState->getNumMatls(); m++) {

        if (!dw->exists(label, m, patch)) {
          std::ostringstream mesg;
          mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i] 
               << ") because it does not exist in DW.\n"
               << "            Patch is: " << *patch << "\n";
          if( handleError( 5, mesg.str(), d_trackingVars[i] ) ) {
            std::cout << "         DW contains (material: " << m << ")\n";
            dw->print();
          }
          continue;
        }
        if (!(start.x() < end.x() && start.y() < end.y() && start.z() < end.z())) {
          std::ostringstream mesg;
          mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i] 
               << ") because the start is greater than the end location:\n"
               << "start: " << start << "\n"
               << "end: " << start << "\n";
          handleError( 6, mesg.str(), d_trackingVars[i] );
          continue;
        }
        if (td->getSubType()->getType() != TypeDescription::double_type &&
            td->getSubType()->getType() != TypeDescription::Vector) {
          // only allow *Variable<double> and *Variable<Vector> for now
          std::ostringstream mesg;
          mesg << "WARNING: VarTracker: Not printing requested variable (" << d_trackingVars[i]
               << ") because its type is not supported:\n"
               << "             " << td->getName() << "\n";
          handleError( 7, mesg.str(), d_trackingVars[i] );
          continue;
        }

        // pending the task that allocates the var, we may not have allocated it yet
        GridVariableBase* v;
        switch (td->getType()) {
        case TypeDescription::CCVariable:
        case TypeDescription::NCVariable:
        case TypeDescription::SFCXVariable:
        case TypeDescription::SFCYVariable:
        case TypeDescription::SFCZVariable: 
          v = dynamic_cast<GridVariableBase*>(dw->d_varDB.get(label, m, patch));
          break;
        default: 
          throw InternalError("Cannot track var type of non-grid-type",__FILE__,__LINE__); break;
        }
        
        start = Max(start, v->getLow());
        end = Min(end, v->getHigh());
        if (!(start.x() < end.x() && start.y() < end.y() && start.z() < end.z())) 
          continue;

        if (!printedHeader) {
          std::string location;
          switch (when) {
          case PRINT_BEFORE_COMM: location = " before communication of "; break;
          case PRINT_BEFORE_EXEC: location = " before execution of "; break;
          case PRINT_AFTER_EXEC: location = " after execution of "; break;
          }
          std::cout << d_myworld->myrank() << location << *dt << endl;
          printedHeader = true;
        }
        
        if (!printedVarName) {
          std::cout << d_myworld->myrank() << "  Variable: " << d_trackingVars[i] << ", DW " << dw->getID() << ", Patch " << patch->getID() << ", Matl " << m << endl;
          if (d_trackingVars[i] == "rho_CC")
            std::cout << "  RHO: " << dw->getID() << " original input " << d_trackingDWs[i] << endl;
        }
            

        switch (td->getSubType()->getType()) {
        case TypeDescription::double_type: 
        {
          GridVariable<double>* var = dynamic_cast<GridVariable<double>*>(v);
          
          for (int z = start.z(); z < end.z(); z++) {
            for (int y = start.y(); y < end.y(); y++) {
              std::cout << d_myworld->myrank() << "  ";
              for (int x = start.x(); x < end.x(); x++) {
                IntVector c(x,y,z);
                std::cout << " " << c << ": " << (*var)[c];
              }
              std::cout << endl;
            }
            std::cout << endl;
          }
        }
        break;
        case TypeDescription::Vector: 
        {
          GridVariable<Vector>* var = dynamic_cast<GridVariable<Vector>*>(v);
          
          for (int z = start.z(); z < end.z(); z++) {
            for (int y = start.y(); y < end.y(); y++) {
              std::cout << d_myworld->myrank() << "  ";
              for (int x = start.x(); x < end.x(); x++) {
                IntVector c(x,y,z);
                std::cout << " " << c << ": " << (*var)[c];
              }
              std::cout << endl;
            }
            std::cout << endl;
          }
        }
        break;
        default: break;
        } // end case variable type
      } // end for m(aterials)
    } // end for p(atches)
  } // end for i : trackingVars.size()
} // end printTrackedVars()

//______________________________________________________________________
//
LoadBalancer*
SchedulerCommon::getLoadBalancer()
{
  UintahParallelPort* lbp = getPort("load balancer");
  LoadBalancer* lb = dynamic_cast<LoadBalancer*>(lbp);
  return lb;
}

//______________________________________________________________________
//
void
SchedulerCommon::addTaskGraph(Scheduler::tgType type)
{
  MALLOC_TRACE_TAG_SCOPE("SchedulerCommon::addTaskGraph");
  TaskGraph* tg = scinew TaskGraph(this, d_myworld, type);
  tg->initialize();
  d_graphs.push_back(tg);
}

//______________________________________________________________________
//
void
SchedulerCommon::addTask(Task* task, const PatchSet* patches,
                         const MaterialSet* matls)
{
  MALLOC_TRACE_TAG_SCOPE("SchedulerCommon::addTask");

  // Save the DW map
  task->setMapping(d_dwmap);
  dbg << "adding Task: " << task->getName() << ", # patches: ";
  if( patches ) dbg << patches->size();
  else          dbg << "0";
  dbg << ", # matls: " ;
  if( matls ) dbg << matls->size();
  else          dbg << "0";
  dbg << "\n";

  d_graphs[d_graphs.size()-1]->addTask(task, patches, matls);
  d_numTasks++;

  if (task->maxGhostCells > d_maxGhost) d_maxGhost = task->maxGhostCells;
  if (task->maxLevelOffset > d_maxLevelOffset) d_maxLevelOffset = task->maxLevelOffset;

  // add to init-requires.  These are the vars which require from the OldDW that we'll
  // need for checkpointing, switching, and the like.
  // In the case of treatAsOld Vars, we handle them because something external to the taskgraph
  // needs it that way (i.e., Regridding on a restart requires checkpointed refineFlags).
  for( const Task::Dependency* dep = task->getRequires(); dep != 0; dep = dep->next ) {
    if(isOldDW(dep->mapDataWarehouse()) || d_treatAsOldVars.find(dep->var->getName()) != d_treatAsOldVars.end()) {
      d_initRequires.push_back(dep);
      d_initRequiredVars.insert(dep->var);
    }
  }

  // for the treat-as-old vars, go through the computes and add them.
  // we can (probably) safely assume that we'll avoid duplicates, since if they were inserted 
  // in the above, they wouldn't need to be marked as such
  for( const Task::Dependency* dep = task->getComputes(); dep != 0; dep = dep->next ) {
    d_computedVars.insert(dep->var);
    if(d_treatAsOldVars.find(dep->var->getName()) != d_treatAsOldVars.end()) {
      d_initRequires.push_back(dep);
      d_initRequiredVars.insert(dep->var);
    }
  }

  //__________________________________
  // create reduction task if computes included one or more reduction vars
  for( const Task::Dependency* dep = task->getComputes(); dep != 0; dep = dep->next ) {

    if( dep->var->typeDescription()->isReductionVariable() ) {
      int levelidx = dep->reductionLevel ? dep->reductionLevel->getIndex() : -1;
      int dw = dep->mapDataWarehouse();

      if (dep->var->allowsMultipleComputes()) {
        if (dbg.active()) {
          dbg << d_myworld->myrank() << " Skipping Reduction task for multi compute variable: "
              << dep->var->getName() << " on level " << levelidx << ", DW " << dw << '\n';
        }
        continue;
      }

      if (dbg.active()) {
        dbg << d_myworld->myrank() << " Creating Reduction task for variable: " << dep->var->getName()
            << " on level " << levelidx << ", DW " << dw << '\n';
      }

      std::ostringstream taskname;
      taskname << "Reduction: " << dep->var->getName() << ", level " << levelidx << ", dw " << dw;

      Task* newtask = scinew Task(taskname.str(), Task::Reduction);

      int d_dwmap[Task::TotalDWs];

      for (int i = 0; i < Task::TotalDWs; i++) {
        d_dwmap[i] = Task::InvalidDW;
      }

      d_dwmap[Task::OldDW] = Task::NoDW;
      d_dwmap[Task::NewDW] = dw;
      newtask->setMapping(d_dwmap);

      if (dep->matls != 0) {
        newtask->modifies(dep->var, dep->reductionLevel, dep->matls, Task::OutOfDomain);
        for (int i = 0; i < dep->matls->size(); i++) {
          int maltIdx = dep->matls->get(i);
          VarLabelMatl<Level> key(dep->var, maltIdx, dep->reductionLevel);
          d_reductionTasks[key] = newtask;
        }
      }
      else {
        for (int m = 0; m < task->getMaterialSet()->size(); m++) {
          newtask->modifies(dep->var, dep->reductionLevel, task->getMaterialSet()->getSubset(m), Task::OutOfDomain);
          for (int i = 0; i < task->getMaterialSet()->getSubset(m)->size(); i++) {
            int maltIdx = task->getMaterialSet()->getSubset(m)->get(i);
            VarLabelMatl<Level> key(dep->var, maltIdx, dep->reductionLevel);
            d_reductionTasks[key] = newtask;
          }
        }
      }

      d_graphs[d_graphs.size() - 1]->addTask(newtask, 0, 0);
      d_numTasks++;
    }
  }
}

//______________________________________________________________________
//
void
SchedulerCommon::releaseLoadBalancer()
{
  releasePort("load balancer");
}

//______________________________________________________________________
//
void
SchedulerCommon::initialize(int numOldDW /* =1 */, int numNewDW /* =1 */)
{

  // doesn't really do anything except initialize/clear the taskgraph
  //   if the default parameter values are used
  int numDW = numOldDW+numNewDW;
  int oldnum = static_cast<int>(d_dws.size());

  // in AMR cases we will often need to move from many new DWs to one.  In those cases, move the last NewDW to be the next new one.
  if (oldnum - d_numOldDWs > 1) {
    d_dws[numDW-1] = d_dws[oldnum-1];
  }

  // Clear out the data warehouse so that memory will be freed
  for(int i=numDW;i<oldnum;i++)
    d_dws[i]=0;
  d_dws.resize(numDW);
  for(;oldnum < numDW; oldnum++)
    d_dws[oldnum] = 0;
  d_numOldDWs = numOldDW;

  // clear the taskd_graphs, and set the first one
  for (unsigned i = 0; i < d_graphs.size(); i++) {
    delete d_graphs[i];
  }

  d_numParticleGhostCells = 0;

  d_graphs.clear();

  d_initRequires.clear();
  d_initRequiredVars.clear();
  d_computedVars.clear();
  d_numTasks = 0;

  d_maxGhost = 0;
  d_maxLevelOffset = 0;

  d_reductionTasks.clear();
  addTaskGraph(NormalTaskGraph);

}

void SchedulerCommon::setParentDWs(DataWarehouse* parent_old_dw, DataWarehouse* parent_new_dw)
{
  OnDemandDataWarehouse* pold = dynamic_cast<OnDemandDataWarehouse*>(parent_old_dw);
  OnDemandDataWarehouse* pnew = dynamic_cast<OnDemandDataWarehouse*>(parent_new_dw);
  if(parent_old_dw && parent_new_dw){
    ASSERT(pold != 0);
    ASSERT(pnew != 0);
    ASSERT(d_numOldDWs > 2);
    d_dws[0]=pold;
    d_dws[1]=pnew;
  }
}

void SchedulerCommon::clearMappings()
{
  for(int i=0;i<Task::TotalDWs;i++)
    d_dwmap[i]=-1;
}

void SchedulerCommon::mapDataWarehouse(Task::WhichDW which, int dwTag)
{
  ASSERTRANGE(which, 0, Task::TotalDWs);
  ASSERTRANGE(dwTag, 0, static_cast<int>(d_dws.size()));
  d_dwmap[which]=dwTag;
}

//______________________________________________________________________
//
DataWarehouse*
SchedulerCommon::get_dw(int idx)
{
  ASSERTRANGE(idx, 0, static_cast<int>(d_dws.size()));
  return d_dws[idx].get_rep();
}

//______________________________________________________________________
//
DataWarehouse*
SchedulerCommon::getLastDW(void)
{
  return get_dw(static_cast<int>(d_dws.size()) - 1);
}

void 
SchedulerCommon::advanceDataWarehouse(const GridP& grid, bool initialization /*=false*/)
{
  dbg << "advanceDataWarehouse, numDWs = " << d_dws.size() << '\n';
  ASSERT(d_dws.size() >= 2);
  // The last becomes last old, and the rest are new
  d_dws[d_numOldDWs-1] = d_dws[d_dws.size()-1];
  if (d_dws.size() == 2 && d_dws[0] == 0) {
    // first datawarehouse -- indicate that it is the "initialization" dw.
    int generation = d_generation++;
    d_dws[1] = scinew OnDemandDataWarehouse(d_myworld, this, generation, grid,
                                          true /* initialization dw */);
  } else {
    for(int i=d_numOldDWs;i< static_cast<int>(d_dws.size());i++) {
      // in AMR initial cases, you can still be in initialization when you advance again
      replaceDataWarehouse(i, grid, initialization);
    }
  }
}

void SchedulerCommon::fillDataWarehouses(const GridP& grid)
{
  MALLOC_TRACE_TAG_SCOPE("SchedulerCommon::fillDatawarehouses");
  for(int i=d_numOldDWs;i< static_cast<int>(d_dws.size());i++)
    if(!d_dws[i])
      replaceDataWarehouse(i, grid);
}

void SchedulerCommon::replaceDataWarehouse(int index, const GridP& grid, bool initialization /*=false*/)
{
  d_dws[index] = scinew OnDemandDataWarehouse(d_myworld, this, d_generation++, grid, initialization);
  if (initialization) {
    return;
  }
  for (unsigned i = 0; i < d_graphs.size(); i++) {
    DetailedTasks* dts = d_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->copyoutDWKeyDatabase(d_dws[index]);
    }
  }
  d_dws[index]->doReserve();
}

void SchedulerCommon::setRestartable(bool restartable)
{
  this->d_restartable = restartable;
}

const std::vector<const Patch*>* SchedulerCommon::
getSuperPatchExtents(const VarLabel* label, int matlIndex, const Patch* patch,
                     Ghost::GhostType requestedGType, int requestedNumGCells,
                     IntVector& requiredLow, IntVector& requiredHigh,
                     IntVector& requestedLow, IntVector& requestedHigh) const
{
  const SuperPatch* connectedPatchGroup =
    m_locallyComputedPatchVarMap->getConnectedPatchGroup(patch);
  if (connectedPatchGroup == 0)
    return 0;

  SuperPatch::Region requestedExtents = connectedPatchGroup->getRegion();
  SuperPatch::Region requiredExtents = connectedPatchGroup->getRegion();  
  
  // expand to cover the entire connected patch group
  [[maybe_unused]] bool containsGivenPatch = false;
  for (unsigned int i = 0; i < connectedPatchGroup->getBoxes().size(); i++) {
    // get the minimum extents containing both the expected ghost cells
    // to be needed and the given ghost cells.
    const Patch* memberPatch = connectedPatchGroup->getBoxes()[i];

    Patch::VariableBasis basis = Patch::translateTypeToBasis(label->typeDescription()->getType(), true);
    
    IntVector lowOffset=IntVector(0,0,0), highOffset=IntVector(0,0,0);
    
    //set requiredLow and requiredHigh as extents without ghost cells
    memberPatch->computeExtents(basis, label->getBoundaryLayer(), lowOffset, highOffset,
                                requiredLow, requiredHigh);

    //compute ghost cell offsets
    Patch::getGhostOffsets(basis, requestedGType, requestedNumGCells,
                           lowOffset, highOffset);
    
    //set requestedLow and requestedHigh as extents with ghost cells
    memberPatch->computeExtents(basis, label->getBoundaryLayer(), lowOffset, highOffset,
                                requestedLow, requestedHigh);

    SuperPatch::Region requiredRegion =
      SuperPatch::Region(requiredLow, requiredHigh);
    requiredExtents = requiredExtents.enclosingRegion(requiredRegion);
    SuperPatch::Region requestedRegion =
      SuperPatch::Region(requestedLow, requestedHigh);
    requestedExtents = requestedExtents.enclosingRegion(requestedRegion);
    if (memberPatch == patch)
      containsGivenPatch = true;
  }

  ASSERT(containsGivenPatch);
  
  requiredLow = requiredExtents.low_;
  requiredHigh = requiredExtents.high_;
  requestedLow = requestedExtents.low_;
  requestedHigh = requestedExtents.high_;

  // requested extents must enclose the required extents at lesst.
  ASSERTEQ(Min(requiredLow, requestedLow), requestedLow);
  ASSERTEQ(Max(requiredHigh, requestedHigh), requestedHigh);
  
  return &connectedPatchGroup->getBoxes();
}

//______________________________________________________________________
//
void
SchedulerCommon::logMemoryUse()
{
  if(!d_memlogfile){
    std::ostringstream fname;
    fname << "uintah_memuse.log.p" << std::setw(5) << std::setfill('0') << d_myworld->myrank() << "." << d_myworld->size();
    d_memlogfile = scinew std::ofstream(fname.str().c_str());
    if(!*d_memlogfile){
      std::cerr << "Error opening file: " << fname.str() << '\n';
    }
  }
  *d_memlogfile << '\n';
  unsigned long total = 0;
  for(int i=0;i< static_cast<int>(d_dws.size());i++){
    char* name;
    if(i==0)
      name=const_cast<char*>("OldDW");
    else if(i== static_cast<int>(d_dws.size())-1)
      name=const_cast<char*>("NewDW");
    else
      name=const_cast<char*>("IntermediateDW");
    if (d_dws[i])
      d_dws[i]->logMemoryUse(*d_memlogfile, total, name);
  }
  
  for (unsigned i = 0; i < d_graphs.size(); i++) {
    DetailedTasks* dts = d_graphs[i]->getDetailedTasks();
    if (dts) {
      dts->logMemoryUse(*d_memlogfile, total, "Taskgraph");
    }
  }
  *d_memlogfile << "Total: " << total << '\n';
  d_memlogfile->flush();
}

//______________________________________________________________________
//
// Makes and returns a map that maps strings to VarLabels of
// that name and a list of material indices for which that
// variable is valid (according to d_allcomps in graph).
Scheduler::VarLabelMaterialMap* SchedulerCommon::makeVarLabelMaterialMap()
{
  VarLabelMaterialMap* result = scinew VarLabelMaterialMap;
  for (unsigned i = 0; i < d_graphs.size(); i++) {
    d_graphs[i]->makeVarLabelMaterialMap(result);
  }
  return result;
}
     
void SchedulerCommon::doEmitTaskGraphDocs()
{
  d_emit_taskgraph=true;
}

void SchedulerCommon::compile()
{
  GridP grid = const_cast<Grid*>(getLastDW()->getGrid());
  GridP oldGrid;

  if (d_dws[0]) {
    oldGrid = const_cast<Grid*>(get_dw(0)->getGrid());
  }

  if(d_numTasks > 0){

    dbg << d_myworld->myrank() << " SchedulerCommon starting compile\n";
    
    // pass the first to the rest, so we can share the scrubcountTable
    DetailedTasks* first = 0;
    for (unsigned i = 0; i < d_graphs.size(); i++) {
      if (d_graphs.size() > 1) {
        dbg << d_myworld->myrank() << "  Compiling graph#" << i << " of " 
            << d_graphs.size() << endl;
      }

      DetailedTasks* dts = 
        d_graphs[i]->createDetailedTasks(useInternalDeps(), first, grid, oldGrid);

      if (!first) {
        first = dts;
      }
    }
    verifyChecksum();
    dbg << d_myworld->myrank() << " SchedulerCommon finished compile\n";
  } else {
    // NOTE: this was added with scheduleRestartInititalize() support (for empty TGs)
    //  Even when numTasks_ <= 0, the code below executed and did nothing worthwhile... seemingly
    //  Shouldn't need to do this work without tasks though -APH 03/12/15
    return; // no tasks and nothing to do
  }

  m_locallyComputedPatchVarMap->reset();

#if 1
  for (int i = 0; i < grid->numLevels(); i++) {
    const PatchSubset* patches = 
      getLoadBalancer()->getPerProcessorPatchSet(grid->getLevel(i))->getSubset(d_myworld->myrank());
    if (patches->size() > 0) {
      m_locallyComputedPatchVarMap->addComputedPatchSet(patches);
    }
  }
#else
  for (unsigned i = 0; i < d_graphs.size(); i++) { 
    DetailedTasks* dts = d_graphs[i]->getDetailedTasks();
    
    if (dts != 0) {    
      // figure out the locally computed patches for each variable.
      for (int i = 0; i < dts->numLocalTasks(); i++) {
        const DetailedTask* dt = dts->localTask(i);
        for(const Task::Dependency* comp = dt->getTask()->getComputes();
            comp != 0; comp = comp->next){
          if (comp->var->typeDescription()->getType() != TypeDescription::ReductionVariable) {
            constHandle<PatchSubset> patches =
              comp->getPatchesUnderDomain(dt->getPatches());
            m_locallyComputedPatchVarMap->addComputedPatchSet(patches.get_rep());
          }
        }
      }
    }
  }
#endif
  for(unsigned int dw=0;dw<d_dws.size();dw++) {
    if (d_dws[dw].get_rep()) {
      for (unsigned i = 0; i < d_graphs.size(); i++) { 
        DetailedTasks* dts = d_graphs[i]->getDetailedTasks();
        dts->copyoutDWKeyDatabase(d_dws[dw]);
      }
      d_dws[dw]->doReserve();
    }
  }
  m_locallyComputedPatchVarMap->makeGroups();
}

bool SchedulerCommon::isOldDW(int idx) const
{
  ASSERTRANGE(idx, 0, static_cast<int>(d_dws.size()));
  return idx < d_numOldDWs;
}

bool SchedulerCommon::isNewDW(int idx) const
{
  ASSERTRANGE(idx, 0, static_cast<int>(d_dws.size()));
  return idx >= d_numOldDWs;
}

//______________________________________________________________________
//
void
SchedulerCommon::finalizeTimestep()
{
  finalizeNodes(d_myworld->myrank());
  for(unsigned int i=d_numOldDWs;i<d_dws.size();i++)
    d_dws[i]->finalize();
}

//______________________________________________________________________
//
void
SchedulerCommon::scheduleAndDoDataCopy(const GridP& grid, SimulationInterface* sim)
{
  double start = Time::currentSeconds();
  // TODO - use the current initReqs and push them back, instead of doing this...
  // clear the old list of vars and matls
  for (unsigned i = 0; i < d_label_matls.size(); i++)
    for ( label_matl_map::iterator iter = d_label_matls[i].begin(); iter != d_label_matls[i].end(); iter++)
      if (iter->second->removeReference())
        delete iter->second;
  
  d_label_matls.clear();
  d_label_matls.resize(grid->numLevels());

  // produce a map from all tasks' requires from the Old DW.  Store the varlabel and matls
  // TODO - only do this ONCE.
  for (unsigned t = 0; t < d_graphs.size(); t++) {
    TaskGraph* tg = d_graphs[t];
    for (int i = 0; i < tg->getNumTasks(); i++) {
      Task* task = tg->getTask(i);
      if (task->getType() == Task::Output)
        continue;  
      for( const Task::Dependency* dep = task->getRequires(); dep != 0; dep=dep->next ){
        bool copyThisVar = dep->whichdw == Task::OldDW;
        // override to manually copy a var
        if (!copyThisVar){
         
          if (d_copyDataVars.find(dep->var->getName()) != d_copyDataVars.end()) {
            copyThisVar = true;
          }
        }

        // Overide the logic above.  There are PerPatch variables that cannot/shouldn't be copied to the new grid,
        // for example PerPatch<FileInfo>.
        if (d_notCopyDataVars.count(dep->var->getName()) > 0  ){
          copyThisVar = false;
        }

        if (copyThisVar) {
          if (dep->var->typeDescription()->getType() == TypeDescription::ReductionVariable)
            // we will take care of reduction variables in a different section
            continue;
          
          // check the level on the case where variables are only computed on certain levels
          const PatchSet* ps = task->getPatchSet();
          int level = -1;
          if (dep->patches) // just in case the task is over multiple levels...
            level = getLevel(dep->patches)->getIndex();
          else if (ps)
            level = getLevel(ps)->getIndex();
          
          // we don't want data with an invalid level, or requiring from a different level (remember, we are
          // using an old task graph).  That willbe copied later (and chances are, it's to modify anyway).
          if (level == -1 || level > grid->numLevels()-1 || dep->patches_dom == Task::CoarseLevel || 
              dep->patches_dom == Task::FineLevel)
            continue;
          
          const MaterialSubset* matSubset = (dep->matls != 0) ?
            dep->matls : dep->task->getMaterialSet()->getUnion();
          
          
          // if var was already found, make a union of the materials
          MaterialSubset* matls = scinew MaterialSubset(matSubset->getVector());
          matls->addReference();
          
          MaterialSubset* union_matls;
          union_matls = d_label_matls[level][dep->var];
          
          if (union_matls) {
            for (int i = 0; i < union_matls->size(); i++) {
              if (!matls->contains(union_matls->get(i))) {
                matls->add(union_matls->get(i)); 
              } 
            }
            if (union_matls->removeReference()) {
              delete union_matls;
            }
          }
          matls->sort();
          d_label_matls[level][dep->var] = matls;
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
  
  DataWarehouse* oldDataWarehouse = this->get_dw(0);
  DataWarehouse* newDataWarehouse = this->getLastDW();

  oldDataWarehouse->setScrubbing(DataWarehouse::ScrubNone);
  newDataWarehouse->setScrubbing(DataWarehouse::ScrubNone);
  const Grid* oldGrid = oldDataWarehouse->getGrid();
  std::vector<Task*> dataTasks;
  std::vector<Handle<PatchSet> > refinePatchSets(grid->numLevels(), (PatchSet*)0);
  std::vector<Handle<PatchSet> > copyPatchSets(grid->numLevels(), (PatchSet*)0);
  SchedulerP sched(dynamic_cast<Scheduler*>(this));

  d_sharedState->setCopyDataTimestep(true);

  for (int i = 0; i < grid->numLevels(); i++) {
    LevelP newLevel = grid->getLevel(i);
    //const PatchSubset  *patches = getLoadBalancer()->getPerProcessorPatchSet(newLevel)->getSubset(d_myworld->myrank());
    if (i > 0) {
      if (i >= oldGrid->numLevels()) {
        // new level - refine everywhere
        refinePatchSets[i] = const_cast<PatchSet*>(newLevel->eachPatch());
        copyPatchSets[i] = scinew PatchSet;
      }
      // find patches with new space - but temporarily, refine everywhere... 
      else if (i < oldGrid->numLevels()) {
        refinePatchSets[i] = scinew PatchSet;
        copyPatchSets[i] = scinew PatchSet;

        std::vector<int> myPatchIDs;
        LevelP oldLevel = oldDataWarehouse->getGrid()->getLevel(i);
        
        // go through the patches, and find if there are patches that weren't entirely 
        // covered by patches on the old grid, and interpolate them.  
        // then after, copy the data, and if necessary, overwrite interpolated data
        const PatchSubset *ps = 
          getLoadBalancer()->getPerProcessorPatchSet(newLevel)->getSubset(d_myworld->myrank());

        // for each patch I own
        for (int p = 0; p < ps->size(); p++) {
          const Patch *newPatch = ps->get(p);
        
          // get the low/high for what we'll need to get
          IntVector lowIndex, highIndex;
          //newPatch->computeVariableExtents(Patch::CellBased, IntVector(0,0,0), Ghost::None, 0, lowIndex, highIndex);
          lowIndex = newPatch->getCellLowIndex();
          highIndex = newPatch->getCellHighIndex();
          
          // find if area on the new patch was not covered by the old patches
          IntVector dist = highIndex-lowIndex;
          int totalCells = dist.x()*dist.y()*dist.z();
          int sum = 0;
          Patch::selectType oldPatches;
          oldLevel->selectPatches(lowIndex, highIndex, oldPatches);
         
          //compute volume of overlapping regions
          for (int old = 0; old < oldPatches.size(); old++) {

            const Patch* oldPatch = oldPatches[old];
            IntVector oldLow = oldPatch->getCellLowIndex();
            IntVector oldHigh = oldPatch->getCellHighIndex();

            IntVector low = Max(oldLow, lowIndex);
            IntVector high = Min(oldHigh, highIndex);
            IntVector dist = high-low;
            sum += dist.x()*dist.y()*dist.z();
          }  // for oldPatches
          if (sum != totalCells) {
            if (Uintah::Parallel::usingMPI()) {
              myPatchIDs.push_back(newPatch->getID());
            } else {
              refinePatchSets[i]->add(newPatch);
            }
          } else {
            if (!Uintah::Parallel::usingMPI()) {
              copyPatchSets[i]->add(newPatch);
            }
          }
        } // for patch

        if (Uintah::Parallel::usingMPI()) {

          //Gather size from all processors
          int mycount = myPatchIDs.size();
          std::vector<int> counts(d_myworld->size());
          MPI_Allgather(&mycount, 1, MPI_INT, &counts[0], 1, MPI_INT, d_myworld->getComm());

          //compute recieve array offset and size
          std::vector<int> displs(d_myworld->size());
          int pos = 0;

          for (int p = 0; p < d_myworld->size(); p++) {
            displs[p] = pos;
            pos += counts[p];
          }

          std::vector<int> allPatchIDs(pos);  //receive array;
          MPI_Allgatherv(&myPatchIDs[0], counts[d_myworld->myrank()], 
                         MPI_INT, &allPatchIDs[0], &counts[0], &displs[0], MPI_INT,
                         d_myworld->getComm());
          //make refinePatchSets from patch ids
          std::set<int> allPatchIDset(allPatchIDs.begin(), allPatchIDs.end());

          for (auto iter = newLevel->patchesBegin(); iter != newLevel->patchesEnd(); ++iter) {
            Patch* newPatch = *iter;
            if (allPatchIDset.find(newPatch->getID()) != allPatchIDset.end()) {
              refinePatchSets[i]->add(newPatch);
            }
            else {
              copyPatchSets[i]->add(newPatch);
            }
          }
        }  // using MPI
      }
      if (refinePatchSets[i]->size() > 0) {
        dbg << d_myworld->myrank() << "  Calling scheduleRefine for patches " << *refinePatchSets[i].get_rep() << endl;
        sim->scheduleRefine(refinePatchSets[i].get_rep(), sched);
      }
    } else {
      refinePatchSets[i] = scinew PatchSet;
      copyPatchSets[i] = const_cast<PatchSet*>(newLevel->eachPatch());
    }

    //__________________________________
    //  Scheduling for copyDataToNewGrid
    if (copyPatchSets[i]->size() > 0) {
      dataTasks.push_back(scinew Task("SchedulerCommon::copyDataToNewGrid", this, 
                                      &SchedulerCommon::copyDataToNewGrid));

      for (auto iter = d_label_matls[i].begin(); iter != d_label_matls[i].end(); iter++) {
        const VarLabel* var = iter->first;
        MaterialSubset* matls = iter->second;

        dataTasks.back()->requires(Task::OldDW, var, 0, Task::OtherGridDomain, matls, 
                                   Task::NormalDomain, Ghost::None, 0);
        dbg << "  Scheduling copy for var " << *var << " matl " << *matls << " Copies: "
            << *copyPatchSets[i].get_rep() << endl;
        dataTasks.back()->computes(var, matls);
      }
      addTask(dataTasks.back(), copyPatchSets[i].get_rep(), d_sharedState->allMaterials());
    }

    //__________________________________
    //  Scheduling for modifyDataOnNewGrid
    if (refinePatchSets[i]->size() > 0) {
      dataTasks.push_back(scinew Task("SchedulerCommon::modifyDataOnNewGrid", this, 
                                      &SchedulerCommon::copyDataToNewGrid));

      for (auto iter = d_label_matls[i].begin(); iter != d_label_matls[i].end(); iter++) {
        const VarLabel* var = iter->first;
        MaterialSubset* matls = iter->second;

        dataTasks.back()->requires(Task::OldDW, var, 0, Task::OtherGridDomain, matls, 
                                   Task::NormalDomain, Ghost::None, 0);
        dbg << "  Scheduling modify for var " << *var << " matl " << *matls << " Modifies: "
            << *refinePatchSets[i].get_rep() << endl;
        dataTasks.back()->modifies(var, matls);
      }
      addTask(dataTasks.back(), refinePatchSets[i].get_rep(), d_sharedState->allMaterials());
    }

    if (i > 0) {
      sim->scheduleRefineInterface(newLevel, sched, 0, 1);
    }
  }

  // set so the load balancer will make an adequate neighborhood, as the default
  // neighborhood isn't good enough for the copy data timestep
  d_sharedState->setCopyDataTimestep(true);  //-- do we still need this?  - BJW
#if !defined( _WIN32 ) && !defined( DISABLE_SCI_MALLOC )
  const char* tag = AllocatorSetDefaultTag("DoDataCopy");
#endif
  this->compile(); 
  d_sharedState->regriddingCompilationTime += Time::currentSeconds() - start;
  
  // save these and restore them, since the next execute will append the scheduler's, and we don't want to.
  double executeTime = d_sharedState->taskExecTime;
  double globalCommTime = d_sharedState->taskGlobalCommTime;
  double localCommTime = d_sharedState->taskLocalCommTime;
  
  start = Time::currentSeconds();
  this->execute();
#if !defined( _WIN32 ) && !defined( DISABLE_SCI_MALLOC )
  AllocatorSetDefaultTag(tag);
#endif

  // copy reduction variables to the new dw
  std::vector<VarLabelMatl<Level> > levelVariableInfo;
  oldDataWarehouse->getVarLabelMatlLevelTriples(levelVariableInfo);
  
  newDataWarehouse->unfinalize();
  for ( unsigned int i = 0; i < levelVariableInfo.size(); i++ ) {
    VarLabelMatl<Level> currentReductionVar = levelVariableInfo[i];

    if (currentReductionVar.label_->typeDescription()->isReductionVariable()) {

      // std::cout << "REDUNCTION:  Label(" << std::setw(15) << currentReductionVar.label_->getName() << "): Patch(" << reinterpret_cast<int>(currentReductionVar.level_) << "): Material(" << currentReductionVar.matlIndex_ << ")" << endl; 
      const Level* oldLevel = currentReductionVar.domain_;
      const Level* newLevel = NULL;
      if (oldLevel && oldLevel->getIndex() < grid->numLevels() ) {
        if (oldLevel->getIndex() >= grid->numLevels()) {
          // the new grid no longer has this level
          continue;
        }
        newLevel = (newDataWarehouse->getGrid()->getLevel( oldLevel->getIndex() )).get_rep();
      }
   
      //Either both levels need to be null or both need to exist (null levels mean global data)
      if(!oldLevel || newLevel)
      {
        ReductionVariableBase* v = dynamic_cast<ReductionVariableBase*>(currentReductionVar.label_->typeDescription()->createInstance());
        oldDataWarehouse->get(*v, currentReductionVar.label_, currentReductionVar.domain_, currentReductionVar.matlIndex_);
        newDataWarehouse->put(*v, currentReductionVar.label_, newLevel, currentReductionVar.matlIndex_);
        delete v; // copied on the put command
      }
    }
  }

  newDataWarehouse->refinalize();

  d_sharedState->regriddingCopyDataTime += Time::currentSeconds() - start;
  d_sharedState->taskExecTime = executeTime;
  d_sharedState->taskGlobalCommTime = globalCommTime;
  d_sharedState->taskLocalCommTime = localCommTime;
  d_sharedState->setCopyDataTimestep(false);
}


void
SchedulerCommon::copyDataToNewGrid(const ProcessorGroup*, const PatchSubset* patches,
                                   const MaterialSubset* matls, DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  dbg << "SchedulerCommon::copyDataToNewGrid() BGN on patches " << *patches  << endl;

  OnDemandDataWarehouse* oldDataWarehouse = dynamic_cast<OnDemandDataWarehouse*>(old_dw);
  OnDemandDataWarehouse* newDataWarehouse = dynamic_cast<OnDemandDataWarehouse*>(new_dw);

  // For each patch in the patch subset which contains patches in the new grid
  for ( int p = 0; p < patches->size(); p++ ) {
    const Patch* newPatch = patches->get(p);
    const Level* newLevel = newPatch->getLevel();

    // to create once per matl instead of once per matl-var
    std::vector<ParticleSubset*> oldsubsets(d_sharedState->getNumMatls()), newsubsets(d_sharedState->getNumMatls());

    // If there is a level that didn't exist, we don't need to copy it
    if ( newLevel->getIndex() >= oldDataWarehouse->getGrid()->numLevels() ) {
      continue;
    }
    
    // find old patches associated with this patch
    LevelP oldLevel = oldDataWarehouse->getGrid()->getLevel( newLevel->getIndex() );
    
    for ( label_matl_map::iterator iter = d_label_matls[oldLevel->getIndex()].begin(); 
          iter != d_label_matls[oldLevel->getIndex()].end(); iter++) {
      const VarLabel* label = iter->first;
      MaterialSubset* var_matls = iter->second;


      // get the low/high for what we'll need to get
      Patch::VariableBasis basis = Patch::translateTypeToBasis(label->typeDescription()->getType(), true);
      IntVector newLowIndex, newHighIndex;
      newPatch->computeVariableExtents(basis, IntVector(0,0,0), Ghost::None, 0, newLowIndex, newHighIndex);

      //__________________________________
      //  Loop over materials
      for (int m = 0; m < var_matls->size(); m++) {
        int matl = var_matls->get(m);

        if (!matls->contains(matl)) {
          //std::cout << "We are skipping material " << currentVar.matlIndex_ << endl;
          continue;
        }

        //__________________________________
        //  Grid Variables
        if (label->typeDescription()->getType() != TypeDescription::ParticleVariable) {
          Patch::selectType oldPatches;
          oldLevel->selectPatches(newLowIndex, newHighIndex, oldPatches);
          
          for ( int oldIdx = 0;  oldIdx < oldPatches.size(); oldIdx++) {
            const Patch* oldPatch = oldPatches[oldIdx];
            
            if(!oldDataWarehouse->exists(label, matl, oldPatch))
              continue; // see comment about oldPatchToTest in ScheduleAndDoDataCopy
            IntVector oldLowIndex;
            IntVector oldHighIndex;
            
            if (newLevel->getIndex() > 0) {
              oldLowIndex = oldPatch->getLowIndexWithDomainLayer(basis);
              oldHighIndex = oldPatch->getHighIndexWithDomainLayer(basis);
            }
            else {
              oldLowIndex = oldPatch->getExtraLowIndex(basis, label->getBoundaryLayer());
              oldHighIndex = oldPatch->getExtraHighIndex(basis, label->getBoundaryLayer());
            }
            
            IntVector copyLowIndex = Max(newLowIndex, oldLowIndex);
            IntVector copyHighIndex = Min(newHighIndex, oldHighIndex);
            
            // based on the selectPatches above, we might have patches we don't want to use, so prune them here.
            if (copyLowIndex.x() >= copyHighIndex.x() || copyLowIndex.y() >= copyHighIndex.y() || copyLowIndex.z() >= copyHighIndex.z())
              continue;
            
            
            switch(label->typeDescription()->getType()){
            case TypeDescription::NCVariable:
            case TypeDescription::CCVariable:
            case TypeDescription::SFCXVariable:
            case TypeDescription::SFCYVariable:
            case TypeDescription::SFCZVariable:
            {
              if(!oldDataWarehouse->exists(label, matl, oldPatch))
                SCI_THROW(UnknownVariable(label->getName(), oldDataWarehouse->getID(), oldPatch, matl,
                                          "in copyDataTo GridVariableBase", __FILE__, __LINE__));
              std::vector<Variable *> varlist;
              oldDataWarehouse->d_varDB.getlist(label, matl, oldPatch, varlist);
              GridVariableBase* v=NULL;

              IntVector srclow = copyLowIndex;
              IntVector srchigh = copyHighIndex;

              for (unsigned int i = 0 ; i < varlist.size(); ++i) {
                v = dynamic_cast<GridVariableBase*>(varlist[i]);

                ASSERT(v->getBasePointer()!=0);

                //restrict copy to data range
                srclow =  Max(copyLowIndex, v->getLow());
                srchigh = Min(copyHighIndex, v->getHigh());
                if (srclow.x() >= srchigh.x() || srclow.y() >= srchigh.y() || srclow.z() >= srchigh.z()) continue;

                if ( !newDataWarehouse->exists(label, matl, newPatch) ) {
                  GridVariableBase* newVariable = v->cloneType();
                  newVariable->rewindow( newLowIndex, newHighIndex );
                  newVariable->copyPatch( v, srclow, srchigh);
                  newDataWarehouse->d_varDB.put(label, matl, newPatch, newVariable, isCopyDataTimestep(), false);
                } else {
                  GridVariableBase* newVariable = 
                    dynamic_cast<GridVariableBase*>(newDataWarehouse->d_varDB.get(label, matl, newPatch ));
                  // make sure it exists in the right region (it might be ghost data)
                  newVariable->rewindow(newLowIndex, newHighIndex);
                  if (oldPatch->isVirtual()) {
                    // it can happen where the old patch was virtual and this is not
                    GridVariableBase* tmpVar = newVariable->cloneType();
                    tmpVar->copyPointer(*v);
                    tmpVar->offset(oldPatch->getVirtualOffset());
                    newVariable->copyPatch( tmpVar, srclow, srchigh );
                    delete tmpVar;
                  }
                  else
                    newVariable->copyPatch( v, srclow, srchigh );
                }
              }
            }
            break;
            case TypeDescription::PerPatch:
            {
            }
            break;
            default:
              SCI_THROW(InternalError("Unknown variable type in copyData: "+label->getName(), __FILE__, __LINE__));
            } // end switch
          } // end oldPatches
        }
        else {
          //__________________________________
          //  Particle Variables
          ParticleSubset* oldsub = oldsubsets[matl];
          if (!oldsub) {
            // collect the particles from the range encompassing this patch.  Use interior cells since
            // extracells aren't collected across processors in the data copy, and they don't matter
            // for particles anyhow (but we will have to reset the bounds to copy the data)
            oldsub = oldDataWarehouse->getParticleSubset(matl, newPatch->getLowIndexWithDomainLayer(Patch::CellBased),
                                                         newPatch->getHighIndexWithDomainLayer(Patch::CellBased), 
                                                         newPatch, d_reloc_new_posLabel, oldLevel.get_rep());
            oldsubsets[matl] = oldsub;
            oldsub->addReference();
          }


          ParticleSubset* newsub = newsubsets[matl];
          // it might have been created in Refine
          if (!newsub) {
            if (!newDataWarehouse->haveParticleSubset(matl, newPatch))
              newsub = newDataWarehouse->createParticleSubset(oldsub->numParticles(), matl, newPatch);
            else {
              newsub = newDataWarehouse->getParticleSubset(matl, newPatch);
              ASSERT(newsub->numParticles() == 0);
              newsub->addParticles(oldsub->numParticles());
            }
            newsubsets[matl] = newsub;
          }

          ParticleVariableBase* newv = dynamic_cast<ParticleVariableBase*>(label->typeDescription()->createInstance());
          newv->allocate(newsub);
          // don't get and copy if there were no old patches
          if (oldsub->getNeighbors().size() > 0) {

            constParticleVariableBase* var = newv->cloneConstType();
            oldDataWarehouse->get(*var, label, oldsub);

            // reset the bounds of the old var's data so copyData doesn't complain
            ParticleSubset* tempset = scinew ParticleSubset(oldsub->numParticles(), matl, newPatch,
                                                            newPatch->getExtraCellLowIndex(), newPatch->getExtraCellHighIndex());
            const_cast<ParticleVariableBase*>(&var->getBaseRep())->setParticleSubset(tempset);
            newv->copyData(&var->getBaseRep());
            delete var; //pset and tempset are deleted with it.
          }
          newDataWarehouse->put(*newv, label, true);
          delete newv; // the container is copied
        }
      } // end matls
    } // end label_matls
    for (unsigned i = 0; i < oldsubsets.size(); i++)
      if (oldsubsets[i] && oldsubsets[i]->removeReference())
        delete oldsubsets[i];
  } // end patches

    // d_lock.writeUnlock(); Do we need this?

  dbg << "SchedulerCommon::copyDataToNewGrid() END" << endl;
}

//______________________________________________________________________
//
void
SchedulerCommon::scheduleParticleRelocation(const LevelP& level,
                                            const VarLabel* old_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& old_labels,
                                            const VarLabel* new_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& new_labels,
                                            const VarLabel* particleIDLabel,
                                            const MaterialSet* matls, int which)
{
  dbg << "Inside scheduleParticleRelocation:"<< __FILE__<<":"<<__LINE__<<std::endl;
  if (which == 1) {
    if (d_reloc_new_posLabel) {
      ASSERTEQ(d_reloc_new_posLabel, new_posLabel);
    }
    d_reloc_new_posLabel = new_posLabel;
    UintahParallelPort* lbp = getPort("load balancer");
    LoadBalancer* lb = dynamic_cast<LoadBalancer*>(lbp);

    dbg << "Calling Relocate::scheduleParticleRelocation::d_reloc1"<< __FILE__<<":"<<__LINE__<<std::endl;
    d_reloc1.scheduleParticleRelocation( this, d_myworld, lb, level,
                                        old_posLabel, old_labels,
                                        new_posLabel, new_labels,
                                        particleIDLabel, matls );
    releasePort("load balancer");
  }

  if(which == 2){
    if (d_reloc_new_posLabel)
      ASSERTEQ(d_reloc_new_posLabel, new_posLabel);
    d_reloc_new_posLabel = new_posLabel;
    UintahParallelPort* lbp = getPort("load balancer");
    LoadBalancer* lb = dynamic_cast<LoadBalancer*>(lbp);
    dbg << "Calling Relocate::scheduleParticleRelocation::d_reloc2"<< __FILE__<<":"<<__LINE__<<std::endl;
    d_reloc2.scheduleParticleRelocation( this, d_myworld, lb, level,
                                        old_posLabel, old_labels,
                                        new_posLabel, new_labels,
                                        particleIDLabel, matls );
    releasePort("load balancer");
  }
}

//______________________________________________________________________
//
void
SchedulerCommon::scheduleParticleRelocation(const LevelP& coarsestLevelwithParticles,
                                            const VarLabel* old_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& old_labels,
                                            const VarLabel* new_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& new_labels,
                                            const VarLabel* particleIDLabel,
                                            const MaterialSet* matls)
{
  if (d_reloc_new_posLabel)
    ASSERTEQ(d_reloc_new_posLabel, new_posLabel);
  d_reloc_new_posLabel = new_posLabel;
  UintahParallelPort* lbp = getPort("load balancer");
  LoadBalancer* lb = dynamic_cast<LoadBalancer*>(lbp);
  d_reloc1.scheduleParticleRelocation( this, d_myworld, lb, coarsestLevelwithParticles,
                                      old_posLabel, old_labels,
                                      new_posLabel, new_labels,
                                      particleIDLabel, matls );
  releasePort("load balancer");
}

//______________________________________________________________________
//
void
SchedulerCommon::scheduleParticleRelocation( const LevelP&                           coarsestLevelwithParticles,
                                             const VarLabel*                         posLabel,
                                             const std::vector<std::vector<const VarLabel*> >& otherLabels,
                                             const MaterialSet*                      matls )
{
  d_reloc_new_posLabel = posLabel;
  UintahParallelPort* lbp = getPort("load balancer");
  LoadBalancer* lb = dynamic_cast<LoadBalancer*>(lbp);
  d_reloc1.scheduleParticleRelocation(this, d_myworld, lb, coarsestLevelwithParticles, posLabel, otherLabels, matls);
  releasePort("load balancer");
}

void SchedulerCommon::overrideVariableBehavior(const std::string& var, bool treatAsOld, 
                                               bool copyData, bool noScrub,
                                               bool notCopyData, bool noCheckpoint)
{
  // treat variable as an "old" var - will be checkpointed, copied, and only scrubbed from an OldDW
  if (treatAsOld) {
    d_treatAsOldVars.insert(var);
  }
  
  // manually copy variable between AMR levels
  if (copyData) {
    d_copyDataVars.insert(var);
    d_noScrubVars.insert(var);
  }
  
  // ignore copying this variable between AMR levels
  if (notCopyData) {
    d_notCopyDataVars.insert(var);
  }
    
  // set variable not to scrub (normally when needed between a normal taskgraph
  // and the regridding phase)
  if (noScrub) {
    d_noScrubVars.insert(var);
  }
  
  // do not checkpoint this variable.
  if ( noCheckpoint ){
    d_notCheckpointVars.insert(var);
  }

}

//______________________________________________________________________
// output the task name and the level it's executing on.
// and each of the patches
void
SchedulerCommon::printTask( ostream& out, DetailedTask* task )
{
  out << std::left;
  out.width(70);
  out << task->getTask()->getName();
  if(task->getPatches()){

    out << " \t on patches ";
    const PatchSubset* patches = task->getPatches();
    for(int p=0;p<patches->size();p++){
      if(p != 0)
        out << ", ";
      out << patches->get(p)->getID();
    }
    
    if (task->getTask()->getType() != Task::OncePerProc) {
      const Level* level = getLevel(patches);
      out << "\t  L-"<< level->getIndex();
    }
  }
}


//______________________________________________________________________
//  Output the task name and the level it's executing on
//  only first patch of that level
void 
SchedulerCommon::printTaskLevels( const ProcessorGroup* d_myworld, 
                                  DebugStream & out, 
                                  DetailedTask* task )
{
  if(out.active() ){
    if(task->getPatches()){

      if (task->getTask()->getType() != Task::OncePerProc) {
      
        const PatchSubset* taskPatches = task->getPatches();
        
        const Level* level = getLevel(taskPatches);
        const Patch* firstPatch = level->getPatch(0);
        
        if(taskPatches->contains(firstPatch)){
        
          out << d_myworld->myrank() << "   ";
          out << std::left;
          out.width(70);
          out << task->getTask()->getName();
          out << " \t  Patch " << firstPatch->getGridIndex() << "\t L-"<< level->getIndex() << "\n";
        }
      }
    }
  }  //debugstream active
} 
