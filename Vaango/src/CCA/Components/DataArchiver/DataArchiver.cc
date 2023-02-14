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

#include <CCA/Components/DataArchiver/DataArchiver.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/ModelInterface.h>
#include <CCA/Ports/OutputContext.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/StringUtil.h>

#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Endian.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/Timers/Timers.hpp>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <strings.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include <libxml/xmlwriter.h>

// TODO - BJW - if multilevel reduction doesn't work, fix all
//        getMaterialSet(0)

#define PADSIZE 1024L
#define ALL_LEVELS 99

#define OUTPUT 0
#define CHECKPOINT 1
#define CHECKPOINT_GLOBAL 2

#define XML_TEXTWRITER 1
#undef XML_TEXTWRITER

namespace {

#ifdef HAVE_PIDX
Uintah::DebugStream dbgPIDX("DataArchiverPIDX",
                            "DataArchiver",
                            "Data archiver PIDX debug stream",
                            false);
#endif

} // namespace

namespace Uintah {

//  Used externally, keep them visible outside this unit
Dout g_DA_dbg("DataArchiver_DBG",
              "DataArchiver",
              "general debug information",
              false);

DataArchiver::DataArchiver(const ProcessorGroup* myworld, int udaSuffix)
  : UintahParallelComponent(myworld)
  , d_udaSuffix(udaSuffix)
{
  d_sync_io_label =
    VarLabel::create("sync_io_vl", CCVariable<float>::getTypeDescription());

  d_tmpMatSubset = scinew MaterialSubset();
  d_tmpMatSubset->add(-1);
  d_tmpMatSubset->addReference();
}

DataArchiver::~DataArchiver()
{
  DOUTR(g_DA_dbg, "DataArchiver::~DataArchiver()");
  VarLabel::destroy(d_sync_io_label);

  if (d_tmpMatSubset && d_tmpMatSubset->removeReference()) {
    delete d_tmpMatSubset;
  }
}

void
DataArchiver::getComponents()
{
  d_simulator = dynamic_cast<SimulationInterface*>(getPort("simulator"));
  if (!d_simulator) {
    throw InternalError(
      "dynamic_cast of 'd_simulator' failed!", __FILE__, __LINE__);
  }

  d_load_balancer = dynamic_cast<LoadBalancer*>(getPort("load balancer"));
  if (!d_load_balancer) {
    throw InternalError(
      "dynamic_cast of 'd_load_balancer' failed!", __FILE__, __LINE__);
  }
}

void
DataArchiver::releaseComponents()
{
  releasePort("simulator");
  releasePort("load balancer");

  d_simulator     = nullptr;
  d_load_balancer = nullptr;

  d_mat_manager = nullptr;
}

void
DataArchiver::problemSetup(const ProblemSpecP& params,
                           const ProblemSpecP& restart_ps,
                           const MaterialManagerP& mat_manager)
{
  DOUTR(g_DA_dbg, "DataArchiver::ProblemSetup");

  d_mat_manager = mat_manager;
  d_upsFile     = params;

  ProblemSpecP p = params->findBlock("DataArchiver");

  if (restart_ps) {

    ProblemSpecP insitu_ps = restart_ps->findBlock("InSitu");

    if (insitu_ps != nullptr) {

      bool haveModifiedVars;
      insitu_ps = insitu_ps->get("haveModifiedVars", haveModifiedVars);

      d_simulator->haveModifiedVars(haveModifiedVars);

      if (haveModifiedVars) {
        std::stringstream tmpstr;
        tmpstr << "DataArchiver found previously modified variables that "
               << "have not been merged into the checkpoint restart "
               << "input.xml file from the from index.xml file.\n"
               << "The modified variables can be found in the "
               << "index.xml file under the 'InSitu' block.\n"
               << "Once merged, change the variable 'haveModifiedVars' in "
               << "the 'InSitu' block in the checkpoint restart timestep.xml "
               << "file to 'false'";

        throw ProblemSetupException(tmpstr.str(), __FILE__, __LINE__);
      } else {
        proc0cout << "DataArchiver found previously modified vars. "
                  << "Assuming the checkpoint restart input.xml file "
                  << "has been updated.\n";
      }
    }
  }

  // PIDX related
  string type;
  p->getAttribute("type", type);
#ifdef HAVE_PIDX
  if (type == "pidx" || type == "PIDX") {
    d_outputFileFormat = PIDX;
    d_PIDX_flags.problemSetup(p);

    // Debug:
    d_PIDX_flags.print();
  }
#endif

  d_outputDoubleAsFloat = p->findBlock("outputDoubleAsFloat") != nullptr;

  // For outputing the sim time and/or time step with the global vars
  p->get("timeStep", d_outputGlobalVarsTimeStep); // default false
  p->get("simTime", d_outputGlobalVarsSimTime);   // default true

  // For modulating the output frequency global vars. By default
  // they are output every time step. Note: Frequency > OnTimeStep
  p->get("frequency", d_outputGlobalVarsFrequency);   // default 1
  p->get("onTimeStep", d_outputGlobalVarsOnTimeStep); // default 0

  if (d_outputGlobalVarsOnTimeStep >= d_outputGlobalVarsFrequency) {
    proc0cout << "Error: the frequency of outputing the global vars "
              << d_outputGlobalVarsFrequency << " "
              << "is less than or equal to the time step ordinality "
              << d_outputGlobalVarsOnTimeStep << " "
              << ". Resetting the ordinality to ";

    if (d_outputGlobalVarsFrequency > 1) {
      d_outputGlobalVarsOnTimeStep = 1;
    } else {
      d_outputGlobalVarsOnTimeStep = 0;
    }
    proc0cout << d_outputGlobalVarsOnTimeStep << std::endl;
  }

  // set to false if restartSetup is called - we can't do it there
  // as the first timestep doesn't have any tasks
  d_outputInitTimeStep = p->findBlock("outputInitTimeStep") != 0;

  if (d_filebase == "") {
    p->require("filebase", d_filebase);
  }

  // get output timestep or time interval info
  d_outputInterval = 0;
  if (!p->get("outputTimeStepInterval", d_outputTimeStepInterval)) {
    d_outputTimeStepInterval = 0;
  }

  if (!p->get("outputInterval", d_outputInterval) &&
      d_outputTimeStepInterval == 0) {
    d_outputInterval = 0.0; // default
  }

  if (d_outputInterval != 0.0 && d_outputTimeStepInterval != 0) {
    throw ProblemSetupException(
      "Use <outputInterval> or <outputTimeStepInterval>, not both",
      __FILE__,
      __LINE__);
  }

  // set default compression mode - can be "gzip" or ""
  string defaultCompressionMode = "";
  if (p->get("compression", defaultCompressionMode)) {
    VarLabel::setDefaultCompressionMode(defaultCompressionMode);
  }

  if (params->findBlock("ParticlePosition")) {
    params->findBlock("ParticlePosition")
      ->getAttribute("label", d_particlePositionName);
  }

  // get the variables to save
  d_saveLabelNames.clear(); // we can problemSetup multiple times on a component
                            // Switch, clear the old ones.
  std::map<std::string, std::string> attributes;
  SaveNameItem saveItem;
  ProblemSpecP save = p->findBlock("save");

  if (save == nullptr) {
    proc0cout << "\nWARNING: No data will be saved as none was specified to be "
                 "saved in the .ups file!\n\n";
  }

  while (save != nullptr) {
    attributes.clear();
    save->getAttributes(attributes);
    saveItem.labelName       = attributes["label"];
    saveItem.compressionMode = attributes["compression"];
    try {
      saveItem.matls = ConsecutiveRangeSet(attributes["material"]);
    } catch (const ConsecutiveRangeSetException& e) {
      throw ProblemSetupException("'" + attributes["material"] + "'" +
                                    " cannot be parsed as a set of material" +
                                    " indices for saving '" +
                                    saveItem.labelName + "'",
                                  __FILE__,
                                  __LINE__);
    }

    if (saveItem.matls.size() == 0) {
      // if materials aren't specified, all valid materials will be saved
      saveItem.matls = ConsecutiveRangeSet::all;
    }

    try {
      // level values are normal for exact level values (i.e.,
      //   0, 2-4 will save those levels.
      //   Leave blank for all levels
      //   Also handle negative values to be relative to the finest levels
      //   I.e., -3--1 would be the top three levels.
      saveItem.levels = ConsecutiveRangeSet(attributes["levels"]);
    } catch (const ConsecutiveRangeSetException& e) {
      throw ProblemSetupException("'" + attributes["levels"] + "'" +
                                    " cannot be parsed as a set of levels" +
                                    " for saving '" + saveItem.labelName + "'",
                                  __FILE__,
                                  __LINE__);
    }

    if (saveItem.levels.size() == 0) {
      // if materials aren't specified, all valid materials will be saved
      saveItem.levels = ConsecutiveRangeSet(ALL_LEVELS, ALL_LEVELS);
    }

    //  bullet proofing: must save p.x
    //  in addition to other particle variables "p.*"
    if (saveItem.labelName == d_particlePositionName ||
        saveItem.labelName == "p.xx") {
      d_saveP_x = true;
    }

    string::size_type pos = saveItem.labelName.find("p.");
    if (pos != string::npos && saveItem.labelName != d_particlePositionName) {
      d_saveParticleVariables = true;
    }

    d_saveLabelNames.push_back(saveItem);

    save = save->findNextBlock("save");
  }

  if (d_saveP_x == false && d_saveParticleVariables == true) {
    throw ProblemSetupException(" You must save " + d_particlePositionName +
                                  " when saving other particle variables",
                                __FILE__,
                                __LINE__);
  }

  // get checkpoint information
  d_checkpointInterval         = 0.0;
  d_checkpointTimeStepInterval = 0;
  d_checkpointWallTimeStart    = 0;
  d_checkpointWallTimeInterval = 0;
  d_checkpointCycle            = 2; /* 2 is the smallest number that is safe
                                       (always keeping an older copy for backup) */
  d_checkpointLastTimeStep = false;

  ProblemSpecP checkpoint = p->findBlock("checkpoint");
  if (checkpoint != nullptr) {

    attributes.clear();
    checkpoint->getAttributes(attributes);

    std::string interval              = attributes["interval"];
    std::string timestepInterval      = attributes["timestepInterval"];
    std::string walltimeStart         = attributes["walltimeStart"];
    std::string walltimeInterval      = attributes["walltimeInterval"];
    std::string walltimeStartHours    = attributes["walltimeStartHours"];
    std::string walltimeIntervalHours = attributes["walltimeIntervalHours"];
    std::string cycle                 = attributes["cycle"];
    std::string lastTimeStep          = attributes["lastTimestep"];

    if (interval != "") {
      d_checkpointInterval = atof(interval.c_str());
    }
    if (timestepInterval != "") {
      d_checkpointTimeStepInterval = atoi(timestepInterval.c_str());
    }
    if (walltimeStart != "") {
      d_checkpointWallTimeStart = atof(walltimeStart.c_str());
    }
    if (walltimeInterval != "") {
      d_checkpointWallTimeInterval = atof(walltimeInterval.c_str());
    }
    if (walltimeStartHours != "") {
      d_checkpointWallTimeStart = atof(walltimeStartHours.c_str()) * 3600.0;
    }
    if (walltimeIntervalHours != "") {
      d_checkpointWallTimeInterval =
        atof(walltimeIntervalHours.c_str()) * 3600.0;
    }
    if (cycle != "") {
      d_checkpointCycle = atoi(cycle.c_str());
    }
    if (lastTimeStep == "true") {
      d_checkpointLastTimeStep = true;
    }

    // Verify that an interval was specified:
    if (interval == "" && timestepInterval == "" && walltimeInterval == "" &&
        walltimeIntervalHours == "") {
      throw ProblemSetupException(
        "ERROR: \n  <checkpoint> must specify either interval, "
        "timestepInterval, walltimeInterval",
        __FILE__,
        __LINE__);
    }
  }

  // Can't use both checkpointInterval and checkpointTimeStepInterval and
  // checkpointWallTimeInterval.
  if (((int)(d_checkpointInterval > 0.0) +
       (int)(d_checkpointTimeStepInterval > 0) +
       (int)(d_checkpointWallTimeInterval > 0)) > 2) {
    throw ProblemSetupException("Use <checkpoint interval=...> or <checkpoint "
                                "timestepInterval=...>, not both",
                                __FILE__,
                                __LINE__);
  }
  // Can't have a WallTimeStart without a WallTimeInterval.
  if (d_checkpointWallTimeStart > 0.0 && d_checkpointWallTimeInterval == 0) {
    throw ProblemSetupException(
      "<checkpoint walltimeStart must have a corresponding walltimeInterval",
      __FILE__,
      __LINE__);
  }

  d_lastTimeStepLocation = "invalid";
  d_isOutputTimeStep     = false;

  // set up the next output and checkpoint time
  d_nextOutputTime     = 0.0;
  d_nextOutputTimeStep = d_outputInitTimeStep ? 0 : 1;

  d_nextCheckpointTime     = d_checkpointInterval;
  d_nextCheckpointTimeStep = d_checkpointTimeStepInterval + 1;
  d_nextCheckpointWallTime = d_checkpointWallTimeStart;

  if (d_checkpointInterval > 0) {
    proc0cout << "Checkpointing:" << std::setw(16) << " Every "
              << d_checkpointInterval << " physical seconds.\n";
  }
  if (d_checkpointTimeStepInterval > 0) {
    proc0cout << "Checkpointing:" << std::setw(16) << " Every "
              << d_checkpointTimeStepInterval << " timesteps.\n";
  }
  if (d_checkpointWallTimeInterval > 0) {
    proc0cout << "Checkpointing:" << std::setw(16) << " Every "
              << d_checkpointWallTimeInterval << " wall clock seconds,"
              << " starting after " << d_checkpointWallTimeStart
              << " seconds.\n";
  }
}

void
DataArchiver::outputProblemSpec(ProblemSpecP& root_ps)
{
  DOUTR(g_DA_dbg, "Doing: DataArchiver::outputProblemSpec");

  // When running in situ the user can modify start variables. These
  // varibales and values are recorded in the index.xml BUT are NOT
  // incorporated into the checkpoint input.xml file. The user must
  // merge them by hand. This flag reminds them to perform this step.
  if (d_simulator->haveModifiedVars()) {

    ProblemSpecP root = root_ps->getRootNode();

    ProblemSpecP is_ps = root->findBlockWithOutAttribute("InSitu");

    if (is_ps == nullptr) {
      is_ps = root->appendChild("InSitu");
    }

    is_ps->appendElement("haveModifiedVars", d_simulator->haveModifiedVars());
  }
}

void
DataArchiver::initializeOutput(const ProblemSpecP& params,
                               [[maybe_unused]] const GridP& grid)
{
  if (d_outputInterval == 0.0 && d_outputTimeStepInterval == 0 &&
      d_checkpointInterval == 0.0 && d_checkpointTimeStepInterval == 0 &&
      d_checkpointWallTimeInterval == 0) {
    return;
  }

  if (getUseLocalFileSystems()) {
    setupLocalFileSystems();
  } else {
    setupSharedFileSystem();
  }

  // Wait for all ranks to finish verifying shared file system....
  Uintah::MPI::Barrier(d_myworld->getComm());

  if (d_writeMeta) {

    saveSVNinfo();

    // Create index.xml:
    string inputname = d_outputDir.getName() + "/input.xml";
    params->output(inputname.c_str());

    std::cout << "Saving original .ups file in UDA...\n";
    Dir ups_location(pathname(params->getFile()));
    ups_location.copy(basename(params->getFile()), d_outputDir);

    createIndexXML(d_outputDir);

    // create checkpoints/index.xml (if we are saving checkpoints)
    if (d_checkpointInterval > 0.0 || d_checkpointTimeStepInterval > 0 ||
        d_checkpointWallTimeInterval > 0) {
      d_checkpointsDir = d_outputDir.createSubdir("checkpoints");
      createIndexXML(d_checkpointsDir);
    }
  } else {
    d_checkpointsDir = d_outputDir.getSubdir("checkpoints");
  }

  // Sync up before every rank can use the base dir.
  Uintah::MPI::Barrier(d_myworld->getComm());

#ifdef HAVE_PIDX
  // StandAlone/restart_merger calls initializeOutput but has no grid.
  if (grid == nullptr) {
    throw InternalError("No grid - can not use PIDX", __FILE__, __LINE__);
  }

  // Setup for PIDX
  if (savingAsPIDX()) {
    if (d_pidx_requested_nth_rank == -1) {
      d_pidx_requested_nth_rank = d_load_balancer->getNthRank();

      if (d_pidx_requested_nth_rank > 1) {
        proc0cout << "Input file requests output to be saved by every "
                  << d_pidx_requested_nth_rank << "th processor.\n"
                  << "  - However, setting output to every processor "
                  << "until a checkpoint is reached." << std::endl;
        d_load_balancer->setNthRank(1);
        d_load_balancer->possiblyDynamicallyReallocate(grid,
                                                       LoadBalancer::REGRID_LB);
        setSaveAsPIDX();
      }
    }
  }
#endif

} // end initializeOutput()

// to be called after problemSetup and initializeOutput get called
void
DataArchiver::restartSetup(Dir& restartFromDir,
                           int startTimeStep,
                           int timestep,
                           double time,
                           bool fromScratch,
                           bool removeOldDir)
{
  d_outputInitTimeStep = false;
  if (d_writeMeta && !fromScratch) {
    // partial copy of dat files
    copyDatFiles(
      restartFromDir, d_outputDir, startTimeStep, timestep, removeOldDir);

    copySection(restartFromDir, d_outputDir, "index.xml", "restarts");
    copySection(restartFromDir, d_outputDir, "index.xml", "variables");
    copySection(restartFromDir, d_outputDir, "index.xml", "globals");

    // partial copy of index.xml and timestep directories and
    // similarly for checkpoints
    copyTimeSteps(
      restartFromDir, d_outputDir, startTimeStep, timestep, removeOldDir);

    Dir checkpointsFromDir = restartFromDir.getSubdir("checkpoints");
    bool areCheckpoints    = true;

    if (time > 0) {
      // the restart_merger doesn't need checkpoints, and calls this with
      // time=0.
      copyTimeSteps(checkpointsFromDir,
                    d_checkpointsDir,
                    startTimeStep,
                    timestep,
                    removeOldDir,
                    areCheckpoints);
      copySection(
        checkpointsFromDir, d_checkpointsDir, "index.xml", "variables");
      copySection(checkpointsFromDir, d_checkpointsDir, "index.xml", "globals");
    }

    if (removeOldDir) {
      // Try to remove the old dir...
      if (!Dir::removeDir(restartFromDir.getName().c_str())) {

        std::cout
          << "WARNING! In DataArchiver.cc::restartSetup(), removeDir() "
             "failed to remove an old checkpoint directory... Running file "
             "system check now.\n";

        // Something strange happened... let's test the filesystem...
        std::stringstream error_stream;
        if (!testFilesystem(
              restartFromDir.getName(), error_stream, Parallel::getMPIRank())) {

          std::cout << error_stream.str();
          std::cout.flush();

          // The file system just gave us some problems...
          printf("WARNING: Filesystem check failed on processor %d\n",
                 Parallel::getMPIRank());
        }

        // Verify that "system works"
        int code = system("echo how_are_you");
        if (code != 0) {
          printf("WARNING: test of system call failed\n");
        }
      }
    }
  } else if (d_writeMeta) { // Just add <restart from = ".." timestep = "..">
                            // tag.
    copySection(restartFromDir, d_outputDir, "index.xml", "restarts");

    string iname          = d_outputDir.getName() + "/index.xml";
    ProblemSpecP indexDoc = loadDocument(iname);

    if (timestep >= 0) {
      addRestartStamp(indexDoc, restartFromDir, timestep);
    }

    indexDoc->output(iname.c_str());
    // indexDoc->releaseDocument();
  }
}

// This is called after problemSetup. It will copy the dat & checkpoint files to
// the new directory. This also removes the global (dat) variables from the
// saveLabels variables
void
DataArchiver::postProcessUdaSetup(Dir& fromDir)
{
  DOUTR(g_DA_dbg, " Doing: DataArchiver::postprocessUdaSetup");

  // copy files
  if (d_writeMeta) {
    d_fromDir = fromDir;
    copyDatFiles(fromDir, d_outputDir, 0, -1, false);
    copySection(fromDir, d_outputDir, "index.xml", "globals");
    proc0cout << "*** Copied dat files to:   " << d_outputDir.getName()
              << std::endl;

    // copy checkpoints
    Dir checkpointsFromDir = fromDir.getSubdir("checkpoints");
    Dir checkpointsToDir   = d_outputDir.getSubdir("checkpoints");
    std::string me         = checkpointsFromDir.getName();
    if (validDir(me)) {
      checkpointsToDir.remove(
        "index.xml",
        false); // this file is created upstream when it shouldn't have
      checkpointsFromDir.copy(d_outputDir);
      proc0cout << "\n*** Copied checkpoints to: " << d_checkpointsDir.getName()
                << std::endl;
      proc0cout << "    Only using 1 processor to copy so this will be slow "
                   "for large checkpoint directories\n"
                << std::endl;
    }

    // copy input.xml.orig if it exists
    std::string there = d_outputDir.getName();
    std::string here  = fromDir.getName() + "/input.xml.orig";
    if (validFile(here)) {
      fromDir.copy(
        "input.xml.orig",
        d_outputDir); // use OS independent copy functions, needed by mira
      proc0cout << "*** Copied input.xml.orig to: " << there << std::endl;
    }

    // copy the original ups file if it exists
    std::vector<std::string> ups;
    fromDir.getFilenamesBySuffix("ups", ups);

    if (ups.size() != 0) {
      fromDir.copy(
        ups[0],
        d_outputDir); // use OS independent copy functions, needed by mira
      proc0cout << "*** Copied ups file (" << ups[0] << ") to: " << there
                << std::endl;
    }
    proc0cout << "\n" << std::endl;
  }

  // removed the global (dat) variables from the saveLabels
  std::string iname     = fromDir.getName() + "/index.xml";
  ProblemSpecP indexDoc = loadDocument(iname);

  ProblemSpecP globals = indexDoc->findBlock("globals");
  if (globals != nullptr) {

    ProblemSpecP variable = globals->findBlock("variable");
    while (variable != nullptr) {
      std::string varname;
      if (!variable->getAttribute("name", varname)) {
        throw InternalError(
          "global variable name attribute not found", __FILE__, __LINE__);
      }

      auto it = d_saveLabelNames.begin();
      while (it != d_saveLabelNames.end()) {
        if ((*it).labelName == varname) {
          it = d_saveLabelNames.erase(it);
        } else {
          it++;
        }
      }
      variable = variable->findNextBlock("variable");
    }
  }

  //  Read in the timestep indices from the restart uda and store them
  //  Use the indices when creating the timestep directories
  ProblemSpecP ts_ps = indexDoc->findBlock("timesteps");
  ProblemSpecP ts    = ts_ps->findBlock("timestep");
  int timestep       = -9;
  int count          = 1;

  while (ts != nullptr) {
    ts->get(timestep);
    d_restartTimeStepIndices[count] = timestep;

    ts = ts->findNextBlock("timestep");
    count++;
  }

  d_restartTimeStepIndices[0] = d_restartTimeStepIndices[1];

  // Set checkpoint outputIntervals
  d_checkpointInterval         = 0.0;
  d_checkpointTimeStepInterval = 0;
  d_checkpointWallTimeInterval = 0;
  d_nextCheckpointTimeStep     = SHRT_MAX;

  // output every timestep -- each timestep is transferring data
  d_outputInitTimeStep     = true;
  d_outputInterval         = 0.0;
  d_outputTimeStepInterval = 1;
  d_doPostProcessUda       = true;
}

void
DataArchiver::copySection(Dir& fromDir,
                          Dir& toDir,
                          string filename,
                          string section)
{
  // copy chunk labeled section between index.xml files
  std::string iname     = fromDir.getName() + "/" + filename;
  ProblemSpecP indexDoc = loadDocument(iname);

  iname                   = toDir.getName() + "/" + filename;
  ProblemSpecP myIndexDoc = loadDocument(iname);

  ProblemSpecP sectionNode = indexDoc->findBlock(section);
  if (sectionNode != nullptr) {
    ProblemSpecP newNode = myIndexDoc->importNode(sectionNode, true);

    // replace whatever was in the section previously
    ProblemSpecP mySectionNode = myIndexDoc->findBlock(section);
    if (mySectionNode != nullptr) {
      myIndexDoc->replaceChild(newNode, mySectionNode);
    } else {
      myIndexDoc->appendChild(newNode);
    }
  }

  myIndexDoc->output(iname.c_str());
}

void
DataArchiver::addRestartStamp(ProblemSpecP indexDoc, Dir& fromDir, int timestep)
{
  // add restart history to restarts section
  ProblemSpecP restarts = indexDoc->findBlock("restarts");
  if (restarts == nullptr) {
    restarts = indexDoc->appendChild("restarts");
  }

  // restart from <dir> at timestep
  ProblemSpecP restartInfo = restarts->appendChild("restart");
  restartInfo->setAttribute("from", fromDir.getName().c_str());

  std::ostringstream timestep_str;
  timestep_str << timestep;

  restartInfo->setAttribute("timestep", timestep_str.str().c_str());
}

void
DataArchiver::copyTimeSteps(Dir& fromDir,
                            Dir& toDir,
                            int startTimeStep,
                            int maxTimeStep,
                            bool removeOld,
                            bool areCheckpoints /*=false*/)
{
  string old_iname         = fromDir.getName() + "/index.xml";
  ProblemSpecP oldIndexDoc = loadDocument(old_iname);
  string iname             = toDir.getName() + "/index.xml";
  ProblemSpecP indexDoc    = loadDocument(iname);

  ProblemSpecP oldTimeSteps = oldIndexDoc->findBlock("timesteps");

  ProblemSpecP ts;
  if (oldTimeSteps != nullptr) {
    ts = oldTimeSteps->findBlock("timestep");
  }

  // while we're at it, add restart information to index.xml
  if (maxTimeStep >= 0) {
    addRestartStamp(indexDoc, fromDir, maxTimeStep);
  }

  // create timesteps element if necessary
  ProblemSpecP timesteps = indexDoc->findBlock("timesteps");
  if (timesteps == nullptr) {
    timesteps = indexDoc->appendChild("timesteps");
  }

  // copy each timestep
  int timestep;
  while (ts != nullptr) {
    ts->get(timestep);
    if (timestep >= startTimeStep &&
        (timestep <= maxTimeStep || maxTimeStep < 0)) {
      // copy the timestep directory over
      std::map<std::string, std::string> attributes;
      ts->getAttributes(attributes);

      std::string hrefNode = attributes["href"];
      if (hrefNode == "") {
        throw InternalError(
          "timestep href attribute not found", __FILE__, __LINE__);
      }

      std::string::size_type href_pos = hrefNode.find_first_of("/");

      std::string href = hrefNode;
      if (href_pos != std::string::npos) {
        href = hrefNode.substr(0, href_pos);
      }

      // copy timestep directory
      Dir timestepDir = fromDir.getSubdir(href);
      if (removeOld) {
        timestepDir.move(toDir);
      } else {
        timestepDir.copy(toDir);
      }

      if (areCheckpoints) {
        d_checkpointTimeStepDirs.push_back(toDir.getSubdir(href).getName());
      }

      // add the timestep to the index.xml
      std::ostringstream timestep_str;
      timestep_str << timestep;
      ProblemSpecP newTS =
        timesteps->appendElement("timestep", timestep_str.str().c_str());

      for (auto [node, attr] : attributes) {
        newTS->setAttribute(node, attr);
      }
    }
    ts = ts->findNextBlock("timestep");
  }

  // re-output index.xml
  indexDoc->output(iname.c_str());
  // indexDoc->releaseDocument();

  // we don't need the old document anymore...
  // oldIndexDoc->releaseDocument();
}

void
DataArchiver::copyDatFiles(Dir& fromDir,
                           Dir& toDir,
                           int startTimeStep,
                           int maxTimeStep,
                           bool removeOld)
{
  char buffer[1000];

  // find the dat file via the globals block in index.xml
  std::string iname     = fromDir.getName() + "/index.xml";
  ProblemSpecP indexDoc = loadDocument(iname);

  ProblemSpecP globals = indexDoc->findBlock("globals");
  if (globals != nullptr) {
    ProblemSpecP variable = globals->findBlock("variable");
    // copy data file associated with each variable
    while (variable != nullptr) {
      std::map<std::string, std::string> attributes;
      variable->getAttributes(attributes);

      std::string hrefNode = attributes["href"];

      if (hrefNode == "") {
        throw InternalError(
          "global variable href attribute not found", __FILE__, __LINE__);
      }
      const char* href = hrefNode.c_str();

      std::ifstream datFile((fromDir.getName() + "/" + href).c_str());
      if (!datFile) {
        throw InternalError("DataArchiver::copyDatFiles(): The file \"" +
                              (fromDir.getName() + "/" + href) +
                              "\" could not be opened for reading!",
                            __FILE__,
                            __LINE__);
      }

      std::ofstream copyDatFile((toDir.getName() + "/" + href).c_str(),
                                std::ios::app);
      if (!copyDatFile) {
        throw InternalError("DataArchiver::copyDatFiles(): The file \"" +
                              (toDir.getName() + "/" + href) +
                              "\" could not be opened for writing!",
                            __FILE__,
                            __LINE__);
      }

      // copy up to maxTimeStep lines of the old dat file to the copy
      int timestep = startTimeStep;
      while (datFile.getline(buffer, 1000) &&
             (timestep < maxTimeStep || maxTimeStep < 0)) {
        copyDatFile << buffer << std::endl;
        timestep++;
      }
      datFile.close();

      if (removeOld) {
        fromDir.remove(href, false);
      }

      variable = variable->findNextBlock("variable");
    }
  }
  // indexDoc->releaseDocument();
}

void
DataArchiver::createIndexXML(Dir& dir)
{
  ProblemSpecP rootElem = ProblemSpec::createDocument("Uintah_DataArchive");

  rootElem->appendElement("numberOfProcessors", d_myworld->nRanks());
  rootElem->appendElement("ParticlePosition", d_particlePositionName);

  std::string format = "UDA";
  if (d_outputFileFormat == PIDX) {
    format = "PIDX";
  }
  rootElem->appendElement("outputFormat", format);

  ProblemSpecP metaElem = rootElem->appendChild("Meta");

  // Some systems dont supply a logname
  const char* logname = getenv("LOGNAME");
  if (logname) {
    metaElem->appendElement("username", logname);
  } else {
    metaElem->appendElement("username", "unknown");
  }

  time_t t = time(nullptr);

  // Chop the newline character off the time string so that the Date
  // field will appear properly in the XML
  std::string time_string(ctime(&t));
  Uintah::collapse(time_string);

  metaElem->appendElement("date", time_string.c_str());
  metaElem->appendElement("endianness", endianness().c_str());
  metaElem->appendElement("nBits", (int)sizeof(unsigned long) * 8);

  string iname = dir.getName() + "/index.xml";
  rootElem->output(iname.c_str());
  // rootElem->releaseDocument();
}

void
DataArchiver::finalizeTimeStep(const GridP& grid,
                               SchedulerP& sched,
                               [[maybe_unused]] bool recompile /*=false*/,
                               int addMaterial /*=0*/)
{
  // this function should get called exactly once per timestep
  const double delT = d_simulator->getDelT();

  DOUTR(g_DA_dbg,
        "DataArchiver::finalizeTimeStep,"
          << std::setw(10) << " timestep:" << d_simulator->getTimeStep()
          << ", sim time:" << d_simulator->getSimTime() << ", delT: " << delT);

  beginOutputTimeStep(grid);

  if (((delT != 0 || d_outputInitTimeStep) &&
       !s_wereSavesAndCheckpointsInitialized) ||
      addMaterial != 0 || d_switchState ||
      grid->numLevels() != d_numLevelsInOutput) {

    if (d_outputInterval != 0.0 || d_outputTimeStepInterval != 0) {
      initSaveLabels(sched, delT == 0);

      if (!s_wereSavesAndCheckpointsInitialized && delT != 0) {
        indexAddGlobals(); /* add saved global (reduction)
                              variables to index.xml */
      }
    }

    // This assumes that the TaskGraph doesn't change after the second
    // timestep and will need to change if the TaskGraph becomes dynamic.
    //   We also need to do this again if this is the init timestep
    if (delT != 0) {
      s_wereSavesAndCheckpointsInitialized = true;

      // can't do checkpoints on init timestep....
      if (d_checkpointInterval != 0.0 || d_checkpointTimeStepInterval != 0 ||
          d_checkpointWallTimeInterval != 0) {

        initCheckpoints(sched);
      }
    }
  }

  d_numLevelsInOutput = grid->numLevels();

#if SCI_ASSERTION_LEVEL >= 2
  d_outputCalled.clear();
  d_outputCalled.resize(d_numLevelsInOutput, false);
  d_checkpointCalled.clear();
  d_checkpointCalled.resize(d_numLevelsInOutput, false);
  d_checkpointGlobalCalled = false;
#endif
}

//  Schedule output tasks for the grid variables, particle variables and
//  reduction variables
void
DataArchiver::sched_allOutputTasks(const GridP& grid,
                                   SchedulerP& sched,
                                   bool recompile /*=false*/)
{
  printSchedule(nullptr, g_DA_dbg, "DataArchiver::sched_allOutputTasks");

  // we don't want to schedule more tasks unless we're recompiling
  if (!recompile) {
    return;
  }

  const double delT = d_simulator->getDelT();

  //  Schedule output for non-checkpoint, grid, particle, reduction/sole
  //  variables
  if ((d_outputInterval > 0.0 || d_outputTimeStepInterval > 0) &&
      (delT != 0 || d_outputInitTimeStep)) {

    Task* task = scinew Task(
      "DataArchiver::outputGlobalVars", this, &DataArchiver::outputGlobalVars);

    for (auto& saveItem : d_saveGlobalLabels) {
      const MaterialSubset* mss = saveItem.getMaterialSubset(nullptr);
      task->requires(Task::NewDW, saveItem.label, mss, Task::SearchTG::OldTG);
    }

    task->setType(Task::OutputGlobalVars);
    sched->addTask(task, nullptr, nullptr);

    // schedule task for non-global variables
    sched_outputVariables(d_saveLabels, grid, sched, false);
  }

  // Schedule Checkpoint (reduction variables)
  if (delT != 0 && d_checkpointCycle > 0 &&
      (d_checkpointInterval > 0 || d_checkpointTimeStepInterval > 0 ||
       d_checkpointWallTimeInterval > 0)) {

    printSchedule(nullptr,
                  g_DA_dbg,
                  "DataArchiver::sched_outputVariables (CheckpointGlobal)");

    // output checkpoint timestep
    Task* task =
      scinew Task("DataArchiver::outputVariables (CheckpointReduction)",
                  this,
                  &DataArchiver::outputVariables,
                  CHECKPOINT_GLOBAL);

    for (auto& saveItem : d_checkpointReductionLabels) {
      const MaterialSubset* mss = saveItem.getMaterialSubset(nullptr);
      task->requires(Task::NewDW, saveItem.label, mss, Task::SearchTG::OldTG);
    }

    sched->addTask(task, nullptr, nullptr);

    sched_outputVariables(d_checkpointLabels, grid, sched, true);
  }

#if HAVE_PIDX
  if (d_outputFileFormat == PIDX) {
    /*  Create PIDX communicators (one communicator per AMR level) */

    // FIXME: It appears that this is called 3 or more times before
    //        timestep 0... why is this the case?  Doesn't hurt
    //        anything, but we should figure out why...

    // ARS - It can get called multiple times when regridding.

    DOUTR(g_DA_dbg,
          "  Creating communicatore per AMR level (required for PIDX)\n");
    createPIDXCommunicator(d_checkpointLabels, grid, sched, true);
  }
#endif
}

void
DataArchiver::setOutputTimeStep(bool val, const GridP& grid)
{
  if (d_isOutputTimeStep != val) {
    d_isOutputTimeStep = val;

    // Create the output timestep directories
    if (d_isOutputTimeStep && d_outputFileFormat != PIDX) {
      makeTimeStepDirs(
        d_outputDir, d_saveLabels, grid, &d_lastTimeStepLocation);

      d_outputTimeStepDirs.push_back(d_lastTimeStepLocation);
    }
  }
}

void
DataArchiver::setCheckpointTimeStep(bool val, const GridP& grid)
{
  if (d_isCheckpointTimeStep != val) {
    d_isCheckpointTimeStep = val;

    // Create the output checkpoint directories
    if (d_isCheckpointTimeStep) {
      std::string timestepDir;
      makeTimeStepDirs(
        d_checkpointsDir, d_checkpointLabels, grid, &timestepDir);
      d_checkpointTimeStepDirs.push_back(timestepDir);

      std::string iname = d_checkpointsDir.getName() + "/index.xml";

      ProblemSpecP index;

      if (d_writeMeta) {
        index = loadDocument(iname);

        // store a back up in case it dies while writing index.xml
        std::string ibackup_name =
          d_checkpointsDir.getName() + "/index_backup.xml";
        index->output(ibackup_name.c_str());
      }

      if (d_checkpointCycle > 0 &&
          (int)d_checkpointTimeStepDirs.size() > d_checkpointCycle) {

        if (d_writeMeta) {
          // Remove reference to outdated checkpoint directory from the
          // checkpoint index.
          ProblemSpecP ts   = index->findBlock("timesteps");
          ProblemSpecP temp = ts->getFirstChild();
          ts->removeChild(temp);

          index->output(iname.c_str());

          // remove out-dated checkpoint directory
          Dir expiredDir(d_checkpointTimeStepDirs.front());

          // Try to remove the expired checkpoint directory...
          if (!Dir::removeDir(expiredDir.getName().c_str())) {

            // Something strange happened... let's test the filesystem...
            std::cout << "\nWarning! removeDir() Failed for '"
                      << expiredDir.getName()
                      << "' in DataArchiver.cc::beginOutputTimeStep()\n\n";
            std::stringstream error_stream;

            if (!testFilesystem(
                  expiredDir.getName(), error_stream, Parallel::getMPIRank())) {
              std::cout << error_stream.str();
              std::cout.flush();
              // The file system just gave us some problems...
              printf("WARNING: Filesystem check failed on processor %d\n",
                     Parallel::getMPIRank());
            }
          }
        }
        d_checkpointTimeStepDirs.pop_front();
      }

      if (d_myworld->myRank() == 0) {
        if (d_checkpointCycle == 0) {
          DOUTR(true,
                "WARNING the checkpoint cycle is set to zero. "
                "No checkpoints will be deleted. This may cause unacceptable "
                "disk usage.");
        }

        if ((int)d_checkpointTimeStepDirs.size() > 10) {
          DOUTR(
            true,
            "WARNING there are currently checkpoint "
              << d_checkpointTimeStepDirs.size() << " files. "
              << "This may be excessive and cause unacceptable disk usage.");
        }
      }

      // if (d_writeMeta)
      // index->releaseDocument();
    }
  }
}

void
DataArchiver::beginOutputTimeStep(const GridP& grid)
{
  const int timestep   = d_simulator->getTimeStep();
  const double simTime = d_simulator->getSimTime();
  const double delT    = d_simulator->getDelT();

  d_isOutputTimeStep     = false;
  d_isCheckpointTimeStep = false;

  // do *not* update d_nextOutputTime or others here.  We need the original
  // values to compare if there is a timestep restart.  See
  // reEvaluateOutputTimeStep

  // Check for an output
  bool isOutputTimeStep =
    // Output based on the simulation time.
    (((d_outputInterval > 0.0 && (delT != 0.0 || d_outputInitTimeStep)) &&
      (simTime + delT >= d_nextOutputTime)) ||

     // Output based on the timestep interval.
     ((d_outputTimeStepInterval > 0 && (delT != 0.0 || d_outputInitTimeStep)) &&
      (timestep >= d_nextOutputTimeStep)) ||

     // Output based on the being the last timestep.
     (d_outputLastTimeStep && d_maybeLastTimeStep));

  setOutputTimeStep(isOutputTimeStep, grid);

  // Check for a checkpoint.
  bool isCheckpointTimeStep =
    // Checkpoint based on the simulation time.
    ((d_checkpointInterval > 0.0 && (simTime + delT) >= d_nextCheckpointTime) ||

     // Checkpoint based on the timestep interval.
     (d_checkpointTimeStepInterval > 0 &&
      timestep >= d_nextCheckpointTimeStep) ||

     // Checkpoint based on the being the last timestep.
     (d_checkpointLastTimeStep && d_maybeLastTimeStep));

  // Checkpoint based on the wall time.
  if (d_checkpointWallTimeInterval > 0) {

    if (d_elapsedWallTime >= d_nextCheckpointWallTime) {
      isCheckpointTimeStep = true;
    }
  }

  setCheckpointTimeStep(isCheckpointTimeStep, grid);

  DOUTR(g_DA_dbg,
        "DataArchiver::beginOutputTimeStep isOutputTimeStep: "
          << isOutputTimeStep << " "
          << " isCheckpointTimeStep: " << isCheckpointTimeStep);
} // end beginOutputTimeStep

void
DataArchiver::makeTimeStepDirs(
  Dir& baseDir,
  [[maybe_unused]] std::vector<DataArchiver::SaveItem>& saveLabels,
  const GridP& grid,
  string* pTimeStepDir /* passed back */)
{
  int numLevels    = grid->numLevels();
  int dir_timestep = getTimeStepTopLevel();

  DOUTR(g_DA_dbg,
        "DataArchiver::makeTimeStepDirs timestep: "
          << d_simulator->getTimeStep() << " dir_timestep: " << dir_timestep);

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0') << dir_timestep;
  *pTimeStepDir = baseDir.getName() + "/" + tname.str();

  // Create the directory for this timestep, if necessary
  // It is not gurantteed that the rank holding d_writeMeta will call
  // outputTimstep to create dir before another rank begin to output data.
  // A race condition happens when a rank executes output task and
  // the rank holding d_writeMeta is still compiling task graph.
  // So every rank should try to create dir.
  Dir tdir = baseDir.createSubdirPlus(tname.str());

  // Create the directory for this level, if necessary
  for (int l = 0; l < numLevels; l++) {
    std::ostringstream lname;
    lname << "l" << l;
    Dir ldir = tdir.createSubdirPlus(lname.str());
  }
}

void
DataArchiver::findNext_OutputCheckPointTimeStep(
  const bool restart,
  [[maybe_unused]] const GridP& grid)
{
  const int timeStep   = d_simulator->getTimeStep();
  const double simTime = d_simulator->getSimTime();

#ifdef HAVE_PIDX
  const double delT = d_simulator->getNextDelT();
#endif

  if (restart) {

    // Output based on the simulaiton time.
    if (d_outputInterval > 0.0) {
      d_nextOutputTime =
        std::ceil(simTime / d_outputInterval) * d_outputInterval;
    }
    // Output based on the time step.
    else if (d_outputTimeStepInterval > 0) {
      d_nextOutputTimeStep =
        (timeStep / d_outputTimeStepInterval) * d_outputTimeStepInterval + 1;

      while (d_nextOutputTimeStep <= timeStep) {
        d_nextOutputTimeStep += d_outputTimeStepInterval;
      }
    }

    // Checkpoint based on the simulaiton time.
    if (d_checkpointInterval > 0.0) {
      d_nextCheckpointTime =
        std::ceil(simTime / d_checkpointInterval) * d_checkpointInterval;
    }
    // Checkpoint based on the time step.
    else if (d_checkpointTimeStepInterval > 0) {
      d_nextCheckpointTimeStep = (timeStep / d_checkpointTimeStepInterval) *
                                   d_checkpointTimeStepInterval +
                                 1;
      while (d_nextCheckpointTimeStep <= timeStep) {
        d_nextCheckpointTimeStep += d_checkpointTimeStepInterval;
      }
    }
    // Checkpoint based on the wall time.
    else if (d_checkpointWallTimeInterval > 0) {
      d_nextCheckpointWallTime =
        std::ceil(d_elapsedWallTime / d_checkpointWallTimeInterval) *
        d_checkpointWallTimeInterval;
    }
  }

  // If this timestep was an output/checkpoint timestep, determine
  // when the next one will be.

  // Do *not* do this step in beginOutputTimeStep because the original
  // values are needed to compare with if there is a time step recompute.
  // See reEvaluateOutputTimeStep

  // When outputing/checkpointing using the simulation or wall time
  // check to see if the simulation or wall time went past more than
  // one interval. If so adjust accordingly.

  // Note - it is not clear why but when outputing/checkpointing using
  // time steps the mod function must also be used. This does not
  // affect most simulations except when there are multiple UPS files
  // such as when components are switched. For example:
  // StandAlone/inputs/UCF/Switcher/switchExample3.ups

  else if (d_isOutputTimeStep) {

    // Output based on the simulaiton time.
    if (d_outputInterval > 0.0) {
      if (simTime >= d_nextOutputTime) {
        d_nextOutputTime +=
          std::floor((simTime - d_nextOutputTime) / d_outputInterval) *
            d_outputInterval +
          d_outputInterval;
      }
    }
    // Output based on the time step.
    else if (d_outputTimeStepInterval > 0) {
      if (timeStep >= d_nextOutputTimeStep) {
        d_nextOutputTimeStep +=
          ((timeStep - d_nextOutputTimeStep) / d_outputTimeStepInterval) *
            d_outputTimeStepInterval +
          d_outputTimeStepInterval;
      }
    }
  }

  if (d_isCheckpointTimeStep) {
    // Checkpoint based on the simulaiton time.
    if (d_checkpointInterval > 0.0) {
      if (simTime >= d_nextCheckpointTime) {
        d_nextCheckpointTime +=
          std::floor((simTime - d_nextCheckpointTime) / d_checkpointInterval) *
            d_checkpointInterval +
          d_checkpointInterval;
      }
    }
    // Checkpoint based on the time step.
    else if (d_checkpointTimeStepInterval > 0) {
      if (timeStep >= d_nextCheckpointTimeStep) {
        d_nextCheckpointTimeStep += ((timeStep - d_nextCheckpointTimeStep) /
                                     d_checkpointTimeStepInterval) *
                                      d_checkpointTimeStepInterval +
                                    d_checkpointTimeStepInterval;
      }
    }

    // Checkpoint based on the wall time.
    else if (d_checkpointWallTimeInterval > 0) {

      if (d_elapsedWallTime >= d_nextCheckpointWallTime) {
        d_nextCheckpointWallTime +=
          std::floor((d_elapsedWallTime - d_nextCheckpointWallTime) /
                     d_checkpointWallTimeInterval) *
            d_checkpointWallTimeInterval +
          d_checkpointWallTimeInterval;
      }
    }
  }

  // When saving via PIDX one needs to predict if one is going to do a
  // checkpoint. These are the exact same checks as done in
  // beginOutputTimeStep
#ifdef HAVE_PIDX
  if (savingAsPIDX()) {

    // Check for a checkpoint.
    d_pidx_checkpointing =
      // Checkpoint based on the simulation time.
      ((d_checkpointInterval > 0.0 &&
        (simTime + delT) >= d_nextCheckpointTime) ||

       // Checkpoint based on the timestep interval.
       (d_checkpointTimeStepInterval > 0 &&
        timeStep >= d_nextCheckpointTimeStep) ||

       // Checkpoint based on the being the last timestep.
       (d_checkpointLastTimeStep && d_maybeLastTimeStep));

    // Checkpoint based on the wall time.
    if (d_checkpointWallTimeInterval > 0) {

      if (d_elapsedWallTime >= d_nextCheckpointWallTime) {
        d_pidx_checkpointing = true;
      }
    }

    // Checkpointing
    if (d_pidx_checkpointing) {

      if (d_pidx_requested_nth_rank > 1) {
        proc0cout << "This is a checkpoint time step (" << timeStep
                  << ") - need to recompile with nth proc set to: "
                  << d_pidx_requested_nth_rank << "\n";

        d_load_balancer->setNthRank(d_pidx_requested_nth_rank);
        d_load_balancer->possiblyDynamicallyReallocate(grid,
                                                       LoadBalancer::REGRID_LB);
        setSaveAsUDA();
        d_pidx_need_to_recompile = true;
      }
    }

    // Check for an output which may need to be postponed

    // Output based on the simulation time.
    if (((d_outputInterval > 0.0 && (delT != 0.0 || d_outputInitTimeStep)) &&
         (simTime + delT >= d_nextOutputTime)) ||

        // Output based on the timestep interval.
        ((d_outputTimeStepInterval > 0 &&
          (delT != 0.0 || d_outputInitTimeStep)) &&
         (timeStep >= d_nextOutputTimeStep)) ||

        // Output based on the being the last timestep.
        (d_outputLastTimeStep && d_maybeLastTimeStep)) {

      proc0cout << "This is an output time step: " << timeStep << ".  ";

      // If this is also a checkpoint time step postpone the output
      if (d_pidx_need_to_recompile) {
        postponeNextOutputTimeStep();

        proc0cout << "   Postposing as it is also a checkpoint time step.";
      }

      proc0cout << "\n";
    }
  }
#endif

  DOUTR(g_DA_dbg,
        "DataArchiver::findNext_OutputCheckPointTimeStep  "
          << std::setw(20) << "isOutputTimeStep: " << d_isOutputTimeStep
          << std::setw(20) << ", outputInterval: " << d_outputInterval
          << std::setw(20) << ", nextOutputTime: " << d_nextOutputTime
          << std::setw(20)
          << ", outputTimeStepInterval: " << d_outputTimeStepInterval
          << std::setw(20) << ", nextOutputTimeStep: " << d_nextOutputTimeStep);

  DOUTR(g_DA_dbg,
        "DataArchiver::findNext_OutputCheckPointTimeStep  "
          << std::setw(20) << "isCheckpointTimeStep: " << d_isCheckpointTimeStep
          << std::setw(20) << ", checkpointInterval: " << d_checkpointInterval
          << std::setw(20) << ", nextCheckpointTime: " << d_nextCheckpointTime
          << std::setw(20)
          << ", checkpointTimeStepInterval: " << d_checkpointTimeStepInterval
          << ", nextCheckpointTimeStep: " << d_nextCheckpointTimeStep
          << ", nextCheckpointWallTime: " << d_nextCheckpointWallTime);

} // end findNext_OutputCheckPoint_TimeStep()

void
DataArchiver::recompute_OutputCheckPointTimeStep()
{
  const double simTime = d_simulator->getSimTime();
  const double delT    = d_simulator->getDelT();

  // Called after a timestep recompute. If the new deltT goes
  // beneath the threshold, cancel the output and/or checkpoint
  // timestep.
  if (d_isOutputTimeStep && d_outputInterval > 0.0) {
    if (simTime + delT < d_nextOutputTime) {
      d_isOutputTimeStep = false;
    }
  }

  if (d_isCheckpointTimeStep && d_checkpointInterval > 0.0) {
    if (simTime + delT < d_nextCheckpointTime) {
      d_isCheckpointTimeStep = false;
      d_checkpointTimeStepDirs.pop_back();
    }
  }

#if SCI_ASSERTION_LEVEL >= 2
  d_outputCalled.clear();
  d_outputCalled.resize(d_numLevelsInOutput, false);
  d_checkpointCalled.clear();
  d_checkpointCalled.resize(d_numLevelsInOutput, false);
  d_checkpointGlobalCalled = false;
#endif

  DOUTR(g_DA_dbg,
        "DataArchiver::recompute_OutputCheckPointTimeStep() "
          << " isOutputTimeStep: " << d_isOutputTimeStep
          << " isCheckpointTimeStep: " << d_isCheckpointTimeStep);
}

//  Update the xml files: index.xml, timestep.xml
void
DataArchiver::writeto_xml_files(const GridP& grid)
{
  if (!d_isCheckpointTimeStep && !d_isOutputTimeStep) {
    DOUTR(
      g_DA_dbg,
      "DataArchiver::writeto_xml_files:  Not an output or checkpoint timestep");
    return;
  }

  double simTime = d_simulator->getSimTime();
  double delT    = d_simulator->getDelT();

  if (d_outputPreviousTimeStep || d_checkpointPreviousTimeStep) {
    if (d_simulator->getLastRegridTimeStep() > getTimeStepTopLevel()) {
      DOUT(d_myworld->myRank() == 0,
           "WARNING : Requesting the previous output/checkpoint for time step "
             << getTimeStepTopLevel() << " but the grid changed on time step "
             << d_simulator->getLastRegridTimeStep()
             << ". Not writing the associated XML files.");
      return;
    }

    simTime -= delT;
  }

  Timers::Simple timer;
  timer.start();

  DOUTR(g_DA_dbg, "DataArchiver::writeto_xml_files:V1");

  //__________________________________
  //  Writeto XML files
  // to check for output nth proc
  int dir_timestep = getTimeStepTopLevel();

  // start dumping files to disk
  std::vector<Dir*> baseDirs;
  if (d_isOutputTimeStep) {
    baseDirs.push_back(&d_outputDir);
  }
  if (d_isCheckpointTimeStep) {
    baseDirs.push_back(&d_checkpointsDir);
  }

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0') << dir_timestep;

  // Loop over the directories we are going to save data in.  baseDirs
  // contains either the base vis dump directory, or the base
  // checkpoint directory, or both (if we doing both a checkpoint and
  // a vis dump on the same timestep).
  for (auto& baseDir : baseDirs) {

    // to save the list of vars. up to 2, since in checkpoints, there are two
    // types of vars
    std::vector<std::vector<SaveItem>*> savelist;

    // Reference this timestep in index.xml
    if (d_writeMeta) {

      bool hasGlobals                = false;
      bool dumpingCheckpoint         = false;
      bool save_io_timestep_xml_file = false;

      if (baseDir == &d_outputDir) {
        savelist.push_back(&d_saveLabels);

        if (d_outputFileFormat == PIDX) {
#ifdef HAVE_PIDX
          // When using PIDX, only save "timestep.xml" if the grid had
          // changed... (otherwise we will just point (symlink) to the
          // last one saved.)
          save_io_timestep_xml_file =
            (d_lastOutputOfTimeStepXML < d_simulator->getLastRegridTimeStep());
#else
          throw InternalError(
            "DataArchiver::writeto_xml_files(): PIDX not configured!",
            __FILE__,
            __LINE__);
#endif
        } else {
          save_io_timestep_xml_file = true; // Always save for a legacy UDA
        }
      } else if (baseDir == &d_checkpointsDir) {
        dumpingCheckpoint = true;
        hasGlobals        = (d_checkpointReductionLabels.size() > 0);
        savelist.push_back(&d_checkpointLabels);
        savelist.push_back(&d_checkpointReductionLabels);
      } else {
        throw "DataArchiver::writeto_xml_files(): Unknown directory!";
      }

      std::string iname     = baseDir->getName() + "/index.xml";
      ProblemSpecP indexDoc = loadDocument(iname);

      // if this timestep isn't already in index.xml, add it in
      if (indexDoc == nullptr) {
        continue; // output timestep but no variables scheduled to be saved.
      }
      ASSERT(indexDoc != nullptr);

      // output data pointers
      for (auto& j : savelist) {
        string variableSection =
          j == &d_checkpointReductionLabels ? "globals" : "variables";

        ProblemSpecP vs = indexDoc->findBlock(variableSection);
        if (vs == nullptr) {
          vs = indexDoc->appendChild(variableSection.c_str());
        }

        for (size_t k = 0; k < j->size(); k++) {
          const VarLabel* var = (*j)[k].label;
          bool found          = false;

          for (ProblemSpecP n = vs->getFirstChild(); n != nullptr;
               n              = n->getNextSibling()) {
            if (n->getNodeName() == "variable") {
              std::map<std::string, std::string> attributes;
              n->getAttributes(attributes);
              std::string varname = attributes["name"];
              if (varname == "") {
                throw InternalError("varname not found", __FILE__, __LINE__);
              }
              if (varname == var->getName()) {
                found = true;
                break;
              }
            }
          }

          if (!found) {
            ProblemSpecP newElem = vs->appendChild("variable");
            newElem->setAttribute(
              "type",
              TranslateVariableType(var->typeDescription()->getName(),
                                    baseDir != &d_outputDir));
            newElem->setAttribute("name", var->getName());

            // Save number of materials associated with this variable...
            SaveItem& saveItem          = (*j)[k];
            const MaterialSubset* matls = saveItem.getMaterialSubset(nullptr);

            newElem->setAttribute("numMaterials",
                                  std::to_string(matls->size()));
          }
        }
      }

      // Check if it's the first checkpoint timestep by checking if the
      // "timesteps" field is in checkpoints/index.xml.  If it is then there
      // exists a timestep.xml file already. Use this below to change
      // information in input.xml...
      bool firstCheckpointTimeStep = false;

      ProblemSpecP ts = indexDoc->findBlock("timesteps");
      if (ts == nullptr) {
        ts                      = indexDoc->appendChild("timesteps");
        firstCheckpointTimeStep = (&d_checkpointsDir == baseDir);
      }

      bool found = false;
      for (ProblemSpecP n = ts->getFirstChild(); n != 0;
           n              = n->getNextSibling()) {
        if (n->getNodeName() == "timestep") {
          int readtimestep;
          if (!n->get(readtimestep)) {
            throw InternalError(
              "Error parsing timestep number", __FILE__, __LINE__);
          }
          if (readtimestep == dir_timestep) {
            found = true;
            break;
          }
        }
      }

      if (!found) {
        // add timestep info
        std::string timestepindex = tname.str() + "/timestep.xml";

        std::ostringstream value, timeVal, deltVal;
        value << dir_timestep;
        ProblemSpecP newElem =
          ts->appendElement("timestep", value.str().c_str());
        newElem->setAttribute("href", timestepindex.c_str());
        timeVal << std::setprecision(17) << simTime;
        newElem->setAttribute("time", timeVal.str());
        deltVal << std::setprecision(17) << delT;
        newElem->setAttribute("oldDelt", deltVal.str());
      }

      indexDoc->output(iname.c_str());
      // indexDoc->releaseDocument();

      // Make a timestep.xml file for this time step we need to do it
      // here in case there is a time steps rescompute. Break out the
      // <Grid> and <Data> section of the DOM tree into a separate
      // grid.xml file which can be created quickly and use less
      // memory using the xmlTextWriter functions (streaming output)

      // If grid has changed (or this is a checkpoint), save the
      // timestep.xml data, otherwise just point to the previous
      // version.
      if (save_io_timestep_xml_file || dumpingCheckpoint) {

        // Create the timestep.xml file.
        ProblemSpecP rootElem = ProblemSpec::createDocument("Uintah_timestep");

        // Create a metadata element to store the per-timestep endianness
        ProblemSpecP metaElem = rootElem->appendChild("Meta");

        metaElem->appendElement("endianness", endianness().c_str());
        metaElem->appendElement("nBits", (int)sizeof(unsigned long) * 8);
        metaElem->appendElement("numProcs", d_myworld->nRanks());

        // TimeStep information
        ProblemSpecP timeElem = rootElem->appendChild("Time");
        timeElem->appendElement("timestepNumber", dir_timestep);
        timeElem->appendElement("currentTime", simTime); // + delT);
        timeElem->appendElement("oldDelt", delT);

        // Output grid section:
        //
        // With AMR, we're not guaranteed that a rank has work on a
        // given level.  Quick check to see that, so we don't create a
        // node that points to no data.
        std::string grid_path = baseDir->getName() + "/" + tname.str() + "/";

#if XML_TEXTWRITER
        writeGridTextWriter(hasGlobals, grid_path, grid);
#else
        // Original version:
        writeGridOriginal(hasGlobals, grid, rootElem);

        // Binary Grid version:
        // writeGridBinary( hasGlobals, grid_path, grid );
#endif

        // Add the <Materials> section to the timestep.xml
        GeometryPieceFactory::resetGeometryPiecesOutput();

        // output each components output Problem spec
        d_simulator->outputProblemSpec(rootElem);

        outputProblemSpec(rootElem);

        std::string name =
          baseDir->getName() + "/" + tname.str() + "/timestep.xml";
        rootElem->output(name.c_str());

        // output input.xml & input.xml.orig

        // a small convenience to the user who wants to change things when he
        // restarts let him know that some information to change will need to be
        // done in the timestep.xml file instead of the input.xml file.  Only do
        // this once, though.
        if (firstCheckpointTimeStep) {
          // loop over the blocks in timestep.xml and remove them from
          // input.xml, with some exceptions.
          std::string inputname = d_outputDir.getName() + "/input.xml";
          ProblemSpecP inputDoc = loadDocument(inputname);
          inputDoc->output((inputname + ".orig").c_str());

          for (ProblemSpecP ps = rootElem->getFirstChild(); ps != nullptr;
               ps              = ps->getNextSibling()) {
            std::string nodeName = ps->getNodeName();
            if (nodeName == "Meta" || nodeName == "Time" ||
                nodeName == "Grid" || nodeName == "Data") {
              continue;
            }
            ProblemSpecP removeNode = inputDoc->findBlock(nodeName);
            if (removeNode != 0) {
              std::string comment =
                "The node " + nodeName +
                " has been removed.  Its original values are\n"
                "    in input.xml.orig, and values used for "
                "restarts are in checkpoints/t*****/timestep.xml";
              ProblemSpecP commentNode = inputDoc->makeComment(comment);
              inputDoc->replaceChild(removeNode, commentNode);
            }
          }
          inputDoc->output(inputname.c_str());
          // inputDoc->releaseDocument();
        }
        // rootElem->releaseDocument();

        // copy the component sections of timestep.xml.
        if (d_doPostProcessUda) {
          copy_outputProblemSpec(d_fromDir, d_outputDir);
        }
      } else {
        // This is an IO timestep dump and the grid has not changed,
        // so we do not need to create a new timestep.xml file...  We
        // will just point (symlink) to the most recently dumped
        // "timestep.xml" file.
        std::ostringstream ts_with_xml_dirname;
        ts_with_xml_dirname << "../t" << std::setw(5) << std::setfill('0')
                            << d_lastOutputOfTimeStepXML << "/timestep.xml";

        std::ostringstream ts_without_xml_dirname;
        ts_without_xml_dirname << d_outputDir.getName() << "/t" << std::setw(5)
                               << std::setfill('0') << dir_timestep
                               << "/timestep.xml";

        [[maybe_unused]] int res =
          symlink(ts_with_xml_dirname.str().c_str(),
                  ts_without_xml_dirname.str().c_str());
      }

      if (save_io_timestep_xml_file) {
        // Record the fact that we just wrote the "timestep.xml" file
        // for timestep #: dir_timestep.
        d_lastOutputOfTimeStepXML = dir_timestep;
      }
    } // end if m_writeMeta
  }   // end loop over baseDirs

  double myTime = timer().seconds();
  (*d_runtimeStats)[XMLIOTime] += myTime;
  (*d_runtimeStats)[TotalIOTime] += myTime;
}

void
DataArchiver::writeGridBinary(const bool hasGlobals,
                              const string& grid_path,
                              const GridP& grid)
{
  // Originally the <Grid> was saved in XML.  While this makes it very
  // human readable, for large runs (10K+ patches), the (original)
  // timestep.xml file, and now the grid.xml file become huge, taking
  // a significant amount of disk space, and time to generate/save.
  // To fix this problem, we will now write the data out in binary.
  // The original <Grid> data looks like this:
  //
  // <Grid>
  //   <numLevels>1</numLevels>
  //   <Level>
  //    <numPatches>4</numPatches>
  //    <totalCells>8820</totalCells>
  //    <extraCells>[1,1,1]</extraCells>
  //    <anchor>[0,0,0]</anchor>
  //    <id>0</id>
  //    <cellspacing>[0.025000000000000001,0.025000000000000001,0.049999999999999996]</cellspacing>
  //    <Patch>
  //     <id>0</id>
  //     <proc>0</proc>
  //     <lowIndex>[-1,-1,-1]</lowIndex>
  //     <highIndex>[20,20,4]</highIndex>
  //     <interiorLowIndex>[0,0,0]</interiorLowIndex>
  //     <interiorHighIndex>[20,20,3]</interiorHighIndex>
  //     <nnodes>2646</nnodes>
  //     <lower>[-0.025000000000000001,-0.025000000000000001,-0.049999999999999996]</lower>
  //     <upper>[0.5,0.5,0.19999999999999998]</upper>
  //     <totalCells>2205</totalCells>
  //    </Patch>
  //   </Level>
  // </Grid>

  FILE* fp;
  int marker = DataArchive::GRID_MAGIC_NUMBER;

  std::string grid_filename = grid_path + "grid.xml";

  fp = fopen(grid_filename.c_str(), "wb");
  fwrite(&marker, sizeof(int), 1, fp);

  // NUmber of Levels
  int numLevels = grid->numLevels();
  fwrite(&numLevels, sizeof(int), 1, fp);

  std::vector<std::vector<bool>> procOnLevel(numLevels);

  for (int lev = 0; lev < numLevels; lev++) {

    LevelP level = grid->getLevel(lev);

    int num_patches  = level->numPatches();
    long num_cells   = level->totalCells();
    IntVector eciv   = level->getExtraCells();
    int* extra_cells = eciv.get_pointer();
    double anchor[3];
    Point anchor_pt = level->getAnchor();
    int id          = level->getID();
    double cell_spacing[3];
    Vector cell_spacing_vec = level->dCell();
    IntVector pbiv          = level->getPeriodicBoundaries();
    int* period             = pbiv.get_pointer();

    anchor[0] = anchor_pt.x();
    anchor[1] = anchor_pt.y();
    anchor[2] = anchor_pt.z();

    cell_spacing[0] = cell_spacing_vec.x();
    cell_spacing[1] = cell_spacing_vec.y();
    cell_spacing[2] = cell_spacing_vec.z();

    fwrite(&num_patches, sizeof(int), 1, fp);    // Number of Patches -  100
    fwrite(&num_cells, sizeof(long), 1, fp);     // Number of Cells   - 8000
    fwrite(extra_cells, sizeof(int), 3, fp);     // Extra Cell Info - [1,1,1]
    fwrite(anchor, sizeof(double), 3, fp);       // Anchor Info     - [0,0,0]
    fwrite(period, sizeof(int), 3, fp);          //
    fwrite(&id, sizeof(int), 1, fp);             // ID of Level     - 0
    fwrite(cell_spacing, sizeof(double), 3, fp); // Cell Spacing - [0.1,0.1,0.1]

    procOnLevel[lev].resize(d_myworld->nRanks());

    // Iterate over patches.
    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         ++iter) {
      const Patch* patch = *iter;

      int patch_id = patch->getID();
      int rank_id  = d_load_balancer->getOutputProc(patch);

      int proc               = d_load_balancer->getOutputProc(patch);
      procOnLevel[lev][proc] = true;

      IntVector ecliiv = patch->getExtraCellLowIndex();
      IntVector echiiv = patch->getExtraCellHighIndex();
      int* low_index   = ecliiv.get_pointer();
      int* high_index  = echiiv.get_pointer();

      // Interior indices
      IntVector cliiv   = patch->getCellLowIndex();
      IntVector chiiv   = patch->getCellHighIndex();
      int* i_low_index  = cliiv.get_pointer();
      int* i_high_index = chiiv.get_pointer();

      int num_nodes = patch->getNumExtraNodes();

      Box box = patch->getExtraBox();

      double lower[3];
      const Point& lower_pt = box.lower();
      lower[0]              = lower_pt.x();
      lower[1]              = lower_pt.y();
      lower[2]              = lower_pt.z();

      double upper[3];
      const Point& upper_pt = box.upper();
      upper[0]              = upper_pt.x();
      upper[1]              = upper_pt.y();
      upper[2]              = upper_pt.z();

      int num_cells = patch->getNumExtraCells();

      //                                                     Name: Example:
      fwrite(&patch_id, sizeof(int), 1, fp); // Patch ID              - 0
      fwrite(&rank_id, sizeof(int), 1, fp);  // Process ID            - 0
      fwrite(low_index,
             sizeof(int),
             3,
             fp); // Low Index             - [-1, -1, -1]
      fwrite(high_index,
             sizeof(int),
             3,
             fp); // High Index            - [20, 20,  4]
      fwrite(i_low_index,
             sizeof(int),
             3,
             fp); // Interior Low Index    - [ 0,  0,  0]
      fwrite(i_high_index,
             sizeof(int),
             3,
             fp); // Interior High Index   - [20, 20,  3]
      fwrite(&num_nodes, sizeof(int), 1, fp); // Number of Extra Nodes - 2646
      fwrite(lower,
             sizeof(double),
             3,
             fp); // Lower                 - [-0.025, -0.025,  -0.05]
      fwrite(upper,
             sizeof(double),
             3,
             fp); // Upper                 - [ 0.5, 0.5,  0.2]
      fwrite(&num_cells, sizeof(int), 1, fp); // Total number of cells - 2205
    }
  }

  // Write an end of file marker...
  fwrite(&marker, sizeof(int), 1, fp);

  fclose(fp);

  writeDataTextWriter(hasGlobals, grid_path, grid, procOnLevel);

} // end writeGridBinary()

// Creates the <Grid> section of the XML DOM for the output file.
// This original approach places the <Grid> inside the timestep.xml
// file.
//
void
DataArchiver::writeGridOriginal(const bool hasGlobals,
                                const GridP& grid,
                                ProblemSpecP rootElem)
{
  // With AMR, we're not guaranteed that a proc do work on a given
  // level.  Quick check to see that, so we don't create a node that
  // points to no data
  int numLevels = grid->numLevels();
  std::vector<std::vector<bool>> procOnLevel(numLevels);

  // Break out the <Grid> and <Data> sections and write those to a
  // "grid.xml" section using libxml2's TextWriter which is a
  // streaming output format which doesn't use a DOM tree.

  ProblemSpecP gridElem = rootElem->appendChild("Grid");

  //__________________________________
  //  output level information
  gridElem->appendElement("numLevels", numLevels);

  for (int l = 0; l < numLevels; ++l) {
    LevelP level           = grid->getLevel(l);
    ProblemSpecP levelElem = gridElem->appendChild("Level");

    if (level->getPeriodicBoundaries() != IntVector(0, 0, 0)) {
      levelElem->appendElement("periodic", level->getPeriodicBoundaries());
    }
    levelElem->appendElement("numPatches", level->numPatches());
    levelElem->appendElement("totalCells", level->totalCells());

    if (level->getExtraCells() != IntVector(0, 0, 0)) {
      levelElem->appendElement("extraCells", level->getExtraCells());
    }

    levelElem->appendElement("anchor", level->getAnchor());
    levelElem->appendElement("id", level->getID());
    levelElem->appendElement("cellspacing", level->dCell());

    //__________________________________
    //  Output patch information:

    procOnLevel[l].resize(d_myworld->nRanks());

    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         ++iter) {
      const Patch* patch = *iter;

      IntVector lo    = patch->getCellLowIndex(); // for readability
      IntVector hi    = patch->getCellHighIndex();
      IntVector lo_EC = patch->getExtraCellLowIndex();
      IntVector hi_EC = patch->getExtraCellHighIndex();

      int proc             = d_load_balancer->getOutputProc(patch);
      procOnLevel[l][proc] = true;

      Box box                = patch->getExtraBox();
      ProblemSpecP patchElem = levelElem->appendChild("Patch");

      patchElem->appendElement("id", patch->getID());
      patchElem->appendElement("proc", proc);
      patchElem->appendElement("lowIndex", patch->getExtraCellLowIndex());
      patchElem->appendElement("highIndex", patch->getExtraCellHighIndex());

      if (patch->getExtraCellLowIndex() != patch->getCellLowIndex()) {
        patchElem->appendElement("interiorLowIndex", patch->getCellLowIndex());
      }

      if (patch->getExtraCellHighIndex() != patch->getCellHighIndex()) {
        patchElem->appendElement("interiorHighIndex",
                                 patch->getCellHighIndex());
      }

      patchElem->appendElement("nnodes", patch->getNumExtraNodes());
      patchElem->appendElement("lower", box.lower());
      patchElem->appendElement("upper", box.upper());
      patchElem->appendElement("totalCells", patch->getNumExtraCells());
    }
  }

  ProblemSpecP dataElem = rootElem->appendChild("Data");

  for (int l = 0; l < numLevels; l++) {
    std::ostringstream lname;
    lname << "l" << l;

    // Create a pxxxxx.xml file for each proc doing the outputting.

    for (int i = 0; i < d_myworld->nRanks(); i++) {
      if ((i % d_load_balancer->getNthRank()) != 0 || !procOnLevel[l][i]) {
        continue;
      }

      std::ostringstream pname;
      pname << lname.str() << "/p" << std::setw(5) << std::setfill('0') << i
            << ".xml";

      std::ostringstream procID;
      procID << i;

      ProblemSpecP df = dataElem->appendChild("Datafile");

      df->setAttribute("href", pname.str());
      df->setAttribute("proc", procID.str());
    }
  }

  if (hasGlobals) {
    ProblemSpecP df = dataElem->appendChild("Datafile");
    df->setAttribute("href", "global.xml");
  }

} // end writeGridOriginal()

void
DataArchiver::writeGridTextWriter(const bool hasGlobals,
                                  const string& grid_path,
                                  const GridP& grid)
{
  // With AMR, we're not guaranteed that a proc will perform work on a given
  // level.  procOnLevel tests for that.  We don't want to create a node that
  // points to no data
  int numLevels = grid->numLevels();
  std::vector<std::vector<bool>> procOnLevel(numLevels);

  // Break out the <Grid> and <Data> sections and write those to
  // grid.xml and data.xml files using libxml2's TextWriter which is a
  // streaming output format which doesn't use a DOM tree.

  string name_grid = grid_path + "grid.xml";

  // Create a new XmlWriter for uri, with no compression.
  xmlTextWriterPtr writer_grid;
  writer_grid = xmlNewTextWriterFilename(name_grid.c_str(), 0);
  xmlTextWriterSetIndent(writer_grid, 2);

  const std::string MY_ENCODING = "UTF-8";
  xmlTextWriterStartDocument(
    writer_grid, nullptr, MY_ENCODING.c_str(), nullptr);

  xmlTextWriterStartElement(writer_grid, BAD_CAST "Grid");

  //__________________________________
  //  output level information
  xmlTextWriterWriteFormatElement(
    writer_grid, BAD_CAST "numLevels", "%d", numLevels);

  for (int l = 0; l < numLevels; ++l) {
    LevelP level = grid->getLevel(l);

    xmlTextWriterStartElement(writer_grid, BAD_CAST "Level");

    if (level->getPeriodicBoundaries() != IntVector(0, 0, 0)) {

      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "periodic",
                                      "[%d,%d,%d]",
                                      level->getPeriodicBoundaries().x(),
                                      level->getPeriodicBoundaries().y(),
                                      level->getPeriodicBoundaries().z());
    }

    xmlTextWriterWriteFormatElement(
      writer_grid, BAD_CAST "numPatches", "%d", level->numPatches());
    xmlTextWriterWriteFormatElement(
      writer_grid, BAD_CAST "totalCells", "%ld", level->totalCells());
    xmlTextWriterWriteFormatElement(
      writer_grid, BAD_CAST "totalCells", "%ld", level->totalCells());

    if (level->getExtraCells() != IntVector(0, 0, 0)) {
      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "extraCells",
                                      "[%d,%d,%d]",
                                      level->getExtraCells().x(),
                                      level->getExtraCells().y(),
                                      level->getExtraCells().z());
    }
    xmlTextWriterWriteFormatElement(writer_grid,
                                    BAD_CAST "anchor",
                                    "[%g,%g,%g]",
                                    level->getAnchor().x(),
                                    level->getAnchor().y(),
                                    level->getAnchor().z());
    xmlTextWriterWriteFormatElement(
      writer_grid, BAD_CAST "id", "%d", level->getID());

    xmlTextWriterWriteFormatElement(writer_grid,
                                    BAD_CAST "cellspacing",
                                    "[%.17g,%.17g,%.17g]",
                                    level->dCell().x(),
                                    level->dCell().y(),
                                    level->dCell().z());

    //__________________________________
    //  output patch information
    procOnLevel[l].resize(d_myworld->nRanks());

    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         ++iter) {
      const Patch* patch = *iter;

      IntVector lo    = patch->getCellLowIndex(); // for readability
      IntVector hi    = patch->getCellHighIndex();
      IntVector lo_EC = patch->getExtraCellLowIndex();
      IntVector hi_EC = patch->getExtraCellHighIndex();

      int proc             = d_load_balancer->getOutputProc(patch);
      procOnLevel[l][proc] = true;

      Box box = patch->getExtraBox();

      xmlTextWriterStartElement(writer_grid, BAD_CAST "Patch");

      xmlTextWriterWriteFormatElement(
        writer_grid, BAD_CAST "id", "%d", patch->getID());
      xmlTextWriterWriteFormatElement(writer_grid, BAD_CAST "proc", "%d", proc);
      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "lowIndex",
                                      "[%d,%d,%d]",
                                      patch->getExtraCellLowIndex().x(),
                                      patch->getExtraCellLowIndex().y(),
                                      patch->getExtraCellLowIndex().z());
      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "highIndex",
                                      "[%d,%d,%d]",
                                      patch->getExtraCellHighIndex().x(),
                                      patch->getExtraCellHighIndex().y(),
                                      patch->getExtraCellHighIndex().z());
      if (patch->getExtraCellLowIndex() != patch->getCellLowIndex()) {
        xmlTextWriterWriteFormatElement(writer_grid,
                                        BAD_CAST "interiorLowIndex",
                                        "[%d,%d,%d]",
                                        patch->getCellLowIndex().x(),
                                        patch->getCellLowIndex().y(),
                                        patch->getCellLowIndex().z());
      }
      if (patch->getExtraCellHighIndex() != patch->getCellHighIndex()) {
        xmlTextWriterWriteFormatElement(writer_grid,
                                        BAD_CAST "interiorHighIndex",
                                        "[%d,%d,%d]",
                                        patch->getCellHighIndex().x(),
                                        patch->getCellHighIndex().y(),
                                        patch->getCellHighIndex().z());
      }
      xmlTextWriterWriteFormatElement(
        writer_grid, BAD_CAST "nnodes", "%d", patch->getNumExtraNodes());
      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "lower",
                                      "[%.17g,%.17g,%.17g]",
                                      box.lower().x(),
                                      box.lower().y(),
                                      box.lower().z());
      xmlTextWriterWriteFormatElement(writer_grid,
                                      BAD_CAST "upper",
                                      "[%.17g,%.17g,%.17g]",
                                      box.upper().x(),
                                      box.upper().y(),
                                      box.upper().z());
      xmlTextWriterWriteFormatElement(
        writer_grid, BAD_CAST "totalCells", "%d", patch->getNumExtraCells());
      xmlTextWriterEndElement(writer_grid); // Close Patch
    }
    xmlTextWriterEndElement(writer_grid); // Close Level
  }
  xmlTextWriterEndElement(writer_grid); // Close Grid

  xmlTextWriterEndDocument(
    writer_grid); // Writes output to the timestep.xml file
  xmlFreeTextWriter(writer_grid);

  writeDataTextWriter(hasGlobals, grid_path, grid, procOnLevel);

} // end writeGridTextWriter()

void
DataArchiver::writeDataTextWriter(
  const bool hasGlobals,
  const string& data_path,
  const GridP& grid,
  const std::vector<std::vector<bool>>& procOnLevel)
{
  int numLevels        = grid->numLevels();
  string data_filename = data_path + "data.xml";
  xmlTextWriterPtr data_writer =
    xmlNewTextWriterFilename(data_filename.c_str(), 0);

  xmlTextWriterSetIndent(data_writer, 2);
  xmlTextWriterStartElement(data_writer, BAD_CAST "Data");

  for (int l = 0; l < numLevels; l++) {
    std::ostringstream lname;
    lname << "l" << l;

    // create a pxxxxx.xml file for each proc doing the outputting
    for (int i = 0; i < d_myworld->nRanks(); i++) {

      if (i % d_load_balancer->getNthRank() != 0 || !procOnLevel[l][i]) {
        continue;
      }

      std::ostringstream pname;
      std::ostringstream procID;

      pname << lname.str() << "/p" << std::setw(5) << std::setfill('0') << i
            << ".xml";
      procID << i;

      xmlTextWriterStartElement(data_writer,
                                BAD_CAST "Datafile"); // Open <Datafile>

      xmlTextWriterWriteAttribute(
        data_writer, BAD_CAST "href", BAD_CAST pname.str().c_str());
      xmlTextWriterWriteAttribute(
        data_writer, BAD_CAST "proc", BAD_CAST procID.str().c_str());

      xmlTextWriterEndElement(data_writer); // Close <Datafile>
    }
  }

  if (hasGlobals) {
    xmlTextWriterStartElement(data_writer,
                              BAD_CAST "Datafile"); // Open <Datafile>
    xmlTextWriterWriteAttribute(
      data_writer, BAD_CAST "href", BAD_CAST "global.xml");
    xmlTextWriterEndElement(data_writer); // Close <Datafile>
  }

  xmlTextWriterEndElement(data_writer); // Close <Data>
  xmlTextWriterEndDocument(
    data_writer); // Writes output to the timestep.xml file
  xmlFreeTextWriter(data_writer);

} // end writeDataTextWriter()

void
DataArchiver::sched_outputVariables(
  std::vector<DataArchiver::SaveItem>& saveLabels,
  const GridP& grid,
  SchedulerP& sched,
  bool isThisCheckpoint)
{
  int var_cnt = 0;

  for (int i = 0; i < grid->numLevels(); i++) {

    const LevelP& level = grid->getLevel(i);
    const PatchSet* patches =
      d_load_balancer->getOutputPerProcessorPatchSet(level);

    string taskName = "DataArchiver::outputVariables";
    if (isThisCheckpoint) {
      taskName += "(checkpoint)";
    }

    Task* task = scinew Task(taskName,
                             this,
                             &DataArchiver::outputVariables,
                             isThisCheckpoint ? CHECKPOINT : OUTPUT);

    ///// hacking: add in a dummy requirement to make the checkpoint
    //             task dependent on the IO task so that they run in a
    //             deterministic order
    if (isThisCheckpoint) {
      Ghost::GhostType gn = Ghost::None;
      task->requires(Task::NewDW, d_sync_io_label, gn, 0);
    } else {
      task->computes(d_sync_io_label);
    }

    for (auto& save_item : saveLabels) {
      const MaterialSubset* matlSubset =
        save_item.getMaterialSubset(level.get_rep());

      if (matlSubset != nullptr) {
        task->requires(Task::NewDW,
                       save_item.label,
                       matlSubset,
                       Task::OutOfDomain,
                       Ghost::None,
                       0,
                       Task::SearchTG::OldTG);

        // Do not scrub any variables that are saved so they can be
        // accessed at any time after all of the tasks are finished.
        // This is needed when saving the old data warehouse as well
        // as the new data warehouse after the tasks are finished.
        sched->overrideVariableBehavior(save_item.label->getName(),
                                        false,
                                        false,
                                        !d_scrubSavedVariables,
                                        false,
                                        false);
        var_cnt++;
      }
    }

    std::string msg;
    msg = "DataArchiver::sched_outputVariables for (" +
          Uintah::to_string(var_cnt) + ") variables on";
    printSchedule(level, g_DA_dbg, msg);

    task->setType(Task::Output);
    sched->addTask(task, patches, d_mat_manager->allMaterials());
  }
}

#if HAVE_PIDX
void
DataArchiver::createPIDXCommunicator(std::vector<SaveItem>& saveLabels,
                                     const GridP& grid,
                                     SchedulerP& sched,
                                     bool isThisCheckpoint)
{
  int proc = d_myworld->myRank();

  // Resize the comms back to 0...
  d_pidxComms.clear();

  // Create new MPI Comms
  d_pidxComms.resize(grid->numLevels());

  for (int i = 0; i < grid->numLevels(); i++) {

    const LevelP& level = grid->getLevel(i);
    std::vector<SaveItem>::iterator saveIter;
    const PatchSet* patches =
      d_load_balancer->getOutputPerProcessorPatchSet(level);

    /*
      int color = 0;
      if (patches[d_myworld->myRank()].size() != 0)
        color = 1;
      MPI_Comm_split(d_myworld->getComm(), color, d_myworld->myRank(),
      &(pidxComms[i]));
   */

    int color                      = 0;
    const PatchSubset* patchsubset = patches->getSubset(proc);
    if (patchsubset->empty() == true) {
      color = 0;
    } else {
      color = 1;
    }

    MPI_Comm_split(
      d_myworld->getComm(), color, d_myworld->myRank(), &(d_pidxComms[i]));

    // if (color == 1) {
    //   int nsize;
    //   MPI_Comm_size(d_pidxComms[i], &nsize);
    //   cout << "NewComm Size = " <<  nsize << "\n";
    // }
  }
}
#endif

// be sure to call releaseDocument on the value returned
auto
DataArchiver::loadDocument(string xmlName) -> ProblemSpecP
{
  return ProblemSpecReader().readInputFile(xmlName);
}

auto
DataArchiver::getOutputLocation() const -> const std::string
{
  return d_outputDir.getName();
}

auto
DataArchiver::doesOutputDirExist() const -> bool
{
  return d_outputDir.exists();
}

void
DataArchiver::indexAddGlobals()
{
  DOUTR(g_DA_dbg, "DataArchiver::indexAddGlobals()");

  // add info to index.xml about each global (reduction) var
  // assume for now that global variables that get computed will not
  // change from timestep to timestep
  static bool wereGlobalsAdded = false;
  if (d_writeMeta && !wereGlobalsAdded && d_saveGlobalLabels.empty()) {
    wereGlobalsAdded = true;

    // add saved global (reduction) variables to index.xml
    std::string iname     = d_outputDir.getName() + "/index.xml";
    ProblemSpecP indexDoc = loadDocument(iname);

    ProblemSpecP globals = indexDoc->appendChild("globals");

    for (auto& saveItem : d_saveGlobalLabels) {
      const VarLabel* var = saveItem.label;

      // FIX - multi-level query
      const MaterialSubset* matls =
        saveItem.getMaterialSet(ALL_LEVELS)->getUnion();

      for (int m = 0; m < matls->size(); m++) {
        int matlIndex = matls->get(m);
        std::ostringstream href;
        href << var->getName();

        if (matlIndex < 0) {
          href << ".dat\0";
        } else {
          href << "_" << matlIndex << ".dat\0";
        }

        ProblemSpecP newElem = globals->appendChild("variable");
        newElem->setAttribute("href", href.str());
        newElem->setAttribute(
          "type",
          TranslateVariableType(var->typeDescription()->getName().c_str(),
                                false));
        newElem->setAttribute("name", var->getName().c_str());
      }
    }

    indexDoc->output(iname.c_str());
    // indexDoc->releaseDocument();
  }
} // end indexAddGlobals()

void
DataArchiver::outputGlobalVars(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* /*matls*/,
                               DataWarehouse* /*old_dw*/,
                               DataWarehouse* new_dw)
{
  // When recomputing a time step outputing global vars is scheduled
  // but should be skipped.
  if (d_simulator->getReductionVariable(recomputeTimeStep_name) ||
      d_saveGlobalLabels.empty()) {
    return;
  }

  bool outputGlobalVars = false;

  // If the frequency is greater than 1 check to see if output is
  // needed.
  if (d_outputGlobalVarsFrequency == 1) {
    outputGlobalVars = true;
  }
  // Note: this check is split up so to be assured that the call to
  // isLastTimeStep is last and called only when needed. Unfortunately,
  // the greater the frequency the more often it will be called.
  else if (d_simulator->getTimeStep() % d_outputGlobalVarsFrequency ==
           d_outputGlobalVarsOnTimeStep) {
    outputGlobalVars = true;
  } else {
    // Get the wall time if is needed, otherwise ignore it.
    double walltime(0);

    // The wall time is not available here
    // if (d_simulator->getWallTimeMax() > 0) {
    //   walltime = m_wall_timers.GetWallTime();
    // }
    // else {
    //   walltime = 0;
    // }

    outputGlobalVars = d_simulator->isLastTimeStep(walltime);
  }

  if (!outputGlobalVars) {
    return;
  }

  printTask(patches, g_DA_dbg, " DataArchiver::outputGlobalVars, begin");

  Timers::Simple timer;
  timer.start();

  const int timeStep   = getTimeStepTopLevel();
  const double simTime = d_simulator->getSimTime();
  const double delT    = d_simulator->getDelT();

  // Dump the stuff in the global saveset into files in the uda
  // at every timestep

  for (auto& saveItem : d_saveGlobalLabels) {
    const VarLabel* var = saveItem.label;

    // FIX, see above
    const MaterialSubset* matls =
      saveItem.getMaterialSet(ALL_LEVELS)->getUnion();

    for (int m = 0; m < matls->size(); m++) {
      int matlIndex = matls->get(m);
      DOUTR(g_DA_dbg,
            "    Reduction " << var->getName() << " matl: " << matlIndex);

      std::ostringstream filename;
      filename << d_outputDir.getName() << "/" << var->getName();
      if (matlIndex < 0) {
        filename << ".dat\0";
      } else {
        filename << "_" << matlIndex << ".dat\0";
      }

#ifdef __GNUG__
      std::ofstream out(filename.str().c_str(), std::ios::app);
#else
      std::ofstream out(filename.str().c_str(), std::ios_base::app);
#endif
      if (!out) {
        throw ErrnoException("DataArchiver::outputReduction(): The file \"" +
                               filename.str() +
                               "\" could not be opened for writing!",
                             errno,
                             __FILE__,
                             __LINE__);
      }

      out << std::setprecision(17);

      // For outputing the sim time and/or time step with the global vars
      if (d_outputGlobalVarsTimeStep) { // default false
        out << std::setw(10) << timeStep << "\t";
      }

      if (d_outputGlobalVarsSimTime) { // default true
        out << std::setprecision(17) << simTime + delT << "\t";
      }

      // Output the global var for this material index.
      new_dw->print(out, var, nullptr, matlIndex);
      out << std::endl;
    }
  }

  double myTime = timer().seconds();
  (*d_runtimeStats)[OutputGlobalIOTime] += myTime;
  (*d_runtimeStats)[TotalIOTime] += myTime;

  printTask(patches, g_DA_dbg, " DataArchiver::outputGlobalVars end");
}

void
DataArchiver::outputVariables([[maybe_unused]] const ProcessorGroup* pg,
                              const PatchSubset* patches,
                              const MaterialSubset* /*matls*/,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw,
                              int type)
{
  // IMPORTANT - this function should only be called once per processor per
  // level per type (files will be opened and closed, and those operations are
  // heavy on parallel file systems)

  // return if not an outpoint/checkpoint timestep
  if ((!d_isOutputTimeStep && type == OUTPUT) ||
      (!d_isCheckpointTimeStep &&
       (type == CHECKPOINT || type == CHECKPOINT_GLOBAL))) {
    return;
  }

  if ((d_outputPreviousTimeStep || d_checkpointPreviousTimeStep) &&
      d_simulator->getLastRegridTimeStep() > getTimeStepTopLevel()) {

    DOUT(d_myworld->myRank() == 0,
         "WARNING : Requesting the previous output/checkpoint for time step "
           << getTimeStepTopLevel() << " but the grid changed on time step "
           << d_simulator->getLastRegridTimeStep()
           << ". Not writing an output/checkpoint file.");
    return;
  }

  printTask(patches, g_DA_dbg, "DataArchiver::outputVariables  begin");

#if SCI_ASSERTION_LEVEL >= 2
  // double-check to make sure only called once per level
  int levelid = type != CHECKPOINT_GLOBAL ? getLevel(patches)->getIndex() : -1;

  if (type == OUTPUT) {
    ASSERT(d_outputCalled[levelid] == false);
    d_outputCalled[levelid] = true;
  } else if (type == CHECKPOINT) {
    ASSERT(d_checkpointCalled[levelid] == false);
    d_checkpointCalled[levelid] = true;
  } else {
    ASSERT(d_checkpointGlobalCalled == false);
    d_checkpointGlobalCalled = true;
  }
#endif

  std::vector<SaveItem>& saveLabels =
    (type == OUTPUT       ? d_saveLabels
     : type == CHECKPOINT ? d_checkpointLabels
                          : d_checkpointReductionLabels);

  DataWarehouse* dw;

  if (d_outputPreviousTimeStep || d_checkpointPreviousTimeStep) {
    dw = old_dw;
  } else {
    dw = new_dw;
  }

  //__________________________________
  // debugging output
  // this task should be called once per variable (per patch/matl subset).
  if (g_DA_dbg) {
    std::ostringstream msg;
    if (type == CHECKPOINT_GLOBAL) {
      msg << "    global";
    } else /* if (type == OUTPUT || type == CHECKPOINT) */ {
      if (type == CHECKPOINT) {
        msg << "    checkpoint ";
      }
      msg << "    patches: ";
      for (int p = 0; p < patches->size(); p++) {
        if (p != 0) {
          msg << ", ";
        }
        if (patches->get(p) == 0) {
          msg << -1;
        } else {
          msg << patches->get(p)->getID();
        }
      }
    }

    msg << " on timestep: " << getTimeStepTopLevel();
    printTask(patches, g_DA_dbg, msg.str());
  }

  //__________________________________
  Dir dir;
  if (type == OUTPUT) {
    dir = d_outputDir;
  } else {
    dir = d_checkpointsDir;
  }

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0')
        << getTimeStepTopLevel(); // could be modified by reduceUda

  Dir tdir = dir.getSubdir(tname.str());
  Dir ldir;

  std::string xmlFilename;
  std::string dataFilebase;
  std::string dataFilename;
  const Level* level = nullptr;

  // find the xml filename and data filename that we will write to
  // Normal reductions will be handled by outputReduction, but checkpoint
  // reductions call this function, and we handle them differently.
  if (type == OUTPUT || type == CHECKPOINT) {

    // find the level and level number associated with this patch
    ASSERT(patches->size() != 0);
    ASSERT(patches->get(0) != 0);

    level = patches->get(0)->getLevel();

#if SCI_ASSERTION_LEVEL >= 1
    for (int i = 0; i < patches->size(); i++) {
      ASSERT(patches->get(i)->getLevel() == level);
    }
#endif

    std::ostringstream lname;
    lname << "l" << level->getIndex(); // Hard coded - steve
    ldir = tdir.getSubdir(lname.str());

    std::ostringstream pname;
    pname << "p" << std::setw(5) << std::setfill('0') << d_myworld->myRank();
    xmlFilename  = ldir.getName() + "/" + pname.str() + ".xml";
    dataFilebase = pname.str() + ".data";
    dataFilename = ldir.getName() + "/" + dataFilebase;
  } else {
    xmlFilename  = tdir.getName() + "/global.xml";
    dataFilebase = "global.data";
    dataFilename = tdir.getName() + "/" + dataFilebase;
  }

  Timers::Simple timer;
  timer.start();

  size_t totalBytes = 0; // total bytes saved over all variables

  if (d_outputFileFormat == UDA || type == CHECKPOINT_GLOBAL) {
    d_outputLock.lock();

    {
      ProblemSpecP doc = ProblemSpec::createDocument("Uintah_Output");
      ASSERT(doc != nullptr);

      // Find the end of the file
      ProblemSpecP n = doc->findBlock("Variable");

      long cur = 0;
      while (n != nullptr) {
        ProblemSpecP endNode = n->findBlock("end");
        ASSERT(endNode != nullptr);

        long end = atol(endNode->getNodeValue().c_str());
        if (end > cur) {
          cur = end;
        }
        n = n->findNextBlock("Variable");
      }

      // Open the data file
      //
      // Note: At least one time on a BGQ machine (Vulcan@LLNL), with 160K
      // patches, a single checkpoint file failed to open, and it 'crashed' the
      // simulation.  As the other processes on the node successfully opened
      // their file, it is possible that a second open call would have
      // succeeded. (The original error no was 71.)  Therefore I am using a
      // while loop and counting the 'tries'.
      int tries            = 1;
      int flags            = O_WRONLY | O_CREAT | O_TRUNC; // file-opening flags
      const char* filename = dataFilename.c_str();
      int fd               = open(filename, flags, 0666);

      while (fd == -1) {

        if (tries >= 50) {
          std::ostringstream msg;
          msg << "DataArchiver::output(): Failed to open file '" << dataFilename
              << "' (after 50 tries).";
          std::cerr << msg.str() << "\n";
          throw ErrnoException(msg.str(), errno, __FILE__, __LINE__);
        }

        fd = open(filename, flags, 0666);
        tries++;
      }

      if (tries > 1) {
        proc0cout << "WARNING: There was a glitch in trying to open the "
                     "checkpoint file: "
                  << dataFilename << ". It took " << tries
                  << " tries to successfully open it.";
      }

      // loop over variables
      for (auto& save_item : saveLabels) {
        const VarLabel* var             = save_item.label;
        const MaterialSubset* var_matls = save_item.getMaterialSubset(level);
        if (var_matls == nullptr) {
          continue;
        }

        //  debugging output
        if (g_DA_dbg) {
          std::ostringstream msg;
          msg << "    " << var->getName() << ", materials: ";
          for (int m = 0; m < var_matls->size(); m++) {
            if (m != 0) {
              msg << ", ";
            }
            msg << var_matls->get(m);
          }
          printTask(patches, g_DA_dbg, msg.str());
        }

        // loop through patches and materials
        for (int p = 0; p < (type == CHECKPOINT_GLOBAL ? 1 : patches->size());
             p++) {
          const Patch* patch;
          int patchID;

          if (type == CHECKPOINT_GLOBAL) {
            // to consolidate into this function, force patch = 0
            patch   = nullptr;
            patchID = -1;
          } else {
            patch   = patches->get(p);
            patchID = patch->getID();
          }

          // Write info for this variable to current index file
          for (int m = 0; m < var_matls->size(); m++) {

            // add info for this variable to the current xml file
            int matlIndex = var_matls->get(m);

            // Variables may not exist when we get here due to something whacky
            // with weird AMR stuff...
            ProblemSpecP pdElem = doc->appendChild("Variable");

            pdElem->appendElement("variable", var->getName());
            pdElem->appendElement("index", matlIndex);
            pdElem->appendElement("patch", patchID);
            pdElem->setAttribute(
              "type",
              TranslateVariableType(var->typeDescription()->getName().c_str(),
                                    type != OUTPUT));

            if (var->getBoundaryLayer() != IntVector(0, 0, 0)) {
              pdElem->appendElement("boundaryLayer", var->getBoundaryLayer());
            }

            // Pad appropriately
            if (cur % PADSIZE != 0) {
              long pad   = PADSIZE - cur % PADSIZE;
              char* zero = scinew char[pad];
              memset(zero, 0, pad);
              int err = (int)write(fd, zero, pad);
              if (err != pad) {
                std::cerr << "Error writing to file: " << filename
                          << ", errno=" << errno << '\n';
                SCI_THROW(ErrnoException("DataArchiver::output (write call)",
                                         errno,
                                         __FILE__,
                                         __LINE__));
              }
              cur += pad;
              delete[] zero;
            }
            ASSERTEQ(cur % PADSIZE, 0);
            pdElem->appendElement("start", cur);

            // output data to data file
            OutputContext oc(fd,
                             filename,
                             cur,
                             pdElem,
                             d_outputDoubleAsFloat && type != CHECKPOINT);
            totalBytes += dw->emit(oc, var, matlIndex, patch);

            pdElem->appendElement("end", oc.cur);
            pdElem->appendElement("filename", dataFilebase.c_str());

#if SCI_ASSERTION_LEVEL >= 1
            struct stat st;
            int s = fstat(fd, &st);

            if (s == -1) {
              std::cerr << "fstat error - file: " << filename
                        << ", errno=" << errno << '\n';
              throw ErrnoException(
                "DataArchiver::output (stat call)", errno, __FILE__, __LINE__);
            }
            ASSERTEQ(oc.cur, st.st_size);
#endif

            cur = oc.cur;
          } // matls
        }   // patches
      }     // save items

      // close files and handles
      int s = close(fd);
      if (s == -1) {
        std::cerr << "Error closing file: " << filename << ", errno=" << errno
                  << '\n';
        throw ErrnoException(
          "DataArchiver::output (close call)", errno, __FILE__, __LINE__);
      }

      doc->output(xmlFilename.c_str());
      // doc->releaseDocument();
    } // end of locked section

    d_outputLock.unlock();
  } // end if uda or global

#if HAVE_PIDX
  //______________________________________________________________________
  //
  //  ToDo
  //      Multiple patches per core. Needed for outputNthProc (Sidharth)
  //      Turn off output from inside of PIDX (Sidharth)
  //      Disable need for MPI in PIDX (Sidharth)
  //      Fix ints issue in PIDX (Sidharth)
  //      Do we need the memset calls?
  //      Is Variable::emitPIDX() and Variable::readPIDX() efficient?
  //      Should we be using calloc() instead of malloc+memset
  //
  if (d_outputFileFormat == PIDX && type != CHECKPOINT_GLOBAL) {

    //__________________________________
    // create the xml dom for this variable
    ProblemSpecP doc = ProblemSpec::createDocument("Uintah_Output-PIDX");
    ASSERT(doc != nullptr);
    ProblemSpecP n = doc->findBlock("Variable");
    while (n != nullptr) {
      n = n->findNextBlock("Variable");
    }

    PIDXOutputContext pidx;
    std::vector<TypeDescription::Type> GridVarTypes =
      pidx.getSupportedVariableTypes();

    // Bulletproofing
    isVarTypeSupported(saveLabels, GridVarTypes);

    // Loop over the grid variable types.
    for (auto& TD : GridVarTypes) {

      // Find all of the variables of this type only.
      std::vector<SaveItem> saveTheseLabels =
        findAllVariablesWithType(saveLabels, TD);

      if (saveTheseLabels.size() > 0) {
        std::string dirName = pidx.getDirectoryName(TD);

        Dir myDir = ldir.getSubdir(dirName);

        totalBytes += saveLabels_PIDX(
          pg, patches, dw, type, saveTheseLabels, TD, ldir, dirName, doc);
      }
    }

    // write the xml

    // The following line was commented out... but with it on, we can
    // restart with PIDX.  However, we need to be able to restart
    // without this.  Uncomment the following line to make progress on
    // restart of PIDX - see line above.
    // doc->output( xmlFilename.c_str() );
  }
#endif

  double myTime   = timer().seconds();
  double byteToMB = 1024 * 1024;

  if (type == OUTPUT) {
    (*d_runtimeStats)[OutputIOTime] += myTime;
    (*d_runtimeStats)[OutputIORate] += (double)totalBytes / (byteToMB * myTime);
  } else if (type == CHECKPOINT) {
    (*d_runtimeStats)[CheckpointIOTime] += myTime;
    (*d_runtimeStats)[CheckpointIORate] +=
      (double)totalBytes / (byteToMB * myTime);
  } else { // type == CHECKPOINT_GLOBAL
    (*d_runtimeStats)[CheckpointGlobalIOTime] += myTime;
    (*d_runtimeStats)[CheckpointGlobalIORate] +=
      (double)totalBytes / (byteToMB * myTime);
  }

  (*d_runtimeStats)[TotalIOTime] += myTime;

  printTask(patches, g_DA_dbg, "DataArchiver::outputVariables  end");

} // end output()

//  output only the savedLabels of a specified type description in PIDX format.
auto
DataArchiver::saveLabels_PIDX(
  [[maybe_unused]] const ProcessorGroup* pg,
  [[maybe_unused]] const PatchSubset* patches,
  [[maybe_unused]] DataWarehouse* new_dw,
  [[maybe_unused]] int type,
  [[maybe_unused]] std::vector<SaveItem>& saveLabels,
  [[maybe_unused]] const TypeDescription::Type TD,
  [[maybe_unused]] Dir ldir,                   // uda/timestep/levelIndex
  [[maybe_unused]] const std::string& dirName, // CCVars, SFC*Vars
  [[maybe_unused]] ProblemSpecP& doc) -> size_t
{

  size_t totalBytesSaved = 0;
#ifdef HAVE_PIDX
  if (dbgPIDX.active()) {
    dbgPIDX << "saveLabels_PIDX()\n";
  }

  const int timeStep = d_simulator->getTimeStep();

  const int levelid  = getLevel(patches)->getIndex();
  const Level* level = getLevel(patches);

  int nSaveItems = saveLabels.size();
  std::vector<int> nSaveItemMatls(nSaveItems);

  int rank = pg->myRank();
  int rc   = -9; // PIDX return code

  //__________________________________
  // Count up the number of variables that will
  // be output. Each variable can have a different number of
  // materials and live on a different number of levels

  int count                      = 0;
  int actual_number_of_variables = 0;

  for (std::vector<SaveItem>::iterator saveIter = saveLabels.begin();
       saveIter != saveLabels.end();
       ++saveIter) {
    const MaterialSubset* var_matls = saveIter->getMaterialSubset(level);

    if (var_matls == nullptr) {
      continue;
    }

    nSaveItemMatls[count] = var_matls->size();

    ++count;
    actual_number_of_variables += var_matls->size();
  }

  if (actual_number_of_variables == 0) {
    // Don't actually have any variables of the specified type on this level.
    return totalBytesSaved;
  }

  // must use this format or else file won't be written
  // inside of uda
  Dir myDir = ldir.getSubdir(dirName); // uda/timestep/level/<CCVars, SFC*Vars>
  std::string idxFilename(dirName + ".idx"); // <CCVars, SFC*Vars....>.idx
  std::string full_idxFilename(
    myDir.getName() + ".idx"); // uda/timestep/level/<CCVars, SFC*Vars...>.idx

  PIDXOutputContext pidx;

  pidx.setOutputDoubleAsFloat((d_outputDoubleAsFloat && type == OUTPUT));

  //__________________________________
  // define the level extents for this variable type
  IntVector lo;
  IntVector hi;
  level->computeVariableExtents(TD, lo, hi);

  PIDX_point level_size;
  pidx.setLevelExtents("DataArchiver::saveLabels_PIDX", lo, hi, level_size);

  // Can this be run in serial without doing a MPI initialize

  d_PIDX_flags.print();

  if (TD != Uintah::TypeDescription::Type::ParticleVariable) {

    pidx.initialize(full_idxFilename,
                    timeStep,
                    /*d_myworld->getComm()*/ d_pidxComms[levelid],
                    d_PIDX_flags,
                    level_size,
                    type);

  } else {
    pidx.initializeParticles(
      full_idxFilename, timeStep, d_pidxComms[levelid], level_size, type);

    PIDX_physical_point physical_global_size;
    IntVector zlo = { 0, 0, 0 };
    IntVector ohi = { 1, 1, 1 };
    BBox b;
    level->getSpatialRange(b);

    PIDX_set_physical_point(physical_global_size,
                            b.max().x() - b.min().x(),
                            b.max().y() - b.min().y(),
                            b.max().z() - b.min().z());

    PIDX_set_physical_dims(pidx.file, physical_global_size);
  }

  //__________________________________
  // allocate memory for pidx variable descriptor array

  rc = PIDX_set_variable_count(pidx.file, actual_number_of_variables);
  pidx.checkReturnCode(
    rc,
    "DataArchiver::saveLabels_PIDX -PIDX_set_variable_count failure",
    __FILE__,
    __LINE__);

  size_t pidxVarSize = sizeof(PIDX_variable*);
  pidx.varDesc       = (PIDX_variable**)malloc(pidxVarSize * nSaveItems);
  memset(pidx.varDesc, 0, pidxVarSize * nSaveItems);

  for (int i = 0; i < nSaveItems; i++) {
    pidx.varDesc[i] = (PIDX_variable*)malloc(pidxVarSize * nSaveItemMatls[i]);
    memset(pidx.varDesc[i], 0, pidxVarSize * nSaveItemMatls[i]);
  }

  //__________________________________
  //  PIDX Diagnostics
  if (rank == 0 && dbgPIDX.active()) {
    printf("[PIDX] IDX file name = %s\n", (char*)idxFilename.c_str());
    printf("[PIDX] levelExtents = %d %d %d\n",
           (hi.x() - lo.x()),
           (hi.y() - lo.y()),
           (hi.z() - lo.z()));
    printf("[PIDX] Total number of variable = %d\n", nSaveItems);
  }

  //__________________________________
  //
  unsigned char*** patch_buffer;
  patch_buffer = (unsigned char***)malloc(sizeof(unsigned char**) *
                                          actual_number_of_variables);

  int vc  = 0; // Variable Counter
  int vcm = 0; // Variable Counter Material

  for (std::vector<SaveItem>::iterator saveIter = saveLabels.begin();
       saveIter != saveLabels.end();
       ++saveIter) {
    const VarLabel* label = saveIter->label;

    const MaterialSubset* var_matls = saveIter->getMaterialSubset(level);
    if (var_matls == nullptr) {
      continue;
    }

    //__________________________________
    //  determine data subtype sizes
    const TypeDescription* td      = label->typeDescription();
    const TypeDescription* subtype = td->getSubType();

    // set values depending on the variable's subtype
    char data_type[512];
    int sample_per_variable = 1;
    size_t varSubType_size  = -9;

    string type_name = "float64";
    int the_size     = sizeof(double);

    switch (subtype->getType()) {

      case Uintah::TypeDescription::Type::long64_type:
        sample_per_variable = 1;
        the_size            = sizeof(
          double); // FIXME: what is the Uintah version of a 64 bit integer?
        type_name = "int64";
        break;
      case Uintah::TypeDescription::Type::Matrix3:
        sample_per_variable = 9;
        break;
      case Uintah::TypeDescription::Type::Point:
        sample_per_variable = 3;
        break;
      case Uintah::TypeDescription::Type::Stencil4:
        sample_per_variable = 4;
        break;
      case Uintah::TypeDescription::Type::Stencil7:
        sample_per_variable = 7;
        break;
      case Uintah::TypeDescription::Type::Vector:
        sample_per_variable = 3;
        break;

      case Uintah::TypeDescription::Type::IntVector:
        sample_per_variable = 3;
        the_size            = sizeof(int);
        type_name           = "int32";
        break;
      case Uintah::TypeDescription::Type::int_type:
        the_size  = sizeof(int);
        type_name = "int32";
        break;
      case Uintah::TypeDescription::Type::double_type:
        if (pidx.isOutputDoubleAsFloat()) { // Take into account saving doubles
                                            // as floats
          the_size  = sizeof(float);
          type_name = "float32";
        }
        break;
      case Uintah::TypeDescription::Type::float_type:
        the_size  = sizeof(float);
        type_name = "float32";
        break;

      default:
        std::ostringstream warn;
        warn << "DataArchiver::saveLabels_PIDX:: (" << label->getName() << " "
             << td->getName() << " ) has not been implemented\n";
        throw InternalError(warn.str(), __FILE__, __LINE__);
    }

    varSubType_size = sample_per_variable * the_size;
    sprintf(data_type, "%d*%s", sample_per_variable, type_name.c_str());

    //__________________________________
    //  materials loop
    for (int m = 0; m < var_matls->size(); m++) {
      int matlIndex = var_matls->get(m);
      std::string var_mat_name;

      std::ostringstream s;
      s << m;
      var_mat_name = label->getName() + "_m" + s.str();

      bool isParticle = (label->typeDescription()->getType() ==
                         Uintah::TypeDescription::Type::ParticleVariable);

      rc = PIDX_variable_create((char*)var_mat_name.c_str(),
                                /* isParticle, <- Add this to PIDX spec*/
                                varSubType_size * 8,
                                data_type,
                                &(pidx.varDesc[vc][m]));

      pidx.checkReturnCode(
        rc,
        "DataArchiver::saveLabels_PIDX - PIDX_variable_create failure",
        __FILE__,
        __LINE__);

      patch_buffer[vcm] =
        (unsigned char**)malloc(sizeof(unsigned char*) * patches->size());

      //__________________________________
      //  patch Loop
      for (int p = 0; p < (type == CHECKPOINT_GLOBAL ? 1 : patches->size());
           ++p) {
        const Patch* patch;

        if (type == CHECKPOINT_GLOBAL) {
          patch = 0;
        } else {
          patch = patches->get(p);
          PIDX_point patchOffset;
          PIDX_point patchSize;
          PIDXOutputContext::patchExtents patchExts;

          pidx.setPatchExtents("DataArchiver::saveLabels_PIDX",
                               patch,
                               level,
                               label->getBoundaryLayer(),
                               td,
                               patchExts,
                               patchOffset,
                               patchSize);

          //__________________________________
          // debugging
          if (dbgPIDX.active() && isProc0_macro) {
            proc0cout << rank << " taskType: " << type
                      << "  PIDX:  " << std::setw(15) << label->getName()
                      << "  " << td->getName() << " Patch: " << patch->getID()
                      << " L-" << level->getIndex()
                      << ",  sample_per_variable: " << sample_per_variable
                      << " varSubType_size: " << varSubType_size
                      << " dataType: " << data_type << "\n";
            patchExts.print(cout);
          }

          //__________________________________
          //  Read in Array3 data to t-buffer
          size_t arraySize;
          if (td->getType() ==
              Uintah::TypeDescription::Type::ParticleVariable) {

            if (label->getName() == d_particlePositionName) {
              int var_index;
              PIDX_get_current_variable_index(pidx.file, &var_index);
              PIDX_set_particles_position_variable_index(pidx.file, var_index);
            }

            ParticleVariableBase* pv =
              new_dw->getParticleVariable(label, matlIndex, patch);
            void* ptr;    // Pointer to data? but we don't use it so not digging
                          // into what it is for.
            string elems; // We don't use this either.
            pv->getSizeInfo(elems, arraySize, ptr);

            int num_particles = arraySize / (sample_per_variable * the_size);

            patch_buffer[vcm][p] = (unsigned char*)malloc(arraySize);

            // Figure out the # of particles in order to create the patch_buffer
            // of the correct size... Or pass in a vector and figure it out
            // after this call somewhere...
            //
            // if we don't have to do the above comments, then we can
            // consolidate this code with the code below.
            //

            new_dw->emitPIDX(
              pidx, label, matlIndex, patch, patch_buffer[vcm][p], arraySize);

            PIDX_physical_point physical_local_offset, physical_local_size;
            PIDX_set_physical_point(
              physical_local_size,
              patch->getBox().upper().x() - patch->getBox().lower().x(),
              patch->getBox().upper().y() - patch->getBox().lower().y(),
              patch->getBox().upper().z() - patch->getBox().lower().z());
            PIDX_set_physical_point(physical_local_offset,
                                    patch->getBox().lower().x(),
                                    patch->getBox().lower().y(),
                                    patch->getBox().lower().z());

            rc = PIDX_variable_write_particle_data_physical_layout(
              pidx.varDesc[vc][m],
              physical_local_offset,
              physical_local_size,
              patch_buffer[vcm][p],
              num_particles,
              PIDX_row_major);
          } else {

            //__________________________________
            // allocate memory for the grid variables
            arraySize = varSubType_size * patchExts.totalCells_EC;

            patch_buffer[vcm][p] = (unsigned char*)malloc(arraySize);
            memset(patch_buffer[vcm][p], 0, arraySize);

            new_dw->emitPIDX(
              pidx, label, matlIndex, patch, patch_buffer[vcm][p], arraySize);

            IntVector extra_cells = patch->getExtraCells();

            patchOffset[0] = patchOffset[0] - lo.x() - extra_cells[0];
            patchOffset[1] = patchOffset[1] - lo.y() - extra_cells[1];
            patchOffset[2] = patchOffset[2] - lo.z() - extra_cells[2];

            rc = PIDX_variable_write_data_layout(pidx.varDesc[vc][m],
                                                 patchOffset,
                                                 patchSize,
                                                 patch_buffer[vcm][p],
                                                 PIDX_row_major);

            pidx.checkReturnCode(rc,
                                 "DataArchiver::saveLabels_PIDX - "
                                 "PIDX_variable_write_data_layout failure",
                                 __FILE__,
                                 __LINE__);
          }

#if 0 // to hardwire buffer values for debugging.
          pidx.hardWireBufferValues(t_buffer, patchExts, arraySize, sample_per_variable );
#endif

          //__________________________________
          //  debugging
          if (dbgPIDX.active()) {
            pidx.printBufferWrap("DataArchiver::saveLabels_PIDX    BEFORE  "
                                 "PIDX_variable_write_data_layout",
                                 subtype->getType(),
                                 sample_per_variable,
                                 patchExts.lo_EC,
                                 patchExts.hi_EC,
                                 patch_buffer[vcm][p],
                                 arraySize);
          }

          totalBytesSaved += arraySize;

#if 0
          //__________________________________
          //  populate the xml dom This layout allows us to highjack
          //  all of the existing data structures in DataArchive
          ProblemSpecP var_ps = doc->appendChild("Variable");
          var_ps->setAttribute( "type",      TranslateVariableType( td->getName().c_str(), type != OUTPUT ) );
          var_ps->appendElement("variable",  label->getName());
          var_ps->appendElement("index",     matlIndex);
          var_ps->appendElement("patch",     patch->getID());
          var_ps->appendElement("start",     vcm);
          var_ps->appendElement("end",       vcm);          // redundant
          var_ps->appendElement("filename",  idxFilename);

          if (label->getBoundaryLayer() != IntVector(0,0,0)) {
            var_ps->appendElement("boundaryLayer", label->getBoundaryLayer());
          }
#endif
        } // is checkpoint?
      }   //  Patches

      rc = PIDX_append_and_write_variable(pidx.file, pidx.varDesc[vc][m]);

      int compression_type;
      PIDX_get_compression_type(pidx.file, &compression_type);

      // TODO uncomment when we will use latest PIDX version
      // if(compression_type == PIDX_CHUNKING_ZFP)
      //  PIDX_set_lossy_compression_bit_rate(pidx.file, pidx.varDesc[vc][m],
      //  m_PIDX_flags.d_checkpointFlags.compressionBitrate);

      pidx.checkReturnCode(rc,
                           "DataArchiver::saveLabels_PIDX - "
                           "PIDX_append_and_write_variable failure",
                           __FILE__,
                           __LINE__);
      vcm++;
    } //  Materials

    vc++;
  } //  Variables

  rc = PIDX_close(pidx.file);
  pidx.checkReturnCode(rc,
                       "DataArchiver::saveLabels_PIDX - PIDX_close failure",
                       __FILE__,
                       __LINE__);

  //__________________________________
  // free buffers
  for (int i = 0; i < actual_number_of_variables; i++) {
    for (int p = 0; p < (type == CHECKPOINT_GLOBAL ? 1 : patches->size());
         p++) {
      free(patch_buffer[i][p]);
      patch_buffer[i][p] = 0;
    }
    free(patch_buffer[i]);
    patch_buffer[i] = 0;
  }

  free(patch_buffer);
  patch_buffer = 0;

  //__________________________________
  //  free memory
  for (int i = 0; i < nSaveItems; i++) {
    free(pidx.varDesc[i]);
  }

  free(pidx.varDesc);
  pidx.varDesc = 0;

  if (dbgPIDX.active()) {
    dbgPIDX << "end saveLabels_PIDX()\n";
  }

#endif

  return totalBytesSaved;
} // end saveLabels_PIDX();

//  Return a vector of saveItems with a common typeDescription
auto
DataArchiver::findAllVariablesWithType(const std::vector<SaveItem>& saveLabels,
                                       const TypeDescription::Type type)
  -> std::vector<DataArchiver::SaveItem>
{
  std::vector<SaveItem> myItems;
  for (const auto& save_item : saveLabels) {
    const VarLabel* label        = save_item.label;
    TypeDescription::Type myType = label->typeDescription()->getType();
    if (myType == type) {
      myItems.push_back(save_item);
    }
  }
  return myItems;
}

//  throw exception if saveItems type description is NOT supported
void
DataArchiver::isVarTypeSupported(
  const std::vector<SaveItem>& saveLabels,
  const std::vector<TypeDescription::Type>& pidxVarTypes)
{
  for (const auto& save_item : saveLabels) {
    const VarLabel* label         = save_item.label;
    const TypeDescription* myType = label->typeDescription();

    bool found = false;
    for (const auto& type : pidxVarTypes) {
      if (myType->getType() == type) {
        found = true;
        continue;
      }
    }

    // throw exception if this type isn't supported
    if (found == false) {
      std::ostringstream warn;
      warn << "DataArchiver::saveLabels_PIDX:: ( " << label->getName() << ",  "
           << myType->getName() << " ) has not been implemented";
      throw InternalError(warn.str(), __FILE__, __LINE__);
    }
  }
}

void
DataArchiver::makeVersionedDir()
{
  unsigned int dirMin = 0;
  unsigned int dirNum = 0;
  unsigned int dirMax = 0;

  bool dirCreated = false;

  // My (Dd) understanding of this while loop is as follows: We try to
  // make a new directory starting at 0 (000), then 1 (001), and then
  // doubling the number each time that the directory already exists.
  // Once we find a directory number that does not exist, we start
  // back at the minimum number that we had already tried and work
  // forward to find an actual directory number that we can use.
  //
  // Eg: If 001 and 002 already exist, then the algorithm tries to
  // create 001, fails (it already exists), tries to created 002,
  // fails again, tries to create 004, succeeds, deletes 004 because
  // we aren't sure that it is the smallest new number, creates 003,
  // succeeds, and since 003 is the smallest (dirMin) the algorithm
  // stops.
  //
  // This routine has been re-written to not use the Dir class because
  // the Dir class throws exceptions and exceptions (under SGI CC) add
  // a huge memory penalty in one (same penalty if more than) are
  // thrown.  This causes memory usage to change if the very first
  // directory (000) can be created (because no exception is thrown
  // and thus no memory is allocated for exceptions).
  //
  // If there is a real error, then we can throw an exception because
  // we don't care about the memory penalty.

  string dirName;

  // first check to see if the suffix passed in on the command line
  // makes a valid uda dir

  if (d_udaSuffix != -1) {
    std::ostringstream name;
    name << d_filebase << "." << std::setw(3) << std::setfill('0')
         << d_udaSuffix;
    dirName = name.str();

    int code = MKDIR(dirName.c_str(), 0777);
    if (code == 0) { // Created the directory successfully
      dirCreated = true;
    } else if (errno != EEXIST) {
      std::cerr << "makeVersionedDir: Error making directory: " << name.str()
                << "\n";
      throw ErrnoException("DataArchiver.cc: mkdir failed for some "
                           "reason besides dir already exists",
                           errno,
                           __FILE__,
                           __LINE__);
    }
  }

  // if that didn't work, go ahead with the real algorithm
  while (!dirCreated) {
    std::ostringstream name;
    name << d_filebase << "." << std::setw(3) << std::setfill('0') << dirNum;
    dirName = name.str();

    int code = MKDIR(dirName.c_str(), 0777);
    if (code == 0) { // Created the directory successfully
      dirMax = dirNum;
      if (dirMax == dirMin) {
        dirCreated = true;
      } else {
        int code = rmdir(dirName.c_str());
        if (code != 0) {
          throw ErrnoException(
            "DataArchiver.cc: rmdir failed", errno, __FILE__, __LINE__);
        }
      }
    } else {
      if (errno != EEXIST) {
        std::cerr << "makeVersionedDir: Error making directory: " << name.str()
                  << "\n";
        throw ErrnoException("DataArchiver.cc: mkdir failed for some "
                             "reason besides dir already exists",
                             errno,
                             __FILE__,
                             __LINE__);
      }
      dirMin = dirNum + 1;
    }

    if (!dirCreated) {
      if (dirMax == 0) {
        if (dirNum == 0) {
          dirNum = 1;
        } else {
          dirNum *= 2;
        }
      } else {
        dirNum = dirMin + ((dirMax - dirMin) / 2);
      }
    }
  }
  // Move the symbolic link to point to the new version. We need to be careful
  // not to destroy data, so we only create the link if no file/dir by that
  // name existed or if it's already a link.
  bool make_link = false;
  struct stat sb;
  int rc = LSTAT(d_filebase.c_str(), &sb);
  if ((rc != 0) && (errno == ENOENT)) {
    make_link = true;
  } else if ((rc == 0) && (S_ISLNK(sb.st_mode))) {
    unlink(d_filebase.c_str());
    make_link = true;
  }
  if (make_link) {
    [[maybe_unused]] auto status = symlink(dirName.c_str(), d_filebase.c_str());
  }

  std::cout << "DataArchiver created " << dirName << std::endl;
  d_outputDir = Dir(dirName);

} // end makeVersionedDir()

//______________________________________________________________________
//  Determine which labels will be saved.
void
DataArchiver::initSaveLabels(SchedulerP& sched, bool initTimeStep)
{
  // if this is the initTimeStep, then don't complain about saving all the vars,
  // just save the ones you can.  They'll most likely be around on the next
  // timestep.
  SaveItem saveItem;
  d_saveGlobalLabels.clear();
  d_saveLabels.clear();
  d_saveLabels.reserve(d_saveLabelNames.size());

  auto pLabelMatlMap = sched->makeVarLabelMaterialMap();

  for (const auto& saveLabel : d_saveLabelNames) {

    // go through each of the saveLabelNames we created in problemSetup
    //   see if that variable has been created, set the compression mode
    //   make sure that the scheduler shows that that it has been scheduled
    //   to be computed.  Then save it to saveItems.
    VarLabel* var = VarLabel::find(saveLabel.labelName);
    if (var == nullptr) {
      if (initTimeStep) {
        continue;
      } else {
        throw ProblemSetupException(saveLabel.labelName +
                                      " variable not found to save.",
                                    __FILE__,
                                    __LINE__);
      }
    }
    if (saveLabel.compressionMode != "") {
      var->setCompressionMode(saveLabel.compressionMode);
    }

    auto found = pLabelMatlMap->find(var->getName());

    if (found == pLabelMatlMap->end()) {
      if (initTimeStep) {
        // ignore this on the init timestep, cuz lots of vars aren't computed on
        // the init timestep
        DOUTR(g_DA_dbg,
              "Ignoring var " << saveLabel.labelName
                              << " on initialization timestep\n");
        continue;
      } else {
        throw ProblemSetupException(saveLabel.labelName +
                                      " variable not computed for saving.",
                                    __FILE__,
                                    __LINE__);
      }
    }

    saveItem.label = var;
    saveItem.matlSet.clear();
    for (ConsecutiveRangeSet::iterator iter = saveLabel.levels.begin();
         iter != saveLabel.levels.end();
         iter++) {

      ConsecutiveRangeSet matlsToSave =
        (ConsecutiveRangeSet((*found).second)).intersected(saveLabel.matls);
      saveItem.setMaterials(*iter, matlsToSave, d_prevMatls, d_prevMatlSet);

      if ((saveLabel.matls != ConsecutiveRangeSet::all) &&
          (saveLabel.matls != matlsToSave)) {
        throw ProblemSetupException(
          saveLabel.labelName +
            " variable not computed for all materials specified to save.",
          __FILE__,
          __LINE__);
      }
    }

    const TypeDescription::Type varType =
      saveItem.label->typeDescription()->getType();

    if (varType == TypeDescription::Type::ReductionVariable ||
        varType == TypeDescription::Type::SoleVariable) {
      d_saveGlobalLabels.push_back(saveItem);
    } else {
      d_saveLabels.push_back(saveItem);
    }
  }

  // d_saveLabelNames.clear();
  DOUTR(g_DA_dbg, "end of initSaveLabels\n");
}

void
DataArchiver::initCheckpoints(SchedulerP& sched)
{
  DOUTR(g_DA_dbg, "initCheckpoints called\n");

  using DepVector           = std::vector<const Task::Dependency*>;
  const DepVector& initreqs = sched->getInitialRequires();

  // special variables to not checkpoint
  const std::set<string>& notCheckPointVars = sched->getNotCheckPointVars();

  SaveItem saveItem;
  d_checkpointReductionLabels.clear();
  d_checkpointLabels.clear();

  // label -> level -> matls
  // we can't store them in two different maps, since we need to know
  // which levels depend on which materials
  using LabelType = std::map<string, std::map<int, ConsecutiveRangeSet>>;
  LabelType label_map;

  for (const auto& dep : initreqs) {

    // define patchset
    const PatchSubset* patchSubset = (dep->m_patches != nullptr)
                                       ? dep->m_patches
                                       : dep->m_task->getPatchSet()->getUnion();

    // adjust the patchSubset if the dependency requires coarse or fine level
    // patches
    constHandle<PatchSubset> patches;
    if (dep->m_patches_dom == Task::CoarseLevel ||
        dep->m_patches_dom == Task::FineLevel) {
      patches     = dep->getPatchesUnderDomain(patchSubset);
      patchSubset = patches.get_rep();
    }

    // Define the Levels
    ConsecutiveRangeSet levels;
    for (int i = 0; i < patchSubset->size(); i++) {
      const Patch* patch = patchSubset->get(i);
      levels.addInOrder(patch->getLevel()->getIndex());
    }

    // **********************************************************************
    // NOTE: Sole vars are not quite get handled as expected. Unlike
    // reduction variables they do NOT have a material dependency
    // (dep->m_matls) which would return a material of -1. As such,
    // carve off a special material index that is also -1. So to
    // match what is in the data warehouse.

    // When putting a sole variable into the data warehouse it wants
    // a material index of -1. Which is similar to a reduction
    // variable. The bottomline is that sole variable should be
    // handled like reduction variables through out are not.
    // **********************************************************************
    const MaterialSubset* matSubset;

    if (dep->m_matls) {
      matSubset = dep->m_matls;
    }
    // Special case (hack) so sole variables have a material index of -1.
    else if (dep->m_var->typeDescription()->getType() ==
             TypeDescription::Type::SoleVariable) {
      matSubset = d_tmpMatSubset;
    } else {
      matSubset = dep->m_task->getMaterialSet()->getUnion();
    }

    // The matSubset is assumed to be in ascending order or
    // addInOrder will throw an exception.
    ConsecutiveRangeSet matls;
    matls.addInOrder(matSubset->getVector().begin(),
                     matSubset->getVector().end());

    for (int level : levels) {
      ConsecutiveRangeSet& unionedVarMatls =
        label_map[dep->m_var->getName()][level];
      unionedVarMatls = unionedVarMatls.unioned(matls);
    }
  }

  d_checkpointLabels.reserve(label_map.size());
  bool hasDelT = false;
  for (auto& [label_name, label_set_map] : label_map) {
    VarLabel* var = VarLabel::find(label_name);
    if (var == nullptr) {
      throw ProblemSetupException(
        label_name + " variable not found to checkpoint.", __FILE__, __LINE__);
    }

    saveItem.label = var;
    saveItem.matlSet.clear();

    for (auto& [id, range_set] : label_set_map) {
      saveItem.setMaterials(id, range_set, d_prevMatls, d_prevMatlSet);

      if (string(var->getName()) == "delT") {
        hasDelT = true;
      }
    }

    // Skip this variable if the default behavior of variable has been
    // overwritten. For example ignore checkpointing PerPatch<FileInfo> variable
    bool skipVar = (notCheckPointVars.count(saveItem.label->getName()) > 0);

    if (!skipVar) {
      if (saveItem.label->typeDescription()->isReductionVariable()) {
        d_checkpointReductionLabels.push_back(saveItem);
      } else {
        d_checkpointLabels.push_back(saveItem);
      }
    }
  }

  if (!hasDelT) {
    VarLabel* var = VarLabel::find(delT_name);
    if (var == nullptr) {
      throw ProblemSetupException(
        "delT variable not found to checkpoint.", __FILE__, __LINE__);
    }
    saveItem.label = var;
    saveItem.matlSet.clear();

    ConsecutiveRangeSet globalMatl("-1");
    saveItem.setMaterials(-1, globalMatl, d_prevMatls, d_prevMatlSet);

    ASSERT(saveItem.label->typeDescription()->getType() ==
             TypeDescription::Type::ReductionVariable ||
           saveItem.label->typeDescription()->getType() ==
             TypeDescription::Type::SoleVariable);

    d_checkpointReductionLabels.push_back(saveItem);
  }
}

void
DataArchiver::SaveItem::setMaterials(int level_id,
                                     const ConsecutiveRangeSet& matls,
                                     ConsecutiveRangeSet& prevMatls,
                                     MaterialSetP& prevMatlSet)
{
  // reuse material sets when the same set of materials is used for different
  // SaveItems in a row -- easier than finding all reusable material set, but
  // effective in many common cases.
  if ((prevMatlSet != nullptr) && (matls == prevMatls)) {
    matlSet[level_id] = prevMatlSet;
  } else {
    MaterialSetP& m = matlSet[level_id];
    m               = scinew MaterialSet();

    std::vector<int> matlVec;
    matlVec.reserve(matls.size());

    for (int matl : matls) {
      matlVec.push_back(matl);
    }

    m->addAll(matlVec);
    prevMatlSet = m;
    prevMatls   = matls;
  }
}

//  Find the materials to output on this level for this saveItem
auto
DataArchiver::SaveItem::getMaterialSubset(const Level* level) const
  -> const MaterialSubset*
{
  // search done by absolute level, or relative to end of levels (-1
  // finest, -2 second finest,...)
  auto iter                       = matlSet.end();
  const MaterialSubset* var_matls = nullptr;

  if (level) {
    int L_index   = level->getIndex();
    int maxLevels = level->getGrid()->numLevels();

    iter = matlSet.find(L_index);

    if (iter == matlSet.end()) {
      iter = matlSet.find(L_index - maxLevels);
    }

    if (iter == matlSet.end()) {
      iter = matlSet.find(ALL_LEVELS);
    }

    if (iter != matlSet.end()) {
      var_matls = iter->second.get_rep()->getUnion();
    }
  } else { // Globals variables that are level independent:

    for (auto& [mat_id, mat_set] : matlSet) {
      var_matls = getMaterialSet(mat_id)->getUnion();
      break;
    }
  }
  return var_matls;
}

auto
DataArchiver::needRecompile(const GridP& /*grid*/) -> bool
{
  bool retVal = false;

#ifdef HAVE_PIDX
  retVal = (retVal || d_pidx_need_to_recompile || d_pidx_restore_nth_rank);
#endif

  return retVal;
  return false;
}

void
DataArchiver::recompile([[maybe_unused]] const GridP& grid)
{

#ifdef HAVE_PIDX
  if (d_pidx_requested_nth_rank > 1) {
    if (d_pidx_restore_nth_rank) {
      proc0cout << "This is the time step following a checkpoint - "
                << "need to put the task graph back with a recompile - "
                << "setting nth output to 1\n";
      d_load_balancer->setNthRank(1);
      d_load_balancer->possiblyDynamicallyReallocate(grid,
                                                     LoadBalancer::REGRID_LB);
      setSaveAsPIDX();
      d_pidx_restore_nth_rank = false;
    }

    if (d_pidx_need_to_recompile) {
      // Don't need to recompile on the next time step as it will
      // happen on this one.  However, the nth rank value will
      // need to be restored after this time step, so set
      // pidx_restore_nth_rank to true.
      d_pidx_need_to_recompile = false;
      d_pidx_restore_nth_rank  = true;
    }
  }
#endif
}

auto
DataArchiver::TranslateVariableType(string type, bool isThisCheckpoint)
  -> string
{
  if (d_outputDoubleAsFloat && !isThisCheckpoint) {
    if (type == "CCVariable<double>") {
      return "CCVariable<float>";
    } else if (type == "NCVariable<double>") {
      return "NCVariable<float>";
    } else if (type == "SFCXVariable<double>") {
      return "SFCXVariable<float>";
    } else if (type == "SFCYVariable<double>") {
      return "SFCYVariable<float>";
    } else if (type == "SFCZVariable<double>") {
      return "SFCZVariable<float>";
    } else if (type == "ParticleVariable<double>") {
      return "ParticleVariable<float>";
    }
  }
  return type;
}

auto
DataArchiver::isLabelSaved(const string& label) const -> bool
{
  if (d_outputInterval == 0.0 && d_outputTimeStepInterval == 0) {
    return false;
  }

  for (const auto& save_label : d_saveLabelNames) {
    if (save_label.labelName == label) {
      return true;
    }
  }
  return false;
}

// Allow the component to set the output interval
void
DataArchiver::setOutputInterval(double newinv)
{
  if (d_outputInterval != newinv) {
    d_outputInterval         = newinv;
    d_outputTimeStepInterval = 0;

    d_nextOutputTime     = 0.0;
    d_nextOutputTimeStep = 0;
  }
}

// Allow the component to set the output timestep interval
void
DataArchiver::setOutputTimeStepInterval(int newinv)
{
  if (d_outputTimeStepInterval != newinv) {

    d_outputInterval         = 0;
    d_outputTimeStepInterval = newinv;

    // Reset the output so force one to happen.
    d_nextOutputTime     = 0.0;
    d_nextOutputTimeStep = 0;
  }
}

// Allow the component to set the checkpoint interval
void
DataArchiver::setCheckpointInterval(double newinv)
{
  if (d_checkpointInterval != newinv) {

    d_checkpointInterval         = newinv;
    d_checkpointTimeStepInterval = 0;
    d_checkpointWallTimeInterval = 0;

    // Reset the check point so force one to happen.
    d_nextCheckpointTime     = 0.0;
    d_nextCheckpointTimeStep = 0.0;
    d_nextCheckpointWallTime = 0.0;

    // If needed create checkpoints/index.xml
    if (!d_checkpointsDir.exists()) {

      if (d_myworld->myRank() == 0) {
        d_checkpointsDir = d_outputDir.createSubdir("checkpoints");
        createIndexXML(d_checkpointsDir);
      }
    }

    // Sync up before every rank can use the checkpoints dir
    Uintah::MPI::Barrier(d_myworld->getComm());
  }
}

// Allow the component to set the checkpoint timestep interval
void
DataArchiver::setCheckpointTimeStepInterval(int newinv)
{
  if (d_checkpointTimeStepInterval != newinv) {

    d_checkpointInterval         = 0;
    d_checkpointTimeStepInterval = newinv;
    d_checkpointWallTimeInterval = 0;

    // Reset the check point so force one to happen.
    d_nextCheckpointTime     = 0.0;
    d_nextCheckpointTimeStep = 0.0;
    d_nextCheckpointWallTime = 0.0;

    // If needed create checkpoints/index.xml
    if (!d_checkpointsDir.exists()) {

      if (d_myworld->myRank() == 0) {
        d_checkpointsDir = d_outputDir.createSubdir("checkpoints");
        createIndexXML(d_checkpointsDir);
      }
    }

    // Sync up before every rank can use the checkpoints dir
    Uintah::MPI::Barrier(d_myworld->getComm());
  }
}

// Allow the component to set the checkpoint wall time interval
void
DataArchiver::setCheckpointWallTimeInterval(int newinv)
{
  if (d_checkpointWallTimeInterval != newinv) {

    d_checkpointInterval         = 0;
    d_checkpointTimeStepInterval = 0;
    d_checkpointWallTimeInterval = newinv;

    // Reset the check point so force one to happen.
    d_nextCheckpointTime     = 0.0;
    d_nextCheckpointTimeStep = 0.0;
    d_nextCheckpointWallTime = 0.0;

    // If needed create checkpoints/index.xml
    if (!d_checkpointsDir.exists()) {

      if (d_myworld->myRank() == 0) {
        d_checkpointsDir = d_outputDir.createSubdir("checkpoints");
        createIndexXML(d_checkpointsDir);
      }
    }

    // Sync up before every rank can use the checkpoints dir
    Uintah::MPI::Barrier(d_myworld->getComm());
  }
}

//  This will copy the portions of the timestep.xml from the old uda
//  to the new uda.  Specifically, the sections related to the components.
void
DataArchiver::copy_outputProblemSpec(Dir& fromDir, Dir& toDir)
{
  int dir_timestep = getTimeStepTopLevel(); // could be modified by reduceUda

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0') << dir_timestep;

  // define the from/to directories & files
  std::string fromPath = fromDir.getName() + "/" + tname.str();
  Dir myFromDir        = Dir(fromPath);
  std::string fromFile = fromPath + "/timestep.xml";

  std::string toPath = toDir.getName() + "/" + tname.str();
  Dir myToDir        = Dir(toPath);
  std::string toFile = toPath + "/timestep.xml";

  //  loop over the blocks in timestep.xml
  //  and copy the component related nodes
  ProblemSpecP inputDoc = loadDocument(fromFile);

  for (ProblemSpecP ps = inputDoc->getFirstChild(); ps != nullptr;
       ps              = ps->getNextSibling()) {
    string nodeName = ps->getNodeName();

    if (nodeName == "Meta" || nodeName == "Time" || nodeName == "Grid" ||
        nodeName == "Data") {
      continue;
    }
    std::cout << "   Now copying the XML node (" << std::setw(20) << nodeName
              << ")"
              << " from: " << fromFile << " to: " << toFile << std::endl;
    copySection(myFromDir, myToDir, "timestep.xml", nodeName);
  }
}

// Check to see if a output directory exists for this time step
auto
DataArchiver::outputTimeStepExists(unsigned int dir_timestep) -> bool
{
  if (dir_timestep == (unsigned int)-1) {
    return true;
  }

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0') << dir_timestep;
  std::string timestepDir = d_outputDir.getName() + "/" + tname.str();

  for (auto& var : d_outputTimeStepDirs) {

    if (var == timestepDir) {
      return true;
    }
  }

  return false;
}

// Check to see if a checkpoint exists for this time step
auto
DataArchiver::checkpointTimeStepExists(unsigned int dir_timestep) -> bool
{
  if (dir_timestep == (unsigned int)-1) {
    return true;
  }

  std::ostringstream tname;
  tname << "t" << std::setw(5) << std::setfill('0') << dir_timestep;
  std::string timestepDir = d_checkpointsDir.getName() + "/" + tname.str();

  for (auto& var : d_checkpointTimeStepDirs) {

    if (var == timestepDir) {
      return true;
    }
  }

  return false;
}

// If your using reduceUda then use use a mapping that's defined in
// reduceUdaSetup()
auto
DataArchiver::getTimeStepTopLevel() -> int
{
  const int timestep =
    d_simulator->getTimeStep() -
    int(d_outputPreviousTimeStep || d_checkpointPreviousTimeStep);

  if (d_doPostProcessUda) {
    return d_restartTimeStepIndices[timestep];
  } else {
    return timestep;
  }
}

// Called to set the elapsed wall time
void
DataArchiver::setElapsedWallTime(double val)
{
  // Set the elasped wall time if the checkpoint is based on the wall time.
  if (d_checkpointWallTimeInterval > 0) {

    // When using the wall time, rank 0 determines the wall time and
    // sends it to all other ranks.
    Uintah::MPI::Bcast(&val, 1, MPI_DOUBLE, 0, d_myworld->getComm());

    d_elapsedWallTime = val;
  }
}

// Called to set the cycle, i.e. how many checkpoints to save
void
DataArchiver::setCheckpointCycle(int val)
{
  d_checkpointCycle = val;
}

//    This is a good idea and should be looked into with git  -Todd
void
DataArchiver::saveSVNinfo()
{
#if 0
  string svn_diff_file =
    string( sci_getenv("SCIRUN_OBJDIR") ) + "/svn_diff.txt";

  if( !validFile( svn_diff_file ) ) {
    cout << "\n"
         << "WARNING: 'svn diff' file '" << svn_diff_file
         << "' does not appear to exist!\n"
         << "\n";
  }
  else {
    string svn_diff_out = d_outputDir.getName() + "/svn_diff.txt";
    string svn_diff_on = string( sci_getenv("SCIRUN_OBJDIR") ) + "/.do_svn_diff";
    if( !validFile( svn_diff_on ) ) {
      cout << "\n"
           << "WARNING: Adding 'svn diff' file to UDA, "
           << "but AUTO DIFF TEXT CREATION is OFF!\n"
           << "         svn_diff.txt may be out of date!  "
           << "Saving as 'possible_svn_diff.txt'.\n"
           << "\n";

      svn_diff_out = d_outputDir.getName() + "/possible_svn_diff.txt";
    }

    copyFile( svn_diff_file, svn_diff_out );
  }
#endif
}

// Verifies that all processes can see the same file system (as rank 0).
void
DataArchiver::setupSharedFileSystem()
{
  Timers::Simple timer;
  timer.start();

  // Verify that all MPI processes can see the common file system (with rank 0).
  std::string fs_test_file_name;

  if (d_myworld->myRank() == 0) {

    d_writeMeta = true;
    // Create a unique file name, using hostname + pid
    char hostname[MAXHOSTNAMELEN];
    if (gethostname(hostname, MAXHOSTNAMELEN) != 0) {
      strcpy(hostname, "hostname-unknown");
    }

    std::ostringstream test_filename_stream;
    test_filename_stream << "vaango_filesystem_test-" << hostname << "-"
                         << getpid();

    // Create the test file...
    FILE* tmpout = fopen(test_filename_stream.str().c_str(), "w");
    if (!tmpout) {
      throw ErrnoException("fopen failed for " + test_filename_stream.str(),
                           errno,
                           __FILE__,
                           __LINE__);
    }
    fprintf(tmpout, "\n"); // Test writing to file...

    if (fflush(tmpout) != 0) {
      throw ErrnoException("fflush", errno, __FILE__, __LINE__);
    }

    if (fsync(fileno(tmpout)) != 0) { // Test syncing a file.
      throw ErrnoException("fsync", errno, __FILE__, __LINE__);
    }

    if (fclose(tmpout) != 0) { // Test closing the file.
      throw ErrnoException("fclose", errno, __FILE__, __LINE__);
    }

    // While the following has never before been necessary, it turns out that
    // the "str()" operator on an ostringstream creates a temporary buffer
    // that can be deleted at any time and so using ".c_str()" on it may return
    // garbage.  To avoid this, we need to copy the ".str()" output into our
    // own string, and then use the ".c_str()" on that non-temporary string.
    const std::string temp_string = test_filename_stream.str();

    const char* outbuf = temp_string.c_str();
    int outlen         = (int)strlen(outbuf);

    // Broadcast test filename length, and then broadcast the actual name.
    Uintah::MPI::Bcast(&outlen, 1, MPI_INT, 0, d_myworld->getComm());
    Uintah::MPI::Bcast(
      const_cast<char*>(outbuf), outlen, MPI_CHAR, 0, d_myworld->getComm());
    fs_test_file_name = test_filename_stream.str();
  } else {
    d_writeMeta = false; // Only rank 0 will emit meta data...

    // All other ranks receive from rank 0 (code above) the length,
    // and then name of the file that we are going to look for...
    int inlen;
    Uintah::MPI::Bcast(&inlen, 1, MPI_INT, 0, d_myworld->getComm());
    char* inbuf = scinew char[inlen + 1];
    Uintah::MPI::Bcast(inbuf, inlen, MPI_CHAR, 0, d_myworld->getComm());
    inbuf[inlen]      = '\0';
    fs_test_file_name = inbuf;

    delete[] inbuf;
  }

  if (d_myworld->myRank() !=
      0) { // Make sure everyone else can see the temp file...

    struct stat st;
    int s = stat(fs_test_file_name.c_str(), &st);

    if ((s != 0) || !S_ISREG(st.st_mode)) {
      std::cerr << "Stat'ing of file: " << fs_test_file_name
                << " failed with errno = " << errno << "\n";
      throw ErrnoException("stat", errno, __FILE__, __LINE__);
    }
  }

  Uintah::MPI::Barrier(d_myworld->getComm()); // Wait until everyone has check
                                              // for the file before proceeding.

  if (d_writeMeta) {

    int s = unlink(fs_test_file_name.c_str()); // Remove the tmp file...
    if (s != 0) {
      std::cerr << "Cannot unlink file: " << fs_test_file_name << '\n';
      throw ErrnoException("unlink", errno, __FILE__, __LINE__);
    }

    makeVersionedDir();
    // Send UDA name to all other ranks.
    std::string udadirname = d_outputDir.getName();

    // Broadcast uda dir name length, and then broadcast the actual name.
    const char* outbuf = udadirname.c_str();
    int outlen         = (int)strlen(outbuf);

    Uintah::MPI::Bcast(&outlen, 1, MPI_INT, 0, d_myworld->getComm());
    Uintah::MPI::Bcast(
      const_cast<char*>(outbuf), outlen, MPI_CHAR, 0, d_myworld->getComm());
  } else {

    // Receive the name of the UDA from rank 0...
    int inlen;
    Uintah::MPI::Bcast(&inlen, 1, MPI_INT, 0, d_myworld->getComm());
    char* inbuf = scinew char[inlen + 1];
    Uintah::MPI::Bcast(inbuf, inlen, MPI_CHAR, 0, d_myworld->getComm());
    inbuf[inlen] = '\0';

    d_outputDir = Dir(inbuf);
    delete[] inbuf;
  }

  if (d_myworld->myRank() == 0) {
    std::cerr << "Verified shared file system in " << timer().seconds()
              << " seconds.\n";
  }

} // end setupSharedFileSystem()

//______________________________________________________________________
//
// setupLocalFileSystems()
//
// This is the old method of checking for shared vs local file systems
// and determining which node(s) will write the UDA meta data.  Rank 0
// creates a file name, sends it to all other ranks.  All other ranks
// use that basename to create subsequent files ("basename-rank").
// Then all ranks (except 0) start looking for these files in order
// (starting one past their rank).  If a file is found, it means that
// the rank is on a node that has another (lower) rank - and that
// lower rank will do the writing.
//
void
DataArchiver::setupLocalFileSystems()
{
  Timers::Simple timer;
  timer.start();

  // See how many shared filesystems that we have
  string basename;
  if (d_myworld->myRank() == 0) {
    // Create a unique string, using hostname+pid
    char* base = strdup(d_filebase.c_str());
    char* p    = base + strlen(base);

    for (; p >= base; p--) {
      if (*p == '/') {
        *++p = 0; // keep trailing slash
        break;
      }
    }

    if (*base) {
      free(base);
      base = strdup(".");
    }

    char hostname[MAXHOSTNAMELEN];
    if (gethostname(hostname, MAXHOSTNAMELEN) != 0) {
      strcpy(hostname, "unknown???");
    }

    std::ostringstream ts;
    ts << base << "-" << hostname << "-" << getpid();
    if (*base) {
      free(base);
    }

    std::string test_string = ts.str();
    const char* outbuf      = test_string.c_str();
    int outlen              = (int)strlen(outbuf);

    Uintah::MPI::Bcast(&outlen, 1, MPI_INT, 0, d_myworld->getComm());
    Uintah::MPI::Bcast(
      const_cast<char*>(outbuf), outlen, MPI_CHAR, 0, d_myworld->getComm());
    basename = test_string;
  } else {
    int inlen;
    Uintah::MPI::Bcast(&inlen, 1, MPI_INT, 0, d_myworld->getComm());
    char* inbuf = scinew char[inlen + 1];
    Uintah::MPI::Bcast(inbuf, inlen, MPI_CHAR, 0, d_myworld->getComm());
    inbuf[inlen] = '\0';
    basename     = inbuf;
    delete[] inbuf;
  }

  // Create a file, of the name p0_hostname-p0_pid-processor_number
  std::ostringstream myname;
  myname << basename << "-" << d_myworld->myRank() << ".tmp";
  std::string fname = myname.str();

  // This will be an empty file, everything is encoded in the name anyway
  FILE* tmpout = fopen(fname.c_str(), "w");
  if (!tmpout) {
    throw ErrnoException(
      "fopen failed for " + fname, errno, __FILE__, __LINE__);
  }
  fprintf(tmpout, "\n");
  if (fflush(tmpout) != 0) {
    throw ErrnoException("fflush", errno, __FILE__, __LINE__);
  }

  if (fsync(fileno(tmpout)) != 0) {
    throw ErrnoException("fsync", errno, __FILE__, __LINE__);
  }

  if (fclose(tmpout) != 0) {
    throw ErrnoException("fclose", errno, __FILE__, __LINE__);
  }

  Uintah::MPI::Barrier(d_myworld->getComm());

  // See who else we can see
  d_writeMeta = true;
  int i;

  for (i = 0; i < d_myworld->myRank(); i++) {
    std::ostringstream name;
    name << basename << "-" << i << ".tmp";
    struct stat st;
    int s = stat(name.str().c_str(), &st);
    if (s == 0 && S_ISREG(st.st_mode)) {
      // File exists, we do NOT need to emit metadata
      d_writeMeta = false;
      break;
    } else if (errno != ENOENT) {
      std::cerr << "Cannot stat file: " << name.str() << ", errno=" << errno
                << '\n';
      throw ErrnoException("stat", errno, __FILE__, __LINE__);
    }
  }

  Uintah::MPI::Barrier(d_myworld->getComm());

  if (d_writeMeta) {
    makeVersionedDir();
    std::string fname = myname.str();
    FILE* tmpout      = fopen(fname.c_str(), "w");
    if (!tmpout) {
      throw ErrnoException("fopen", errno, __FILE__, __LINE__);
    }
    std::string dirname = d_outputDir.getName();
    fprintf(tmpout, "%s\n", dirname.c_str());
    if (fflush(tmpout) != 0) {
      throw ErrnoException("fflush", errno, __FILE__, __LINE__);
    }
#if defined(__APPLE__)
    if (fsync(fileno(tmpout)) != 0) {
      throw ErrnoException("fsync", errno, __FILE__, __LINE__);
    }
#elif !defined(__bgq__) // __bgq__ is defined on Blue Gene Q computers...
    if (fdatasync(fileno(tmpout)) != 0) {
      throw ErrnoException("fdatasync", errno, __FILE__, __LINE__);
    }
#endif
    if (fclose(tmpout) != 0) {
      throw ErrnoException("fclose", errno, __FILE__, __LINE__);
    }
  }

  Uintah::MPI::Barrier(d_myworld->getComm());

  if (!d_writeMeta) {
    std::ostringstream name;
    name << basename << "-" << i << ".tmp";
    std::ifstream in(name.str().c_str());
    if (!in) {
      throw InternalError(
        "DataArchiver::initializeOutput(): The file \"" + name.str() +
          "\" not found on second pass for filesystem discovery!",
        __FILE__,
        __LINE__);
    }
    std::string dirname;
    in >> dirname;
    d_outputDir = Dir(dirname);
  }

  int count = d_writeMeta ? 1 : 0;
  int nunique;
  // This is an AllReduce, not a reduce.  This is necessary to
  // ensure that all processors wait before they remove the tmp files
  Uintah::MPI::Allreduce(
    &count, &nunique, 1, MPI_INT, MPI_SUM, d_myworld->getComm());
  if (d_myworld->myRank() == 0) {
    std::cerr << "Discovered " << nunique << " unique filesystems in "
              << timer().seconds() << " seconds\n";
  }
  // Remove the tmp files...
  int s = unlink(myname.str().c_str());
  if (s != 0) {
    std::cerr << "Cannot unlink file: " << myname.str() << '\n';
    throw ErrnoException("unlink", errno, __FILE__, __LINE__);
  }
}

} // namespace Uintah