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

#ifndef VAANGO_CCA_COMPONENTS_DATAARCHIVER_DataArchiver_H
#define VAANGO_CCA_COMPONENTS_DATAARCHIVER_DataArchiver_H

#include <CCA/Ports/Output.h>
#include <CCA/Ports/PIDXOutputContext.h>

#include <CCA/Components/Schedulers/RuntimeStatsEnum.h>

#include <Core/Containers/ConsecutiveRangeSet.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/MaterialSetP.h>
#include <Core/OS/Dir.h>
#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Util/Assert.h>

namespace Uintah {

class DataWarehouse;
class SimulationInterface;
class LoadBalancer;
class Level;

//! Handles outputting the data.
class DataArchiver
  : public Output
  , public UintahParallelComponent
{
public:
  inline static bool s_wereSavesAndCheckpointsInitialized{ false };

public:
  DataArchiver(const ProcessorGroup* myworld, int udaSuffix = -1);

  ~DataArchiver() override;

  // Disallow copy and move
  DataArchiver(const DataArchiver&) = delete;
  DataArchiver(DataArchiver&&)      = delete;
  auto
  operator=(const DataArchiver&) -> DataArchiver& = delete;
  auto
  operator=(DataArchiver&&) -> DataArchiver& = delete;

  // Methods for managing the components attached via the ports.
  void
  setComponents([[maybe_unused]] UintahParallelComponent* parent) override
  {
  }

  void
  getComponents() override;

  void
  releaseComponents() override;

  //! Sets up when the DataArchiver will output and what data, according
  //! to params.  Also stores state to keep track of time and timesteps
  //! in the simulation.  (If you only need to use DataArchiver to copy
  //! data, then you can pass a nullptr MaterialManager
  void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_ps,
               const MaterialManagerP& mat_manager) override;

  virtual void
  outputProblemSpec(ProblemSpecP& root_ps);

  //! This function will set up the output for the simulation.  As part
  //! of this it will output the input.xml and index.xml in the uda
  //! directory.  Call after calling all problemSetups.
  void
  initializeOutput(const ProblemSpecP& params, const GridP& grid) override;

  //! Call this when restarting from a checkpoint after calling
  //! problemSetup.  This will copy timestep directories and dat
  //! files up to the specified timestep from restartFromDir if
  //! fromScratch is false and will set time and timestep variables
  //! appropriately to continue smoothly from that timestep.
  //! If timestep is negative, then all timesteps will getg copied
  //! if they are to be copied at all (fromScratch is false).
  void
  restartSetup(Dir& restartFromDir,
               int startTimestep,
               int timestep,
               double time,
               bool fromScratch,
               bool removeOldDir) override;

  //! Call this after calling problemSetup.  It will copy the data
  //! and checkpoint files over and make it ignore
  //! dumping reduction variables.
  void
  postProcessUdaSetup(Dir& fromDir) override;

  //! Copy a section between udas' index.xml.
  void
  copySection(Dir& fromDir, Dir& toDir, std::string file, std::string section);

  //! Copy a section from another uda's to our index.xml.
  void
  copySection(Dir& fromDir, std::string section)
  {
    copySection(fromDir, d_outputDir, "index.xml", section);
  }

  //! Checks to see if this is an output timestep.
  //! If it is, setup directories and xml files that we need to output.
  //! Call once per timestep, and if recompiling,
  //! after all the other tasks are scheduled.
  void
  finalizeTimestep(const GridP&,
                   SchedulerP&,
                   bool recompile  = false,
                   int addMaterial = 0) override;

  //! schedule the output tasks if we are recompiling the taskgraph.
  void
  sched_allOutputTasks(const GridP&,
                       SchedulerP&,
                       bool recompile = false) override;

  //! Find the next times to output
  //! Call after timestep has completed.
  void
  findNext_OutputCheckPointTimestep(bool restart, const GridP&) override;

  // Called after a time step recompute where delta t is adjusted
  // to make sure an output and/or checkpoint time step is needed.
  void
  recompute_OutputCheckPointTimestep() override;

  //! Write metadata to xml files.
  //! Call after timestep has completed.
  void
  writeto_xml_files(const GridP& grid) override;

  //! Returns as a string the name of the top of the output directory.
  const std::string
  getOutputLocation() const override;

  // Returns bool, does the outputdir exist
  bool
  doesOutputDirExist() const override;

  // Normally saved vars are scrubbed if not needed for the next
  // time step. Bypass scubbing when running in situ or if wanting
  // to save the previous time step.
  void
  setScrubSavedVariables(bool val) override
  {
    d_scrubSavedVariables = val;
  };

  //! Asks if we need to recompile the task graph.
  bool
  needRecompile(const GridP& grid) override;

  void
  recompile(const GridP& grid) override;

  //! The task that handles the outputting.  Scheduled in finalizeTimestep.
  //! Handles outputs and checkpoints and differentiates between them in the
  //! last argument.  Outputs as binary the data acquired from VarLabel in
  //! p_dir.
  void
  outputVariables(const ProcessorGroup*,
                  const PatchSubset* patch,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw,
                  int type);

  //! Task that handles outputting non-checkpoint variables.
  //! Scheduled in finalizeTimestep.
  void
  outputGlobalVars(const ProcessorGroup*,
                   const PatchSubset* patch,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  // Get the time the next output will occur
  double
  getNextOutputTime() const override
  {
    return d_nextOutputTime;
  }

  // Get the time step the next output will occur
  int
  getNextOutputTimestep() const override
  {
    return d_nextOutputTimestep;
  }

  // Pushes output back by one time step.
  void
  postponeNextOutputTimestep() override
  {
    ++d_nextOutputTimestep;
  }

  // Get the time/time step/wall time of the next checkpoint will occur
  double
  getNextCheckpointTime() const override
  {
    return d_nextCheckpointTime;
  }

  int
  getNextCheckpointTimestep() const override
  {
    return d_nextCheckpointTimestep;
  }

  int
  getNextCheckpointWallTime() const override
  {
    return d_nextCheckpointWallTime;
  }

  // Returns true if data will be output this time step
  void
  setOutputTimestep(bool val, const GridP& grid) override;

  bool
  isOutputTimestep() const override
  {
    return d_isOutputTimestep;
  }

  // Returns true if data will be checkpointed this time step
  void
  setCheckpointTimestep(bool val, const GridP& grid) override;

  bool
  isCheckpointTimestep() const override
  {
    return d_isCheckpointTimestep;
  }

  //! Get the directory of the current time step for outputting info.
  const std::string&
  getLastTimestepOutputLocation() const override
  {
    return d_lastTimestepLocation;
  }

  bool
  isLabelSaved(const std::string& label) const override;

  //! Allow a component to define the output and checkpoint interval on the fly.
  void
  setOutputInterval(double inv) override;

  double
  getOutputInterval() const override
  {
    return d_outputInterval;
  }

  void
  setOutputTimestepInterval(int inv) override;

  int
  getOutputTimestepInterval() const override
  {
    return d_outputTimestepInterval;
  }

  void
  setCheckpointInterval(double inv) override;

  double
  getCheckpointInterval() const override
  {
    return d_checkpointInterval;
  }

  void
  setCheckpointTimestepInterval(int inv) override;

  int
  getCheckpointTimestepInterval() const override
  {
    return d_checkpointTimestepInterval;
  }

  void
  setCheckpointWallTimeInterval(int inv) override;

  int
  getCheckpointWallTimeInterval() const override
  {
    return d_checkpointWallTimeInterval;
  }

  bool
  savingAsPIDX() const override
  {
    return (d_outputFileFormat == PIDX);
  }

  // Instructs the DataArchive to save data using the original UDA format or
  // using PIDX.
  void
  setSaveAsUDA() override
  {
    d_outputFileFormat = UDA;
  }

  void
  setSaveAsPIDX() override
  {
    d_outputFileFormat = PIDX;
  }

  void
  maybeLastTimestep(bool val) override
  {
    d_maybeLastTimestep = val;
  };

  bool
  maybeLastTimestep() override
  {
    return d_maybeLastTimestep;
  };

  void
  setSwitchState(bool val) override
  {
    d_switchState = val;
  }

  bool
  getSwitchState() const override
  {
    return d_switchState;
  }

  void
  setElapsedWallTime(double val) override;

  double
  getElapsedWallTime() const override
  {
    return d_elapsedWallTime;
  };

  void
  setCheckpointCycle(int val) override;

  double
  getCheckpointCycle() const override
  {
    return d_checkpointCycle;
  };

  void
  setUseLocalFileSystems(bool val) override
  {
    d_useLocalFileSystems = val;
  };

  bool
  getUseLocalFileSystems() const override
  {
    return d_useLocalFileSystems;
  };

  void
  setRuntimeStats(
    ReductionInfoMapper<RuntimeStatsEnum, double>* runtimeStats) override
  {
    d_runtimeStats = runtimeStats;
  };

  // Returns true if an output or checkpoint exists for the time step
  bool
  outputTimestepExists(unsigned int ts) override;

  bool
  checkpointTimestepExists(unsigned int ts) override;

public:
  //! problemSetup parses the ups file into a list of these
  //! (d_saveLabelNames)
  struct SaveNameItem
  {
    std::string labelName;
    std::string compressionMode;
    ConsecutiveRangeSet matls;
    ConsecutiveRangeSet levels;
  };

  class SaveItem
  {
  public:
    void
    setMaterials(int level_id,
                 const ConsecutiveRangeSet& matls,
                 ConsecutiveRangeSet& prevMatls,
                 MaterialSetP& prevMatlSet);

    [[nodiscard]] auto
    getMaterialSet(int level_id) const -> const MaterialSet*
    {
      return matlSet.at(level_id).get_rep();
    }

    auto
    getMaterialSubset(const Level* level) const -> const MaterialSubset*;

    const VarLabel* label;

    std::map<int, MaterialSetP> matlSet;
  };

private:
  //         PIDX related
#if HAVE_PIDX
  PIDXOutputContext::PIDX_flags d_PIDX_flags; // Contains the knobs & switches

  std::vector<MPI_Comm>
    m_pidxComms; // Array of MPI Communicators for PIDX usage...

  // creates communicator every AMR level required for PIDX
  void
  createPIDXCommunicator(std::vector<SaveItem>& saveLabels,
                         const GridP& grid,
                         SchedulerP& sched,
                         bool isThisACheckpoint);

  // Timestep # of the last time we saved "timestep.xml". -1 == not
  // yet saved. Only save timestep.xml as needed (ie, when a
  // regrid occurs), otherwise a given timestep will refer (symlink)
  // to the last time it was saved.  Note, this is in reference to
  // IO timesteps.  We always generate and save timestep.xml for
  // Checkpoint output.
#endif

  // output the all of the saveLabels in PIDX format
  auto
  saveLabels_PIDX(const ProcessorGroup* pg,
                  const PatchSubset* patches,
                  DataWarehouse* new_dw,
                  int type,
                  std::vector<SaveItem>& saveLabels,
                  const TypeDescription::Type TD,
                  Dir ldir,                   // uda/timeStep/levelIndex
                  const std::string& dirName, // CCVars, SFC*Vars
                  ProblemSpecP& doc) -> size_t;

  // Searches through "saveLabels" and returns all the SaveItems that are of the
  // same "type".
  auto
  findAllVariablesWithType(const std::vector<SaveItem>& saveLabels,
                           const TypeDescription::Type type)
    -> std::vector<DataArchiver::SaveItem>;

  // bulletproofing so user can't save unsupported var type
  void
  isVarTypeSupported(const std::vector<SaveItem>& saveLabels,
                     const std::vector<TypeDescription::Type>& pidxVarTypes);

  // Writes out the <Grid> and <Data> sections into the
  // timestep.xml file by creating a DOM and then writing it out.
  void
  writeGridOriginal(const bool hasGlobals,
                    const GridP& grid,
                    ProblemSpecP rootElem);

  // Writes out the <Grid> and <Data> sections (respectively) to separate files
  // (that are associated with timestep.xml) using a XML streamer.
  void
  writeGridTextWriter(const bool hasGlobals,
                      const std::string& grid_path,
                      const GridP& grid);

  void
  writeDataTextWriter(const bool hasGlobals,
                      const std::string& data_path,
                      const GridP& grid,
                      const std::vector<std::vector<bool>>& procOnLevel);

  // Writes out the <Grid> section (associated with timestep.xml) to a separate
  // binary file.
  void
  writeGridBinary(const bool hasGlobals,
                  const std::string& grid_path,
                  const GridP& grid);

  //     Non PIDX
  //! returns a ProblemSpecP reading the xml file xmlName.
  //! You will need to that you need to call ProblemSpec::releaseDocument
  auto
  loadDocument(std::string xmlName) -> ProblemSpecP;

  //! creates the uda directory with a trailing version suffix
  void
  makeVersionedDir();

  void
  initSaveLabels(SchedulerP& sched, bool initTimestep);

  void
  initCheckpoints(SchedulerP& sched);

  //! helper for beginOutputTimestep - creates and writes
  //! the necessary directories and xml files to begin the
  //! output timestep.
  void
  makeTimestepDirs(Dir& dir,
                   std::vector<SaveItem>& saveLabels,
                   const GridP& grid,
                   std::string* pTimestepDir /* passed back */);

  //! helper for finalizeTimestep - schedules a task for each var's output
  void
  sched_outputVariables(std::vector<SaveItem>& saveLabels,
                        const GridP& grid,
                        SchedulerP& sched,
                        bool isThisCheckpoint);

  //! Helper for finalizeTimestep - determines if, based on the current
  //! time and timestep, this will be an output or checkpoint timestep.
  void
  beginOutputTimestep(const GridP& grid);

  //! helper for initializeOutput - writes the initial index.xml file,
  //! both setting the d_indexDoc var and writing it to disk.
  void
  createIndexXML(Dir& dir);

  //! helper for restartSetup - adds the restart field to index.xml
  void
  addRestartStamp(ProblemSpecP indexDoc, Dir& fromDir, int timestep);

  //! helper for restartSetup - copies the timestep directories AND
  //! timestep entries in index.xml
  void
  copyTimesteps(Dir& fromDir,
                Dir& toDir,
                int startTimestep,
                int maxTimestep,
                bool removeOld,
                bool areCheckpoints = false);

  //! helper for restartSetup - copies the reduction dat files to
  //! new uda dir (from startTimestep to maxTimestep)
  void
  copyDatFiles(Dir& fromDir,
               Dir& toDir,
               int startTimestep,
               int maxTimestep,
               bool removeOld);

  //! add saved global (reduction) variables to index.xml
  void
  indexAddGlobals();

  // setupLocalFileSystems() and setupSharedFileSystem() are used to
  // create the UDA (versioned) directory.  setupLocalFileSystems()
  // is old method of determining which ranks should output UDA
  // metadata and handles the case when each node has its own local
  // file system (as opposed to a shared file system across all
  // nodes). setupLocalFileSystems() will only be used if
  // specifically turned on via a command line arg to sus when
  // running using MPI.
  void
  setupLocalFileSystems();

  void
  setupSharedFileSystem(); // Verifies that all ranks see a shared FS.

  void
  saveSVNinfo();

private:
  enum OutputFileFormat
  {
    UDA,
    PIDX
  };

  //! This is if you want to pass in the uda extension on the command line
  int d_udaSuffix{ -1 };

  bool d_isOutputTimestep{ false };     //!< set if this is an output timestep
  bool d_isCheckpointTimestep{ false }; //!< set if a checkpoint timestep

  //! Wheter or not p.x is saved
  bool d_saveP_x{ false };
  std::string d_particlePositionName{ "p.x" };

  // For postprocessUda
  bool d_doPostProcessUda{ false };

  // Output file format
  OutputFileFormat d_outputFileFormat{ UDA };

  //! index.xml and checkpoints/index.xml
  ProblemSpecP d_XMLIndexDoc{ nullptr };
  ProblemSpecP d_CheckpointXMLIndexDoc{ nullptr };

  // If the <DataArchiver> section of the .ups file contains:
  //   <outputDoubleAsFloat />
  // Then we will set the d_OutputDoubleAsFloat boolean to true
  // and we will try to output floats instead of doubles.
  //
  // NOTE: This does not affect checkpoints as they will
  //       always be outputting doubles for accuracy.
  bool d_outputDoubleAsFloat{ false };

  // This is the number of times the DataArchiver will retry
  // a file system operation before it gives up and throws
  // an exception.
  int d_fileSystemRetrys{ 10 };

  //! The number of levels the DA knows about.  If this changes, we need to
  //! redo output and Checkpoint tasks.
  int d_numLevelsInOutput{ 0 };

  //! Represents whether this proc will output non-processor-specific
  //! files
  bool d_writeMeta{ false };

  // Hacky variable to ensure that PIDX checkpoint and IO tasks that
  // happen to fall on the same time step run in a serialized manner
  // (as it appears that PIDX is not thread safe).  If there was a
  // better way to synchronize tasks, we should do that...
  VarLabel* d_sync_io_label;

  MaterialSubset* d_tmpMatSubset{ nullptr };

private:
  //! i.e., filebase.000
  std::string d_filebase{ "" };

  SimulationInterface* d_simulator{ nullptr };
  LoadBalancer* d_load_balancer{ nullptr };

  //! pointer to simulation state, to get timestep and time info
  MaterialManagerP d_mat_manager;

  // Only one of these should be non-zero.  The value is read from the .ups
  // file.
  double d_outputInterval{ 0 };      // In seconds.
  int d_outputTimestepInterval{ 0 }; // Number of time steps.

  double d_nextOutputTime{ 0 };  // used when d_outputInterval != 0
  int d_nextOutputTimestep{ 0 }; // used when d_outputTimestepInterval != 0

  // Output the last time step.
  bool d_outputLastTimestep{ false };

  Dir d_outputDir; //!< top of uda dir

  //! Whether or not to save the initialization timestep
  bool d_outputInitTimestep{ false };

  //! last timestep dir (filebase.000/t#)
  std::string d_lastTimestepLocation{ "invalid" };

  int d_lastOutputOfTimestepXML{ -1 };

  //! string for uda dir (actual dir will have postpended numbers
  // List of current output dirs
  std::list<std::string> d_outputTimestepDirs;

  double d_elapsedWallTime{ 0 };
  bool d_maybeLastTimestep{ false };

  //! Whether or not particle vars are saved
  //! Requires p.x to be set
  bool d_saveParticleVariables{ false };

  bool d_switchState{ false };

  // Tells the data archiver that we are running with each MPI node
  // having a separate file system.  (Simulation defaults to running
  // on a shared file system.)
  bool d_useLocalFileSystems{ false };

  //! d_saveLabelNames is a temporary list containing VarLabel
  //! names to be saved and the materials to save them for.  The
  //! information will be basically transferred to d_saveLabels or
  //! d_saveReductionLabels after mapping VarLabel names to their
  //! actual VarLabel*'s.
  std::list<SaveNameItem> d_saveLabelNames;
  std::vector<SaveItem> d_saveLabels;
  std::vector<SaveItem> d_saveGlobalLabels;

  // for efficiency of SaveItem's
  ConsecutiveRangeSet d_prevMatls;
  MaterialSetP d_prevMatlSet{ nullptr };

  //! d_checkpointLabelNames is a temporary list containing
  //! the names of labels to save when checkpointing
  std::vector<SaveItem> d_checkpointLabels;
  std::vector<SaveItem> d_checkpointReductionLabels;

  // Only one of these should be non-zero.
  double d_checkpointInterval{ 0 };      // In seconds.
  int d_checkpointTimestepInterval{ 0 }; // In seconds.

  // How much real time (in seconds) to wait for checkpoint can be
  // used with or without one of the above two.  WalltimeStart
  // cannot be used without walltimeInterval.
  int d_checkpointWallTimeStart{ 0 };    // Amount of (real) time to wait before
                                         // first checkpoint.
  int d_checkpointWallTimeInterval{ 0 }; // Amount of (real) time to between
                                         // checkpoints.

  //! How many checkpoint dirs to keep around
  int d_checkpointCycle{ 2 };

  //! Top of checkpoints dir
  Dir d_checkpointsDir{ "" };

  //! List of current checkpoint dirs
  std::list<std::string> d_checkpointTimestepDirs;
  double d_nextCheckpointTime{ 0 }; //!< used when d_checkpointInterval != 0

  //!< used when d_checkpointTimestepInterval != 0
  int d_nextCheckpointTimestep{ 0 };

  //!< used when d_checkpointWalltimeInterval != 0
  int d_nextCheckpointWallTime{ 0 };

  bool d_outputPreviousTimestep{ false };
  bool d_checkpointPreviousTimestep{ false };
  bool d_checkpointLastTimestep{ false };

  //-----------------------------------------------------------
  // RNJ -
  //
  // In order to avoid having to open and close index.xml,
  // p<xxxxx>.xml, and p<xxxxx>.data when we want to update
  // each variable, we will keep track of some XML docs and
  // file handles and only open and close them once per
  // timestep if they are needed.
  //-----------------------------------------------------------

  // We need to have two separate XML Index Docs
  // because it is possible to do an output
  // and a checkpoint at the same time.

  ProblemSpecP d_upsFile{ nullptr };

  // Each level needs it's own data file handle
  // and if we are outputting and checkpointing
  // at the same time we need two different sets.
  // Also store the filename for error-tracking purposes.
  std::map<int, std::pair<int, char*>> d_DataFileHandles;
  std::map<int, std::pair<int, char*>> d_CheckpointDataFileHandles;

  // Each level needs it's own XML Data Doc
  // and if we are outputting and checkpointing
  // at the same time we need two different sets.
  std::map<int, ProblemSpecP> d_XMLDataDocs;
  std::map<int, ProblemSpecP> d_CheckpointXMLDataDocs;

  //  used for migrating timestep directories
  std::map<int, int> d_restartTimestepIndices;

  Dir d_fromDir{ "" }; // keep track of the original uda

  void
  copy_outputProblemSpec(Dir& fromDir, Dir& toDir);

  // returns either the top level timestep or if reduceUda is used
  // a value from the index.xml file
  auto
  getTimestepTopLevel() -> int;

  // Normally saved vars are scrubbed if not needed for the next
  // time step. By pass scubbing when running in situ or if wanting
  // to save the previous time step.
  bool d_scrubSavedVariables{ true };

  // The following four variables affect the global var output only.
  // For outputing the sim time and/or time step with the global vars
  bool d_outputGlobalVarsTimestep{ false };
  bool d_outputGlobalVarsSimTime{ true };

  // For modulating the output frequency global vars. By default
  // they are output every time step. Note: Frequency > OnTimestep
  unsigned int d_outputGlobalVarsFrequency{ 1 };
  unsigned int d_outputGlobalVarsOnTimestep{ 0 };

  auto
  TranslateVariableType(std::string type, bool isThisCheckpoint) -> std::string;
  ReductionInfoMapper<RuntimeStatsEnum, double>* d_runtimeStats;

#ifdef HAVE_PIDX
  bool d_pidx_need_to_recompile{ false };
  bool d_pidx_restore_nth_rank{ false };
  int d_pidx_requested_nth_rank{ -1 };
  bool d_pidx_checkpointing{ false };
#endif

#if SCI_ASSERTION_LEVEL >= 2
  // double-check to make sure that DA::output is only called once per level per
  // processor per type
  std::vector<bool> d_outputCalled;
  std::vector<bool> d_checkpointCalled;
  bool d_checkpointGlobalCalled{ false };
#endif

  Uintah::MasterLock d_outputLock;
};

} // End namespace Uintah

#endif
