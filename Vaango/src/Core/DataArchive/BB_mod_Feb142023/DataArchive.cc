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

#include <Core/DataArchive/DataArchive.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/InputContext.h>
#include <CCA/Ports/LoadBalancer.h>

#include <Core/Containers/OffsetArray1.h>

#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/UnknownVariable.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <Core/Math/MiscMath.h>

#include <Core/Util/Assert.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Timers/Timers.hpp>
#include <Core/Util/XMLUtils.h>

#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <sys/param.h>
#include <unistd.h>

namespace {
// Usage: export SCI_DEBUG="BC_dbg:+"
// There can be multiple DataArchives so must be static.
static Uintah::Dout DA_dbg{ "DA_dbg",
                            "DataArchive",
                            "DataArchive debug info",
                            false };
} // namespace

namespace Uintah {

DataArchive::DataArchive(const std::string& filebase,
                         int processor /* = 0 */,
                         int numProcessors /* = 1 */,
                         bool verbose /* = true */)
  : d_file_base(filebase)
  , d_processor(processor)
  , d_num_processors(numProcessors)
{
#ifdef STATIC_BUILD
  if (!s_types_initialized) {
    s_types_initialized = true;
    // For static builds, sometimes the Uintah types (CCVariable, etc) do not
    // get automatically registered to the Uintah type system... this call
    // forces that to happen. proc0cout << "Loading Uintah var types into type
    // system (static build).\n";
    instantiateVariableTypes();
  }
#endif

  if (d_file_base == "") {
    throw InternalError(
      "DataArchive::DataArchive 'filebase' cannot be empty (\"\").",
      __FILE__,
      __LINE__);
  }

  while (d_file_base[d_file_base.length() - 1] == '/') {
    // Remove '/' from the end of the filebase (if there is one).
    d_file_base = d_file_base.substr(0, filebase.length() - 1);
  }

  string index = d_file_base + "/index.xml";
  if (verbose && processor == 0) {
    // proc0cout << "Parsing " << index << "\n";
  }

  d_index_file = fopen(index.c_str(), "r");
  if (d_index_file == nullptr) {
    throw InternalError(
      "DataArchive::DataArchive() failed to open index xml file.",
      __FILE__,
      __LINE__);
  }

  d_global_endianness = "";
  d_global_num_bits   = -1;
  queryEndiannessAndBits(d_index_file, d_global_endianness, d_global_num_bits);

  queryAndSetParticlePositionName(d_index_file);

  d_file_format = NOT_SPECIFIED;
  queryAndSetFileFormat(d_index_file);
}

DataArchive::~DataArchive()
{
  // The d_created_var_labels member variable, is used to keep track of
  // the VarLabels (for each of the data fields found in the data
  // archive we are reading data out of) for which a varLabel does not
  // already exist.  Now that we have read in the data, we no longer
  // need these temporary VarLabels, so delete them to avoid a memory
  // leak.  Note, most of these VarLabels will be 're-created' by
  // individual components when they go to access their data.

  for ([[maybe_unused]] auto& [name, label_p] : d_created_var_labels) {
    VarLabel::destroy(label_p);
  }
}

void
DataArchive::queryAndSetParticlePositionName(FILE* doc)
{
  rewind(doc);

  while (true) {

    std::string line = UintahXML::getLine(doc);
    if (line == "") {
      return;
    } else if (line.compare(0, 18, "<ParticlePosition>") == 0) {
      std::vector<std::string> pieces = UintahXML::splitXMLtag(line);
      d_particle_position_name        = pieces[1];
      return;
    }
  }
}

void
DataArchive::queryAndSetFileFormat(FILE* doc)
{
  rewind(doc);

  while (true) {

    std::string line = UintahXML::getLine(doc);
    if (line == "") { // ??? end of file???

      // bulletproofing
      if (d_file_format == NOT_SPECIFIED) {
        // proc0cout << "Warning: Reading in an UDA that is missing the
        // <outputFormat> tag... "
        //           << "defaulting to type old 'UDA'.\n";
        d_file_format = UDA;
      }
      return;
    } else if (line.compare(0, 14, "<outputFormat>") == 0) {
      std::vector<std::string> pieces = UintahXML::splitXMLtag(line);

      const std::string format = Uintah::string_toupper(pieces[1]);

      if (format == "PIDX") {
        d_file_format = PIDX;
      } else {
        d_file_format = UDA;
      }
      return;
    }
  }
}

// Static, so can be called from either DataArchive or TimeData.
void
DataArchive::queryEndiannessAndBits(FILE* doc, string& endianness, int& numBits)
{
  rewind(doc);
  bool found = ProblemSpec::findBlock("<Meta>", doc);

  if (found) {
    while (true) {

      string line = UintahXML::getLine(doc);
      if (line == "" || line == "</Meta>") {
        return;
      } else {
        std::vector<std::string> pieces = UintahXML::splitXMLtag(line);
        if (pieces[0] == "<endianness>") {
          endianness = pieces[1];
        } else if (pieces[0] == "<nBits>") {
          numBits = atoi(pieces[1].c_str());
        }
      }
    }
  }
}

void
DataArchive::queryProcessors(unsigned int& nProcs)
{
  rewind(d_index_file); // Start looking from the top of the file.

  bool found = ProblemSpec::findBlock("<Uintah_DataArchive>", d_index_file);

  if (!found) {
    throw InternalError("DataArchive::queryProcessors 'Uintah_DataArchive' "
                        "node not found in index.xml",
                        __FILE__,
                        __LINE__);
  }

  while (true) {

    std::string line = UintahXML::getLine(d_index_file);
    if (line == "" || line == "</Uintah_DataArchive>") {
      return;
    } else {
      std::vector<std::string> pieces = UintahXML::splitXMLtag(line);
      if (pieces[0] == "<numberOfProcessors>") {
        nProcs = atoi(pieces[1].c_str());
        return;
      }
    }
  }
}

void
DataArchive::queryTimesteps(std::vector<int>& index,
                            std::vector<double>& times,
                            std::vector<double>& oldDelTs)
{
  Timers::Simple timer;
  timer.start();
  if (d_time_data.size() == 0) {
    d_lock.lock();
    if (d_time_data.size() == 0) {

      rewind(d_index_file);
      bool found = ProblemSpec::findBlock("<timesteps>", d_index_file);

      if (!found) {
        std::ostringstream out;
        out << "DataArchive::queryTimesteps 'timesteps' node not found in "
               "index.xml";
        throw InternalError(out.str(), __FILE__, __LINE__);
      }

      while (true) {

        std::string line = UintahXML::getLine(d_index_file);
        if (line == "" || line == "</timesteps>") {
          break;
        } else if (line.compare(0, 10, "<timestep ") == 0) {

          ProblemSpec ts_doc(line);

          std::map<std::string, std::string> attributes;
          ts_doc.getAttributes(attributes);
          std::string tsfile = attributes["href"];
          if (tsfile == "") {
            throw InternalError(
              "DataArchive::queryTimesteps:timestep href not found",
              __FILE__,
              __LINE__);
          }

          int timestepNumber;
          double currentTime;
          double delT;

          // Usually '.../timestep.xml'
          std::string ts_path_and_filename = d_file_base + "/" + tsfile;
          ProblemSpecP timestepDoc         = 0;

          std::string::size_type deliminator_index = tsfile.find("/");
          std::string tnumber(tsfile, 0, deliminator_index);

          // Usually '.../grid.xml'
          std::string grid_path_and_filename =
            d_file_base + "/" + tnumber + "/" + "grid.xml";

          if (attributes["time"] == "") {
            // This block if for earlier versions of the index.xml file that do
            // not contain time information as attributes of the timestep field.
            throw InternalError(
              "DataArchive::queryTimesteps:Cannot find Time block",
              __FILE__,
              __LINE__);
          } else {
            // This block will read delt and time info from the index.xml file
            // instead of opening every single timestep.xml file to get this
            // information
            std::istringstream timeVal(attributes["time"]);
            std::istringstream delTVal(attributes["oldDelt"]);
            std::istringstream timestepVal(ts_doc.getNodeValue());

            timeVal >> currentTime;
            delTVal >> delT;
            timestepVal >> timestepNumber;

            if (!timeVal || !timestepVal) {
              std::cout
                << "WARNING: DataArchive.cc:  std::stringstream failed...\n";
            }
          }

          d_ts_indices.push_back(timestepNumber);
          d_ts_times.push_back(currentTime);
          d_ts_old_delTs.push_back(delT);
          d_time_data.push_back(TimeData(this, ts_path_and_filename));
        }
      } // end while
    }
    d_lock.unlock();
  }

  index    = d_ts_indices;
  times    = d_ts_times;
  oldDelTs = d_ts_old_delTs;

  DOUT(DA_dbg,
       "DataArchive::queryTimesteps completed in " << timer().seconds()
                                                   << " seconds");
}

DataArchive::TimeData&
DataArchive::getTimeData(int index)
{

  ASSERTRANGE(index, 0, (int)d_time_data.size());
  TimeData& td = d_time_data[index];
  if (!td.d_initialized) {
    td.init();
  }

  std::list<int>::iterator is_cached =
    std::find(last_N_timesteps.begin(), last_N_timesteps.end(), index);

  if (is_cached != last_N_timesteps.end()) {
    // It's in the list, so yank it in preperation for putting it at
    // the top of the list.
    DOUT(DA_dbg, "Already cached, putting at top of list.\n");
    last_N_timesteps.erase(is_cached);
  } else {
    DOUT(DA_dbg, "Not in list.\n");
    // Not in the list.  If the list is maxed out, purge the cache
    // of the last item by removing it from the list.  If
    // timestep_cache_size is <= 0, there is an unlimited size to
    // the cache, so don't purge.
    DOUT(DA_dbg,
         "timestep_cache_size = " << timestep_cache_size
                                  << ", last_N_timesteps.size() = "
                                  << last_N_timesteps.size() << "\n");
    if (timestep_cache_size > 0 &&
        (int)(last_N_timesteps.size()) >= timestep_cache_size) {
      int cacheTimestep = last_N_timesteps.back();
      last_N_timesteps.pop_back();
      DOUT(DA_dbg, "Making room.  Purging index " << cacheTimestep << "\n");
      d_time_data[cacheTimestep].purgeCache();
    }
  }

  // Finally insert our new candidate at the top of the list.
  last_N_timesteps.push_front(index);
  return td;
}

int
DataArchive::queryPatchwiseProcessor(const Patch* patch, int index)
{
  d_lock.lock();
  TimeData& timedata = getTimeData(index);

  int proc =
    timedata.d_patchInfo[patch->getLevel()->getIndex()][patch->getLevelIndex()]
      .proc;
  d_lock.unlock();
  return proc;
}

GridP
DataArchive::queryGrid(int index, const ProblemSpecP& ups, bool assignBCs)
{
  Timers::Simple timer;
  timer.start();

  d_lock.lock();

  TimeData& timedata = getTimeData(index);

  // Based on the timestep path and file name (eg: .../timestep.xml), we need
  // to cut off the associated path so that we can find the path to grid.xml.
  std::string::size_type path_length =
    timedata.d_ts_path_and_filename.rfind("/") + 1;
  std::string path(timedata.d_ts_path_and_filename, 0, path_length);
  std::string grid_filename = path + "grid.xml";

  FILE* fp_grid = fopen(grid_filename.c_str(), "r");

  // Check if the grid.xml is present, and use that, if it isn't, then use the
  // grid information that is stored in timestep.xml.
  bool grid_xml_is_binary = false;
  if (fp_grid == nullptr) {
    // Could not open grid.xml, just go with timestep.xml.
    fp_grid = fopen(timedata.d_ts_path_and_filename.c_str(), "r");
  } else {
    // Determine if the grid is written in ASCII xml, or in binary.
    unsigned int marker = -1;
    [[maybe_unused]] size_t dat = fread(&marker, sizeof(marker), 1, fp_grid);

    if (marker == GRID_MAGIC_NUMBER) {
      grid_xml_is_binary = true;
    } else {
      // FIXME: do we need to reset the file pointer here?
    }
  }

  if (fp_grid == nullptr) {
    throw InternalError("DataArchive::queryGrid() failed to open input file.\n",
                        __FILE__,
                        __LINE__);
  }

  GridP grid = scinew Grid();
  std::vector<std::vector<int>> procMap; // One vector<int> per level.

  if (grid_xml_is_binary) {
    grid->readLevelsFromFileBinary(fp_grid, procMap);
  } else {
    grid->readLevelsFromFile(fp_grid, procMap);
  }

  fclose(fp_grid);

  // Check to see if the grid has already been reconstructed and that
  // the cell scaling has not changed. Cell scale check can be removed
  // if Uintah is no longer used for visualization
  if (timedata.d_grid != nullptr) {
    d_lock.unlock();
    return timedata.d_grid;
  }

  if (ups && assignBCs) { // 'ups' is non-null only for restarts.

    ProblemSpecP grid_ps = ups->findBlock("Grid");
    grid->assignBCS(grid_ps, nullptr);
  }

  timedata.d_patchInfo.clear();
  timedata.d_matlInfo.clear();

  for (int levelIndex = 0; levelIndex < grid->numLevels(); levelIndex++) {

    // Initialize timedata with empty vectors for this level:
    timedata.d_patchInfo.push_back(std::vector<PatchData>());
    timedata.d_matlInfo.push_back(std::vector<bool>());

    // Now pull out the patch processor information that we got during the grid
    // creation and put it in the timedata struct.
    std::vector<int>& procMapForLevel = procMap[levelIndex];

    for (auto proc_map : procMapForLevel) {
      PatchData pi;
      pi.proc = proc_map;
      timedata.d_patchInfo[levelIndex].push_back(pi);
    }
  }

  timedata.d_grid = grid;

  d_lock.unlock();
  grid->performConsistencyCheck();

  timedata.d_grid = grid;

  DOUT(DA_dbg,
       "DataArchive::queryGrid completed in " << timer().seconds()
                                              << " seconds");
  return grid;
}

void
DataArchive::queryLifetime(double& /*min*/, double& /*max*/, particleId /*id*/)
{
  std::cerr << "DataArchive::lifetime not finished\n";
}

void
DataArchive::queryLifetime(double& /*min*/,
                           double& /*max*/,
                           const Patch* /*patch*/)
{
  std::cerr << "DataArchive::lifetime not finished\n";
}

// PIDX hack: fix me... how to store this info correctly?
std::unordered_map<std::string, ConsecutiveRangeSet> var_materials;

void
DataArchive::queryVariables(std::vector<std::string>& names,
                            std::vector<int>& num_matls,
                            std::vector<const Uintah::TypeDescription*>& types)
{
  Timers::Simple timer;
  timer.start();

  d_lock.lock();

  rewind(d_index_file); // Start at beginning of file.
  bool found = ProblemSpec::findBlock("<variables>", d_index_file);
  if (!found) {
    throw InternalError(
      "DataArchive::queryVariables:variables section not found\n",
      __FILE__,
      __LINE__);
  }

  queryVariables(d_index_file, names, num_matls, types);

  // PIDX hack:
  if (var_materials.size() == 0) {

    for (unsigned int i = 0; i < num_matls.size(); i++) {

      ConsecutiveRangeSet set;
      for (int j = 0; j < num_matls[i]; j++) {
        set.addInOrder(j);
      }
      var_materials[names[i]] = set;
    }
  }
  // end PIDX hack.

  d_lock.unlock();
  DOUT(DA_dbg,
       "DataArchive::queryVariables completed in " << timer().seconds()
                                                   << " seconds\n");
}

void
DataArchive::queryGlobals(std::vector<std::string>& names,
                          std::vector<int>& num_matls,
                          std::vector<const Uintah::TypeDescription*>& types)
{
  Timers::Simple timer;
  timer.start();

  d_lock.lock();

  rewind(d_index_file); // Start looking from the top of the file.

  bool result = ProblemSpec::findBlock("<globals>", d_index_file);
  if (!result) {
    d_lock.unlock();
    return;
  }

  queryVariables(d_index_file, names, num_matls, types, true);

  d_lock.unlock();

  DOUT(DA_dbg,
       "DataArchive::queryGlobals completed in " << timer().seconds()
                                                 << " seconds\n");
}

void
DataArchive::queryVariables([[maybe_unused]] FILE* fp,
                            std::vector<std::string>& names,
                            std::vector<int>& num_matls,
                            std::vector<const Uintah::TypeDescription*>& types,
                            bool globals /* = false */)
{
  // Assuming that fp points to the line following "<variables>"...
  std::string end_block;
  if (globals) {
    end_block = "</globals>";
  } else {
    end_block = "</variables>";
  }

  while (true) {

    std::string line = UintahXML::getLine(d_index_file);
    if (line == "" || line == end_block) {
      break;
    } else if (line.compare(0, 10, "<variable ") == 0) {

      ProblemSpec ts_doc(line);

      std::map<std::string, std::string> attributes;
      ts_doc.getAttributes(attributes);
      std::string the_type = attributes["type"];
      if (the_type == "") {
        throw InternalError("DataArchive::queryVariables() - 'type' not found",
                            __FILE__,
                            __LINE__);
      }

      const TypeDescription* td = TypeDescription::lookupType(the_type);

      if (!td) {
        static TypeDescription* unknown_type = 0;
        if (!unknown_type) {
          unknown_type = scinew TypeDescription(TypeDescription::Type::Unknown,
                                                "-- unknown type --",
                                                false,
                                                MPI_Datatype(-1));
        }
        td = unknown_type;
      }

      types.push_back(td);
      std::string name = attributes["name"];
      if (name == "") {
        throw InternalError("DataArchive::queryVariables() - 'name' not found",
                            __FILE__,
                            __LINE__);
      }
      names.push_back(name);

      std::string num_matls_string = attributes["numMaterials"];
      if (num_matls_string == "") {
        num_matls.push_back(
          -1); // Either a global var (thus only one material), or an old UDA
               // that did not save this info...
      } else {
        num_matls.push_back(atoi(num_matls_string.c_str()));
      }

    } else {
      throw InternalError(
        "DataArchive::queryVariables() - bad data in variables block.",
        __FILE__,
        __LINE__);
    }
  }
}

#if HAVE_PIDX
// Notes: You must free() the returned "dataPIDX" when you are done with it.

bool
DataArchive::setupQueryPIDX(PIDX_access& access,
                            PIDX_file& idxFile,
                            PIDX_variable& varDesc,
                            const LevelP& level,
                            const TypeDescription* type,
                            const string& name,
                            const int matlIndex,
                            const int timeIndex)
{
  TimeData& timedata = getTimeData(timeIndex);

  //__________________________________
  //  Creating access
  PIDX_create_access(&access);

  if (Parallel::usingMPI()) {
    if (d_pidx_comms.size() == 0) {
      // cout << "No pidx comms, using main comm on this: " << this << "\n";
      // this is a hack because createPIDXCommunicator() was not called in the
      // case of running under (at least) compare_uda... need to figure out the
      // right way to do this.
      MPI_Comm comm = Parallel::getRootProcessorGroup()->getComm();
      PIDX_set_mpi_access(access, comm);
    } else {
      // cout << Uintah::Parallel::getMPIRank() << ": setting up pidx comm\n";
      PIDX_set_mpi_access(access, d_pidx_comms[level->getIndex()]);
      // cout << Uintah::Parallel::getMPIRank() << ": done setting up pidx
      // comm\n";
    }
  }
  //__________________________________
  //  Open idx file
  std::ostringstream levelPath;
  levelPath << timedata.d_ts_directory << "l" << level->getIndex()
            << "/"; // uda/timestep/level/

  std::string typeStr;
  if (type->getType() == TypeDescription::Type::CCVariable) {
    typeStr = "CCVars";
  } else if (type->getType() == TypeDescription::Type::SFCXVariable) {
    typeStr = "SFCXVars";
  } else if (type->getType() == TypeDescription::Type::SFCYVariable) {
    typeStr = "SFCYVars";
  } else if (type->getType() == TypeDescription::Type::SFCZVariable) {
    typeStr = "SFCZVars";
  } else if (type->getType() == TypeDescription::Type::NCVariable) {
    typeStr = "NCVars";
  } else if (type->getType() == TypeDescription::Type::ParticleVariable) {
    typeStr = "ParticleVars";
  } else {
    typeStr = "NOT_IMPLEMENTED";
  }

  std::string idxFilename = levelPath.str() + typeStr + ".idx";

  PIDX_point global_size;

  // cout << Uintah::Parallel::getMPIRank() << ": open pidx file: " <<
  // idxFilename << ", looking for var: " << name << "\n";

  int ret = PIDX_file_open(
    idxFilename.c_str(), PIDX_MODE_RDONLY, access, global_size, &idxFile);
  // PIDXOutputContext::checkReturnCode( ret,"DataArchive::setupQueryPIDX() -
  // PIDX_file_open failure", __FILE__, __LINE__ );
  if (ret != PIDX_success) {
    // Is this a good idea?  Or are we potentially masking true problems.
    cout << Uintah::Parallel::getMPIRank()
         << ": pidx file not found so ending pidx portion of QUERY() ON PROC\n";
    return false;
  }

  //__________________________________
  //  Extra Calls that _MAY_ be needed
  // PIDX_point global_size;
  // ret = PIDX_get_dims(idxFile, global_size);          // returns the
  // levelSize  Is this needed? PIDXOutputContext::checkReturnCode(
  // ret,"DataArchive::setupQueryPIDX() - PIDX_get_dims failure", __FILE__,
  // __LINE__);

  // cout << Uintah::Parallel::getMPIRank() << ": pidx get var count\n";

  int variable_count = 0; ///< Number of fields in PIDX file
  ret                = PIDX_get_variable_count(idxFile, &variable_count);
  PIDXOutputContext::checkReturnCode(
    ret,
    "DataArchive::setupQueryPIDX() - PIDX_get_variable_count failure",
    __FILE__,
    __LINE__);

  // cout << Uintah::Parallel::getMPIRank() << ": pidx var count is " <<
  // variable_count << "\n";

  // int me;
  // PIDX_get_current_time_step(idxFile, &me);
  // cout << " PIDX file has currentl timestep: " << me << std::endl;

  //__________________________________
  //  set locations in PIDX file for querying variable
  // cout << Uintah::Parallel::getMPIRank() << ": pidx set time\n";
  int timestep = d_ts_indices[timeIndex];
  ret          = PIDX_set_current_time_step(idxFile, timestep);
  PIDXOutputContext::checkReturnCode(
    ret,
    "DataArchive::setupQueryPIDX() - PIDX_set_current_time_step failure",
    __FILE__,
    __LINE__);

  std::ostringstream mstr;
  mstr << "_m" << matlIndex; // Add _m# to name of variable.
  std::string full_name = name + mstr.str();

  // cout << Uintah::Parallel::getMPIRank() << ": setting var: " << full_name <<
  // "\n";
  ret = PIDX_set_current_variable_by_name(idxFile, full_name.c_str());
  // proc0cout << "ret is " << ret << ", was looking for" << name << "\n";

  if (ret != PIDX_success) {
    // PIDXOutputContext::checkReturnCode( ret, "DataArchive::setupQueryPIDX() -
    // PIDX_set_current_variable_index failure", __FILE__, __LINE__ );
    cout << Uintah::Parallel::getMPIRank()
         << ": variable not found so ending pidx portion of QUERY()\n";
    return false;
  }

  //__________________________________
  // read IDX file for variable desc
  // cout << Uintah::Parallel::getMPIRank() << ": pidx get current var\n";

  ret = PIDX_get_current_variable(idxFile, &varDesc);

  // cout << Uintah::Parallel::getMPIRank() << ": PIDX GET CURRENT VAR done\n";

  PIDXOutputContext::checkReturnCode(
    ret,
    "DataArchive::setupQueryPIDX() - PIDX_get_current_variable failure",
    __FILE__,
    __LINE__);

  // proc0cout << "Read in PIDX variable: " << varDesc->var_name << ", looking
  // for UDA var: " << name << "\n"; // DEBUG PRINTOUT

#if 0
  if( my_var_name != varDesc->var_name ) {
    throw( "Failed sanity check, pidx var name is different from what I thought it should be..." );
  }
#endif
  return true;

} // end setupQueryPIDX()

void
DataArchive::queryPIDX(BufferAndSizeTuple* data,
                       const PIDX_variable& varDesc,
                       const TypeDescription* type,
                       const string& name,
                       const int matlIndex,
                       const Patch* patch,
                       const int timeIndex)
{
  // cout << Uintah::Parallel::getMPIRank()
  //      << ": pidx query called for        VARIABLE: " << name
  //      << ", material index " << matlIndex
  //      << ", Level " << (patch ? patch->getLevel()->getIndex() : -1)
  //      << ", patch " << (patch ? patch->getID() : -1)
  //      << ", time index " << timeIndex << "\n";

  if (d_file_format != PIDX) {
    throw InternalError(
      "queryPIDX() called on non-PIDX data archive", __FILE__, __LINE__);
  }

  //__________________________________
  //  bulletproofing
  if (isPIDXEnabled() == false) {
    std::ostringstream error;
    error << "\nERROR DataArchive::queryPIDX()\n"
          << "The uda you are trying to open was written using the PIDX file "
             "format.\n"
          << "However, you did not configure Uintah to support PIDX.\n"
          << "You must re-configure and compile with PIDX enabled.";
    throw InternalError(error.str(), __FILE__, __LINE__);
  }

  // TimeData    & timedata = getTimeData( timeIndex );
  // const Patch * real_patch = patch->getRealPatch();
  // int           patchid = real_patch->getID();

  //__________________________________
  //  open PIDX
  //  TO DO:
  //    - do we need  calls to PIDX_get_variable_count() PIDX_get_dims()??

  if (!patch) {
    // bullet proofing
    throw InternalError(
      "Error: queryPIDX() requires a non-null patch.", __FILE__, __LINE__);
  }
  // cout << Uintah::Parallel::getMPIRank() << ": starting pidx portion of
  // query()\n";

  PIDXOutputContext pidx;
  const Level* level = patch->getLevel();

  //__________________________________
  // define the level extents for this variable type
  IntVector lo;
  IntVector hi;
  level->findCellIndexRange(lo, hi);
  PIDX_point level_size;
  pidx.setLevelExtents("DataArchive::queryPIDX()", lo, hi, level_size);

  //__________________________________
  // define patch extents
  PIDX_point patchOffset;
  PIDX_point patchSize;
  PIDXOutputContext::patchExtents patchExts;

  const IntVector boundary_layer(0, 0, 0);
  pidx.setPatchExtents("DataArchive::queryPIDX()",
                       patch,
                       level,
                       boundary_layer,
                       type,
                       patchExts,
                       patchOffset,
                       patchSize);

  if (dbg.active() && isProc0_macro) {
    patchExts.print(cout);
  }

  // debugging
  //  if (dbg.active()||true ){
  //    cout << Uintah::Parallel::getMPIRank() << ": Query:  file: " <<
  //    idxFilename
  //         << "    " << name
  //         << " ts: " << timestep
  //         << " mi: " <<  matlIndex
  //         << " pID: "   << patchid
  //         << " lev: "     << level->getIndex() << " ------ " // << "\n"
  //         << "    " << varDesc->var_name
  //         << " type_name: " << varDesc->type_name
  //         << " varIndex: " << varIndex
  //         << " vals_/_samp: " << varDesc->vps
  //         << " bits_/_samp: "<< bits_per_sample
  //         << " arrySz " << arraySize << "\n";
  //  }

  int ret;

  if (type->getType() == Uintah::TypeDescription::Type::ParticleVariable) {

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

    ret = PIDX_variable_read_particle_data_layout(varDesc,
                                                  physical_local_offset,
                                                  physical_local_size,
                                                  (void**)&(data->buffer),
                                                  (uint64_t*)&(data->size),
                                                  PIDX_row_major);
    PIDXOutputContext::checkReturnCode(
      ret,
      "DataArchive::queryPIDX() - PIDX_variable_read_particle_data_layout "
      "failure",
      __FILE__,
      __LINE__);
  } else {

    //__________________________________
    // Allocate memory and read in data from PIDX file  Need to use patch_buffer
    // !!!

    int values_per_sample = varDesc->vps;
    int bits_per_sample   = 0;
    int ret =
      PIDX_default_bits_per_datatype(varDesc->type_name, &bits_per_sample);
    PIDXOutputContext::checkReturnCode(
      ret,
      "DataArchive::queryPIDX() - PIDX_default_bits_per_datatype failure",
      __FILE__,
      __LINE__);

    size_t arraySize =
      (bits_per_sample / 8) * patchExts.totalCells_EC * values_per_sample;
    data->buffer = (unsigned char*)malloc(arraySize);
    data->size   = arraySize;
    memset(data->buffer, 0, arraySize);

    IntVector extra_cells = patch->getExtraCells();

    patchOffset[0] = patchOffset[0] - lo.x() - extra_cells[0];
    patchOffset[1] = patchOffset[1] - lo.y() - extra_cells[1];
    patchOffset[2] = patchOffset[2] - lo.z() - extra_cells[2];

    // cout << Uintah::Parallel::getMPIRank() << ": level: " <<
    // level->getIndex() << ", patchoffset: " << patchOffset[0] << ", " <<
    // patchOffset[1] << ", " << patchOffset[2]
    //      << ", patchsize: " << patchSize[0] << ", " << patchSize[1] << ", "
    //      << patchSize[2] << "\n";

    ret = PIDX_variable_read_data_layout(
      varDesc, patchOffset, patchSize, data->buffer, PIDX_row_major);

    PIDXOutputContext::checkReturnCode(
      ret,
      "DataArchive::queryPIDX() - PIDX_variable_read_data_layout failure",
      __FILE__,
      __LINE__);

    //__________________________________
    // debugging
    if (dbg.active()) {
      pidx.printBufferWrap("DataArchive::query    AFTER  close",
                           type->getSubType()->getType(),
                           varDesc->vps,
                           patchExts.lo_EC,
                           patchExts.hi_EC,
                           data->buffer,
                           arraySize);
    }
  }
  // cout << Uintah::Parallel::getMPIRank() << ": ENDING pidx portion of
  // query()\n";

} // end queryPIDX()

//
// For now, wrapping the PIDX setupQueryPIDX() and queryPIDX() functions in
// this convenience function in order to be able to read in a single variable
// at a time.  This is what is used (indirectly) by puda/compare_uda/etc.
//
bool
DataArchive::queryPIDXSerial(Variable& var,
                             const string& name,
                             const int matlIndex,
                             const Patch* patch,
                             const int timeIndex)
{
  // cout << Uintah::Parallel::getMPIRank()
  //      << ": queryPIDXSerial() called for VARIABLE: " << name
  //      << ", material index " << matlIndex
  //      << ", Level " << (patch ? patch->getLevel()->getIndex() : -1)
  //      << ", patch " << (patch ? patch->getID() : -1)
  //      << ", time index " << timeIndex << "\n";

  if (!isPIDXFormat()) {
    throw InternalError(
      "DataArchive::queryPIDXSerial() called on non-PIDX DataArchive...",
      __FILE__,
      __LINE__);
  }

  PIDX_file idxFile;
  PIDX_variable varDesc;
  PIDX_access access;

  const LevelP& level               = patch->getLevelP();
  const Uintah::TypeDescription* td = var.virtualGetTypeDescription();
  BufferAndSizeTuple* data          = new BufferAndSizeTuple();

  if (td->getType() == TypeDescription::ReductionVariable) {
    // Bulletproofing
    throw InternalError(
      "DataArchive::queryPIDXSerial() - should never get here...",
      __FILE__,
      __LINE__);
  }

  bool found = setupQueryPIDX(
    access, idxFile, varDesc, level, td, name, matlIndex, timeIndex);
  if (!found) {
    return false;
  }

  queryPIDX(data, varDesc, td, name, matlIndex, patch, timeIndex);

  int ret = PIDX_close(idxFile);
  PIDXOutputContext::checkReturnCode(
    ret,
    "DataArchive::queryPIDXSerial() - PIDX_close failure",
    __FILE__,
    __LINE__);

  ret = PIDX_close_access(access);
  PIDXOutputContext::checkReturnCode(
    ret,
    "DataArchive::queryPIDXSerial() - PIDX_close_access failure",
    __FILE__,
    __LINE__);

  //__________________________________
  // Allocate memory for grid or particle variables
  if (td->getType() == TypeDescription::Type::ParticleVariable) {

    int numParticles = data->size;

    if (patch->isVirtual()) {
      throw InternalError(
        "DataArchive::query: Particle query on virtual patches "
        "not finished.  We need to adjust the particle positions to virtual "
        "space...",
        __FILE__,
        __LINE__);
    }

    psetDBType::key_type key(matlIndex, patch);
    ParticleSubset* psubset       = 0;
    psetDBType::iterator psetIter = d_pset_DB.find(key);

    if (psetIter != d_pset_DB.end()) {
      psubset = (*psetIter).second.get_rep();
    }

    if (psubset == 0 || (int)psubset->numParticles() != data->size) {
      psubset = scinew ParticleSubset(numParticles, matlIndex, patch);
      cout << "numParticles: " << numParticles << "\n";

      cout << "d_pset size: " << d_pset_DB.size() << "\n";
      cout << "1. key is: " << key.first << "\n";
      cout << "2. key is: " << key.second << "\n";
      d_pset_DB[key] = psubset;
    }
    (static_cast<ParticleVariableBase*>(&var))->allocate(psubset);
    //      (dynamic_cast<ParticleVariableBase*>(&var))->allocate(psubset);
  } else if (td->getType() == TypeDescription::PerPatch ||
             td->getType() == TypeDescription::SoleVariable ||
             td->getType() == TypeDescription::ReductionVariable) {
  } else { // Grid Var
    const IntVector bl(0, 0, 0);
    var.allocate(patch, bl);
  }

  bool swap_bytes = false; // FIX ME this should not be hard coded!

  var.readPIDX(data->buffer, data->size, swap_bytes);

  free(data->buffer);
  free(data);

  return true;
}
#endif // HAVE_PIDX

bool
DataArchive::query(Variable& var,
                   const std::string& name,
                   int matlIndex,
                   const Patch* patch,
                   int index,
                   DataFileInfo* dfi /* = 0 */)
{
  Timers::Simple timer;
  timer.start();

#if !defined(DISABLE_SCI_MALLOC)
  const char* tag = AllocatorSetDefaultTag("QUERY");
#endif

  d_lock.lock();
  TimeData& timedata = getTimeData(index);
  d_lock.unlock();

  ASSERT(timedata.d_initialized);

  // make sure info for this patch gets parsed from p*****.xml.
  d_lock.lock();
  timedata.parsePatch(patch);
  d_lock.unlock();

  VarData& varinfo = timedata.d_varInfo[name];
  string data_filename;
  int patchid;
  VarType varType = BLANK;

  if (patch) {

    if (d_file_format == PIDX) {
#if HAVE_PIDX
      return queryPIDXSerial(var, name, matlIndex, patch, timeIndex);
#else
      throw InternalError(
        "DataArchive::query() called for PIDX UDA - but PIDX not configured.",
        __FILE__,
        __LINE__);
#endif
    }

    varType = PATCH_VAR;

    // we need to use the real_patch (in case of periodic boundaries) to get the
    // data, but we need the passed in patch to allocate the patch to the proper
    // virtual region... (see var.allocate below)
    const Patch* real_patch         = patch->getRealPatch();
    [[maybe_unused]] int levelIndex = real_patch->getLevel()->getIndex();
    [[maybe_unused]] int patchIndex = real_patch->getLevelIndex();

    PatchData& patchinfo =
      timedata.d_patchInfo[real_patch->getLevel()->getIndex()]
                          [real_patch->getLevelIndex()];
    ASSERT(patchinfo.parsed);

    patchid = real_patch->getID();

    std::ostringstream ostr;
    // append l#/datafilename to the directory
    ostr << timedata.d_ts_directory << "l" << patch->getLevel()->getIndex()
         << "/" << patchinfo.datafilename;
    data_filename = ostr.str();
  } else {
    varType = GLOBAL_VAR;

    // reference reduction file 'global.data' will a null patch
    patchid       = -1;
    data_filename = timedata.d_ts_directory + timedata.d_globaldata;
  }

  // on a call from restartInitialize, we already have the information from the
  // dfi, otherwise get it from the hash table info
  DataFileInfo datafileinfo;
  if (!dfi) {
    // if this is a virtual patch, grab the real patch, but only do that here -
    // in the next query, we want the data to be returned in the virtual
    // coordinate space
    std::vector<VarnameMatlPatch>::iterator iter =
      std::find(timedata.d_datafileInfoIndex.begin(),
                timedata.d_datafileInfoIndex.end(),
                VarnameMatlPatch(name, matlIndex, patchid));

    if (iter == timedata.d_datafileInfoIndex.end()) {
      std::cerr
        << "VARIABLE NOT FOUND: " << name << ", material index " << matlIndex
        << ", patch " << patch->getID() << ", time index " << index
        << "\nPlease make sure the correct material index is specified\n";
      throw InternalError(
        "DataArchive::query:Variable not found", __FILE__, __LINE__);
    }

    int pos = std::distance(timedata.d_datafileInfoIndex.begin(), iter);
    dfi     = &timedata.d_datafileInfoValue[pos];
  }

  const TypeDescription* td = var.virtualGetTypeDescription();
  ASSERT(td->getName() == varinfo.type);

  // Allocate memory for grid or particle variables
  if (td->getType() == TypeDescription::Type::ParticleVariable) {
    if (dfi->numParticles == -1) {
      throw InternalError(
        "DataArchive::query:Cannot get numParticles", __FILE__, __LINE__);
    }

    if (patch->isVirtual()) {
      throw InternalError(
        "DataArchive::query: Particle query on virtual patches "
        "not finished.  We need to adjust the particle positions to virtual "
        "space...",
        __FILE__,
        __LINE__);
    }

    psetDBType::key_type key(matlIndex, patch);
    ParticleSubset* psubset       = 0;
    psetDBType::iterator psetIter = d_pset_DB.find(key);

    if (psetIter != d_pset_DB.end()) {
      psubset = (*psetIter).second.get_rep();
    }

    if (psubset == 0 || (int)psubset->numParticles() != dfi->numParticles) {
      psubset = scinew ParticleSubset(dfi->numParticles, matlIndex, patch);
      //      std::cout << "numParticles: " << dfi->numParticles << "\n";
      //      std::cout << "d_pset size: " << d_pset_DB.size() << "\n";
      //      std::cout << "1. key is: " << key.first << "\n";
      //      std::cout << "2. key is: " << key.second << "\n";
      d_pset_DB[key] = psubset;
    }
    (static_cast<ParticleVariableBase*>(&var))->allocate(psubset);
  } else if (td->getType() != TypeDescription::Type::PerPatch ||
             td->getType() != TypeDescription::Type::SoleVariable ||
             td->getType() != TypeDescription::Type::ReductionVariable) {
  } else {
    var.allocate(patch, varinfo.boundaryLayer);
  }

  // open data file Standard Uda Format
  if (d_file_format == UDA || varType == GLOBAL_VAR) {
    int fd = open(data_filename.c_str(), O_RDONLY);

    if (fd == -1) {
      std::cerr << "Error opening file: " << data_filename.c_str()
                << ", errno=" << errno << '\n';
      throw ErrnoException(
        "DataArchive::query (open call)", errno, __FILE__, __LINE__);
    }

    off_t ls = lseek(fd, dfi->start, SEEK_SET);

    if (ls == -1) {
      std::cerr << "Error lseek - file: " << data_filename.c_str()
                << ", errno=" << errno << '\n';
      throw ErrnoException(
        "DataArchive::query (lseek call)", errno, __FILE__, __LINE__);
    }

    InputContext ic(fd, data_filename.c_str(), dfi->start);

    Timers::Simple read_timer;
    read_timer.start();

    var.read(ic,
             dfi->end,
             timedata.d_swapBytes,
             timedata.d_nBytes,
             varinfo.compression);

    DOUT(DA_dbg,
         "DataArchive::query: time to read raw data: " << read_timer().seconds()
                                                       << " seconds");
    ASSERTEQ(dfi->end, ic.cur);

    int s = close(fd);
    if (s == -1) {
      std::cerr << "Error closing file: " << data_filename.c_str()
                << ", errno=" << errno << '\n';
      throw ErrnoException(
        "DataArchive::query (close call)", errno, __FILE__, __LINE__);
    }
  }

#if !defined(DISABLE_SCI_MALLOC)
  AllocatorSetDefaultTag(tag);
#endif

  DOUT(DA_dbg,
       "DataArchive::query() completed in " << timer().seconds() << " seconds");

  return true;
}

bool
DataArchive::query(Variable& var,
                   const string& name,
                   int matlIndex,
                   const Patch* patch,
                   int timeIndex,
                   Ghost::GhostType gt,
                   int numGhostcells)
{
  if (numGhostcells == 0) {
    return query(var, name, matlIndex, patch, timeIndex, nullptr);
  } else {
    d_lock.lock();
    TimeData& td = getTimeData(timeIndex);
    d_lock.unlock();
    td.parsePatch(patch); // make sure vars is actually populated
    if (td.d_varInfo.find(name) != td.d_varInfo.end()) {
      VarData& varinfo            = td.d_varInfo[name];
      const TypeDescription* type = TypeDescription::lookupType(varinfo.type);
      IntVector low, high;
      patch->computeVariableExtents(
        type->getType(), varinfo.boundaryLayer, gt, numGhostcells, low, high);
      queryRegion(
        var, name, matlIndex, patch->getLevel(), timeIndex, low, high);
    } else {
      std::cerr
        << "VARIABLE NOT FOUND: " << name << ", material index " << matlIndex
        << ", patch " << patch->getID() << ", time index " << timeIndex
        << "\nPlease make sure the correct material index is specified\n";
      throw InternalError(
        "DataArchive::query:Variable not found", __FILE__, __LINE__);
    }
    return true;
  }
}

void
DataArchive::queryRegion(Variable& var,
                         const string& name,
                         int matlIndex,
                         const Level* level,
                         int timeIndex,
                         IntVector low,
                         IntVector high)
{
  // NOTE - this is not going to do error checking like making sure the entire
  // volume is filled.
  //   We'll assume that if there were bad regions, they would have been
  //   caught in the simulation.
  GridVariableBase* gridvar = dynamic_cast<GridVariableBase*>(&var);
  ASSERT(gridvar);
  gridvar->allocate(low, high);

  d_lock.lock();
  TimeData& td = getTimeData(timeIndex);
  d_lock.unlock();
  const TypeDescription* type = 0;
  Patch::VariableBasis basis =
    Patch::NodeBased; // not sure if this is a reasonable default...
  Patch::selectType patches;

  level->selectPatches(low, high, patches);
  for (size_t i = 0; i < patches.size(); i++) {
    const Patch* patch = patches[i];

    if (type == 0) {
      td.parsePatch(patch); // make sure varInfo is loaded
      VarData& varinfo = td.d_varInfo[name];
      type             = TypeDescription::lookupType(varinfo.type);
      basis            = Patch::translateTypeToBasis(type->getType(), false);
    }
    IntVector l, h;

    l = Max(patch->getExtraLowIndex(basis, IntVector(0, 0, 0)), low);
    h = Min(patch->getExtraHighIndex(basis, IntVector(0, 0, 0)), high);
    if (l.x() >= h.x() || l.y() >= h.y() || l.z() >= h.z()) {
      continue;
    }
    GridVariableBase* tmpVar = gridvar->cloneType();
    query(*tmpVar, name, matlIndex, patch, timeIndex);

    if (patch->isVirtual()) {
      // if patch is virtual, it is probable a boundary layer/extra cell that
      // has been requested (from AMR) let Bryan know if this doesn't work. We
      // need to adjust the source but not the dest by the virtual offset
      tmpVar->offset(patch->getVirtualOffset());
    }
    try {
      gridvar->copyPatch(tmpVar, l, h);
    } catch (InternalError& e) {
      std::cout << " Bad range: " << low << " " << high
                << ", patch intersection: " << l << " " << h << " actual patch "
                << patch->getLowIndex(basis) << " "
                << patch->getHighIndex(basis)
                << " var range: " << tmpVar->getLow() << " "
                << tmpVar->getHigh() << std::endl;
      throw e;
    }
    delete tmpVar;
  }
}

void
DataArchive::findPatchAndIndex(GridP grid,
                               Patch*& patch,
                               particleIndex& idx,
                               long64 particleID,
                               int matlIndex,
                               int levelIndex,
                               int index)
{
  Patch* local = patch;
  if (patch != nullptr) {
    ParticleVariable<long64> var;
    query(var, "p.particleID", matlIndex, patch, index);
    //  cerr<<"var["<<idx<<"] = "<<var[idx]<<endl;
    if (idx <
          static_cast<particleIndex>(var.getParticleSubset()->numParticles()) &&
        var[idx] == particleID) {
      return;
    } else {
      ParticleSubset* subset = var.getParticleSubset();
      for (auto p_idx : *subset) {
        if (var[p_idx] == particleID) {
          idx = p_idx;
          return;
        }
      }
    }
  }
  patch = nullptr;
  //   for (int level_nr = 0;
  //        (level_nr < grid->numLevels()) && (patch == nullptr); level_nr++)
  //        {

  //     const LevelP level = grid->getLevel(level_nr);
  const LevelP level = grid->getLevel(levelIndex);

  for (Level::const_patch_iterator iter = level->patchesBegin();
       (iter != level->patchesEnd()) && (patch == nullptr);
       iter++) {
    if (*iter == local) {
      continue;
    }
    ParticleVariable<long64> var;
    query(var, "p.particleID", matlIndex, *iter, index);
    ParticleSubset* subset = var.getParticleSubset();
    for (auto p_idx : *subset) {
      if (var[p_idx] == particleID) {
        patch = *iter;
        idx   = p_idx;
        //      cerr<<"var["<<p_idx<<"] = "<<var[*p_iter]<<endl;
        break;
      }
    }

    if (patch != nullptr) {
      break;
    }
  }
  //  }
}

#if HAVE_PIDX
void
DataArchive::createPIDXCommunicator(const GridP& grid, LoadBalancer* lb)
{
  int rank = Uintah::Parallel::getMPIRank();
  // cout << rank << ": entering createPIDXCommunicator() for this: " << this <<
  // "\n";

  lb->possiblyDynamicallyReallocate(grid, LoadBalancer::RESTART_LB);

  // Resize the comms back to 0...
  d_pidx_comms.clear();

  // Create new MPI Comms
  d_pidx_comms.resize(grid->numLevels());

  // cout << rank << ": number of levels: " << grid->numLevels() << "\n";

  for (int lev = 0; lev < grid->numLevels(); lev++) {

    const LevelP& level     = grid->getLevel(lev);
    const PatchSet* patches = lb->getOutputPerProcessorPatchSet(level);

    int color                      = 0;
    const PatchSubset* patchsubset = patches->getSubset(rank);

    color = !patchsubset->empty();
    // cout << rank << ": createPIDXCommunicator() color is: " << color << "\n";

    MPI_Comm comm = Parallel::getRootProcessorGroup()->getComm();
    MPI_Comm_split(comm, color, rank, &(d_pidx_comms[lev]));

    // Debug:
    if (color == 1) {
      int nsize;
      MPI_Comm_size(d_pidx_comms[lev], &nsize);
      // cout << rank << ": NewComm Size = " <<  nsize << " on level: " << lev
      // << "\n";
    }
  }
  // cout << rank << ": leaving createPIDXCommunicator() with " <<
  // d_pidx_comms.size() << " pidx comms\n";
}
#endif

void
DataArchive::restartInitialize(int index,
                               const GridP& grid,
                               DataWarehouse* dw,
                               LoadBalancer* lb,
                               double* pTime)
{
  std::vector<int> ts_indices;
  std::vector<double> times;
  queryTimesteps(ts_indices, times);

  std::vector<std::string> names;
  std::vector<int> num_matls;
  std::vector<const TypeDescription*> typeDescriptions;
  queryVariables(names, num_matls, typeDescriptions);
  queryGlobals(names, num_matls, typeDescriptions);

  std::map<std::string, VarLabel*> varMap;
  std::map<std::string, int>
    varNameToNumMatlsMap; // FIXEME: has to be a better way to do this...

  for (unsigned i = 0; i < names.size(); i++) {
    VarLabel* vl = VarLabel::find(names[i]);
    if (vl == nullptr) {
      //      proc0cout << "Warning, VarLabel for " << names[i] << " was not
      //      found... attempting to create.\n"
      //          << "However, it is possible that this may cause problems
      //          down the road...\n";
      //***** THIS ASSUMES A SINGLE GHOST CELL ***** BE CAREFUL ********
      // check if we have extracells specified. This affects Wasatch only and
      // should have no impact on other components.
      // const bool hasExtraCells = (grid->getPatchByID(0,0)->getExtraCells()
      // != Uintah::IntVector(0,0,0)); if extracells are specified, then
      // create varlabels that are consistent with Wasatch varlabels.
      vl = VarLabel::create(names[i], typeDescriptions[i], IntVector(0, 0, 0));

      // At the end of this routine, we will need to delete the VarLabels that
      // we create here in order to avoid a memory leak.
      d_created_var_labels[names[i]] = vl;
    }
    varMap[names[i]]               = vl;
    varNameToNumMatlsMap[names[i]] = num_matls[i];
  }

  TimeData& timedata = getTimeData(index);

  *pTime = times[index];

  if (lb) {
    lb->restartInitialize(this, index, timedata.d_ts_path_and_filename, grid);
  }

  // set here instead of the SimCont because we need the DW ID to be set
  // before saving particle subsets
  dw->setID(ts_indices[index]);

  if (d_file_format == UDA) {
    // make sure to load all the data so we can iterate through it
    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);
      for (int p = 0; p < level->numPatches(); p++) {
        const Patch* patch = level->getPatch(p);
        if (lb->getPatchwiseProcessorAssignment(patch) == d_processor) {
          timedata.parsePatch(patch);
        }
      }
    }
  }

  // iterate through all entries in the VarData hash table, and loading the
  // variables if that data belongs on this processor
  if (d_file_format == UDA) {
    std::vector<VarnameMatlPatch>::iterator iter;

    for (iter = timedata.d_datafileInfoIndex.begin();
         iter != timedata.d_datafileInfoIndex.end();
         ++iter) {

      int pos = std::distance(timedata.d_datafileInfoIndex.begin(), iter);
      VarnameMatlPatch& key = *iter;
      DataFileInfo& data    = timedata.d_datafileInfoValue[pos];

      // get the Patch from the Patch ID (ID of -1 = nullptr - for reduction
      // vars)
      const Patch* patch =
        key.patchid_ == -1 ? nullptr : grid->getPatchByID(key.patchid_, 0);
      int matl = key.matlIndex_;

      VarLabel* label = varMap[key.name_];
      if (label == 0) {
        throw UnknownVariable(key.name_,
                              dw->getID(),
                              patch,
                              matl,
                              "on DataArchive::scheduleRestartInitialize",
                              __FILE__,
                              __LINE__);
      }

      if (!patch || !lb ||
          lb->getPatchwiseProcessorAssignment(patch) == d_processor) {
        Variable* var = label->typeDescription()->createInstance();
        query(*var, key.name_, matl, patch, index, &data);

        ParticleVariableBase* particles;
        if ((particles = dynamic_cast<ParticleVariableBase*>(var))) {
          if (!dw->haveParticleSubset(matl, patch)) {
            dw->saveParticleSubset(particles->getParticleSubset(), matl, patch);
          } else {
            ASSERTEQ(dw->getParticleSubset(matl, patch),
                     particles->getParticleSubset());
          }
        }

        dw->put(var, label, matl, patch);
        delete var; // should have been cloned when it was put
      }
    }
  } else { // Reading PIDX UDA
#if HAVE_PIDX
    // cout << "Here\n";
    createPIDXCommunicator(grid, lb);

    // LEVEL LOOP

    for (int lev_num = 0; lev_num < grid->numLevels(); lev_num++) {

      // cout << Uintah::Parallel::getMPIRank() << ":    on level: " << lev_num
      // << "\n";

      LevelP level = grid->getLevel(lev_num);

      // VARIABLE LOOP

      for (map<std::string, VarLabel*>::iterator varMapIter = varMap.begin();
           varMapIter != varMap.end();
           ++varMapIter) {

        const VarLabel* label       = varMapIter->second;
        const TypeDescription* type = label->typeDescription();
        const string var_name       = varMapIter->first;
        int number_of_materials;

        if (type->getType() == TypeDescription::ReductionVariable ||
            type->getType() == TypeDescription::SoleVariable) {
          number_of_materials = 0;

          // read in the reduction and sole var data here

          // FIXME: should this happen here or outside of the level loop?

          Variable* red_var = label->typeDescription()->createInstance();
          // Hard-coded to 0 as only reduction (global) variables are stored in
          // d_datafileInfoValue when using PIDX.
          DataFileInfo& data = timedata.d_datafileInfoValue[0];
          bool found;

          // cout << Uintah::Parallel::getMPIRank() << ": calling query() on: "
          // << var_name << "\n";

          found = query(
            *red_var, varMapIter->first, -1, nullptr, timestep_index, &data);
          if (!found) {
            // cout << "Warning: " << *label << " (reduction var) not found...
            // skipping.\n";
            delete red_var;
            continue;
          }
          // cout << Uintah::Parallel::getMPIRank() << ": putting reduction var
          // " << var_name << " info into dw, patch info: " << lev_num << "\n";
          dw->put(red_var, label, -1, nullptr);
          delete red_var;
        } else {
          number_of_materials = varNameToNumMatlsMap[var_name];
        }

        // cout << Uintah::Parallel::getMPIRank() << ": READING in var: " <<
        // *label << " with number of matls: " << number_of_materials << "\n";
        // MATERIAL LOOP

        for (int matl = 0; matl < number_of_materials;
             matl++) { // FIXME CHANGE "1" TO CORRECT VALUE!!!!!!!!!!!!!!!!

          // cout << Uintah::Parallel::getMPIRank() << ":      looking for matl:
          // " << matl << "\n";

          map<VarnameMatlPatch, BufferAndSizeTuple*> dataBufferMap;

          PIDX_file idxFile;
          PIDX_variable varDesc;
          PIDX_access access;

          // Non-reduction and non-sole variables:
          if (type->getType() != TypeDescription::ReductionVariable &&
              type->getType() != TypeDescription::SoleVariable) {
            // cout << Uintah::Parallel::getMPIRank() << ": calling
            // setupQueryPIDX()\n";
            bool found = setupQueryPIDX(access,
                                        idxFile,
                                        varDesc,
                                        level,
                                        type,
                                        var_name,
                                        matl,
                                        timestep_index);
            // cout << Uintah::Parallel::getMPIRank() << ": done calling
            // setupQueryPIDX()\n";

            if (!found) {
              // Did not find this var/material... skipping...
              // cout << Uintah::Parallel::getMPIRank() << ":         not found,
              // skipping....\n";
              continue;
            }
          }

          // PATCH LOOP

          for (int patch_index = 0; patch_index < level->numPatches();
               patch_index++) {

            // cout << Uintah::Parallel::getMPIRank() << ":         patch_index:
            // " << patch_index << ".\n";

            const Patch* patch = level->getPatch(patch_index);

            if (lb->getPatchwiseProcessorAssignment(patch) == d_processor) {

              // Only read in data that belongs to this processor.

              // cout << Uintah::Parallel::getMPIRank() << ":         is my
              // patch.\n";

              // Non-reduction and non-sole variables:
              if (type->getType() != TypeDescription::ReductionVariable &&
                  type->getType() != TypeDescription::SoleVariable) {

                // cout << Uintah::Parallel::getMPIRank() << ":         create
                // buffer\n";
                BufferAndSizeTuple* data = new BufferAndSizeTuple();

                // cout << Uintah::Parallel::getMPIRank() << ": a) data->buffer
                // is " << (data->buffer == nullptr ? "nullptr" :
                // (void*)(data->buffer))
                //      << " and size: " << data->size << "\n";

                queryPIDX(
                  data, varDesc, type, var_name, matl, patch, timestep_index);

#if 0
                if( data->buffer == nullptr ) {
                  cout << Uintah::Parallel::getMPIRank() << ": b) DATA IS nullptr and size: " << data->size << "\n";
                }
                else {
                  cout << Uintah::Parallel::getMPIRank() << ": b) data is at " << (void*)(data->buffer) << " and size: " << data->size << "\n";
                }
#endif

                //                if( !found ) {
                //                  proc0cout << "Warning: " << *label << " not
                //                  found... skipping.\n";
                //                 // delete var;
                //                  continue;
                //                }

                VarnameMatlPatch vmp(label->getName(), matl, patch->getID());
                // cout << Uintah::Parallel::getMPIRank() << ": inserting data
                // tubple: " << data << " into dataBufferMap\n";

                dataBufferMap[vmp] = data;

                // cout << Uintah::Parallel::getMPIRank() << ": just put var in
                // dw for patch: " << patch->getID() << ", size of map is: " <<
                // dataBufferMap.size()
                //      << " and map is: " << &dataBufferMap << "\n";
              } else {
                throw InternalError(
                  "DataArchive - should never get here", __FILE__, __LINE__);
              }
            } // end if this is my patch
          }   // end for patch

          // Non-reduction and non-sole variables:
          if (type->getType() != TypeDescription::ReductionVariable &&
              type->getType() != TypeDescription::SoleVariable) {

            // Finish up reading in the PIDX data... in the above loop in the
            // queryPIDX() calls, we provided PIDX with a list of all the
            // patches that we need... Now we call PIDX_close() to tell PIDX to
            // actually read that data. Once we have the raw data from PIDX
            // (stored in the buffers pointed to in the dataBufferMap), we copy
            // it into the actual Uintah vars.

            // Close idx file and access.
            int ret = PIDX_close(idxFile);
            PIDXOutputContext::checkReturnCode(
              ret,
              "DataArchive::query() - PIDX_close failure",
              __FILE__,
              __LINE__);

            ret = PIDX_close_access(access);
            PIDXOutputContext::checkReturnCode(
              ret,
              "DataArchive::query() - PIDX_close_access failure",
              __FILE__,
              __LINE__);

            // Iterate over buffer map and put data into actual variables...
            // cout << Uintah::Parallel::getMPIRank() << ": SIZE OF MAP is now:
            // " << dataBufferMap.size()
            //      << " and map is: " << &dataBufferMap << "\n";

            for (map<VarnameMatlPatch, BufferAndSizeTuple*>::iterator iter =
                   dataBufferMap.begin();
                 iter != dataBufferMap.end();
                 ++iter) {

              BufferAndSizeTuple* data = iter->second;

              // cout << Uintah::Parallel::getMPIRank() << ": the data tuple is:
              // " << data << "\n";

              const VarnameMatlPatch& vmp = iter->first;

              // I am assuming that this "getPatchByID" call isn't very
              // expensive... need to verify...
              const Patch* patch = grid->getPatchByID(vmp.patchid_, 0);

              Variable* var = label->typeDescription()->createInstance();
              const IntVector bl(0, 0, 0);

              if (label->typeDescription()->getType() ==
                  TypeDescription::Type::ParticleVariable) {

                int matlIndex = vmp.matlIndex_;

                psetDBType::key_type key(matlIndex, patch);
                ParticleSubset* psubset       = 0;
                psetDBType::iterator psetIter = d_psetDB.find(key);

                if (psetIter != d_psetDB.end()) {
                  psubset = (*psetIter).second.get_rep();
                }

                if (psubset == 0 ||
                    (int)psubset->numParticles() != data->size) {
                  psubset = scinew ParticleSubset(data->size, matlIndex, patch);
                  //      cout << "numParticles: " << dfi->numParticles << "\n";
                  //      cout << "d_pset size: " << d_psetDB.size() << "\n";
                  //      cout << "1. key is: " << key.first << "\n";
                  //      cout << "2. key is: " << key.second << "\n";
                  d_psetDB[key] = psubset;
                }
                (static_cast<ParticleVariableBase*>(var))->allocate(psubset);

              } else { // Non particle var

                var->allocate(patch, bl);
              }

              //__________________________________
              // Now move the dataPIDX buffer into the array3 variable

              // cout << Uintah::Parallel::getMPIRank() << ": c) var: " <<
              // vmp.name_
              //      << ", data->buffer is " << (data->buffer == nullptr ?
              //      "nullptr" : (void*)(data->buffer)) << " and size: " <<
              //      data->size << "\n";

              var->readPIDX(data->buffer, data->size, timedata.d_swapBytes);

              // cout << Uintah::Parallel::getMPIRank() << ": done with
              // readPIDX\n";

              free(data->buffer);
              free(data);

              dw->put(var,
                      label,
                      matl,
                      patch); // fixme, clean up duplicate usage - (ie, fix
                              // handling of matl)
              // proc0cout << "done with put\n";
            } // end for dataBufferMap iteration

          } // end if not reduction var

          // cout << Uintah::Parallel::getMPIRank() << ": done with matl: " <<
          // matl << ".\n";

        } // end for matl

        // cout << Uintah::Parallel::getMPIRank() << ": done with var: " <<
        // var_name << ".\n";

      } // end for varMapIter

      // cout << Uintah::Parallel::getMPIRank() << ": done with level: " <<
      // lev_num << ".\n";

    } // end for lev_num

#else
    throw InternalError(
      "Asked to read a PIDX UDA, but PIDX not configured.", __FILE__, __LINE__);
#endif
  } // end PIDX reading else
} // end restartInitialize

//______________________________________________________________________
//  This method is a specialization of restartInitialize().
//  It's only used by the reduceUda component
void
DataArchive::postprocess_ReadUda(const ProcessorGroup* pg,
                                 const int timeIndex,
                                 const GridP& grid,
                                 const PatchSubset* patches,
                                 DataWarehouse* dw,
                                 LoadBalancer* lb)
{
  std::vector<int> timesteps;
  std::vector<double> times;
  std::vector<std::string> names;
  std::vector<int> num_matls;
  std::vector<const TypeDescription*> typeDescriptions;

  queryTimesteps(timesteps, times);
  queryVariables(names, num_matls, typeDescriptions);
  queryGlobals(names, num_matls, typeDescriptions);

  // create varLabels if they don't already exist
  std::map<std::string, VarLabel*> varMap;
  for (unsigned i = 0; i < names.size(); i++) {
    VarLabel* vl = VarLabel::find(names[i]);

    if (vl == nullptr) {
      vl = VarLabel::create(names[i], typeDescriptions[i], IntVector(0, 0, 0));
      d_created_var_labels[names[i]] = vl;
    }

    varMap[names[i]] = vl;
  }

  dw->setID(timeIndex);

  TimeData& timedata = getTimeData(timeIndex);

  // make sure to load all the data so we can iterate through it
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    timedata.parsePatch(patch);
  }

  //__________________________________
  // Iterate through all entries in the VarData hash table, and load the
  // variables.
  for (std::vector<VarnameMatlPatch>::iterator iter =
         timedata.d_datafileInfoIndex.begin();
       iter != timedata.d_datafileInfoIndex.end();
       ++iter) {

    int pos = std::distance(timedata.d_datafileInfoIndex.begin(), iter);
    VarnameMatlPatch& key = *iter;
    DataFileInfo& data    = timedata.d_datafileInfoValue[pos];

    // get the Patch from the Patch ID (ID of -1 = nullptr - for reduction
    // vars)
    const Patch* patch =
      key.patchid_ == -1 ? nullptr : grid->getPatchByID(key.patchid_, 0);
    int matl = key.matlIndex_;

    VarLabel* label = varMap[key.name_];

    if (label == nullptr) {
      continue;
    }

    // If this proc does not own this patch
    // then ignore the variable
    int proc = lb->getPatchwiseProcessorAssignment(patch);
    if (proc != pg->myRank()) {
      continue;
    }

    // Put the data in the DataWarehouse.
    Variable* var = label->typeDescription()->createInstance();
    query(*var, key.name_, matl, patch, timeIndex, &data);

    ParticleVariableBase* particles;
    if ((particles = dynamic_cast<ParticleVariableBase*>(var))) {
      if (!dw->haveParticleSubset(matl, patch)) {
        dw->saveParticleSubset(particles->getParticleSubset(), matl, patch);
      } else {
        ASSERTEQ(dw->getParticleSubset(matl, patch),
                 particles->getParticleSubset());
      }
    }

    dw->put(var, label, matl, patch);
    delete var;
  }

} // end reduceUda_ReadUda()

bool
DataArchive::queryRestartTimestep(int& timestep)
{
  // FYI: There was a bug that caused the <restarts> in UDAs to look like
  // this:
  //
  //    <restarts/>
  //    <restart from="advect.uda.001" timestep="22"/>
  //
  // I fixed this recently and now they look like this.
  //    <restarts>
  //       <restart from="advect.uda.001" timestep="22"/>
  //    </restarts>
  //
  // However, to handle both cases, instead of looking for the <restarts>
  // block I just search the entire file for "<restart ...".  This is fine
  // because we create and define this file, so there would never be another
  // "<restart
  // ...>" anywhere else... I hope.

  ProblemSpec* restart_ps = nullptr;

  rewind(d_index_file); // Start parsing from top of file.
  while (true) {

    string line = UintahXML::getLine(d_index_file);
    if (line == "") {
      break;
    } else if (line.compare(0, 9, "<restart ") == 0) {

      if (restart_ps) {
        delete restart_ps;
      }
      restart_ps = scinew ProblemSpec(line);
    }
  }

  if (restart_ps != nullptr) {

    // Found (the last) "<restart " node.

    std::map<std::string, std::string> attributes;
    restart_ps->getAttributes(attributes);
    string ts_num_str = attributes["timestep"];
    if (ts_num_str == "") {
      throw InternalError(
        "DataArchive::queryVariables() - 'timestep' not found",
        __FILE__,
        __LINE__);
    }
    timestep = atoi(ts_num_str.c_str());
    delete restart_ps;
    return true;
  } else {
    return false;
  }
}

// We want to cache at least a single timestep, so that we don't have
// to reread the timestep for every patch queried.  This sets the
// cache size to one, so that this condition is held.
void
DataArchive::turnOffXMLCaching()
{
  setTimestepCacheSize(1);
}

// Sets the number of timesteps to cache back to the default_cache_size
void
DataArchive::turnOnXMLCaching()
{
  setTimestepCacheSize(default_cache_size);
}

// Sets the timestep cache size to whatever you want.  This is useful
// if you want to override the default cache size determined by
// TimeHashMaps.
void
DataArchive::setTimestepCacheSize([[maybe_unused]] int new_size)
{
  d_lock.lock();
  // Now we need to reduce the size
  int current_size = (int)last_N_timesteps.size();
  DOUT(DA_dbg, "current_size = " << current_size << "\n");
  if (timestep_cache_size >= current_size) {
    // everything's fine
    d_lock.unlock();
    return;
  }

  int kill_count = current_size - timestep_cache_size;
  DOUT(DA_dbg, "kill_count = " << kill_count << "\n");
  for (int i = 0; i < kill_count; i++) {
    int cacheTimestep = last_N_timesteps.back();
    DOUT(DA_dbg, "Making room.  Purging time index " << cacheTimestep << "\n");

    last_N_timesteps.pop_back();
    d_time_data[cacheTimestep].purgeCache();
  }
  d_lock.unlock();
}

DataArchive::TimeData::TimeData(DataArchive* da,
                                const string& timestepPathAndFilename)
  : d_initialized(false)
  , d_ts_path_and_filename(timestepPathAndFilename)
  , d_parent_da(da)
{
  d_ts_directory = timestepPathAndFilename.substr(
    0, timestepPathAndFilename.find_last_of('/') + 1);
}

DataArchive::TimeData::~TimeData()
{
  purgeCache();
}

void
DataArchive::TimeData::init()
{
  d_initialized = true;
  // Pull the list of data xml files from the timestep.xml file.

  FILE* ts_file = fopen(d_ts_path_and_filename.c_str(), "r");

  if (ts_file == nullptr) {
    // FIXME: add more info to exception.
    throw ProblemSetupException(
      "Failed to open timestep file.", __FILE__, __LINE__);
  }

  // Handle endianness and number of bits
  string endianness = d_parent_da->d_global_endianness;
  int numbits       = d_parent_da->d_global_num_bits;

  DataArchive::queryEndiannessAndBits(ts_file, endianness, numbits);

  if (endianness == "" || numbits == -1) {
    // This will only happen on a very old UDA.
    throw ProblemSetupException(
      "endianness and/or numbits missing", __FILE__, __LINE__);
  }

  d_swapBytes = endianness != string(Uintah::endianness());
  d_nBytes    = numbits / 8;

  bool found = false;

  // Based on the timestep path and file name (eg: .../timestep.xml), we need
  // to cut off the associated path so that we can find the path to data.xml.
  string::size_type path_length = d_ts_path_and_filename.rfind("/") + 1;
  string path(d_ts_path_and_filename, 0, path_length);
  string data_filename = path + "data.xml";

  FILE* data_file = fopen(data_filename.c_str(), "r");

  string looked_in = data_filename;

  if (data_file != nullptr) {
    // If the data.xml file exists, look in it.
    found = ProblemSpec::findBlock("<Data>", data_file);
  } else {
    // Otherwise, look in the original timestep.xml file.
    found     = ProblemSpec::findBlock("<Data>", ts_file);
    looked_in = d_ts_path_and_filename;
  }

  if (!found) {
    throw InternalError(
      "Cannot find <Data> in timestep file", __FILE__, __LINE__);
  }

  bool done = false;
  while (!done) {

    std::string line = "";
    if (data_file != nullptr) {
      line = UintahXML::getLine(data_file);
    } else {
      line = UintahXML::getLine(ts_file);
    }

    if (line == "" || line == "</Data>") {
      done = true;
    } else if (line.compare(0, 10, "<Datafile ") == 0) {

      ProblemSpec ts_doc(line);

      std::map<std::string, std::string> attributes;
      ts_doc.getAttributes(attributes);
      string datafile = attributes["href"];
      if (datafile == "") {
        throw InternalError("DataArchive::TimeData::init() - 'href' not found",
                            __FILE__,
                            __LINE__);
      }
      std::string proc = attributes["proc"];

      // WARNING: QWERTY: READ THIS Dav...

      /* - Remove this check for restarts.  We need to accurately
         determine which patch goes on which proc, and for the moment
         we need to be able to parse all pxxxx.xml files.  --BJW
         if (proc != "") {
         int procnum = atoi(proc.c_str());
         if ((procnum % numProcessors) != processor) {
         continue;
         }
         }
      */
      if (datafile == "global.xml") {
        // Assuming that global.xml will always be small and thus using normal
        // xml lib parsing...
        parseFile(d_ts_directory + datafile, -1, -1);
      } else {

        // Get the level info out of the xml file: should be lX/pxxxxx.xml.
        unsigned level = 0;
        std::string::size_type start =
          datafile.find_first_of("l", 0, datafile.length() - 3);
        std::string::size_type end = datafile.find_first_of("/");
        if (start != std::string::npos && end != std::string::npos &&
            end > start && end - start <= 2) {
          level = atoi(datafile.substr(start + 1, end - start).c_str());
        }

        if (level >= d_xmlFilenames.size()) {
          d_xmlFilenames.resize(level + 1);
          d_xmlParsed.resize(level + 1);
        }

        std::string filename = d_ts_directory + datafile;
        d_xmlFilenames[level].push_back(filename);
        d_xmlParsed[level].push_back(false);
      }
    } else {
      throw InternalError(
        "DataArchive::TimeData::init() - bad line in <Data> block...",
        __FILE__,
        __LINE__);
    }
  } // end while()

  fclose(ts_file);

  if (data_file) {
    fclose(data_file);
  }
}

void
DataArchive::TimeData::purgeCache()
{
  d_grid = 0;

  d_datafileInfoIndex.clear();
  d_datafileInfoValue.clear();

  d_patchInfo.clear();
  d_varInfo.clear();
  d_xmlFilenames.clear();
  d_xmlParsed.clear();
  d_initialized = false;
}

// This is the function that parses the p*****.xml file for a single
// processor.
void
DataArchive::TimeData::parseFile(const string& filename,
                                 int levelNum,
                                 int basePatch)
{
  // parse the file
  ProblemSpecP top = ProblemSpecReader().readInputFile(filename);

  // materials are the same for all patches on a level - don't parse them for
  // more than one file
  bool addMaterials = levelNum >= 0 && d_matlInfo[levelNum].size() == 0;

  for (ProblemSpecP vnode = top->getFirstChild(); vnode != nullptr;
       vnode              = vnode->getNextSibling()) {
    if (vnode->getNodeName() == "Variable") {
      std::string varname;
      if (!vnode->get("variable", varname)) {
        throw InternalError("Cannot get variable name", __FILE__, __LINE__);
      }

      int patchid;
      if (!vnode->get("patch", patchid) && !vnode->get("region", patchid)) {
        throw InternalError("Cannot get patch id", __FILE__, __LINE__);
      }

      int index;
      if (!vnode->get("index", index)) {
        throw InternalError("Cannot get index", __FILE__, __LINE__);
      }

      if (addMaterials) {
        // set the material to existing.  index+1 to use matl -1
        if (index + 1 >= (int)d_matlInfo[levelNum].size()) {
          d_matlInfo[levelNum].resize(index + 2);
        }
        d_matlInfo[levelNum][index] = true;
      }

      std::map<std::string, std::string> attributes;
      vnode->getAttributes(attributes);

      string type = attributes["type"];
      if (type == "") {
        throw InternalError("DataArchive::query:Variable doesn't have a type",
                            __FILE__,
                            __LINE__);
      }
      long start;
      if (!vnode->get("start", start)) {
        throw InternalError(
          "DataArchive::query:Cannot get start", __FILE__, __LINE__);
      }
      long end;
      if (!vnode->get("end", end)) {
        throw InternalError(
          "DataArchive::query:Cannot get end", __FILE__, __LINE__);
      }
      string filename;
      if (!vnode->get("filename", filename)) {
        throw InternalError(
          "DataArchive::query:Cannot get filename", __FILE__, __LINE__);
      }

      // not required
      string compressionMode = "";
      IntVector boundary(0, 0, 0);
      int numParticles = -1;

      vnode->get("compression", compressionMode);
      vnode->get("boundaryLayer", boundary);
      vnode->get("numParticles", numParticles);

      if (d_varInfo.find(varname) == d_varInfo.end()) {
        VarData& varinfo      = d_varInfo[varname];
        varinfo.type          = type;
        varinfo.compression   = compressionMode;
        varinfo.boundaryLayer = boundary;
      } else if (compressionMode != "") {
        // For particles variables of size 0, the uda doesn't say it
        // has a compressionMode...  (FYI, why is this?  Because it is
        // ambiguous... if there is no data, is it compressed?)
        //
        // To the best of my understanding, we only look at the variables
        // stats the first time we encounter it... even if there are multiple
        // materials. So we run into a problem is the variable has 0 data the
        // first time it is looked at... The problem there is that it doesn't
        // mark it as being compressed, and therefore the next time we see
        // that variable (eg, in another material) we (used to) assume it was
        // not compressed... the following lines compenstate for this problem:
        VarData& varinfo    = d_varInfo[varname];
        varinfo.compression = compressionMode;
      }

      if (levelNum == -1) { // global file (reduction vars)
        d_globaldata = filename;
      } else {
        ASSERTRANGE(patchid - basePatch, 0, (int)d_patchInfo[levelNum].size());

        PatchData& patchinfo = d_patchInfo[levelNum][patchid - basePatch];
        if (!patchinfo.parsed) {
          patchinfo.parsed       = true;
          patchinfo.datafilename = filename;
        }
      }
      VarnameMatlPatch vmp(varname, index, patchid);
      DataFileInfo dummy;

      if (std::find(d_datafileInfoIndex.begin(),
                    d_datafileInfoIndex.end(),
                    vmp) != d_datafileInfoIndex.end()) {
        // std::cerr <<  "Duplicate variable name: " << name << std::endl;
      } else {
        DataFileInfo dfi(start, end, numParticles);
        d_datafileInfoIndex.push_back(vmp);
        d_datafileInfoValue.push_back(dfi);
      }
    } else if (vnode->getNodeType() != ProblemSpec::TEXT_NODE) {
      std::cerr << "WARNING: Unknown element in Variables section: "
                << vnode->getNodeName() << '\n';
    }
  }
  // top->releaseDocument();
}

void
DataArchive::TimeData::parsePatch(const Patch* patch)
{
  ASSERT(d_grid != nullptr);
  if (!patch) {
    proc0cout << "parsePatch called with null patch....\n";
    return;
  }

  const Patch* real_patch = patch->getRealPatch();

  // make sure the data for this patch has been processed.
  // Return straightaway if we have parsed this patch
  int levelIndex       = real_patch->getLevel()->getIndex();
  int levelBasePatchID = real_patch->getLevel()->getPatch(0)->getID();
  int patchIndex       = real_patch->getLevelIndex();

  PatchData& patchinfo = d_patchInfo[levelIndex][patchIndex];
  if (patchinfo.parsed) {
    return;
  }

  if (d_parent_da->d_file_format == PIDX) {
    d_xmlParsed[levelIndex][patchIndex] = true;
    patchinfo.parsed                    = true;
    return;
  }

  // If this is a newer uda, the patch info in the grid will store the
  // processor where the data is
  if (patchinfo.proc != -1) {
    std::ostringstream file;
    file << d_ts_directory << "l" << (int)real_patch->getLevel()->getIndex()
         << "/p" << std::setw(5) << std::setfill('0') << (int)patchinfo.proc
         << ".xml";
    parseFile(file.str(), levelIndex, levelBasePatchID);
  } else {
    // Try making a guess as to the processor.  First go is to try
    // the processor of the same index as the patch.  Many datasets
    // have only one patch per processor, so this is a reasonable
    // first attempt.  Future attemps could perhaps be smarter.
    if (!patchinfo.parsed && patchIndex < (int)d_xmlParsed[levelIndex].size() &&
        !d_xmlParsed[levelIndex][patchIndex]) {
      parseFile(
        d_xmlFilenames[levelIndex][patchIndex], levelIndex, levelBasePatchID);
      d_xmlParsed[levelIndex][patchIndex] = true;
    }

    // failed the guess - parse the entire dataset for this level
    if (!patchinfo.parsed) {
      for (unsigned proc = 0; proc < d_xmlFilenames[levelIndex].size();
           proc++) {
        parseFile(
          d_xmlFilenames[levelIndex][proc], levelIndex, levelBasePatchID);
        d_xmlParsed[levelIndex][proc] = true;
      }
    }
  }
}

// Parses the timestep xml file for <oldDelt>
//
double
DataArchive::getOldDelt(int restart_index)
{
  TimeData& timedata = getTimeData(restart_index);
  FILE* fp           = fopen(timedata.d_ts_path_and_filename.c_str(), "r");
  if (fp == nullptr) {
    throw InternalError(
      "DataArchive::setOldDelt() failed open datafile.", __FILE__, __LINE__);
  }
  // Note, old UDAs had a <delt> flag, but that was deprecated long ago in
  // favor of the <oldDelt> flag which is what we are going to look for here.

  while (true) {

    std::string line = UintahXML::getLine(fp);

    if (line == "") {
      fclose(fp);
      throw InternalError("DataArchive::setOldDelt() failed to find <oldDelt>.",
                          __FILE__,
                          __LINE__);
    } else if (line.compare(0, 9, "<oldDelt>") == 0) {
      std::vector<std::string> pieces = UintahXML::splitXMLtag(line);

      fclose(fp);
      return atof(pieces[1].c_str());
    }
  }
}

// Parses the timestep xml file and skips the <Meta>, <Grid>, and <Data>
// sections, returning everything else.  This function assumes that the
// timestep.xml file was created by us and is in the correct order - in other
// words, anything after </Data> is component related, and everything before
// it can be removed.
//
ProblemSpecP
DataArchive::getTimestepDocForComponent(int restart_index)
{
  TimeData& timedata = getTimeData(restart_index);
  FILE* fp           = fopen(timedata.d_ts_path_and_filename.c_str(), "r");

  if (fp == nullptr) {
    throw InternalError(
      "DataArchive::getTimespecDocForComponent() failed open datafile.",
      __FILE__,
      __LINE__);
  }

  string buffer = "<Uintah_timestep>";

  while (true) {

    string line = UintahXML::getLine(fp);

    buffer.append(line);

    if (line == "</Uintah_timestep>") {
      break;
    }
  }

  fclose(fp);

  ProblemSpec* result = new ProblemSpec(buffer);

  return result;
}

ConsecutiveRangeSet
DataArchive::queryMaterials(const string& varname,
                            const Patch* patch,
                            int index)
{
  ConsecutiveRangeSet matls;

  if (d_file_format == PIDX) {
    matls = var_materials[varname];
    return matls;
  }

  Timers::Simple timer;
  timer.start();

  d_lock.lock();

  TimeData& timedata = getTimeData(index);
  timedata.parsePatch(patch);

  for (unsigned i = 0;
       i < timedata.d_matlInfo[patch->getLevel()->getIndex()].size();
       i++) {
    // i-1, since the matlInfo is adjusted to allow -1 as entries
    VarnameMatlPatch vmp(varname, i - 1, patch->getRealPatch()->getID());
    DataFileInfo dummy;

    if (std::find(timedata.d_datafileInfoIndex.begin(),
                  timedata.d_datafileInfoIndex.end(),
                  vmp) != timedata.d_datafileInfoIndex.end()) {
      matls.addInOrder(i - 1);
    }
  }

  d_lock.unlock();
  DOUT(DA_dbg,
       "DataArchive::queryMaterials completed in " << timer().seconds()
                                                   << " seconds\n");

  return matls;
}

int
DataArchive::queryNumMaterials(const Patch* patch, int index)
{
  Timers::Simple timer;
  timer.start();

  d_lock.lock();

  TimeData& timedata = getTimeData(index);

  timedata.parsePatch(patch);

  int numMatls = -1;

  for (unsigned i = 0;
       i < timedata.d_matlInfo[patch->getLevel()->getIndex()].size();
       i++) {
    if (timedata.d_matlInfo[patch->getLevel()->getIndex()][i]) {
      numMatls++;
    }
  }

  DOUT(DA_dbg,
       "DataArchive::queryNumMaterials completed in " << timer().seconds()
                                                      << " seconds\n");

  d_lock.unlock();
  return numMatls;
}

//______________________________________________________________________
//    Does this variable exist on this patch at this timestep
bool
DataArchive::exists(const string& varname,
                    const Patch* patch,
                    const int timeStep)
{
  d_lock.lock();

  TimeData& timedata = getTimeData(timeStep);
  timedata.parsePatch(patch);

  int levelIndex = patch->getLevel()->getIndex();

  for (unsigned i = 0; i < timedata.d_matlInfo[levelIndex].size(); i++) {
    // i-1, since the matlInfo is adjusted to allow -1 as entries
    VarnameMatlPatch vmp(varname, i - 1, patch->getRealPatch()->getID());
    DataFileInfo dummy;

    if (std::find(timedata.d_datafileInfoIndex.begin(),
                  timedata.d_datafileInfoIndex.end(),
                  vmp) != timedata.d_datafileInfoIndex.end()) {
      d_lock.unlock();
      return true;
    }
  }
  d_lock.unlock();

  return false;
}

} // namespace Uintah