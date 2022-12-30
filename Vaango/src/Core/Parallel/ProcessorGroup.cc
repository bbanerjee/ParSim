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

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Parallel/Parallel.h>

#include <iostream>
#include <map>

namespace Uintah {

ProcessorGroup::ProcessorGroup(const ProcessorGroup* parent, MPI_Comm comm,
                               int rank, int size, int threads)
  : d_parent_group(parent)
  , d_comm(comm)
  , d_rank(rank)
  , d_nRanks(size)
  , d_threads(threads)
{
  procName_t procName;
  int procNameLength;

  MPI::Get_processor_name(procName, &procNameLength);
  procName[procNameLength] = '\0';

  procName_t* allProcNames = new procName_t[d_nRanks];
  d_all_proc_indices.resize(d_nRanks);

  // Gather all of the processor names in rank order.
  MPI::Allgather(procName, MPI_MAX_PROCESSOR_NAME + 1, MPI_CHAR, allProcNames,
                 MPI_MAX_PROCESSOR_NAME + 1, MPI_CHAR, d_comm);

  std::map<std::string, unsigned int> name_index;

  // Loop through each rank and map its processor name to a node index.
  for (int i = 0; i < d_nRanks; ++i) {
    std::string proc_name = std::string(allProcNames[i]);

    // If the name is not found add it to the map.
    if (name_index.find(proc_name) == name_index.end()) {
      int nNodes = name_index.size();

      // Add the name to map so to track unique processor names.
      name_index[proc_name] = nNodes;

      // Store the processor name for later look up by node index.
      d_all_proc_names.push_back(proc_name);
    }

    // For each rank save its node index. This index can be
    // used to get the node name.
    d_all_proc_indices[i] = name_index[proc_name];
  }

  delete[] allProcNames;

  // More than one node so create a node based communicator.
  if (name_index.size() > 1) {
    // The node index becomes the "color" while using the world rank
    // which is unique for the ordering.
    MPI::Comm_split(d_comm, d_all_proc_indices[d_rank], d_rank, &d_node_comm);

    // Get the number of ranks for this node and the rank on this node.
    MPI::Comm_size(d_node_comm, &d_node_nRanks);
    MPI::Comm_rank(d_node_comm, &d_node_rank);
  } else {
    // Only one node so no need for a separtate communicator.
    d_node_nRanks = d_nRanks;
    d_node_rank = d_rank;
    d_node_comm = d_comm;
  }
}

ProcessorGroup::~ProcessorGroup() {}

void
ProcessorGroup::setGlobalComm(int num_comms) const
{
  if (d_threads <= 1) {
    return;
  }

  int curr_size = d_global_comms.size();
  if (num_comms <= curr_size) {
    return;
  }

  d_global_comms.resize(num_comms);
  for (int i = curr_size; i < num_comms; i++) {
    if (MPI_Comm_dup(d_comm, &d_global_comms[i]) != MPI_SUCCESS) {
      std::cerr << "Rank: " << d_rank << " - MPI Error in MPI_Comm_dup\n";
      Parallel::exitAll(1);
    }
  }
}

// For this rank get the node name.
std::string
ProcessorGroup::myNodeName() const
{
  return d_all_proc_names[d_all_proc_indices[d_rank]];
}

// For any rank get the node index.
int
ProcessorGroup::getNodeIndexFromRank(unsigned int rank) const
{
  // Make sure the rank is valid.
  if (0 <= rank && rank < d_all_proc_indices.size())
    return d_all_proc_indices[rank];
  else
    return -1;
}

// For any rank get the node name.
std::string
ProcessorGroup::getNodeNameFromRank(unsigned int rank) const
{
  // Make sure the rank is valid.
  if (0 <= rank && rank < d_all_proc_indices.size())
    return d_all_proc_names[d_all_proc_indices[rank]];
  else
    return std::string("");
}

// For any node get the node name.
std::string
ProcessorGroup::getNodeName(unsigned int node) const
{
  // Make sure the node is valid.
  if (0 <= node && node < d_all_proc_names.size())
    return d_all_proc_names[node];
  else
    return std::string("");
}

} // end namespace Uintah