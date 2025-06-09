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

#ifndef __CORE_PARALLEL_PROCESSORGROUP_H__
#define __CORE_PARALLEL_PROCESSORGROUP_H__

#include <Core/Parallel/UintahMPI.h>
#include <vector>
#include <string>

namespace Uintah {

class Parallel;

class ProcessorGroup
{
public:
  ~ProcessorGroup();

  // Eliminate copy, assignment and move
  ProcessorGroup(const ProcessorGroup&) = delete;
  ProcessorGroup& operator=(const ProcessorGroup&) = delete;
  ProcessorGroup(ProcessorGroup&&) = delete;
  ProcessorGroup& operator=(ProcessorGroup&&) = delete;

  std::string myNodeName() const;

  // Returns the total number of MPI nodes in this MPI session.
  int nNodes() const { return d_all_proc_names.size(); }
  int myNode() const { return d_all_proc_indices[d_rank]; }

  // Returns the total number of MPI rank in the node MPI session.
  int myNode_nRanks() const { return d_node_nRanks; }
  int myNode_myRank() const { return d_node_rank; }

  // Returns the total number of MPI ranks in this MPI session.
  int nRanks() const { return d_nRanks; }
  int myRank() const { return d_rank; }

  MPI_Comm getComm() const { return d_comm; }
  MPI_Comm getNodeComm() const { return d_node_comm; }

  MPI_Comm getGlobalComm(int comm_index) const
  {
    if (d_threads <= 1 || comm_index == -1) {
      return d_comm;
    } else {
      return d_global_comms[comm_index];
    }
  }

  void setGlobalComm(int num_comms) const;

  // Utilities for getting node based information.
  std::optional<int> getNodeIndexFromRank(unsigned int rank) const noexcept;
  std::optional<std::string> getNodeNameFromRank(unsigned int rank) const noexcept;
  std::optional<std::string> getNodeName(unsigned int node) const noexcept;

private:
  friend class Parallel;
  ProcessorGroup(const ProcessorGroup* parent, MPI_Comm comm, int rank,
                 int size, int threads);

private:
  const ProcessorGroup* d_parent_group{ nullptr };

  MPI_Comm d_comm{ 0 };
  MPI_Comm d_node_comm{ 0 };
  mutable std::vector<MPI_Comm> d_global_comms;

  int d_node_rank{ -1 };  // MPI rank of this process relative to the node.
  int d_node_nRanks{ 0 }; // Total number of MPI Ranks relative to the node.

  int d_rank{ -1 };  // MPI rank of this process.
  int d_nRanks{ 0 }; // Total number of MPI ranks

  int d_threads{ 0 };

  // For storing all the processor names so to provide a mapping from
  // the name to an index from any rank.
  typedef char procName_t[MPI_MAX_PROCESSOR_NAME + 1];

  // For each rank store an index to its node.
  std::vector<unsigned int> d_all_proc_indices;

  // For each node store its processor name.
  std::vector<std::string> d_all_proc_names;
};

} // End namespace Uintah

#endif //__CORE_PARALLEL_PROCESSORGROUP_H__
