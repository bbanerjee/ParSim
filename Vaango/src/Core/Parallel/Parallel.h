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

#ifndef __CORE_PARALLEL_PARALLEL_H__
#define __CORE_PARALLEL_PARALLEL_H__

#include <iostream>
#include <thread>

// Macro used to reduce output on large parallel runs
//
//   Note, make sure that Uintah::MPI::Init (or Uintah::MPI::Init_thread) is
//   called before using isProc0_macro.
//
#define isProc0_macro                                                          \
  (Uintah::Parallel::getMPIRank() == 0) &&                                     \
    (Uintah::Parallel::getMainThreadID() == std::this_thread::get_id())

#define proc0cout                                                              \
  if (isProc0_macro)                                                           \
  std::cout

#define proc0cerr                                                              \
  if (isProc0_macro)                                                           \
  std::cerr

#define MAX_THREADS 64
#define MAX_HALO_DEPTH 5

namespace Uintah {

class ProcessorGroup;

class Parallel
{
public:
  enum Circumstances
  {
    NormalShutdown,
    Abort
  };

  // Initializes MPI if necessary.
  static void initializeManager(int& argc, char**& arg);

  // Check to see whether initializeManager has been called
  static bool isInitialized();

  // Shut down MPI gracefully
  static void finalizeManager(Circumstances cirumstances = NormalShutdown);

  // Return root context processorgroup
  static ProcessorGroup* getRootProcessorGroup();

  // Return the MPI Rank of this process.  If this is not running
  // under MPI, then 0 is returned.  Rank value is set after call to
  // initializeManager();
  static int getMPIRank();

  // Return the size of MPI_Comm
  static int getMPISize();

  // Return true if this process is to use GPUs, false otherwise
  static bool usingDevice();

  // Set whether or not to use available GPUs
  static void setUsingDevice(bool state);

  // Return the number of threads that a processing element is
  // allowed to use to compute its tasks.
  static int getNumThreads();

  // Return the number of thread partitions that a processing element is
  // allowed to use to compute its tasks.
  static int getNumPartitions();

  // Return the number of threads per partition.
  static int getThreadsPerPartition();

  // Return the ID of the main thread, via std::this_thread::get_id()
  static std::thread::id getMainThreadID();

  // Set the number of task runner threads to the value specified
  static void setNumThreads(int num);

  // Set the number of task runner OMP thread partitions to the value specified
  static void setNumPartitions(int num);

  // Set the number of threads per OMP partition
  static void setThreadsPerPartition(int num);

  // Pass the specified exit code to std::exit()
  static void exitAll(int code);

public:
  Parallel(const Parallel&) = delete;
  Parallel& operator=(const Parallel&) = delete;
  Parallel(Parallel&&) = delete;
  Parallel& operator=(Parallel&&) = delete;

private:
  Parallel();
  ~Parallel();

  static bool s_initialized;
  static bool s_using_device;
  static int s_num_threads;
  static int s_num_partitions;
  static int s_threads_per_partition;
  static int s_world_rank;
  static int s_world_size;
  static std::thread::id s_main_thread_id;
  static ProcessorGroup* s_root_context;
};
} // End namespace Uintah

#endif // __CORE_PARALLEL_PARALLEL_H__
