/*
 * The MIT License
 *
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

#include <StandAlone/Utils/vaango_options.h>
#include <StandAlone/Utils/vaango_utils.h>

#include <Core/Parallel/Parallel.h>

#include <sci_defs/compile_defs.h>

#include <filesystem>
#include <iostream>
#include <string_view>
#include <unistd.h>
#include <vector>

namespace {

static bool s_emit_graphs{ false };
static bool s_local_filesystem{ false };
static bool s_only_validate_ups{ false };
static bool s_postprocess_uda{ false };
static bool s_restart{ false };
static bool s_restart_from_scratch{ true };
static bool s_restart_remove_old_dir{ false };
static bool s_show_config_cmd{ false };
static bool s_show_git_diff{ false };
static bool s_show_git_status{ false };
static bool s_show_version{ false };
static bool s_validate_ups{ true };
static bool s_use_gpu{ false };

static int s_num_partitions{ 0 };
static int s_num_threads{ 0 };
static int s_restart_checkpoint_index{ -1 };
static int s_threads_per_partition{ 0 };
static int s_uda_suffix{ -1 };

static std::string s_uda_dir{ "" };      // for restart
static std::string s_uda_filename{ "" }; // name of the UDA directory
static std::string s_solver_name{ "" };  // empty string defaults to CGSolver

static Uintah::IntVector s_layout{ 1, 1, 1 };

} // end anonymous namespace

namespace Vaango {
namespace Utils {
namespace Options {

void
parse(int argc, char** argv)
{
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if ((arg == "-help") || (arg == "-h")) {
      Vaango::Utils::Options::usage("", "", argv[0]);
    } else if ((arg == "-debug") || (arg == "-d")) {
      Vaango::Utils::print_active_debug_streams();
    } else if (arg == "-nthreads") {
      if (++i == argc) {
        Vaango::Utils::Options::usage(
          "You must provide a number of threads for -nthreads", arg, argv[0]);
      }
      s_num_threads = atoi(argv[i]);
      if (s_num_threads < 1) {
        Vaango::Utils::Options::usage(
          "Number of threads is too small", arg, argv[0]);
      } else if (s_num_threads > MAX_THREADS) {
        Vaango::Utils::Options::usage(
          "Number of threads is out of range. Specify fewer threads, "
          "or increase MAX_THREADS (.../src/Core/Parallel/Parallel.h) and "
          "recompile.",
          arg,
          argv[0]);
      }
    } else if (arg == "-npartitions") {
      if (++i == argc) {
        Vaango::Utils::Options::usage(
          "You must provide a number of thread partitions for -npartitions",
          arg,
          argv[0]);
      }
      s_num_partitions = atoi(argv[i]);
      if (s_num_partitions < 1) {
        Vaango::Utils::Options::usage(
          "Number of thread partitions is too small", arg, argv[0]);
      } else if (s_num_partitions > MAX_THREADS) {
        Vaango::Utils::Options::usage(
          "Number of thread partitions is out of range. Specify fewer thread "
          "partitions, "
          "or increase MAX_THREADS (.../src/Core/Parallel/Parallel.h) and "
          "recompile.",
          arg,
          argv[0]);
      }
    } else if (arg == "-nthreadsperpartition") {
      if (++i == argc) {
        Vaango::Utils::Options::usage(
          "You must provide a number of threads per partition for "
          "-nthreadsperpartition",
          arg,
          argv[0]);
      }
      s_threads_per_partition = atoi(argv[i]);
      if (s_threads_per_partition < 1) {
        Vaango::Utils::Options::usage(
          "Number of threads per partition is too small", arg, argv[0]);
      }
#ifdef _OPENMP
      if (s_threads_per_partition > omp_get_max_threads()) {
        Vaango::Utils::Options::usage(
          "Number of threads per partition must be <= omp_get_max_threads()",
          arg,
          argv[0]);
      }
#endif
    } else if (arg == "-solver") {
      if (++i == argc) {
        Vaango::Utils::Options::usage(
          "You must provide a solver name for -solver", arg, argv[0]);
      }
      s_solver_name = argv[i];
    } else if (arg == "-emit_taskgraphs") {
      s_emit_graphs = true;
    } else if (arg == "-local_filesystem") {
      s_local_filesystem = true;
    } else if (arg == "-restart") {
      s_restart = true;
    } else if (arg == "-uda_suffix") {
      if (i < argc - 1) {
        s_uda_suffix = atoi(argv[++i]);
      } else {
        Vaango::Utils::Options::usage(
          "You must provide a suffix number for -uda_suffix", arg, argv[0]);
      }
    } else if (arg == "-nocopy") { // default anyway, but that's fine
      s_restart_from_scratch = true;
    } else if (arg == "-copy") {
      s_restart_from_scratch   = false;
      s_restart_remove_old_dir = false;
    } else if (arg == "-move") {
      s_restart_from_scratch   = false;
      s_restart_remove_old_dir = true;
    } else if (arg == "-gpucheck") {
      Vaango::Utils::check_gpus();
    }
#ifdef HAVE_CUDA
    else if (arg == "-gpu") {
      s_use_gpu = true;
    }
#endif
    else if (arg == "-t") {
      if (i < argc - 1) {
        s_restart_checkpoint_index = atoi(argv[++i]);
      }
    } else if (arg == "-layout") {
      if (++i == argc) {
        Vaango::Utils::Options::usage(
          "You must provide a vector arg for -layout", arg, argv[0]);
      }
      int ii, jj, kk;
      if (sscanf(argv[i], "%dx%dx%d", &ii, &jj, &kk) != 3) {
        Vaango::Utils::Options::usage(
          "Error parsing -layout", argv[i], argv[0]);
      }
      s_layout = Uintah::IntVector(ii, jj, kk);
    } else if (arg == "-config_cmd") {
      s_show_config_cmd = true;
    } else if (arg == "-git_diff") {
      s_show_git_diff = true;
    } else if (arg == "-git_status") {
      s_show_git_status = true;
    } else if (arg == "-validate") {
      s_only_validate_ups = true;
    } else if (arg == "-do_not_validate") {
      s_validate_ups = false;
    } else if (arg == "-postprocess_uda") {
      s_postprocess_uda = true;
    } else if (arg == "-version" || arg == "-v") {
      s_show_config_cmd = true;
      s_show_git_diff   = true;
      s_show_git_status = true;
      s_show_version    = true;
    } else {
      // A filename was already provided, thus this is an error.
      if (s_uda_filename != "") {
        Vaango::Utils::Options::usage("", arg, argv[0]);
      } else if (argv[i][0] == '-') {
        // Don't allow 'filename' to begin with '-'.
        Vaango::Utils::Options::usage(
          "Error!  It appears that the filename you specified begins with "
          "a '-'.\n"
          "        This is not allowed.  Most likely there is problem with "
          "your\n"
          "        command line.",
          argv[i],
          argv[0]);
      } else {
        s_uda_filename = argv[i];
      }
    }
  }

  if (s_uda_filename == "" && s_show_version == false) {
    Vaango::Utils::Options::usage("No input file specified", "", argv[0]);
  }

  //  bulletproofing
  if (s_restart || s_postprocess_uda) {
    s_uda_dir      = s_uda_filename;
    s_uda_filename = s_uda_filename + "/input.xml";

    // If restarting (etc), make sure that the uda specified is not a
    // symbolic link to an Uda.  This is because the sym link can
    // (will) be updated to point to a new uda, thus creating an
    // inconsistency.  Therefore it is just better not to use the sym
    // link in the first place.
    if (std::filesystem::is_symlink(
          std::filesystem::status(s_uda_dir.c_str()))) {
      std::cout << "\n";
      std::cout
        << "ERROR: " + s_uda_dir +
             " is a symbolic link.  Please use the full name of the UDA.\n";
      std::cout << "\n";
      Vaango::Utils::stop_mpi_and_exit(1);
    }
  }

  if (!s_validate_ups) {
    // Print out warning message here (after Parallel::initializeManager()),
    // so that proc0cout works correctly.
    proc0cout << "\n";
    proc0cout << "WARNING: You have turned OFF .ups file validation... this "
                 "may cause many unforeseen problems\n";
    proc0cout << "         with your simulation run.  It is strongly "
                 "suggested that you leave validation on!\n";
    proc0cout << "\n";
  }

  // Output header
  // helpful for cleaning out old stale udas
  time_t t = time(nullptr);
  std::string time_string(ctime(&t));
  char name[256];
  gethostname(name, 256);

  std::cout << "Date:    " << time_string; // has its own newline
  std::cout << "Machine: " << name << "\n";
  std::cout << "Assertion level: " << SCI_ASSERTION_LEVEL << "\n";
  std::cout << "CFLAGS: " << CFLAGS << "\n";
  std::cout << "CXXFLAGS: " << CXXFLAGS << "\n";

  Vaango::Utils::display_git_info(s_show_git_diff, s_show_git_status);
  Vaango::Utils::display_config_info(s_show_config_cmd);

  if (s_show_version) {
    stop_mpi_and_exit(2);
  }
}

void
usage(const std::string& message,
      const std::string& badarg,
      const std::string& progname)
{
  Vaango::Utils::start_mpi();

  if (Uintah::Parallel::getMPIRank() == 0) {
    std::cerr << "\n";
    if (badarg != "") {
      std::cerr << "Error parsing argument: " << badarg << '\n';
    }
    std::cerr << "\n";
    std::cerr << message << "\n";
    std::cerr << "\n";
    std::cerr << "Usage: " << progname << " [options] <input_file_name>\n\n";
    std::cerr << "Valid options are:\n";
    std::cerr << "-h[elp]              :"
              << " This usage information.\n";
    std::cerr << "-nthreads <#>        :"
              << " Number of threads per MPI process."
              << " Requires the multi-threaded Unified scheduler\n";
    std::cerr << "-layout NxMxO        :"
              << " Eg: 2x1x1.  MxNxO must equal number";
    std::cerr << " of boxes you are using.\n";
    std::cerr << "-emit_taskgraphs     :"
              << " Output taskgraph information\n";
    std::cerr << "-restart             :"
              << " Give the checkpointed uda directory as the input file\n";
    std::cerr << "-postprocess_uda     :"
              << " Passes variables in an uda through"
              << " postprocessing tasks, computing new variables and"
              << " creating a new uda.\n";
    std::cerr << "-uda_suffix <number> :"
              << " Make a new uda dir with <number> as"
              << " the default suffix\n";
    std::cerr << "-t <timestep>        :"
              << " Restart timestep (last checkpoint is"
              << " default. You can use -t 0 for the first checkpoint)\n";
    std::cerr << "-copy                :"
              << " Copy from old uda when restarting\n";
    std::cerr << "-move                :"
              << " Move from old uda when restarting\n";
    std::cerr << "-nocopy              :"
              << " Default: Don't copy or move old uda"
              << " timestep when restarting\n";
    std::cerr << "-d[ebug]             :"
              << " Lists the available debug streams\n";
    std::cerr << "-validate            :"
              << " Verifies the .ups file is valid and quits!\n";
    std::cerr << "-do_not_validate     :"
              << " Skips .ups file validation! Please avoid this flag"
              << " if possible.\n";
    std::cerr << "-cmake_cmd           :"
              << " Display cmake command used to build Vaango\n";
    std::cerr << "-local_filesystem    :"
              << " If using MPI, use this flag if each node has a local disk\n";
#ifdef HAVE_CUDA
    std::cerr << "-gpu                 : "
              << " Use available GPU devices, requires multi-threaded "
              << " Unified scheduler \n";
#endif
    std::cerr << "-gpucheck            : "
              << " Returns 1 if Vaango was compiled with "
              << " CUDA and there is a GPU available. \n";
    std::cerr << "                     : "
              << " Returns 2 if Vaango was not compiled "
              << " with CUDA or there are no GPUs available. \n";
    std::cerr << "-version             :"
              << " Display Vaango and git version\n";
    std::cerr << "-git_diff            :"
              << " Checks for new changes to the code\n";
    std::cerr << "-git_status          :"
              << " Checks for logs of new changes to the code\n";
    std::cerr << "\n\n";
  }
  Vaango::Utils::stop_mpi_and_exit(2);
}

bool
emit_graphs()
{
  return s_emit_graphs;
}

bool
local_filesystem()
{
  return s_local_filesystem;
}

bool
only_validate_ups()
{
  return s_only_validate_ups;
}

bool
postprocess_uda()
{
  return s_postprocess_uda;
}

bool
restart()
{
  return s_restart;
}

bool
restart_from_scratch()
{
  return s_restart_from_scratch;
}

bool
restart_remove_old_dir()
{
  return s_restart_remove_old_dir;
}

bool
show_config_cmd()
{
  return s_show_config_cmd;
}

bool
show_git_diff()
{
  return s_show_git_diff;
}

bool
show_git_status()
{
  return s_show_git_status;
}

bool
show_version()
{
  return s_show_version;
}

bool
validate_ups()
{
  return s_validate_ups;
}

bool
use_gpu()
{
  return s_use_gpu;
}

int
num_partitions()
{
  return s_num_partitions;
}

int
num_threads()
{
  return s_num_threads;
}

int
restart_checkpoint_index()
{
  return s_restart_checkpoint_index;
}

int
threads_per_partition()
{
  return s_threads_per_partition;
}

int
uda_suffix()
{
  return s_uda_suffix;
}

const std::string&
uda_dir()
{
  return s_uda_dir;
}

const std::string&
uda_filename()
{
  return s_uda_filename;
}

const std::string&
solver_name()
{
  return s_solver_name;
}

const Uintah::IntVector&
grid_layout()
{
  return s_layout;
}

} // namespace Options
} // namespace Utils
} // namespace Vaango
