/*
 * The MIT License
 *
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

#include <StandAlone/Utils/vaango_options.h>
#include <StandAlone/Utils/vaango_utils.h>

#include <CCA/Components/Schedulers/KokkosScheduler.h>
#include <Core/Parallel/Parallel.h>

#include <sci_defs/compile_defs.h>
#include <sci_defs/kokkos_defs.h>
#include <sci_defs/cuda_defs.h>

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
static bool s_show_git_diff{ false };
static bool s_show_git_status{ false };
static bool s_show_version{ false };
static bool s_validate_ups{ true };
static bool s_use_gpu{ false };

static int s_num_threads{ 0 };
static int s_restart_checkpoint_index{ -1 };
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
    } else if (arg == "-cpu") {
#if defined(KOKKOS_USING_GPU)
      Uintah::Parallel::setUsingCPU(true);
#endif
    } else if (arg == "-gpucheck") {
#if defined(KOKKOS_USING_GPU)
      if (KokkosScheduler::verifyAnyGpuActive()) {
        std::cout << "At least one GPU detected!" << std::endl;
      } else {
        std::cout << "No GPU detected!" << std::endl;
      }
      Uintah::Parallel::exitAll(1);
#else
      std::cout << "Not compiled for GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-gpu") {
#if defined(KOKKOS_USING_GPU)
      s_use_gpu = true;
      Uintah::Parallel::setUsingDevice(true);
#else
      std::cout << "Not compiled for GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_instances_per_task") {
#if defined(KOKKOS_USING_GPU)
      int kokkos_instances_per_task = 0;
      if (++i == argc) {
        usage("You must provide a number of Kokkos instances per task for "
              "-kokkos_instances_per_task.",
              arg,
              argv[0]);
      }
      kokkos_instances_per_task = atoi(argv[i]);
      if (kokkos_instances_per_task < 1) {
        usage(
          "Number of Kokkos Instances per task is too small.", arg, argv[0]);
        Uintah::Parallel::exitAll(1);
      }
      Uintah::Parallel::setKokkosInstancesPerTask(kokkos_instances_per_task);
#else
      std::cout << "Not compiled for GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_policy") {
#if defined(HAVE_KOKKOS)
      Uintah::Parallel::Kokkos_Policy kokkos_policy = Uintah::Parallel::Kokkos_Team_Policy;
      if (++i == argc) {
        usage("You must provide the policy -kokkos_policy.", arg, argv[0]);
      }

      if (strcmp(argv[i], "team") == 0) {
        kokkos_policy = Uintah::Parallel::Kokkos_Team_Policy;
      } else if (strcmp(argv[i], "range") == 0) {
        kokkos_policy = Uintah::Parallel::Kokkos_Range_Policy;
      } else if (strcmp(argv[i], "mdrange") == 0) {
        kokkos_policy = Uintah::Parallel::Kokkos_MDRange_Policy;
      } else if (strcmp(argv[i], "mdrange_rev") == 0) {
        kokkos_policy = Uintah::Parallel::Kokkos_MDRange_Reverse_Policy;
      } else {
        usage("Unknown Kokkos policy", arg, argv[0]);
        Uintah::Parallel::exitAll(1);
      }

      Uintah::Parallel::setKokkosPolicy(kokkos_policy);
#else
      std::cout << "Not compiled for Kokkos Range Policy support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_leagues_per_loop") {
#if defined(HAVE_KOKKOS)
      int kokkos_leagues_per_loop = 0;
      if (++i == argc) {
        usage("You must provide a number of Kokkos TeamPolicy leagues (work "
              "items) per loop for -kokkos_leagues_per_loop.",
              arg,
              argv[0]);
      }
      kokkos_leagues_per_loop = atoi(argv[i]);
      if (kokkos_leagues_per_loop < 1) {
        usage("Number of Kokkos TeamPolicy leagues (work items) per loop is "
              "too small.",
              arg,
              argv[0]);
        Uintah::Parallel::exitAll(1);
      }
      Uintah::Parallel::setKokkosLeaguesPerLoop(kokkos_leagues_per_loop);
#else
      std::cout << "Not compiled for GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_teams_per_league") {
#if defined(HAVE_KOKKOS)
      int kokkos_teams_per_block = 0;
      if (++i == argc) {
        usage("You must provide a number of Kokkos TeamPolicy teams (threads) "
              "per legaue for -kokkos_teams_per_block.",
              arg,
              argv[0]);
      }
      kokkos_teams_per_block = atoi(argv[i]);
      if (kokkos_teams_per_block < 1) {
        usage("Number of Kokkos TeamPolicy teams (threads) per legaue is too "
              "small.",
              arg,
              argv[0]);
        Uintah::Parallel::exitAll(1);
      }
      Uintah::Parallel::setKokkosTeamsPerLeague(kokkos_teams_per_block);
#else
      std::cout << "Not compiled for GPU support" << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_chunk_size") {
#if defined(HAVE_KOKKOS)
      int kokkos_chunk_size = 0;
      if (++i == argc) {
        usage(
          "You must provide the chunk size -kokkos_chunk_size.", arg, argv[0]);
      }
      kokkos_chunk_size = atoi(argv[i]);
      if (kokkos_chunk_size < 1) {
        usage("The Kokkos chunk size is too small.", arg, argv[0]);
        Uintah::Parallel::exitAll(1);
      }

      if (Uintah::Parallel::getKokkosPolicy() != Uintah::Parallel::Kokkos_Team_Policy &&
          Uintah::Parallel::getKokkosPolicy() != Uintah::Parallel::Kokkos_Range_Policy) {
        Uintah::Parallel::setKokkosPolicy(Uintah::Parallel::Kokkos_Range_Policy);
      }
      Uintah::Parallel::setKokkosChunkSize(kokkos_chunk_size);
#else
      std::cout << "Not compiled for Kokkos GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-kokkos_tile_size") {
#if defined(HAVE_KOKKOS)
      int kokkos_tile_isize = 0;
      if (++i == argc) {
        usage("You must provide the tile i index size -kokkos_tile_size.",
              arg,
              argv[0]);
      }
      kokkos_tile_isize = atoi(argv[i]);
      if (kokkos_tile_isize < 1) {
        usage("The Kokkos tile isize is too small.", arg, argv[0]);
        Uintah::Parallel::exitAll(1);
      }

      int kokkos_tile_jsize = 0;
      if (++i == argc) {
        usage("You must provide the tile j index size -kokkos_tile_size.",
              arg,
              argv[0]);
      }
      kokkos_tile_jsize = atoi(argv[i]);
      if (kokkos_tile_jsize < 1) {
        usage("The Kokkos tile jsize is too small.", arg, argv[0]);
        Uintah::Parallel::exitAll(0);
      }

      int kokkos_tile_ksize = 0;
      if (++i == argc) {
        usage("You must provide the tile k index size -kokkos_tile_size.",
              arg,
              argv[0]);
      }
      kokkos_tile_ksize = atoi(argv[i]);
      if (kokkos_tile_ksize < 1) {
        usage("The Kokkos tile ksize is too small.", arg, argv[0]);
        Uintah::Parallel::exitAll(1);
      }

      if (Uintah::Parallel::getKokkosPolicy() != Uintah::Parallel::Kokkos_MDRange_Policy &&
          Uintah::Parallel::getKokkosPolicy() !=
            Uintah::Parallel::Kokkos_MDRange_Reverse_Policy) {
        Uintah::Parallel::setKokkosPolicy(Uintah::Parallel::Kokkos_MDRange_Policy);
      }

      Uintah::Parallel::setKokkosTileSize(
        kokkos_tile_isize, kokkos_tile_jsize, kokkos_tile_ksize);
#else
      std::cout << "Not compiled for Kokkos GPU support." << std::endl;
      Uintah::Parallel::exitAll(0);
#endif
    } else if (arg == "-taskname_to_time") {
      // A hidden command line option useful for timing GPU tasks by forcing
      // this task name to wait until they can all be launched as a big group.
      // This helps time by avoiding any interleaving of other tasks in the way.
      // This command line option must be paired with two additional arguments.
      // The first being the name of the task
      // The second being the amount of times that task is expected to run in a
      // timestep.
      i++;
      std::string taskName = argv[i];
      i++;
      unsigned int amountTaskNameExpectedToRun = atoi(argv[i]);
      Uintah::Parallel::setTaskNameToTime(taskName);
      Uintah::Parallel::setAmountTaskNameExpectedToRun(
        amountTaskNameExpectedToRun);
    } else if (arg == "-t") {
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

  std::cout << "--------------------------------------------------------\n";
  std::cout << "Date:    " << time_string; // has its own newline
  std::cout << "Machine: " << name << "\n";
  std::cout << "Assertion level: " << SCI_ASSERTION_LEVEL << "\n";
  std::cout << "--------------------------------------------------------\n";

  Vaango::Utils::display_git_info(s_show_git_diff, s_show_git_status);
  std::cout << "--------------------------------------------------------\n";

  auto compile_cmd = Vaango::Utils::get_vaango_compile_command("vaango.cc");
  std::cout << "Compile cmd: " << compile_cmd << "\n";
  std::cout << "--------------------------------------------------------\n";

  if (s_show_version) {
    stop_mpi_and_exit(2);
  }
}

void
usage(const std::string& message,
      const std::string& badarg,
      const std::string& progname)
{
#if defined(KOKKOS_USING_GPU)
  std::string kokkos_str{"KokkosCUDA"};
#endif

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
    std::cerr << "-d[ebug]             :"
              << " Lists the available debug streams\n";
#if defined(KOKKOS_USING_GPU)
    std::cerr << "-cpu                            : "
              << " Use the CPU based MPI or Unified scheduler instead of the " 
              << kokkos_str << ".\n";
    std::cerr << "-gpu                            : "
              << " Use available GPU devices, requires multi-threaded "
              << " Unified scheduler \n";
    std::cerr << "-gpucheck                       : "
              << " Checks if there is a GPU available for the " 
              << kokkos_str << ".\n";
    std::cerr << "-kokkos_instances_per_task <#>  : "
              << " Number of Kokkos instances per task (default 1).\n";
    std::cerr << "-kokkos_policy <policy>         : "
              << " Kokkos Execution Policy - team, range, mdrange (default), "
              << " or mdrange_rev.\n";
    std::cerr << "-kokkos_leagues_per_loop <#>    : "
              << " Kokkos TeamPolicy number of leagues (work items) per loop "
              << " (default 1).\n";
    std::cerr << "-kokkos_teams_per_league <#>    : "
              << " Kokkos TeamPolicy number of teams (threads) per Kokkos "
              << " TeamPolicy league (default 256/16).\n";
    std::cerr << "-kokkos_chunk_size <#>          : "
              << " Kokkos TeamPolicy and RangePolicy chunk size.\n";
    std::cerr << "-kokkos_tile_size <# # #>       : "
              << " Kokkos MDRangePolicy tile size.\n";
#endif
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
    std::cerr << "-validate            :"
              << " Verifies the .ups file is valid and quits!\n";
    std::cerr << "-do_not_validate     :"
              << " Skips .ups file validation! Please avoid this flag"
              << " if possible.\n";
    std::cerr << "-cmake_cmd           :"
              << " Display cmake command used to build Vaango\n";
    std::cerr << "-local_filesystem    :"
              << " If using MPI, use this flag if each node has a local disk\n";
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
