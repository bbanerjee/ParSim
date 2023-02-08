#include <StandAlone/Utils/vaango_utils.h>

#include <CCA/Components/Schedulers/UnifiedScheduler.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Environment.h>

namespace Vaango {
namespace Utils {

void
start_mpi()
{
  int argc    = 0;
  char** argv = nullptr;

  // Initialize MPI so that "usage" is only printed by proc 0.
  // (If we are using MPICH, then Uintah::MPI::Init() has already been called.)
  Uintah::Parallel::initializeManager(argc, argv);
}

void
stop_mpi_and_exit(int flag, const std::string& msg)
{
  if (msg != "") {
    std::cerr << msg << "\n";
  }
  Uintah::Parallel::finalizeManager();
  Uintah::Parallel::exitAll(flag);
}

void
stop_mpi_with_abort()
{
  Uintah::Parallel::finalizeManager(Uintah::Parallel::Abort);
}

void
print_active_debug_streams()
{
  start_mpi();
  if (Uintah::Parallel::getMPIRank() == 0) {
    // report all active Dout debug objects
    std::cout << "\nThe following Douts are known. Active Douts are indicated "
                 "with plus sign."
              << std::endl;
    std::cout << "To activate a Dout, set the environment variable 'setenv "
                 "SCI_DEBUG \"Dout_Name:+\"'"
              << std::endl;
    Uintah::Dout::printAll();

    // report all active DebugStreams
    std::cout << "\nThe following DebugStreams are known. Active streams are "
                 "indicated with plus sign."
              << std::endl;
    std::cout << "To activate a DebugStreams set the environment variable "
                 "'setenv SCI_DEBUG \"Debug_Stream_Name:+\"'"
              << std::endl;
    Uintah::DebugStream::printAll();
  }
  stop_mpi_and_exit(2);
}

void
check_gpus()
{
#ifdef HAVE_CUDA
  int retVal = UnifiedScheduler::verifyAnyGpuActive();
  if (retVal == 1) {
    std::cout << "At least one GPU detected!" << std::endl;
  } else {
    std::cout << "No GPU detected!" << std::endl;
  }
  Uintah::Parallel::exitAll(retVal);
#endif
  std::cout << "No GPU detected!" << std::endl;
  Uintah::Parallel::exitAll(2); // If the above didn't exit with a 1, then we
                                // didn't have a GPU, so exit with a 2.
  std::cout << "This doesn't run" << std::endl;
}

void
display_git_info(bool show_git_diff, bool show_git_status)
{
  // Run git commands Uintah
  std::cout << "git branch:" << GIT_BRANCH << "\n";
  std::cout << "git date:   " << GIT_DATE << "\n";
  std::cout << "git hash:   " << GIT_HASH << "\n";

  if (show_git_diff || show_git_status) {
    std::cout << "GIT::\n";
    std::string sdir = std::string(Uintah::sci_getenv("SCIRUN_SRCDIR"));

    if (show_git_diff) {
      std::string cmd =
        "cd " + sdir + "; git --no-pager diff  --no-color --minimal";
      std::cout << "\n git diff::\n";
      std::system(cmd.c_str());
    }

    if (show_git_status) {
      std::string cmd = "cd " + sdir + "; git status  --branch --short";
      std::cout << "\n git status --branch --short::\n";
      std::system(cmd.c_str());

      cmd = "cd " + sdir + "; git log -1  --format=\"%ad %an %H\" | cat";
      std::cout << "\n git log -1::\n";
      std::system(cmd.c_str());
    }
    std::cout << "::GIT\n";
  }
}

void
display_config_info(bool show_config_cmd)
{
  if (show_config_cmd) {
    std::string odir = std::string(Uintah::sci_getenv("SCIRUN_OBJDIR"));
    std::string cmd  = "cd " + odir + "; sed -n '7'p config.log ";
    std::cout << "Configure Command::\n";
    std::system(cmd.c_str());
    std::cout << "\n\n";
  }
}

void
check_malloc()
{
#if defined(DISABLE_SCI_MALLOC)
  if (getenv("MALLOC_STATS")) {
    printf("\nERROR:\n");
    printf("ERROR: Environment variable MALLOC_STATS set, but  "
           "--enable-sci-malloc was not configured...\n");
    printf("ERROR:\n\n");
    Uintah::Parallel::exitAll(1);
  }
  if (getenv("MALLOC_TRACE")) {
    printf("\nERROR:\n");
    printf("ERROR: Environment variable MALLOC_TRACE set, but  "
           "--enable-sci-malloc was not configured...\n");
    printf("ERROR:\n\n");
    Uintah::Parallel::exitAll(1);
  }
  if (getenv("MALLOC_STRICT")) {
    printf("\nERROR:\n");
    printf("ERROR: Environment variable MALLOC_STRICT set, but "
           "--enable-sci-malloc  was not configured...\n");
    printf("ERROR:\n\n");
    Uintah::Parallel::exitAll(1);
  }
#endif
}

} // namespace Utils
} // namespace Vaango
