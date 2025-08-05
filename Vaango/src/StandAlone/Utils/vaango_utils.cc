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

#include <StandAlone/Utils/vaango_utils.h>

#include <CCA/Components/Schedulers/UnifiedScheduler.h>
#include <CCA/Components/Schedulers/KokkosScheduler.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Environment.h>

#include <expected>
#include <filesystem>
#include <memory>
#include "vaango_utils.h"
#include <submodules/json/single_include/nlohmann/json.hpp>

namespace fs = std::filesystem;
namespace Vaango {
namespace Utils {

void set_input_ups_path(const std::string& path)
{
  s_input_ups_path = fs::path(path).parent_path(); 
}

const std::string& get_input_ups_path() 
{
  return s_input_ups_path;
}

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
  try {
    Uintah::Parallel::finalizeManager();
  } catch (const Uintah::InternalError& e) {
    std::cout << std::string(e.message()) << std::endl;
  } catch (const std::exception& e) {
    std::cout << std::string(e.what()) << std::endl;
  }
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
  int retVal = KokkosScheduler::verifyAnyGpuActive();
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

#ifdef HAVE_CUDA
// This function executes a shell command and returns its output.
// It returns std::nullopt if the command fails to execute.
// You should handle the error case where the command fails by checking if the
// optional has a value.
std::optional<std::string>
executeCommand(const std::string& command)
{
  std::array<char, 128> buffer;
  std::string result;

  // Use a custom deleter to ensure pclose is called.
  auto pipe_deleter = [](FILE* f) {
    if (f) {
      pclose(f);
    }
  };
  std::unique_ptr<FILE, decltype(pipe_deleter)> pipe(
    popen(command.c_str(), "r"), pipe_deleter);

  if (!pipe) {
    std::cerr << "Failed to execute command: " << command << std::endl;
    return std::nullopt; // Use std::nullopt to signify an error.
  }

  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }

  // Remove the trailing newline character if one exists.
  if (!result.empty() && result.back() == '\n') {
    result.pop_back();
  }

  return result;
}
#else
std::expected<std::string, std::string>
executeCommand(const std::string& command)
{
  std::array<char, 128> buffer;
  std::string result;

  // Use popen to execute command and capture output
  auto pipe_deleter = [](FILE* f) {
    if (f) {
      pclose(f);
    }
  };
  auto pipe = std::unique_ptr<FILE, decltype(pipe_deleter)>(
    popen(command.c_str(), "r"), pipe_deleter);

  if (!pipe) {
    return std::unexpected("Failed to execute command: " + command);
  }

  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }

  // Remove trailing newline if present
  if (!result.empty() && result.back() == '\n') {
    result.pop_back();
  }

  return result;
}
#endif

struct LiveGitInfo
{
  std::string diff;
  std::string status;
  std::string error;
};

LiveGitInfo
getLiveGitInfo()
{

#ifdef HAVE_CUDA
  LiveGitInfo info;
  // Get git diff
  if (auto diff_result = executeCommand(
        "git --no-pager diff --no-color --minimal 2>/dev/null")) {
    info.diff = *diff_result;
  } else {
    // If the command fails, the optional will be empty.
    // The previous `error` member variable is not updated here because
    // the error is now handled within `executeCommand`.
    // Instead, we just set a general error message and return.
    info.error = "Failed to get git diff.";
    return info;
  }

  // Get git status
  if (auto status_result =
        executeCommand("git status --branch --short 2>/dev/null")) {
    info.status = *status_result;
  } else {
    // Similarly, handle the failure case for the status command.
    info.error = "Failed to get git status.";
    return info;
  }

  return info;
#else
  LiveGitInfo info;
  // Get git diff
  if (auto diff_result = executeCommand(
        "git --no-pager diff --no-color --minimal 2>/dev/null")) {
    info.diff = diff_result.value();
  } else {
    info.error = "Failed to get git diff: " + diff_result.error();
    return info;
  }

  // Get git status
  if (auto status_result =
        executeCommand("git status --branch --short 2>/dev/null")) {
    info.status = status_result.value();
  } else {
    info.error = "Failed to get git status: " + status_result.error();
    return info;
  }
  return info;
#endif

}

void
display_git_info(bool show_git_diff, bool show_git_status)
{
  // Run git commands Uintah
  //std::cout << "CMAKE_SRC_DIR:" << CMAKE_SRC_DIR << "\n";
  std::cout << "git branch: " << GIT_BRANCH << "\n";
  std::cout << "git date:   " << GIT_DATE << "\n";
  std::cout << "git hash:   " << GIT_HASH << "\n";

  //std::cout << std::boolalpha << " show diff: " << show_git_diff << ", " << show_git_status << "\n";
  if (show_git_diff || show_git_status) {
    auto info = getLiveGitInfo();
    std::cout << "GIT::\n";
    std::cout << "git diff:   " << info.diff << "\n";
    std::cout << "git status:   " << info.status << "\n";
    std::cout << "GIT::\n";
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

std::string
get_vaango_compile_command(const std::string& vaango_filename)
{
  #ifndef COMPILE_COMMANDS_PATH
  #define COMPILE_COMMANDS_PATH "None"
  #endif

  std::string compile_commands_path = COMPILE_COMMANDS_PATH;

  // Check if the file exists
  if (!std::filesystem::exists(compile_commands_path)) {
      return "Error: compile_commands.json not found at: " + compile_commands_path;
  }

  //std::cout << "COMPILE_COMMANDS_PATH: " << COMPILE_COMMANDS_PATH << "\n";

  try {

    // Read and parse the JSON file
    std::ifstream file(compile_commands_path);
    if (!file.is_open()) {
      return "Error: Could not open compile_commands.json";
    }

    nlohmann::json compile_commands;
    file >> compile_commands;

    // Search for the file in the compile commands
    for (const auto& entry : compile_commands) {
      if (entry.contains("file") && entry.contains("command")) {
        std::string file_path = entry["file"];

        // Check if the filename appears in the file path
        if (file_path.find(vaango_filename) != std::string::npos) {
          return entry["command"];
        }
      }
    }

    return "Error: No compile command found for file: " + vaango_filename;

  } catch (const nlohmann::json::exception& e) {
    return "Error parsing JSON: " + std::string(e.what());
  } catch (const std::exception& e) {
    return "Error: " + std::string(e.what());
  }
}

} // namespace Utils
} // namespace Vaango
