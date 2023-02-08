#ifndef __VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__
#define __VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__

#include <string>

#include <sci_defs/git_info.h>

namespace Vaango {
namespace Utils {

void
start_mpi();

void
stop_mpi_and_exit(int flag = 2, const std::string& msg = "");

void
stop_mpi_with_abort();

void
print_active_debug_streams();

void
check_gpus();

void
display_git_info(bool show_git_diff, bool show_git_status);

void
display_config_info(bool show_config_cmd);

void
check_malloc();

} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__
