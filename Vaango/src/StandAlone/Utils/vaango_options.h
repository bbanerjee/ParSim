#ifndef __VAANGO_STANDALONE_UTILS_VAANGO_OPTIONS_H__
#define __VAANGO_STANDALONE_UTILS_VAANGO_OPTIONS_H__

#include <Core/Geometry/IntVector.h>

#include <string>

namespace Vaango {
namespace Utils {
namespace Options {

void
parse(int argc, char** argv);

void
usage(const std::string& message,
      const std::string& badarg,
      const std::string& progname);

bool
emit_graphs();

bool
local_filesystem();

bool
only_validate_ups();

bool
postprocess_uda();

bool
restart();

bool
restart_from_scratch();

bool
restart_remove_old_dir();

bool
show_config_cmd();

bool
show_git_diff();

bool
show_git_status();

bool
show_version();

bool
validate_ups();

bool
use_gpu();

int
num_partitions();

int
num_threads();

int
restart_checkpoint_index();

int
threads_per_partition();

int
uda_suffix();

const std::string&
uda_dir();

const std::string&
uda_filename();

const std::string&
solver_name();

const Uintah::IntVector&
grid_layout();

} // namespace Options
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_VAANGO_OPTIONS_H__
