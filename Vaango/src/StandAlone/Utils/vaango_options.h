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
