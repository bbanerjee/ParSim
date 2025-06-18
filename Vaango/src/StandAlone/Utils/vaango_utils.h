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

#ifndef __VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__
#define __VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__

#include <string>

#include <sci_defs/git_info.h>

namespace Vaango {
namespace Utils {

inline static std::string s_input_ups_path{ "" };

void
set_input_ups_path(const std::string& path);

const std::string&
get_input_ups_path();

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
check_malloc();

std::string
get_vaango_compile_command(const std::string& vaango_filename);

} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_VAANGO_UTILS_H__
