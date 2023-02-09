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

#ifndef __VAANGO_STANDALONE_UTILS_COMPARE_UDA_OPTIONS_H__
#define __VAANGO_STANDALONE_UTILS_COMPARE_UDA_OPTIONS_H__

#include <Core/Geometry/IntVector.h>

#include <string>
#include <vector>

namespace Vaango {
namespace Utils {
namespace Options {

void
compare_uda_parse(int argc, char** argv);

void
compare_uda_usage(const std::string& badarg, const std::string& progname);

bool
tolerance_as_warnings();

bool
tolerance_error();

bool
concise();

bool
sort_variables();

bool
include_extra_cells();

bool
strict_types();

int
uda_levels(int level);

double
abs_tolerance();

double
rel_tolerance();

const std::string&
filebase_1();

const std::string&
filebase_2();

const std::vector<std::string>&
ignore_vars();

const std::vector<std::string>&
compare_vars();

void
abort_uncomparable(std::ostringstream& warn);

void
tolerance_failure();

//  parse user input and create a vector of strings.  Deliminiter is ","
std::vector<std::string>
parseVector(const char* input);


} // namespace Options
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_VAANGO_OPTIONS_H__
