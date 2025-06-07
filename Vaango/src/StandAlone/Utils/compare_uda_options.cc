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

#include <StandAlone/Utils/compare_uda_options.h>

#include <Core/Parallel/Parallel.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/FileUtils.h>

#include <sci_defs/compile_defs.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

static bool s_tolerance_as_warnings{ false };
static bool s_tolerance_error{ false };
static bool s_concise{ false };
static bool s_sort_variables{ true };
static bool s_include_extra_cells{ false };
static bool s_strict_types{ true };

static int s_uda_levels[2] = { -9, -9 };

static double s_abs_tolerance{ 1.0e-6 };
static double s_rel_tolerance{ 1.0e-9 };

static std::string s_filebase_1{ "" };
static std::string s_filebase_2{ "" };

static std::vector<std::string> s_ignore_vars;
static std::vector<std::string> s_compare_vars;

} // end anonymous namespace

namespace Vaango {
namespace Utils {
namespace Options {

void
compare_uda_parse(int argc, char** argv)
{
  for (int i = 1; i < argc; i++) {
    std::string s = argv[i];

    // remove "-" and convert to upper case.
    s.erase(std::remove(s.begin(), s.end(), '-'), s.end());
    s = Uintah::string_toupper(s);

    if (s == "ABS_TOLERANCE") {
      if (++i == argc) {
        compare_uda_usage("-abs_tolerance, no value given", argv[0]);
      } else {
        s_abs_tolerance = atof(argv[i]);
      }
    } else if (s == "REL_TOLERANCE") {
      if (++i == argc) {
        compare_uda_usage("-rel_tolerance, no value given", argv[0]);
      } else {
        s_rel_tolerance = atof(argv[i]);
      }
    } else if (s == "AS_WARNINGS") {
      s_tolerance_as_warnings = true;
    } else if (s == "CONCISE") {
      s_concise = true;
    } else if (s == "DONT_SORT") {
      s_sort_variables = false;
    } else if (s == "INCLUDEEXTRACELLS") {
      s_include_extra_cells = true;
    } else if (s == "SKIP_UNKNOWN_TYPES") {
      s_strict_types = false;
    } else if (s == "EXACT") {
      s_abs_tolerance = 0;
      s_rel_tolerance = 0;
    } else if (s == "LEVELS") {
      s_uda_levels[0] = atoi(argv[++i]);
      std::cout << "udaLevels: " << s_uda_levels[0];
      s_uda_levels[1] = atoi(argv[++i]);
      std::cout << "  " << s_uda_levels[1];
    } else if (s == "IGNOREVARIABLES") {
      if (++i == argc) {
        compare_uda_usage("-ignoreVariables, no variable given", argv[0]);
      } else {
        s_ignore_vars = parseVector(argv[i]);
      }
    } else if (s == "COMPAREVARIABLES") {
      if (++i == argc) {
        compare_uda_usage("-compareVariables, no variable given", argv[0]);
      } else {
        s_compare_vars = parseVector(argv[i]);
      }
    } else if (s == "H") {
      compare_uda_usage("", argv[0]);
    } else {
      if (s_filebase_1 != "") {
        if (s_filebase_2 != "") {
          compare_uda_usage(s, argv[0]);
        } else {
          s_filebase_2 = argv[i];
        }
      } else {
        s_filebase_1 = argv[i];
      }
    }
  }

  //__________________________________
  //  bulletproofing
  if (s_filebase_2 == "") {
    std::cerr << "\nYou must specify two archive directories.\n";
    compare_uda_usage("", argv[0]);
  }

  if (!Uintah::validDir(s_filebase_1)) {
    std::cerr << "\nParameter '" << s_filebase_1
              << "' is not a valid directory.\n";
    compare_uda_usage("", argv[0]);
  }
  if (!Uintah::validDir(s_filebase_2)) {
    std::cerr << "\nParameter '" << s_filebase_2
              << "' is not a valid directory.\n";
    compare_uda_usage("", argv[0]);
  }

  std::cerr << "Using absolute tolerance: " << s_abs_tolerance << std::endl;
  std::cerr << "Using relative tolerance: " << s_rel_tolerance << std::endl;

  if (s_include_extra_cells) {
    std::cerr
      << "Comparision domain:  interior patch cells, including extra cells"
      << std::endl;
  } else {
    std::cerr
      << "Comparision domain:  interior patch cells, excluding extra cells"
      << std::endl;
  }

  if (s_uda_levels[0] != -9) {
    std::cerr << "Comparing uda1: Level(" << s_uda_levels[0]
              << ") against uda2: level(" << s_uda_levels[1] << ")" << std::endl;
  }

  if (s_rel_tolerance < 0) {
    std::cerr << "Must have a non-negative value rel_tolerance.\n";
    Uintah::Parallel::exitAll(1);
  }
}

void
compare_uda_usage(const std::string& badarg, const std::string& progname)
{
  if (badarg != "") {
    std::cerr << "\nError parsing argument: " << badarg << '\n';
  }
  std::cerr
    << "\nUsage: " << progname
    << " [options] <UDA archive directory 1> <UDA archive directory 2>\n\n";
  std::cerr << "Valid options are:\n";
  std::cerr << "  -h[elp]\n";
  std::cerr
    << "  -abs_tolerance [double]          (Allowable absolute difference "
       "of any number, default: 1e-9)\n";
  std::cerr
    << "  -rel_tolerance [double]          (Allowable relative difference "
       "of any number, default: 1e-6)\n";
  std::cerr
    << "  -exact                           (Perform an exact comparison, "
       "absolute/relative tolerance = 0)\n";
  std::cerr
    << "  -levels     [int int]            (Optional:  level index for uda "
       "1 and uda 2)\n";
  std::cerr << "  -includeExtraCells               (include extra cells in the "
               "comparison\n";
  std::cerr << "  -as_warnings                     (Treat tolerance errors as "
               "warnings and continue)\n";
  std::cerr
    << "  -concise                         (With '-as_warnings', only print "
       "first incidence of error per var.)\n";
  std::cerr
    << "  -skip_unknown_types              (Skip variable comparisons of "
       "unknown types without error)\n";
  std::cerr
    << "  -ignoreVariables [var1,var2....] (Skip these variables. Comma "
       "delimited list, no spaces.)\n";
  std::cerr << "  -compareVariables[var1,var2....] (Only these variables are "
               "compared. Comma delimited list, no spaces.)\n";
  std::cerr
    << "  -dont_sort                       (Don't sort the variable names "
       "before comparing them)";
  std::cerr
    << "\nNote: The absolute and relative tolerance tests must both fail\n"
    << "      for a comparison to fail.\n\n";
  std::cerr << "  Exit values:\n";
  std::cerr << "    -1:      Comparison failed, variable type not supported.\n";
  std::cerr << "     0:      Comparison passed.\n";
  std::cerr << "     1:      Error in input parameters.\n";
  std::cerr << "     2:      Comparison failed, tolerances exceeded. \n";
  std::cerr << "     5:      The uda directories may not be compared.\n";

  Uintah::Parallel::exitAll(1);
}

bool
tolerance_as_warnings()
{
  return s_tolerance_as_warnings;
}

bool
tolerance_error()
{
  return s_tolerance_error;
}

bool
concise()
{
  return s_concise;
}

bool
sort_variables()
{
  return s_sort_variables;
}

bool
include_extra_cells()
{
  return s_include_extra_cells;
}

bool
strict_types()
{
  return s_strict_types;
}

int
uda_levels(int level)
{
  return s_uda_levels[level];
}

double
abs_tolerance()
{
  return s_abs_tolerance;
}

double
rel_tolerance()
{
  return s_rel_tolerance;
}

const std::string&
filebase_1()
{
  return s_filebase_1;
}

const std::string&
filebase_2()
{
  return s_filebase_2;
}

const std::vector<std::string>&
ignore_vars()
{
  return s_ignore_vars;
}

const std::vector<std::string>&
compare_vars()
{
  return s_compare_vars;
}

void
abort_uncomparable(std::ostringstream& warn)
{
  std::cerr << "\n_______________ERROR:compare_uda___________________\n";

  int error = -9;
  if (s_tolerance_error) {
    std::cerr << " The differences have exceeded the tolerances to the point that "
            "the udas can no longer be compared.\n";
    std::cerr << " Now exiting (2).\n";
    error = 2;
  } else {

    std::cerr << warn.str();
    std::cerr << "\n  The uda directories may not be compared.\n";
    std::cerr << " Now exiting (5).\n";
    error = 5;
  }
  std::cerr << "\n_______________ERROR:compare_uda___________________\n";
  Uintah::Parallel::exitAll(error);
}

void
tolerance_failure()
{
  if (s_tolerance_as_warnings) {
    s_tolerance_error = true;
    std::cerr << std::endl;
  }
  else
    Uintah::Parallel::exitAll(2);
}

//  parse user input and create a vector of strings.  Deliminiter is ","
std::vector<std::string>
parseVector(const char* input)
{
  std::vector<std::string> result;
  std::stringstream ss(input);

  while (ss.good()) {
    std::string substr;
    getline(ss, substr, ',');
    result.push_back(substr);
  }
  return result;
}


} // namespace Options
} // namespace Utils
} // namespace Vaango
