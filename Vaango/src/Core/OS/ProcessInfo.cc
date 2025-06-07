/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

/*
 *  ProcessInfo.cc:
 *
 *  Written by:
 *   Author: Randy Jones
 *   Department of Computer Science
 *   University of Utah
 *   Date: 2004/02/05
 *
 */

#include <Core/OS/ProcessInfo.h>

#include <cstdio>
#include <cstring>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace Uintah {

bool
ProcessInfo::isSupported(int info_type)
{
#if defined(__linux)

  switch (info_type) {
    case MEM_SIZE:
      return true;
    case MEM_RSS:
      return true;
    default:
      return false;
  }

#else
  return false;
#endif
}

unsigned long
ProcessInfo::getInfo(int info_type)
{

#if defined(__linux)

  char statusFileName[MAXPATHLEN];
  sprintf(statusFileName, "/proc/%d/status", getpid());

  FILE* file = fopen(statusFileName, "r");

  if (file != nullptr) {
    unsigned long tempLong = 0;
    char tempString[1024];
    const char* compareString;

    switch (info_type) {
      case MEM_SIZE:
        compareString = "VmSize:";
        break;
      case MEM_RSS:
        compareString = "VmRSS:";
        break;
      default:
        fclose(file);
        return 0;
    }

    while (!feof(file)) {
      [[maybe_unused]] auto status = fscanf(file, "%s", tempString);
      if (!strcmp(tempString, compareString)) {
        [[maybe_unused]] auto stat = fscanf(file, "%ld", &tempLong);
        fclose(file);
        return tempLong * 1024;
      }
    }
  }
#endif

  return 0;

} // unsigned long ProcessInfo::getInfo ( int info_type )

std::string
ProcessInfo::toHumanUnits(unsigned long value)
{
  char tmp[64];

  sprintf(tmp, "%.2lf MBs", value / 1000000.0);
  return tmp;
}

} // namespace Uintah
