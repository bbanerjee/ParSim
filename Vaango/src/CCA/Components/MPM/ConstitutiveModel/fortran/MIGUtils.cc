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
#include <Core/Parallel/Parallel.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

namespace Uintah {

extern "C"
{

  // These functions are called from fortran and thus need to be
  // compiled using "C" naming conventions for the symbols.

  void
  bombed_(char* mes, int len_mes)
  {
    std::cerr << "Code bombed with the following message:" << std::endl;
    for (int i = 0; i < len_mes; i++) {
      putchar(mes[i]);
    }
    std::cerr << "\n";
    exit(1);
    return;
  }

  void
  logmes_(char* mes, int len_mes)
  {
    if (Uintah::Parallel::getMPIRank() == 0) {
      for (int i = 0; i < len_mes; i++) {
        putchar(mes[i]);
      }
      proc0cout << "\n";
    }

    return;
  }

  void
  faterr_(char* mes1, char* mes2, int len_mes1, int len_mes2)
  {
    std::cerr << "FATAL ERROR DETECTED BY ";
    for (int i = 0; i < len_mes1; i++) {
      putchar(mes1[i]);
    }
    std::cerr << ":\n";

    for (int i = 0; i < len_mes2; i++) {
      putchar(mes2[i]);
    }
    std::cerr << "\n";
    exit(1);

    return;
  }

  void
  log_error_(char* mes, int len_mes)
  {
    std::cerr << "**ERROR** Code bombed with the following message:"
              << std::endl;
    for (int i = 0; i < len_mes; i++) {
      putchar(mes[i]);
    }
    std::cerr << "\n";
    exit(1);
    return;
  }

  void
  log_warning_(char* mes, int len_mes)
  {
    if (Uintah::Parallel::getMPIRank() == 0) {
      for (int i = 0; i < len_mes; i++) {
        putchar(mes[i]);
      }
      proc0cout << "\n";
    }

    return;
  }
  void
  log_message_(char* mes, int len_mes)
  {
    if (Uintah::Parallel::getMPIRank() == 0) {
      for (int i = 0; i < len_mes; i++) {
        putchar(mes[i]);
      }
      proc0cout << "\n";
    }

    return;
  }

} // end extern "C"

} // namespace Uintah