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
 * The MIT License
 *
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
 *  SimpleMath.cc: Simple math benchmark.
 *
 *  Written by:
 *   Randy Jones
 *   Department of Computer Science
 *   University of Utah
 *   July 31, 2003
 *
 */

#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Util/Timers/Timers.hpp>
#include <iostream>

using namespace Uintah;
using namespace std;

const int SIZE_DEFAULT = 80;
const int LOOP_DEFAULT = 1;

void
usage(void)
{
  std::cerr << "Usage: SimpleMath <size> [<loop>]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  <size>  This benchmark will use four CCVariable<double>"
            << std::endl;
  std::cerr << "          variables, each with a size*size*size resolution"
            << std::endl;
  std::cerr << "          and perform the operation 'result = a * x + b'"
            << std::endl;
  std::cerr << "          on each of their elements." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  <loop>  The above operation will be repeated <loop> times."
            << std::endl;
}

void
bench1(int loop,
       const IntVector& low,
       const IntVector& high,
       CCVariable<double>& result,
       double a,
       CCVariable<double>& x,
       CCVariable<double>& b)
{
  for (int i = 0; i < loop; i++) {
    for (CellIterator iter(low, high); !iter.done(); iter++) {
      result[*iter] = a * x[*iter] + b[*iter];
    }
  }
}

void
bench2(int loop,
       const IntVector& low,
       const IntVector& high,
       CCVariable<double>& result,
       double a,
       CCVariable<double>& x,
       CCVariable<double>& b)
{
  for (int i = 0; i < loop; i++) {
    IntVector d(high - low);
    int size         = d.x() * d.y() * d.z();
    double* rr       = &result[low];
    const double* xx = &x[low];
    const double* bb = &b[low];
    for (int i = 0; i < size; i++) {
      rr[i] = a * xx[i] + bb[i];
    }
  }
}

int
main(int argc, char** argv)
{
  int size = SIZE_DEFAULT;
  int loop = LOOP_DEFAULT;

  /*
   * Parse arguments
   */
  if (argc > 1) {
    size = atoi(argv[1]);

    if (size <= 0) {
      usage();
      return EXIT_FAILURE;
    }

    if (argc > 2) {
      loop = atoi(argv[2]);

      if (loop <= 0) {
        usage();
        return EXIT_FAILURE;
      }
    }
  } else {
    usage();
    return EXIT_FAILURE;
  }

  std::cout << "Simple Math Benchmark: " << std::endl;
  std::cout << "Resolution (" << size << ", " << size << ", " << size << ")"
            << std::endl;
  std::cout << "Repeating " << loop << " time(s)." << std::endl;

  IntVector low(0, 0, 0);
  IntVector high(size, size, size);

  CCVariable<double> result, x, b;
  double a = 5;

  result.allocate(low, high);
  //  a.allocate( low, high );
  x.allocate(low, high);
  b.allocate(low, high);

  // a.initialize( 5 );
  x.initialize(6);
  b.initialize(2);

  Timers::Simple timer;

  {

    timer.start();
    bench1(loop, low, high, result, a, x, b);
    timer.stop();

    double megaFlops =
      (loop * size * size * size * 2.0) / 1000000.0 / timer().seconds();

    std::cout << "Completed in " << timer().seconds() << " seconds.";
    std::cout << " (" << megaFlops << " MFLOPS)" << std::endl;
  }
  {
    timer.reset(true);
    bench2(loop, low, high, result, a, x, b);
    timer.stop();

    double megaFlops =
      (loop * size * size * size * 2.0) / 1000000.0 / timer().seconds();

    std::cout << "Completed in " << timer().seconds() << " seconds.";
    std::cout << " (" << megaFlops << " MFLOPS)" << std::endl;
  }

  return EXIT_SUCCESS;
}
