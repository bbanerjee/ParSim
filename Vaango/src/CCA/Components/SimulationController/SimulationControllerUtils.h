#ifndef UINTAH_HOMEBREW_SIMULATIONCONTROLLER_UTILS_H
#define UINTAH_HOMEBREW_SIMULATIONCONTROLLER_UTILS_H

namespace Uintah {
  struct double_int
  {
     double val;
     int loc;
     double_int(double val, int loc): val(val), loc(loc) {}
     double_int(): val(0), loc(-1) {}
  };

  inline double stdDeviation(double sum_of_x, double sum_of_x_squares, int n)
  {
    return sqrt((n*sum_of_x_squares - sum_of_x*sum_of_x)/(n*n));
  }
}

#endif
