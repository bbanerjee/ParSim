#ifndef __EMU2DC_TIME_H__
#define __EMU2DC_TIME_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace Emu2DC {

  class Time 
  {
  public:
    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Time& time);

  public:
    Time();
    Time(const Uintah::ProblemSpecP& ps);
    ~Time();

    void initialize(const Uintah::ProblemSpecP& ps);

    void setDelT(double dt) {d_delT = dt;}

    inline double maxTime() const { return d_max_time;}
    inline double delT() const { return d_delT;}
    inline int maxIter() const { return d_max_iter;}
    inline double timeStepFactor() const {return d_factor;}

    double incrementTime(const double delT) {d_cur_time += delT; return d_cur_time;}
    double currentTime() const {return d_cur_time;}

  private:

    // Simulation time 
    double d_max_time;  // Maximum time for which simulation is to be run
    double d_delT;      // Timestep size
    int d_max_iter;     // Maximum number of iterations
    double d_factor;    // Timestep reduction factor

    double d_cur_time;  // current time

  }; // end class

} // end namespace

#endif

