/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MATITI_MPMTIME_H__
#define __MATITI_MPMTIME_H__

#include <Time.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace BrMPM {

  class MPMTime : public Time
  {
//  public:
//    friend std::ostream& operator<<(std::ostream& out, const Matiti::Time& time);

  public:
    MPMTime();
    MPMTime(const Uintah::ProblemSpecP& ps);
    virtual ~MPMTime();

    void initialize(const Uintah::ProblemSpecP& ps);

//    void setDelT(double dt) {d_delT = dt;}

    inline double initTime() const { return d_init_time;}
    inline double finalTime() const { return d_final_time;}
    inline double CFL() const {return d_cfl;}
//    inline double delT() const { return d_delT;}
//    inline int maxIter() const { return d_max_iter;}
//    inline double timeStepFactor() const {return d_factor;}

//    double incrementTime(const double delT) {d_cur_time += delT; return d_cur_time;}
//    double currentTime() const {return d_cur_time;}

  private:

    // Simulation time 
    double d_init_time;  // Initial time for which simulation is to be run
    double d_final_time; // Final time for which simulation is to be run
    double d_cfl;         // CFL Condition
//    double d_delT;      // Timestep size
    int d_max_iter;     // Maximum number of iterations
    double d_factor;    // Timestep reduction factor
    double d_cur_time;  // current time

//    Domain d_domain;    // Object of the Domain class
//    MaterialUP d_mat;   // Object of the Material class

  }; // end class

} // end namespace

#endif

