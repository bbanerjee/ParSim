#ifndef __MATITI_MPMTIME_H__
#define __MATITI_MPMTIME_H__

#include <Time.h>
#include <Domain.h>
#include <MaterialUP.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace Matiti {

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

