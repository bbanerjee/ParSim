#include <Core/Const/const.h>

namespace dem {

// Pi
const REAL Pi = 3.141592653589;

// numerical EPS (NOT machine epsilon)
const REAL EPS = 1.0E-12;

// random number seed (Not a constant)
long idum = -1;

// output field width and precision
const std::size_t OWID = 15; // output width
const std::size_t OPREC =
  6; // output precision, number of digits after decimal dot

// other global variables
std::ofstream debugInf; // debug info
MPI_File overlapInf;    // contact overlap info, parallel I/O
std::size_t iteration;  // iteration number
REAL timeStep;          // time step
REAL timeAccrued;       // accurued time
REAL plan_gravity;      // used in PeriDEMBond.cpp, rho*l^2*g
}
