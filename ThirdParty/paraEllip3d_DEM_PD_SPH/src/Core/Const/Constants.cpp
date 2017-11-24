#include <Core/Const/Constants.h>

namespace dem {

// Pi
const REAL Pi = 3.141592653589;

// numerical EPS (NOT machine epsilon)
const REAL EPS = 1.0E-12;

// random number seed (Not a constant)
long idum = -1;

// output field width 
const std::size_t OWID = 15; 

// output precision, number of digits after decimal dot
const std::size_t OPREC = 6; 

// other global variables
std::ofstream debugInf; // debug info
MPI_File overlapInf;    // contact overlap info, parallel I/O
std::size_t g_iteration;  // iteration number
REAL g_timeStep;          // time step
REAL g_timeAccrued;       // accrued time
REAL plan_gravity;      // used in PeriDEMBond.cpp, rho*l^2*g
}
