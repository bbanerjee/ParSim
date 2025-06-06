#ifndef CONST_H
#define CONST_H
#include <Core/Types/RealTypes.h>
#include <cstddef>
#include <fstream>
#include <mpi.h>

namespace dem {

// Pi
extern const REAL Pi;

// numerical EPS (NOT machine epsilon)
extern const REAL EPS;

// random number seed (NOT a constant)
extern long idum;

// output field width and precision
extern const std::size_t OWID;
extern const std::size_t OPREC;

// other global variables
extern std::ofstream debugInf;
extern MPI_File overlapInf;
extern std::size_t g_iteration;
extern REAL g_timeStep;
extern REAL g_timeAccrued;
extern REAL plan_gravity;
}
#endif
