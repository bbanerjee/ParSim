#ifndef FLUID_H
#define FLUID_H

#include <Boundary/Boundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/Containers.h>
#include <DiscreteElements/Gradation.h>
#include <DiscreteElements/Particle.h>
#include <InputOutput/Parameter.h>
#include <cstddef>
#include <valarray>

namespace dem {

class Fluid
{
  typedef std::valarray<std::valarray<std::valarray<REAL>>> Array3D;
  typedef std::valarray<std::valarray<std::valarray<std::valarray<REAL>>>>
    Array4D;
  typedef std::valarray<
    std::valarray<std::valarray<std::valarray<std::valarray<REAL>>>>>
    Array5D;

private:
  static constexpr REAL Rs = 287.06; // specific gas constant

  std::size_t nx; // nx = number of total cell centers = parts + two boundary
                  // points in x direction
  std::size_t ny; // ny = number of total cell centers = parts + two boundary
                  // points in y direction
  std::size_t nz; // nz = number of total cell centers = parts + two boundary
                  // points in z direction
  REAL dx;
  REAL dy;
  REAL dz;

  REAL RK;           // Runge-Kutta scheme
  REAL CFL;          // Courant-Friedrichs-Lewy condition
  REAL gamma;        // ration of specific heat
  REAL arrayBC[6];   // boundary condition
  REAL z0;           // initial discontinuity plane in Z direction
  REAL rhoL, uL, pL; // unknown
  REAL rhoR, uR, pR; // known
  REAL mach;         // shock Mach number, known
  REAL shockSpeed;   // unknown
  REAL Cd;           // drag coefficient
  int volFrac;       // use grid volume fraction or not

  std::size_t n_dim, n_var, n_integ, var_den, var_eng, var_prs, var_msk;
  std::size_t var_mom[3], var_vel[3];

  Array4D arrayU;
  Array4D arrayUtmp;
  // 4-dimensional, defined at cell centers
  // nx, ny, nz, n_var
  // (a) fixed:
  // arrayU[i][j][k][0]: var_den
  // arrayU[i][j][k][1]: var_mom[0]
  // arrayU[i][j][k][2]: var_mom[1]
  // arrayU[i][j][k][3]: var_mom[2]
  // arrayU[i][j][k][4]: var_eng
  // arrayU[i][j][k][5]: var_vel[0]
  // arrayU[i][j][k][6]: var_vel[1]
  // arrayU[i][j][k][7]: var_vel[2]
  // arrayU[i][j][k][8]: var_prs
  // (b) extended:
  // arrayU[i][j][k][9]: var_msk

  Array4D arrayGridCoord;
  // fluid grid coordinates, 4-dimensional
  // nx, ny, nz, n_dim
  // arrayGridCoord[i][j][k][0]: coord_x
  // arrayGridCoord[i][j][k][1]: coord_y
  // arrayGridCoord[i][j][k][2]: coord_z

  Array4D arrayPenalForce;
  // fluid grid forces, 4-dimensional
  // nx, ny, nz, n_dim
  // arrayPenalForce[i][j][k][0]: force_x
  // arrayPenalForce[i][j][k][1]: force_y
  // arrayPenalForce[i][j][k][2]: force_z

  Array4D arrayPressureForce;
  // fluid grid forces, 4-dimensional
  // nx, ny, nz, n_dim
  // arrayPressureForce[i][j][k][0]: force_x
  // arrayPressureForce[i][j][k][1]: force_y
  // arrayPressureForce[i][j][k][2]: force_z

  Array4D arrayFlux;
  // 4-dimensional, defined at cell centers
  // nx, ny, nz, n_integ
  // arrayFlux[i][j][k][0]: var_den
  // arrayFlux[i][j][k][1]: var_mom[0]
  // arrayFlux[i][j][k][2]: var_mom[1]
  // arrayFlux[i][j][k][3]: var_mom[2]
  // arrayFlux[i][j][k][4]: var_eng

  Array5D arrayRoeFlux;
  Array5D arrayRoeFluxStep2;
  Array5D arrayRoeFluxStep3;
  // 5-dimensional, defined at cell faces
  // nx-1, ny-1, nz-1, n_integ, n_dim
  // arrayRoeFlux[i][j][k][0]: var_den
  // arrayRoeFlux[i][j][k][1]: var_mom[0]
  // arrayRoeFlux[i][j][k][2]: var_mom[1]
  // arrayRoeFlux[i][j][k][3]: var_mom[2]
  // arrayRoeFlux[i][j][k][4]: var_eng

  Array4D arrayRoeFluxTmp;
  // 4-dimensional
  // nx-1, ny-1, nz-1, n_integ

  Array3D arrayH;
  // Enthalpy, 3-dimensional, total enthalpy (not static enthalpy)
  // nx, ny, nz

  Array3D arraySoundSpeed;
  // speed of sound, 3-dimensional
  // nx, ny, nz

public:
  Fluid() {}

  void initParameter(Box& container, Gradation& gradation);
  void initialize();
  void initialCondition();
  REAL calcTimeStep();
  void RankineHugoniot();
  void addGhostPoints();
  void soundSpeed();
  void enthalpy();
  void flux();
  void RoeFlux(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR,
               std::size_t idim, std::size_t i, std::size_t j, std::size_t k);
  void UtoW(); // U - integrated; W - primitive
  void WtoU();
  void rotateIJK();
  void inteStep1();
  void inteStep2();
  void inteStep3();

  void getParticleInfo(ParticlePArray& ptcls);
  void runOneStep();
  void calcParticleForce(ParticlePArray& ptcls, std::ofstream& ofs);
  void penalize();
  void plot(const char*) const;
};

} // name space dem
#endif
