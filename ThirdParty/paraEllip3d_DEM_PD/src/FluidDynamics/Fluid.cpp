#include <Core/Const/const.h>
#include <FluidDynamics/Fluid.h>
#include <Core/Util/Utility.h>
#include <cmath>

namespace dem {

constexpr REAL Fluid::Rs;

void
Fluid::initParameter(Box& container, Gradation& gradation)
{

  RK = util::getParam<REAL>("RK");
  CFL = util::getParam<REAL>("CFL");
  gamma = util::getParam<REAL>("airGamma");
  arrayBC[0] = util::getParam<REAL>("x1Reflecting");
  arrayBC[1] = util::getParam<REAL>("x2Reflecting");
  arrayBC[2] = util::getParam<REAL>("y1Reflecting");
  arrayBC[3] = util::getParam<REAL>("y2Reflecting");
  arrayBC[4] = util::getParam<REAL>("z1Reflecting");
  arrayBC[5] = util::getParam<REAL>("z2Reflecting");
  rhoR = util::getParam<REAL>("rightDensity");
  pR = util::getParam<REAL>("rightPressure");
  uR = util::getParam<REAL>("rightVelocity");
  mach = util::getParam<REAL>("MachNumber");
  Cd = util::getParam<REAL>("Cd");
  std::size_t ptclGrid = util::getParam<std::size_t>("ptclGrid");
  volFrac = util::getParam<REAL>("volFrac");

  REAL minX = container.getMinCorner().x();
  REAL minY = container.getMinCorner().y();
  REAL minZ = container.getMinCorner().z();
  REAL z1Distance = util::getParam<REAL>("z1Distance");
  minZ -= z1Distance;
  z0 = minZ + z1Distance * util::getParam<REAL>("z1Percent");

  REAL maxX = container.getMaxCorner().x();
  REAL maxY = container.getMaxCorner().y();
  REAL maxZ = container.getMaxCorner().z();
  REAL minR = gradation.getPtclMinRadius();

  dx = (minR * 2) / ptclGrid;
  dy = dx;
  dz = dx;
  nx = static_cast<std::size_t>(ceil((maxX - minX) / dx));
  ny = static_cast<std::size_t>(ceil((maxY - minY) / dy));
  nz = static_cast<std::size_t>(ceil((maxZ - minZ) / dz));

  dx = (maxX - minX) / nx;
  dy = (maxY - minY) / ny;
  dz = (maxZ - minZ) / nz;

  nx += 2;
  ny += 2;
  nz += 2;

  // fixed
  n_dim = 3;
  n_var = 0;
  n_integ = 0;

  var_den = n_var++;
  n_integ++;
  var_mom[0] = n_var++;
  n_integ++;
  var_mom[1] = n_var++;
  n_integ++;
  var_mom[2] = n_var++;
  n_integ++;
  var_eng = n_var++;
  n_integ++;

  var_vel[0] = n_var++;
  var_vel[1] = n_var++;
  var_vel[2] = n_var++;
  var_prs = n_var++;

  // extended
  var_msk = n_var++;

  ///*
  debugInf << std::setw(OWID) << "Runge-Kutta" << std::setw(OWID) << RK
           << std::endl;
  debugInf << std::setw(OWID) << "CFL" << std::setw(OWID) << CFL << std::endl;
  debugInf << std::setw(OWID) << "gamma" << std::setw(OWID) << gamma
           << std::endl;
  debugInf << std::setw(OWID) << "x1Rflecting" << std::setw(OWID) << arrayBC[0]
           << std::endl;
  debugInf << std::setw(OWID) << "x2Rflecting" << std::setw(OWID) << arrayBC[1]
           << std::endl;
  debugInf << std::setw(OWID) << "y1Rflecting" << std::setw(OWID) << arrayBC[2]
           << std::endl;
  debugInf << std::setw(OWID) << "y2Rflecting" << std::setw(OWID) << arrayBC[3]
           << std::endl;
  debugInf << std::setw(OWID) << "z1Rflecting" << std::setw(OWID) << arrayBC[4]
           << std::endl;
  debugInf << std::setw(OWID) << "z2Rflecting" << std::setw(OWID) << arrayBC[5]
           << std::endl;
  debugInf << std::setw(OWID) << "z0" << std::setw(OWID) << z0 << std::endl;
  debugInf << std::setw(OWID) << "rhoR" << std::setw(OWID) << rhoR << std::endl;
  debugInf << std::setw(OWID) << "pR" << std::setw(OWID) << pR << std::endl;
  debugInf << std::setw(OWID) << "uR" << std::setw(OWID) << uR << std::endl;
  debugInf << std::setw(OWID) << "Mach" << std::setw(OWID) << mach << std::endl;
  debugInf << std::setw(OWID) << "Cd" << std::setw(OWID) << Cd << std::endl;
  debugInf << std::setw(OWID) << "ptclGrid" << std::setw(OWID) << ptclGrid
           << std::endl;
  debugInf << std::setw(OWID) << "gridSize" << std::setw(OWID) << dx
           << std::endl;
  debugInf << std::setw(OWID) << "volFrac" << std::setw(OWID) << volFrac
           << std::endl;

  /*
  debugInf << "n_var " << n_var << std::endl;
  debugInf << "n_integ " << n_integ << std::endl;
  debugInf << "var_den " << var_den  << std::endl;
  debugInf << "var_mom[0] " << var_mom[0] << std::endl;
  debugInf << "var_mom[1] " << var_mom[1] << std::endl;
  debugInf << "var_mom[2] " << var_mom[2] << std::endl;
  debugInf << "var_eng " << var_eng  << std::endl;
  debugInf << "var_vel[0] " << var_vel[0] << std::endl;
  debugInf << "var_vel[1] " << var_vel[1] << std::endl;
  debugInf << "var_vel[2] " << var_vel[2] << std::endl;
  debugInf << "var_prs " << var_prs  << std::endl;
  debugInf << "var_msk " << var_msk  << std::endl;
  */

  // nx, ny, nz, n_dim
  arrayGridCoord.resize(nx);
  for (std::size_t i = 0; i < arrayGridCoord.size(); ++i) {
    arrayGridCoord[i].resize(ny);
    for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j) {
      arrayGridCoord[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k)
        arrayGridCoord[i][j][k].resize(n_dim);
    }
  }

  // coordinates
  for (std::size_t i = 0; i < arrayGridCoord.size(); ++i)
    for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
      for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {
        arrayGridCoord[i][j][k][0] = (minX - dx / 2) + i * dx;
        arrayGridCoord[i][j][k][1] = (minY - dy / 2) + j * dy;
        arrayGridCoord[i][j][k][2] = (minZ - dz / 2) + k * dz;
      }

  // nx, ny, nz, n_dim
  arrayPenalForce.resize(nx);
  for (std::size_t i = 0; i < arrayPenalForce.size(); ++i) {
    arrayPenalForce[i].resize(ny);
    for (std::size_t j = 0; j < arrayPenalForce[i].size(); ++j) {
      arrayPenalForce[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayPenalForce[i][j].size(); ++k)
        arrayPenalForce[i][j][k].resize(n_dim);
    }
  }

  // nx, ny, nz, n_dim
  arrayPressureForce.resize(nx);
  for (std::size_t i = 0; i < arrayPressureForce.size(); ++i) {
    arrayPressureForce[i].resize(ny);
    for (std::size_t j = 0; j < arrayPressureForce[i].size(); ++j) {
      arrayPressureForce[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayPressureForce[i][j].size(); ++k)
        arrayPressureForce[i][j][k].resize(n_dim);
    }
  }

  // nx, ny, nz, n_var
  arrayU.resize(nx);
  for (std::size_t i = 0; i < arrayU.size(); ++i) {
    arrayU[i].resize(ny);
    for (std::size_t j = 0; j < arrayU[i].size(); ++j) {
      arrayU[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayU[i][j].size(); ++k)
        arrayU[i][j][k].resize(n_var);
    }
  }

  // nx, ny, nz, n_var
  arrayUtmp.resize(nx);
  for (std::size_t i = 0; i < arrayUtmp.size(); ++i) {
    arrayUtmp[i].resize(ny);
    for (std::size_t j = 0; j < arrayUtmp[i].size(); ++j) {
      arrayUtmp[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayUtmp[i][j].size(); ++k)
        arrayUtmp[i][j][k].resize(n_var);
    }
  }

  // nx, ny, nz, n_integ
  arrayFlux.resize(nx);
  for (std::size_t i = 0; i < arrayFlux.size(); ++i) {
    arrayFlux[i].resize(ny);
    for (std::size_t j = 0; j < arrayFlux[i].size(); ++j) {
      arrayFlux[i][j].resize(nz);
      for (std::size_t k = 0; k < arrayFlux[i][j].size(); ++k)
        arrayFlux[i][j][k].resize(n_integ);
    }
  }

  // nx-1, ny-1, nz-1, n_integ, n_dim
  arrayRoeFlux.resize(nx - 1);
  for (std::size_t i = 0; i < arrayRoeFlux.size(); ++i) {
    arrayRoeFlux[i].resize(ny - 1);
    for (std::size_t j = 0; j < arrayRoeFlux[i].size(); ++j) {
      arrayRoeFlux[i][j].resize(nz - 1);
      for (std::size_t k = 0; k < arrayRoeFlux[i][j].size(); ++k) {
        arrayRoeFlux[i][j][k].resize(n_integ);
        for (std::size_t m = 0; m < arrayRoeFlux[i][j][k].size(); ++m)
          arrayRoeFlux[i][j][k][m].resize(n_dim);
      }
    }
  }

  if (RK >= 1) {
    // nx-1, ny-1, nz-1, n_integ, n_dim
    arrayRoeFluxStep2.resize(nx - 1);
    for (std::size_t i = 0; i < arrayRoeFluxStep2.size(); ++i) {
      arrayRoeFluxStep2[i].resize(ny - 1);
      for (std::size_t j = 0; j < arrayRoeFluxStep2[i].size(); ++j) {
        arrayRoeFluxStep2[i][j].resize(nz - 1);
        for (std::size_t k = 0; k < arrayRoeFluxStep2[i][j].size(); ++k) {
          arrayRoeFluxStep2[i][j][k].resize(n_integ);
          for (std::size_t m = 0; m < arrayRoeFluxStep2[i][j][k].size(); ++m)
            arrayRoeFluxStep2[i][j][k][m].resize(n_dim);
        }
      }
    }
  }

  if (RK == 2) {
    // nx-1, ny-1, nz-1, n_integ, n_dim
    arrayRoeFluxStep3.resize(nx - 1);
    for (std::size_t i = 0; i < arrayRoeFluxStep3.size(); ++i) {
      arrayRoeFluxStep3[i].resize(ny - 1);
      for (std::size_t j = 0; j < arrayRoeFluxStep3[i].size(); ++j) {
        arrayRoeFluxStep3[i][j].resize(nz - 1);
        for (std::size_t k = 0; k < arrayRoeFluxStep3[i][j].size(); ++k) {
          arrayRoeFluxStep3[i][j][k].resize(n_integ);
          for (std::size_t m = 0; m < arrayRoeFluxStep3[i][j][k].size(); ++m)
            arrayRoeFluxStep3[i][j][k][m].resize(n_dim);
        }
      }
    }
  }

  // nx-1, ny-1, nz-1, n_integ
  arrayRoeFluxTmp.resize(nx - 1);
  for (std::size_t i = 0; i < arrayRoeFluxTmp.size(); ++i) {
    arrayRoeFluxTmp[i].resize(ny - 1);
    for (std::size_t j = 0; j < arrayRoeFluxTmp[i].size(); ++j) {
      arrayRoeFluxTmp[i][j].resize(nz - 1);
      for (std::size_t k = 0; k < arrayRoeFluxTmp[i][j].size(); ++k) {
        arrayRoeFluxTmp[i][j][k].resize(n_integ);
      }
    }
  }

  // nx, ny, nz
  arrayH.resize(nx);
  for (std::size_t i = 0; i < arrayH.size(); ++i) {
    arrayH[i].resize(ny);
    for (std::size_t j = 0; j < arrayH[i].size(); ++j)
      arrayH[i][j].resize(nz);
  }

  // nx, ny, nz
  arraySoundSpeed.resize(nx);
  for (std::size_t i = 0; i < arraySoundSpeed.size(); ++i) {
    arraySoundSpeed[i].resize(ny);
    for (std::size_t j = 0; j < arraySoundSpeed[i].size(); ++j)
      arraySoundSpeed[i][j].resize(nz);
  }
}

void
Fluid::initialize()
{
  RankineHugoniot();
  initialCondition();
  soundSpeed(); // for printing mach
}

void
Fluid::runOneStep()
{
  inteStep1();

  if (RK >= 1) {
    arrayRoeFluxStep2 = arrayRoeFlux;
    inteStep2();
    if (RK == 2) {
      arrayRoeFluxStep3 = arrayRoeFlux;
      inteStep3();
    }
  }
}

void
Fluid::inteStep1()
{
  addGhostPoints();
  soundSpeed();
  timeStep = std::min(timeStep, calcTimeStep());
  debugInf << std::setw(OWID) << iteration << std::setw(OWID) << timeStep
           << std::setw(OWID) << timeAccrued << std::endl;
  enthalpy();
  rotateIJK();

  // update conservative variables at the next time step
  for (std::size_t i = 1; i < nx - 1; ++i)
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t k = 1; k < nz - 1; ++k) {
        for (std::size_t m = 0; m < n_integ; ++m)
          arrayU[i][j][k][m] -=
            (timeStep / dx *
               (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i - 1][j][k][m][0]) +
             timeStep / dy *
               (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j - 1][k][m][1]) +
             timeStep / dz *
               (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k - 1][m][2]));
      }

  // calculate primitive after finding conservative variables
  UtoW();
}

void
Fluid::inteStep2()
{
  addGhostPoints();
  soundSpeed();
  enthalpy();
  rotateIJK();

  // update conservative variables at the next time step
  for (std::size_t i = 1; i < nx - 1; ++i)
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t k = 1; k < nz - 1; ++k) {
        for (std::size_t m = 0; m < n_integ; ++m)
          arrayU[i][j][k][m] -=
            (timeStep / (2 * RK * dx) *
               (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i - 1][j][k][m][0] +
                (arrayRoeFluxStep2[i][j][k][m][0] -
                 arrayRoeFluxStep2[i - 1][j][k][m][0])) +
             timeStep / (2 * RK * dy) *
               (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j - 1][k][m][1] +
                (arrayRoeFluxStep2[i][j][k][m][1] -
                 arrayRoeFluxStep2[i][j - 1][k][m][1])) +
             timeStep / (2 * RK * dz) *
               (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k - 1][m][2] +
                (arrayRoeFluxStep2[i][j][k][m][2] -
                 arrayRoeFluxStep2[i][j][k - 1][m][2])));
      }

  // calculate primitive after finding conservative variables
  UtoW();
}

void
Fluid::inteStep3()
{
  addGhostPoints();
  soundSpeed();
  enthalpy();
  rotateIJK();

  // update conservative variables at the next time step
  for (std::size_t i = 1; i < nx - 1; ++i)
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t k = 1; k < nz - 1; ++k) {
        for (std::size_t m = 0; m < n_integ; ++m)
          arrayU[i][j][k][m] -=
            (timeStep / (6 * dx) *
               (arrayRoeFlux[i][j][k][m][0] - arrayRoeFlux[i - 1][j][k][m][0] +
                (arrayRoeFluxStep2[i][j][k][m][0] -
                 arrayRoeFluxStep2[i - 1][j][k][m][0]) +
                4 * (arrayRoeFluxStep3[i][j][k][m][0] -
                     arrayRoeFluxStep3[i - 1][j][k][m][0])) +
             timeStep / (6 * dy) *
               (arrayRoeFlux[i][j][k][m][1] - arrayRoeFlux[i][j - 1][k][m][1] +
                (arrayRoeFluxStep2[i][j][k][m][1] -
                 arrayRoeFluxStep2[i][j - 1][k][m][1]) +
                4 * (arrayRoeFluxStep3[i][j][k][m][1] -
                     arrayRoeFluxStep3[i][j - 1][k][m][1])) +
             timeStep / (6 * dz) *
               (arrayRoeFlux[i][j][k][m][2] - arrayRoeFlux[i][j][k - 1][m][2] +
                (arrayRoeFluxStep2[i][j][k][m][2] -
                 arrayRoeFluxStep2[i][j][k - 1][m][2]) +
                4 * (arrayRoeFluxStep3[i][j][k][m][2] -
                     arrayRoeFluxStep3[i][j][k - 1][m][2])));
      }

  // calculate primitive after finding conservative variables
  UtoW();
}

void
Fluid::rotateIJK()
{
  std::size_t id[3][3] = { { 0, 1, 2 }, { 1, 0, 2 }, { 2, 1, 0 } };

  // for x, y, z directions
  for (std::size_t idim = 0; idim < n_dim; ++idim) {
    arrayUtmp = arrayU; // must have same rank and extent

    // switch components
    for (std::size_t i = 0; i < nx; ++i)
      for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t k = 0; k < nz; ++k)
          for (std::size_t jdim = 0; jdim < n_dim; ++jdim) {
            arrayUtmp[i][j][k][var_mom[jdim]] =
              arrayU[i][j][k][var_mom[id[idim][jdim]]];
            arrayUtmp[i][j][k][var_vel[jdim]] =
              arrayU[i][j][k][var_vel[id[idim][jdim]]];
          }

    flux();

    // for local Riemann problem
    for (std::size_t i = 0; i < nx - 1; ++i) {
      for (std::size_t j = 0; j < ny - 1; ++j) {
        for (std::size_t k = 0; k < nz - 1; ++k) {
          std::size_t IL[3] = { i, j, k };
          std::size_t IR[3] = { i, j, k };
          IR[idim] += 1;
          REAL uL[9], uR[9], FL[5], FR[5], HL, HR; // local variable only
          HL = arrayH[IL[0]][IL[1]][IL[2]];
          HR = arrayH[IR[0]][IR[1]][IR[2]];
          for (std::size_t m = 0; m < n_var; ++m) {
            uL[m] = arrayUtmp[IL[0]][IL[1]][IL[2]][m];
            uR[m] = arrayUtmp[IR[0]][IR[1]][IR[2]][m];
          }
          for (std::size_t m = 0; m < n_integ; ++m) {
            FL[m] = arrayFlux[IL[0]][IL[1]][IL[2]][m];
            FR[m] = arrayFlux[IR[0]][IR[1]][IR[2]][m];
          }
          RoeFlux(uL, uR, FL, FR, HL, HR, idim, i, j, k);
        }
      }
    }

    for (std::size_t i = 0; i < nx - 1; ++i)
      for (std::size_t j = 0; j < ny - 1; ++j)
        for (std::size_t k = 0; k < nz - 1; ++k)
          for (std::size_t m = 0; m < n_integ; ++m)
            arrayRoeFluxTmp[i][j][k][m] = arrayRoeFlux[i][j][k][m][idim];

    // switch components back for consistency with u
    for (std::size_t i = 0; i < nx - 1; ++i)
      for (std::size_t j = 0; j < ny - 1; ++j)
        for (std::size_t k = 0; k < nz - 1; ++k)
          for (std::size_t m = 0; m < n_dim; ++m)
            arrayRoeFlux[i][j][k][var_mom[m]][idim] =
              arrayRoeFluxTmp[i][j][k][var_mom[id[idim][m]]];

  } // end of for x, y, z directions
}

void
Fluid::penalize()
{
  // Brinkman penalization
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t m = 0; m < n_dim; ++m) {
          arrayU[i][j][k][var_mom[m]] -=
            arrayU[i][j][k][var_msk] * arrayPenalForce[i][j][k][m] * timeStep;
          arrayU[i][j][k][var_eng] -= arrayU[i][j][k][var_msk] *
                                      arrayPenalForce[i][j][k][m] *
                                      arrayU[i][j][k][var_vel[m]] * timeStep;
        }
}

void
Fluid::addGhostPoints()
{
  // non-reflecting BCs
  for (std::size_t j = 1; j < ny - 1; ++j)
    for (std::size_t k = 1; k < nz - 1; ++k)
      for (std::size_t m = 0; m < n_var; ++m) {
        arrayU[0][j][k][m] = arrayU[1][j][k][m];
        arrayU[nx - 1][j][k][m] = arrayU[nx - 2][j][k][m];
      }

  for (std::size_t i = 1; i < nx - 1; ++i)
    for (std::size_t k = 1; k < nz - 1; ++k)
      for (std::size_t m = 0; m < n_var; ++m) {
        arrayU[i][0][k][m] = arrayU[i][1][k][m];
        arrayU[i][ny - 1][k][m] = arrayU[i][ny - 2][k][m];
      }

  for (std::size_t i = 1; i < nx - 1; ++i)
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t m = 0; m < n_var; ++m) {
        arrayU[i][j][0][m] = arrayU[i][j][1][m];
        arrayU[i][j][nz - 1][m] = arrayU[i][j][nz - 2][m];
      }

  // reflecting BCs
  bool reflecting = false;
  for (double it : arrayBC) {
    if (it > 0) {
      reflecting = true;
      break;
    }
  }

  if (reflecting) {
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t k = 1; k < nz - 1; ++k)
        for (std::size_t m = 0; m < 1; ++m) { // x-direction
          arrayU[0][j][k][var_mom[m]] *= (1 - 2 * arrayBC[0]);
          arrayU[nx - 1][j][k][var_mom[m]] *= (1 - 2 * arrayBC[1]);
          arrayU[0][j][k][var_vel[m]] *= (1 - 2 * arrayBC[0]);
          arrayU[nx - 1][j][k][var_vel[m]] *= (1 - 2 * arrayBC[1]);
        }

    for (std::size_t i = 1; i < nx - 1; ++i)
      for (std::size_t k = 1; k < nz - 1; ++k)
        for (std::size_t m = 1; m < 2; ++m) { // y-direction
          arrayU[i][0][k][var_mom[m]] *= (1 - 2 * arrayBC[2]);
          arrayU[i][ny - 1][k][var_mom[m]] *= (1 - 2 * arrayBC[3]);
          arrayU[i][0][k][var_vel[m]] *= (1 - 2 * arrayBC[2]);
          arrayU[i][ny - 1][k][var_vel[m]] *= (1 - 2 * arrayBC[3]);
        }

    for (std::size_t i = 1; i < nx - 1; ++i)
      for (std::size_t j = 1; j < ny - 1; ++j)
        for (std::size_t m = 2; m < 3; ++m) { // z-direction
          arrayU[i][j][0][var_mom[m]] *= (1 - 2 * arrayBC[4]);
          arrayU[i][j][nz - 1][var_mom[m]] *= (1 - 2 * arrayBC[5]);
          arrayU[i][j][0][var_vel[m]] *= (1 - 2 * arrayBC[4]);
          arrayU[i][j][nz - 1][var_vel[m]] *= (1 - 2 * arrayBC[5]);
        }
  }
}

REAL
Fluid::calcTimeStep()
{
  std::valarray<REAL> allGrid(nx * ny * nz);
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        allGrid[i + j * nx + k * nx * ny] =
          fabs(arrayU[i][j][k][var_vel[0]]) + arraySoundSpeed[i][j][k];
  REAL sx = allGrid.max();

  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        allGrid[i + j * nx + k * nx * ny] =
          fabs(arrayU[i][j][k][var_vel[1]]) + arraySoundSpeed[i][j][k];
  REAL sy = allGrid.max();

  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        allGrid[i + j * nx + k * nx * ny] =
          fabs(arrayU[i][j][k][var_vel[2]]) + arraySoundSpeed[i][j][k];
  REAL sz = allGrid.max();

  std::valarray<REAL> dtMin(3);
  dtMin[0] = dx / sx;
  dtMin[1] = dy / sy;
  dtMin[2] = dz / sz;

  return CFL * dtMin.min();
}

void
Fluid::soundSpeed()
{
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        arraySoundSpeed[i][j][k] =
          sqrt(gamma * arrayU[i][j][k][var_prs] / arrayU[i][j][k][var_den]);
}

void
Fluid::enthalpy()
{
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        arrayH[i][j][k] =
          (arrayU[i][j][k][var_eng] + arrayU[i][j][k][var_prs]) /
          arrayU[i][j][k][var_den];
}

void
Fluid::initialCondition()
{
  // normal shock
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k) {
        if (arrayGridCoord[i][j][k][2] <= z0) {
          arrayU[i][j][k][var_den] = rhoL;
          arrayU[i][j][k][var_vel[2]] = uL;
          arrayU[i][j][k][var_prs] = pL;
        } else {
          arrayU[i][j][k][var_den] = rhoR;
          arrayU[i][j][k][var_vel[2]] = uR;
          arrayU[i][j][k][var_prs] = pR;
        }
      }

  WtoU();
}

void
Fluid::RankineHugoniot()
{
  shockSpeed = mach * sqrt(gamma * pR / rhoR);
  pL = (pR * (1 - gamma) + 2 * rhoR * pow(shockSpeed - uR, 2)) / (1 + gamma);
  rhoL = (pow(rhoR * (shockSpeed - uR), 2) * (1 + gamma)) /
         (rhoR * pow(shockSpeed - uR, 2) * (gamma - 1) + 2 * pR * gamma);
  uL = (rhoR * (shockSpeed - uR) * (2 * shockSpeed + uR * (gamma - 1)) -
        2 * pR * gamma) /
       (rhoR * (shockSpeed - uR) * (1 + gamma));
  ///*
  debugInf << std::setw(OWID) << "rhoL" << std::setw(OWID) << rhoL << std::endl;
  debugInf << std::setw(OWID) << "uL" << std::setw(OWID) << uL << std::endl;
  debugInf << std::setw(OWID) << "pL" << std::setw(OWID) << pL << std::endl;
  debugInf << std::setw(OWID) << "shockSpeed" << std::setw(OWID) << shockSpeed
           << std::endl
           << std::endl;
  //*/
}

void
Fluid::flux()
{
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k) {
        arrayFlux[i][j][k][var_den] =
          arrayUtmp[i][j][k][var_den] * arrayUtmp[i][j][k][var_vel[0]]; // rho*u
        arrayFlux[i][j][k][var_mom[0]] =
          arrayUtmp[i][j][k][var_den] * pow(arrayUtmp[i][j][k][var_vel[0]], 2) +
          arrayUtmp[i][j][k][var_prs]; // rho*u^2 + p
        arrayFlux[i][j][k][var_mom[1]] =
          arrayUtmp[i][j][k][var_den] * arrayUtmp[i][j][k][var_vel[0]] *
          arrayUtmp[i][j][k][var_vel[1]]; // rho*u*v
        arrayFlux[i][j][k][var_mom[2]] =
          arrayUtmp[i][j][k][var_den] * arrayUtmp[i][j][k][var_vel[0]] *
          arrayUtmp[i][j][k][var_vel[2]]; // rho*u*w
        arrayFlux[i][j][k][var_eng] =
          arrayUtmp[i][j][k][var_vel[0]] *
          (arrayUtmp[i][j][k][var_eng] +
           arrayUtmp[i][j][k][var_prs]); // u*(E + p)
      }
}

void
Fluid::RoeFlux(REAL uL[], REAL uR[], REAL FL[], REAL FR[], REAL HL, REAL HR,
               std::size_t idim, std::size_t it, std::size_t jt, std::size_t kt)
{
  REAL avgRho = sqrt(uL[var_den] * uR[var_den]);
  REAL avgH = (sqrt(uL[var_den]) * HL + sqrt(uR[var_den]) * HR) /
              (sqrt(uL[var_den]) + sqrt(uR[var_den]));
  REAL avgU =
    (sqrt(uL[var_den]) * uL[var_vel[0]] + sqrt(uR[var_den]) * uR[var_vel[0]]) /
    (sqrt(uL[var_den]) + sqrt(uR[var_den]));
  REAL avgV =
    (sqrt(uL[var_den]) * uL[var_vel[1]] + sqrt(uR[var_den]) * uR[var_vel[1]]) /
    (sqrt(uL[var_den]) + sqrt(uR[var_den]));
  REAL avgW =
    (sqrt(uL[var_den]) * uL[var_vel[2]] + sqrt(uR[var_den]) * uR[var_vel[2]]) /
    (sqrt(uL[var_den]) + sqrt(uR[var_den]));
  REAL avgSoundSpeed = sqrt(
    (gamma - 1) * (avgH - 0.5 * (avgU * avgU + avgV * avgV + avgW * avgW)));
  if ((avgH - 0.5 * (avgU * avgU + avgV * avgV + avgW * avgW)) < 0)
    debugInf << std::setw(3) << it << std::setw(3) << jt << std::setw(3) << kt
             << std::setw(OWID) << uL[var_den] << std::setw(OWID) << uR[var_den]
             << std::setw(OWID) << HL << std::setw(OWID) << HR
             << std::setw(OWID) << uL[var_vel[0]] << std::setw(OWID)
             << uR[var_vel[0]] << std::setw(OWID) << uL[var_vel[1]]
             << std::setw(OWID) << uR[var_vel[1]] << std::setw(OWID)
             << uL[var_vel[2]] << std::setw(OWID) << uR[var_vel[2]]
             << std::setw(OWID) << avgH << std::setw(OWID) << avgU
             << std::setw(OWID) << avgV << std::setw(OWID) << avgW
             << std::setw(OWID) << avgSoundSpeed << std::endl;

  REAL eigen[5];
  eigen[var_den] = avgU - avgSoundSpeed;
  eigen[var_mom[0]] = eigen[var_mom[1]] = eigen[var_mom[2]] = avgU;
  eigen[var_eng] = avgU + avgSoundSpeed;

  REAL avgWaveStr[5], du[9];
  for (std::size_t i = 0; i < n_var; ++i)
    du[i] = uR[i] - uL[i];

  avgWaveStr[var_den] =
    (du[var_prs] - avgRho * avgSoundSpeed * du[var_vel[0]]) /
    (2 * avgSoundSpeed * avgSoundSpeed);
  avgWaveStr[var_mom[0]] =
    du[var_den] - du[var_prs] / (avgSoundSpeed * avgSoundSpeed);
  avgWaveStr[var_mom[1]] = avgRho * du[var_vel[1]];
  avgWaveStr[var_mom[2]] = avgRho * du[var_vel[2]];
  avgWaveStr[var_eng] =
    (du[var_prs] + avgRho * avgSoundSpeed * du[var_vel[0]]) /
    (2 * avgSoundSpeed * avgSoundSpeed);

  REAL avgK[5][5]; // right eigenvectors
  avgK[var_den][var_den] = 1;
  avgK[var_mom[0]][var_den] = avgU - avgSoundSpeed;
  avgK[var_mom[1]][var_den] = avgV;
  avgK[var_mom[2]][var_den] = avgW;
  avgK[var_eng][var_den] = avgH - avgU * avgSoundSpeed;

  avgK[var_den][var_mom[0]] = 1;
  avgK[var_mom[0]][var_mom[0]] = avgU;
  avgK[var_mom[1]][var_mom[0]] = avgV;
  avgK[var_mom[2]][var_mom[0]] = avgW;
  avgK[var_eng][var_mom[0]] = 0.5 * (avgU * avgU + avgV * avgV + avgW * avgW);

  avgK[var_den][var_mom[1]] = 0;
  avgK[var_mom[0]][var_mom[1]] = 0;
  avgK[var_mom[1]][var_mom[1]] = 1;
  avgK[var_mom[2]][var_mom[1]] = 0;
  avgK[var_eng][var_mom[1]] = avgV;

  avgK[var_den][var_mom[2]] = 0;
  avgK[var_mom[0]][var_mom[2]] = 0;
  avgK[var_mom[1]][var_mom[2]] = 0;
  avgK[var_mom[2]][var_mom[2]] = 1;
  avgK[var_eng][var_mom[2]] = avgW;

  avgK[var_den][var_eng] = 1;
  avgK[var_mom[0]][var_eng] = avgU + avgSoundSpeed;
  avgK[var_mom[1]][var_eng] = avgV;
  avgK[var_mom[2]][var_eng] = avgW;
  avgK[var_eng][var_eng] = avgH + avgU * avgSoundSpeed;

  REAL RF[5];
  for (std::size_t i = 0; i < n_integ; ++i)
    RF[i] = 0.5 * (FL[i] + FR[i]);

  for (std::size_t ie = 0; ie < n_integ; ++ie) {
    for (std::size_t je = 0; je < n_integ; ++je)
      RF[ie] -= 0.5 * avgWaveStr[je] * fabs(eigen[je]) * avgK[ie][je];
    arrayRoeFlux[it][jt][kt][ie][idim] = RF[ie];
  }
}

void
Fluid::UtoW()
{ // converting conservative variables into primitive
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k) {

        for (std::size_t m = 0; m < n_dim; ++m)
          arrayU[i][j][k][var_vel[m]] =
            arrayU[i][j][k][var_mom[m]] / arrayU[i][j][k][var_den];

        arrayU[i][j][k][var_prs] = 0;
        for (std::size_t m = 0; m < n_dim; ++m)
          arrayU[i][j][k][var_prs] +=
            arrayU[i][j][k][var_den] * pow(arrayU[i][j][k][var_vel[m]], 2) / 2;
        arrayU[i][j][k][var_prs] =
          (arrayU[i][j][k][var_eng] - arrayU[i][j][k][var_prs]) * (gamma - 1);
      }
}

void
Fluid::WtoU()
{ // converting primitive variables into conservative
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t m = 0; m < n_dim; ++m)
          arrayU[i][j][k][var_mom[m]] =
            arrayU[i][j][k][var_den] * arrayU[i][j][k][var_vel[m]];

        arrayU[i][j][k][var_eng] = 0;
        for (std::size_t m = 0; m < n_dim; ++m)
          arrayU[i][j][k][var_eng] +=
            arrayU[i][j][k][var_den] * pow(arrayU[i][j][k][var_vel[m]], 2) / 2;
        arrayU[i][j][k][var_eng] += arrayU[i][j][k][var_prs] / (gamma - 1);
      }
}

void
Fluid::getParticleInfo(ParticlePArray& ptcls)
{
  for (const auto& ptcl : ptcls)
    ptcl->clearFluidGrid();

  // 0 ~ (n-1), including boundaries
  for (std::size_t i = 0; i < arrayGridCoord.size(); ++i)
    for (std::size_t j = 0; j < arrayGridCoord[i].size(); ++j)
      for (std::size_t k = 0; k < arrayGridCoord[i][j].size(); ++k) {

        arrayU[i][j][k][var_msk] = 0;
        REAL coord_x = arrayGridCoord[i][j][k][0];
        REAL coord_y = arrayGridCoord[i][j][k][1];
        REAL coord_z = arrayGridCoord[i][j][k][2];

        if (volFrac == 0) {

          for (auto& ptcl : ptcls)
            if (ptcl->surfaceError(Vec(coord_x, coord_y, coord_z)) <=
                0) { // inside particle surface
              arrayU[i][j][k][var_msk] = 1;
              ptcl->recordFluidGrid(i, j, k, 1.0);
            }

        } else if (volFrac == 1) {

          for (auto& ptcl : ptcls) {
            bool in[8];
            in[0] = ptcl->surfaceError(Vec(coord_x - dx / 2, coord_y - dy / 2,
                                           coord_z - dz / 2)) < 0;
            in[1] = ptcl->surfaceError(Vec(coord_x + dx / 2, coord_y - dy / 2,
                                           coord_z - dz / 2)) < 0;
            in[2] = ptcl->surfaceError(Vec(coord_x - dx / 2, coord_y + dy / 2,
                                           coord_z - dz / 2)) < 0;
            in[3] = ptcl->surfaceError(Vec(coord_x + dx / 2, coord_y + dy / 2,
                                           coord_z - dz / 2)) < 0;
            in[4] = ptcl->surfaceError(Vec(coord_x - dx / 2, coord_y - dy / 2,
                                           coord_z + dz / 2)) < 0;
            in[5] = ptcl->surfaceError(Vec(coord_x + dx / 2, coord_y - dy / 2,
                                           coord_z + dz / 2)) < 0;
            in[6] = ptcl->surfaceError(Vec(coord_x - dx / 2, coord_y + dy / 2,
                                           coord_z + dz / 2)) < 0;
            in[7] = ptcl->surfaceError(Vec(coord_x + dx / 2, coord_y + dy / 2,
                                           coord_z + dz / 2)) < 0;

            if (in[0] || in[1] || in[2] || in[3] || in[4] || in[5] || in[6] ||
                in[7]) { // if any vertex is inside particle surface
              arrayU[i][j][k][var_msk] = 1;

              REAL volFraction = 1;
              if (!(in[0] && in[1] && in[2] && in[3] && in[4] && in[5] &&
                    in[6] &&
                    in[7])) { // if any vertex is outside particle surface
                std::size_t fineGrid = 5;
                std::size_t fineCount = 0;
                for (std::size_t vi = 0; vi < fineGrid; ++vi)
                  for (std::size_t vj = 0; vj < fineGrid; ++vj)
                    for (std::size_t vk = 0; vk < fineGrid; ++vk) {
                      if (ptcl->surfaceError(Vec(
                            coord_x - dx / 2 + (0.5 + vi) * dx / fineGrid,
                            coord_y - dy / 2 + (0.5 + vj) * dy / fineGrid,
                            coord_z - dz / 2 + (0.5 + vk) * dz / fineGrid)) <=
                          0)
                        ++fineCount;
                    }
                volFraction = fineCount / pow(fineGrid, 3);
              }

              ptcl->recordFluidGrid(
                i, j, k, volFraction); // no for break, as multiple particles
                                       // could intrude into the same grid
            }
          }

        } // end of if else
      }
}

void
Fluid::calcParticleForce(ParticlePArray& ptcls, std::ofstream& ofs)
{
  // must clear forces each loop, otherwise Fluid::plot prints wrong values;
  // but Fluid::penalize works OK since it uses masks.
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t m = 0; m < n_dim; ++m) {
          arrayPenalForce[i][j][k][m] = 0;
          arrayPressureForce[i][j][k][m] = 0;
        }

  for (const auto& ptcl : ptcls) {
    REAL etaBx = 8.0 / 3.0 * ptcl->getA() / Cd; // local direction x (i.e. a)
    REAL etaBy = 8.0 / 3.0 * ptcl->getB() / Cd; // local direction y (i.e. b)
    REAL etaBz = 8.0 / 3.0 * ptcl->getC() / Cd; // local direction z (i.e. c)

    Vec penalForce = 0, presForce = 0;
    Vec penalMoment = 0, presMoment = 0;
    std::vector<std::vector<REAL>> fluidGrid = ptcl->getFluidGrid();
    for (auto& iter : fluidGrid) {
      std::size_t i = static_cast<std::size_t>(iter[0]);
      std::size_t j = static_cast<std::size_t>(iter[1]);
      std::size_t k = static_cast<std::size_t>(iter[2]);
      REAL volFraction = iter[3];

      REAL coord_x = arrayGridCoord[i][j][k][0];
      REAL coord_y = arrayGridCoord[i][j][k][1];
      REAL coord_z = arrayGridCoord[i][j][k][2];

      REAL uxFluid = arrayU[i][j][k][var_vel[0]];
      REAL uyFluid = arrayU[i][j][k][var_vel[1]];
      REAL uzFluid = arrayU[i][j][k][var_vel[2]];

      Vec dist = Vec(coord_x, coord_y, coord_z) - ptcl->currentPosition();
      // w X r = omga % dist, where % is overloaded as cross product
      Vec omgar = cross(ptcl->currentOmega(), dist); 

      REAL ux = ptcl->currentVel().x() + omgar.x();
      REAL uy = ptcl->currentVel().y() + omgar.y();
      REAL uz = ptcl->currentVel().z() + omgar.z();

      // principal axis decomposition
      Vec globalDelta = Vec(fabs(uxFluid - ux) * (uxFluid - ux),
                            fabs(uyFluid - uy) * (uyFluid - uy),
                            fabs(uzFluid - uz) * (uzFluid - uz));
      Vec localDelta = ptcl->globalToLocal(globalDelta);
      Vec localPenal, globalPenal;
      // localDelta needs to project in local frame in order to calculate local
      // penalization forces
      localPenal.setX(arrayU[i][j][k][var_den] * localDelta.x() / etaBx);
      localPenal.setY(arrayU[i][j][k][var_den] * localDelta.y() / etaBy);
      localPenal.setZ(arrayU[i][j][k][var_den] * localDelta.z() / etaBz);
      globalPenal = ptcl->localToGlobal(localPenal) * volFraction;
      // one grid could have multiple particles intruded, +=, not =
      arrayPenalForce[i][j][k][0] += globalPenal.x();
      arrayPenalForce[i][j][k][1] += globalPenal.y();
      arrayPenalForce[i][j][k][2] += globalPenal.z();

      // restrict pressure gradient grids
      if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1 && k > 0 &&
          k < nz - 1) { // do not use (i-1) for std::size_t because (i-1) is
                        // postive when i=0
        arrayPressureForce[i][j][k][0] =
          (-(arrayU[i + 1][j][k][var_prs] - arrayU[i - 1][j][k][var_prs]) /
           (2 * dx)) *
          volFraction;
        arrayPressureForce[i][j][k][1] =
          (-(arrayU[i][j + 1][k][var_prs] - arrayU[i][j - 1][k][var_prs]) /
           (2 * dy)) *
          volFraction;
        arrayPressureForce[i][j][k][2] =
          (-(arrayU[i][j][k + 1][var_prs] - arrayU[i][j][k - 1][var_prs]) /
           (2 * dz)) *
          volFraction;
      }

      Vec penalForceInc(arrayPenalForce[i][j][k][0], arrayPenalForce[i][j][k][1],
                        arrayPenalForce[i][j][k][2]);
      Vec pressForceInc(arrayPressureForce[i][j][k][0], arrayPressureForce[i][j][k][1],
                        arrayPressureForce[i][j][k][2]);
      penalForce += penalForceInc;
      presForce += pressForceInc;

      // r X F,  % is overloaded as cross product
      penalMoment += cross(dist, penalForceInc);
      presMoment += cross(dist, pressForceInc);

    } // end of fluidGrid loop

    penalForce *= dx * dy * dz;
    presForce *= dx * dy * dz;
    ptcl->addForce(penalForce);
    ptcl->addForce(presForce);

    penalMoment *= dx * dy * dz;
    presMoment *= dx * dy * dz;
    ptcl->addMoment(penalMoment);
    ptcl->addMoment(presMoment);

    if (ptcl->getId() == 1) {
      ofs << std::setw(OWID) << iteration << std::setw(OWID) << timeAccrued
          << std::setw(OWID) << penalForce.x() << std::setw(OWID)
          << penalForce.y() << std::setw(OWID) << penalForce.z()
          << std::setw(OWID) << presForce.x() << std::setw(OWID)
          << presForce.y() << std::setw(OWID) << presForce.z()
          << std::setw(OWID) << penalMoment.x() << std::setw(OWID)
          << penalMoment.y() << std::setw(OWID) << penalMoment.z()
          << std::setw(OWID) << presMoment.x() << std::setw(OWID)
          << presMoment.y() << std::setw(OWID) << presMoment.z()
          << std::setw(OWID) << ptcl->getAccel().x() << std::setw(OWID)
          << ptcl->getAccel().y() << std::setw(OWID)
          << ptcl->getAccel().z() << std::setw(OWID)
          << ptcl->currentVel().x() << std::setw(OWID)
          << ptcl->currentVel().y() << std::setw(OWID)
          << ptcl->currentVel().z() << std::endl;
    }
  } // end of particle loop
}

void
Fluid::plot(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error: Fluid::plot" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "VARIABLES = \"x\"" << std::setw(OWID) << "\"y\""
      << std::setw(OWID) << "\"z\"" << std::setw(OWID) << "\"mach\""
      << std::setw(OWID) << "\"density\"" << std::setw(OWID) << "\"momentum_x\""
      << std::setw(OWID) << "\"momentum_y\"" << std::setw(OWID)
      << "\"momentum_z\"" << std::setw(OWID) << "\"energy\"" << std::setw(OWID)
      << "\"velocity_x\"" << std::setw(OWID) << "\"velocity_y\""
      << std::setw(OWID) << "\"velocity_z\"" << std::setw(OWID)
      << "\"pressure\"" << std::setw(OWID) << "\"temperature\""
      << std::setw(OWID) << "\"mask\"" << std::setw(OWID) << "\"penalFx\""
      << std::setw(OWID) << "\"penalFy\"" << std::setw(OWID) << "\"penalFz\""
      << std::setw(OWID) << "\"pressureFx\"" << std::setw(OWID)
      << "\"pressureFy\"" << std::setw(OWID) << "\"pressureFz\"" << std::endl;

  ofs << "ZONE I=" << nx - 2 << ", J=" << ny - 2 << ", K=" << nz - 2
      << ", DATAPACKING=POINT" << std::endl;

  for (std::size_t k = 1; k < nz - 1; ++k)
    for (std::size_t j = 1; j < ny - 1; ++j)
      for (std::size_t i = 1; i < nx - 1; ++i) {
        ofs << std::setw(OWID) << arrayGridCoord[i][j][k][0] << std::setw(OWID)
            << arrayGridCoord[i][j][k][1] << std::setw(OWID)
            << arrayGridCoord[i][j][k][2] << std::setw(OWID)
            << vnormL2(Vec(arrayU[i][j][k][var_vel[0]],
                         arrayU[i][j][k][var_vel[1]],
                         arrayU[i][j][k][var_vel[2]])) /
                 arraySoundSpeed[i][j][k]
            << std::setw(OWID) << arrayU[i][j][k][var_den] << std::setw(OWID)
            << arrayU[i][j][k][var_mom[0]] << std::setw(OWID)
            << arrayU[i][j][k][var_mom[1]] << std::setw(OWID)
            << arrayU[i][j][k][var_mom[2]] << std::setw(OWID)
            << arrayU[i][j][k][var_eng] << std::setw(OWID)
            << arrayU[i][j][k][var_vel[0]] << std::setw(OWID)
            << arrayU[i][j][k][var_vel[1]] << std::setw(OWID)
            << arrayU[i][j][k][var_vel[2]] << std::setw(OWID)
            << arrayU[i][j][k][var_prs] << std::setw(OWID)
            << arrayU[i][j][k][var_prs] / (Rs * arrayU[i][j][k][var_den])
            << std::setw(OWID) << arrayU[i][j][k][var_msk] << std::setw(OWID)
            << arrayPenalForce[i][j][k][0] << std::setw(OWID)
            << arrayPenalForce[i][j][k][1] << std::setw(OWID)
            << arrayPenalForce[i][j][k][2] << std::setw(OWID)
            << arrayPressureForce[i][j][k][0] << std::setw(OWID)
            << arrayPressureForce[i][j][k][1] << std::setw(OWID)
            << arrayPressureForce[i][j][k][2] << std::endl;
      }

  ofs.close();
}

} // name space dem
