/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/GradientComputer/VelocityGradientComputer.h>
#include <Core/Exceptions/InvalidValue.h>

using namespace Uintah;

VelocityGradientComputer::VelocityGradientComputer(const MPMFlags* Mflag)
  : GradientComputer(Mflag)
{
}

VelocityGradientComputer::VelocityGradientComputer(
  const VelocityGradientComputer* gc)
  : GradientComputer(gc)
{
}

VelocityGradientComputer*
VelocityGradientComputer::clone()
{
  return scinew VelocityGradientComputer(*this);
}

VelocityGradientComputer::~VelocityGradientComputer() {}

// Actually compute velocity gradient
void
VelocityGradientComputer::computeVelGrad(ParticleInterpolator* interpolator,
                                         const double* oodx,
                                         const short pgFld[],
                                         const Point& px,
                                         const Matrix3& pSize,
                                         const Matrix3& pDefGrad_old,
                                         constNCVariable<Vector>& gVelocity,
                                         constNCVariable<Vector>& GVelocity,
                                         Matrix3& velGrad_new)
{
  double numInfluenceNodes = interpolator->size();
  std::vector<IntVector> ni(numInfluenceNodes);
  std::vector<Vector> d_S(numInfluenceNodes);
  if (!flag->d_axisymmetric) {

    // Get the node indices that surround the cell
    interpolator->findCellAndShapeDerivatives(px, ni, d_S, pSize, pDefGrad_old);

    // Fracture
    if (flag->d_fracture) {
      // Special vel grad for fracture
      computeVelocityGradient(velGrad_new,
                              ni,
                              d_S,
                              oodx,
                              pgFld,
                              gVelocity,
                              GVelocity);
    } else {
      // Standard 3d vel grad computation
      computeGrad(velGrad_new, ni, d_S, oodx, gVelocity);
    }
  } else { // axi-symmetric kinematics
    // Get the node indices that surround the cell
    std::vector<double> S(numInfluenceNodes);
    interpolator->findCellAndWeightsAndShapeDerivatives(px,
                                                        ni,
                                                        S,
                                                        d_S,
                                                        pSize,
                                                        pDefGrad_old);
    // x -> r, y -> z, z -> theta
    computeAxiSymVelocityGradient(velGrad_new, ni, d_S, S, oodx, gVelocity, px);
  } // endif (!flag->d_axisymmetric)

  if (std::isnan(velGrad_new.Norm())) {
    std::ostringstream out;
    out << "**ERROR**: Nan in velocity gradient value." << std::endl;
    out << " velGrad = " << velGrad_new << std::endl;
    // out << " ni = " << ni << " d_S = " << d_S << " oodx = " << oodx <<
    // std::endl;
    for (int k = 0; k < flag->d_8or27; k++) {
      out << " gVelocity [" << ni[k] << " = " << gVelocity[ni[k]] << std::endl;
    }
    out << " pDefGrad_old = " << pDefGrad_old << std::endl;
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  return;
}

//-------------------------------------------------------------------------
// Protected methods
//-------------------------------------------------------------------------
void
VelocityGradientComputer::computeAxiSymVelocityGradient(
  Matrix3& velGrad,
  std::vector<IntVector>& ni,
  std::vector<Vector>& d_S,
  std::vector<double>& S,
  const double* oodx,
  constNCVariable<Vector>& gVelocity,
  const Point& px)
{
  // x -> r, y -> z, z -> theta
  for (int k = 0; k < flag->d_8or27; k++) {
    Vector gvel = gVelocity[ni[k]];
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        velGrad(i, j) += gvel[i] * d_S[k][j] * oodx[j];
      }
    }
    velGrad(2, 2) += gvel.x() * d_S[k].z();
  }
}

void
VelocityGradientComputer::computeVelocityGradient(
  Matrix3& velGrad,
  std::vector<IntVector>& ni,
  std::vector<Vector>& d_S,
  const double* oodx,
  const short pgFld[],
  constNCVariable<Vector>& gVelocity,
  constNCVariable<Vector>& GVelocity)
{
  Vector gvel(0., 0., 0);
  for (int k = 0; k < flag->d_8or27; k++) {
    if (pgFld[k] == 1) {
      gvel = gVelocity[ni[k]];
    }
    if (pgFld[k] == 2) {
      gvel = GVelocity[ni[k]];
    }
    for (int j = 0; j < 3; j++) {
      double d_SXoodx = d_S[k][j] * oodx[j];
      for (int i = 0; i < 3; i++) {
        velGrad(i, j) += gvel[i] * d_SXoodx;
      }
    }
  }
}
