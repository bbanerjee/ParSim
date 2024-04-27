/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2018-2023 Parresia Research Limited, NZ
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

#ifndef PARTICLE_INTERPOLATOR_H
#define PARTICLE_INTERPOLATOR_H

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <array>
#include <vector>

#include <Core/Grid/Variables/NCVariable.h>
namespace Uintah {

class Patch;
class Stencil7;
using std::vector;
using Uintah::IntVector;
using Uintah::Point;
using Uintah::Vector;

using Array2_D = std::array<double, 2>;
using Array8_D = std::array<double, 8>;
using Array8_V = std::array<Vector, 8>;
using Array8_P = std::array<Point, 8>;

enum class IntegrationRule
{
  ONE_POINT,
  EIGHT_POINT
};

class ParticleInterpolator
{

public:
  ParticleInterpolator(){};
  virtual ~ParticleInterpolator(){};

  virtual std::unique_ptr<ParticleInterpolator>
  clone(const Patch*) = 0;

  virtual void
  setLcrit([[maybe_unused]] double d_cpdi_lcrit){};

  virtual void
  findCellAndWeights(const Point& p,
                     std::vector<IntVector>& ni,
                     std::vector<double>& S,
                     const Matrix3& size,
                     const Matrix3& defgrad) = 0;

  virtual void
  findCellAndWeights([[maybe_unused]] const Point& p,
                     [[maybe_unused]] std::vector<IntVector>& ni,
                     [[maybe_unused]] std::vector<double>& S){};

  virtual void
  findCellAndShapeDerivatives(const Point& pos,
                              std::vector<IntVector>& ni,
                              std::vector<Vector>& d_S,
                              const Matrix3& size,
                              const Matrix3& defgrad = Matrix3(0)) = 0;

  virtual void
  findCellAndWeightsAndShapeDerivatives(const Point& pos,
                                        std::vector<IntVector>& ni,
                                        std::vector<double>& S,
                                        std::vector<Vector>& d_S,
                                        const Matrix3& size,
                                        const Matrix3& defgrad) = 0;

  //__________________________________
  //  Needed for AMRMPM
  virtual void
  findCellAndWeights(const Point& pos,
                     std::vector<IntVector>& ni,
                     std::vector<double>& S,
                     constNCVariable<Stencil7>& zoi,
                     constNCVariable<Stencil7>& zoi_fine,
                     const bool& getFiner,
                     int& num_cur,
                     int& num_fine,
                     int& num_coarse,
                     const Vector& size,
                     bool coarse_part,
                     const Patch* patch) = 0;

  virtual void
  findCellAndWeights_CFI(const Point& pos,
                         std::vector<IntVector>& ni,
                         std::vector<double>& S,
                         constNCVariable<Stencil7>& zoi) = 0;

  virtual void
  findCellAndWeightsAndShapeDerivatives_CFI(const Point& pos,
                                            std::vector<IntVector>& CFI_ni,
                                            std::vector<double>& S,
                                            std::vector<Vector>& d_S,
                                            constNCVariable<Stencil7>& zoi) = 0;
  virtual int
  size() = 0;

  /*!
   * Find the volume averaged (over a hexahedral particle domain)
   * shape functions and their derivatives for a set of grid points
   * that influence the particle.
   */
  void
  findCellAndWeightsAndShapeDerivativeAverages(
    const Patch* patch,
    const Point& pPosition,
    std::vector<IntVector>& numInfluenceNodes,
    std::vector<double>& S_ip_av,
    std::vector<Vector>& G_ip_av,
    const Matrix3& pSize,
    const Matrix3& pDefGrad,
    bool derivatives = false);

protected:
  /*!
   * Get the vertex positions of the particle parallelepiped
   */
  Array8_P
  getParticleDomainHex8Normalized(const Point& pPosition,
                                  const Matrix3& pSize,
                                  const Matrix3& pDefGrad) const;

  /*!
   * PURPOSE
   *  Get s-point abscissas and weight for quad product Gauss rule
   * ORIGIAL AUTHOR
   *  C. A. Felippa, Jan 1966 (Fortran IV, UCB 7094)
   * RETURNS
   *  the abscissae and weight factors
   *  of the p-th Gauss-Legendre integration rule (=1,2,3,4) over
   *  the interval xi:(-1,+1).
   * INPUT ARGUMENTS:
   *   numIntegrationPts: Number of points in the integration rule (1 to 4)
   *                      If le 1 assume 1; if gt 4 assume 4.
   *   index:  Index of sample point (1 to P)
   *
   * OUTPUT ARGUMENTS:
   *   xi: Abscissa of sample point (zero of Legendre polynomial)
   *   weight: Weight factor
   */
  void
  lineGaussQuadratureWeights(int numIntegrationPts,
                             int index,
                             double& xi,
                             double& weight);

  /*!
   * PURPOSE :
   *   Returns the hexahedron coordinates of sample points and weights for a
   *   s-point Gauss-product integration rule over an isoparametric hexahedron.
   * ORIGINAL AUTHOR :
   *   C. A. Felippa, June 1967 (Fortran IV, UCB 7094)
   * INPUT :
   *   numGaussPtsXi   : Number of Gauss points in the XI direction
   *   indexXi         : Index of sample point in the XI direction
   *   numGaussPtsEta  : Number of Gauss points in the ETA direction
   *   indexEta        : Index of sample point in the ETA direction
   *   numGaussPtsZeta : Number of Gauss points in the ZETA direction
   *   indexZeta       : Index of sample points in the MU direction
   * OUTPUT :
   *   xi, eta, zeta   : Isoparametric coordinates of sample point
   *   weight          : Weight factor
   */
  void
  getGaussPtsAndWeightsHex8(int numGaussPtsXi,
                            int indexXi,
                            int numGaussPtsEta,
                            int indexEta,
                            int numGaussPtsZeta,
                            int indexZeta,
                            double& xi,
                            double& eta,
                            double& zeta,
                            double& weight) const;
  /*!
   * PURPOSE :
   *   Compute signed volume of an 8-node brick element
   * ORIGINAL AUTHOR :
   *   C. A. Felippa, August 1973 (Fortran IV, LMSC 1108)
   * INPUT :
   *   x, y, z : Global coordinates of brick corners
   *   integration_rule :  Specifies the Gauss integration rule:
   *     ONE_POINT   :  1 point rule at center
   *     EIGHT_POINT :  8 point rule (2 x 2 x 2) which gives the exact volume
   *                    for any 8-node brick
   * OUTPUT :
   *   vol : Signed volume of element. vol<=0 flags an error.
   */
  double
  volumeHex8(const Array8_D& x,
             const Array8_D& y,
             const Array8_D& z,
             IntegrationRule integration_rule) const;
  /*
   * PURPOSE :
   *   Computes natural partial derivatives of trilinear shape functions of
   *    isoparametric 8-node brick
   * ASSUMPTIONS:
   *   Nodes are ordered as x-, y-, z-
   *                        x+, y-, z-
   *                        x+, y+, z-
   *                        x-, y+, z-
   *                        x-, y-, z+
   *                        x+, y-, z+
   *                        x+, y+, z+
   *                        x-, y+, z+
   * ORIGINAL AUTHOR :
   *   C. A. Felippa, August 1973 (Fortran IV, LMSC 1108)
   * INPUT :
   *   xi, eta, zeta
   *     Natural coordinates of point at which derivatives are to be evaluated
   * OUTPUT :
   *   dNxi : Array of 8 partial derivatives wrt xi
   *   dNeta : Array of 8 partial derivatives wrt eta
   *   dNzeta : Array of 8 partial derivatives wrt zeta
   */
  void
  naturalDerivativesHex8(double xi,
                         double eta,
                         double zeta,
                         Array8_D& dNxi,
                         Array8_D& dNeta,
                         Array8_D& dNzeta) const;
  /*!
   * PURPOSE
   *  Get s-point abscissas and weight for quad product Gauss rule
   * ORIGINAL AUTHOR
   *  C. A. Felippa, Jan 1966 (Fortran IV, UCB 7094)
   * RETURNS
   *  the abscissae and weight factors
   *  of the p-th Gauss-Legendre integration rule (=1,2,3,4) over
   *  the interval xi:(-1,+1).
   * INPUT ARGUMENTS:
   *   numIntegrationPts: Number of points in the integration rule (1 to 4)
   *                      If le 1 assume 1; if gt 4 assume 4.
   *   index:  Index of sample point (1 to P)
   * OUTPUT ARGUMENTS:
   *   xi: Abscissa of sample point (zero of Legendre polynomial)
   *   weight: Weight factor
   */
  void
  lineGaussQuadratureWeights(int numIntegrationPts,
                             int index,
                             double& xi,
                             double& weight) const;
};
} // namespace Uintah

#endif
