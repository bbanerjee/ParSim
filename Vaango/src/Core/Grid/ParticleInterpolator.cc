/*
 * The MIT License
 *
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

#include <Core/Grid/ParticleInterpolator.h>
#include <Core/Math/MiscMath.h>

using namespace Uintah;

/*!
 * Get the normalized vertex positions of the particle parallelepiped
 */
Array8_P
ParticleInterpolator::getParticleDomainHex8Normalized(const Point& pNormPosition,
                                                      const Matrix3& pSize,
                                                      const Matrix3& pDefGrad) const
{
  // Compute the deformed particle domain basis vectors
  Matrix3 defBasis = pDefGrad * pSize;
  Vector e1 = {defBasis(0,0) * 0.5, defBasis(1,0) * 0.5, defBasis(2,0) * 0.5};
  Vector e2 = {defBasis(0,1) * 0.5, defBasis(1,1) * 0.5, defBasis(2,1) * 0.5};
  Vector e3 = {defBasis(0,2) * 0.5, defBasis(1,2) * 0.5, defBasis(2,2) * 0.5};

  // Compute the particle domain vertex locations
  Array8_P vertices;
  vertices[0] = (pNormPosition - e1 - e2 - e3);
  vertices[1] = (pNormPosition + e1 - e2 - e3);
  vertices[2] = (pNormPosition + e1 + e2 - e3);
  vertices[3] = (pNormPosition - e1 + e2 - e3);

  vertices[4] = (pNormPosition - e1 - e2 + e3);
  vertices[5] = (pNormPosition + e1 - e2 + e3);
  vertices[6] = (pNormPosition + e1 + e2 + e3);
  vertices[7] = (pNormPosition - e1 + e2 + e3);
  
  return vertices;
}

/*!
 * Find the volume averaged (over a hexahedral particle domain)
 * shape functions and their derivatives for a set of grid points
 * that influence the particle.
 */
void 
ParticleInterpolator::findCellAndWeightsAndShapeDerivativeAverages(
  const Patch* patch,
  const Point& pPosition,
  std::vector<IntVector>& numInfluenceNodes,
  std::vector<double>& S_ip_av, std::vector<Vector>& G_ip_av,
  const Matrix3& pSize, const Matrix3& pDefGrad,
  bool derivatives)
{
  // Get the normalized particle position 
  Point pNormPosition = patch->getLevel()->positionToIndex(pPosition);

  // Identify the indices of the eight nodes that influence the particle
  auto p_ix_lower = Floor(pNormPosition.x());
  auto p_iy_lower = Floor(pNormPosition.y());
  auto p_iz_lower = Floor(pNormPosition.z());
  auto p_ix_upper = p_ix_lower + 1;
  auto p_iy_upper = p_iy_lower + 1;
  auto p_iz_upper = p_iz_lower + 1;
  IntVector p_node1(p_ix_lower, p_iy_lower, p_iz_lower);
  IntVector p_node2(p_ix_upper, p_iy_lower, p_iz_lower);
  IntVector p_node3(p_ix_upper, p_iy_upper, p_iz_lower);
  IntVector p_node4(p_ix_lower, p_iy_upper, p_iz_lower);
  IntVector p_node5(p_ix_lower, p_iy_lower, p_iz_upper);
  IntVector p_node6(p_ix_upper, p_iy_lower, p_iz_upper);
  IntVector p_node7(p_ix_upper, p_iy_upper, p_iz_upper);
  IntVector p_node8(p_ix_lower, p_iy_upper, p_iz_upper);

  // Get the normalized vertices of the particle domain
  Array8_P pNormVertices = getParticleDomainHex8Normalized(pNormPosition, pSize, pDefGrad);

  // For the particle position and domain vertices, find the nodes of influence
  for (auto& vertex : pNormVertices) {

    // Get the node indices closest to the particle domain vertex
    auto v_ix_lower = Floor(vertex.x());
    auto v_iy_lower = Floor(vertex.y());
    auto v_iz_lower = Floor(vertex.z());
    auto v_ix_upper = v_ix_lower + 1;
    auto v_iy_upper = v_iy_lower + 1;
    auto v_iz_upper = v_iz_lower + 1;
    IntVector v_node1(v_ix_lower, v_iy_lower, v_iz_lower);
    IntVector v_node2(v_ix_upper, v_iy_lower, v_iz_lower);
    IntVector v_node3(v_ix_upper, v_iy_upper, v_iz_lower);
    IntVector v_node4(v_ix_lower, v_iy_upper, v_iz_lower);
    IntVector v_node5(v_ix_lower, v_iy_lower, v_iz_upper);
    IntVector v_node6(v_ix_upper, v_iy_lower, v_iz_upper);
    IntVector v_node7(v_ix_upper, v_iy_upper, v_iz_upper);
    IntVector v_node8(v_ix_lower, v_iy_upper, v_iz_upper);


  }

  // Find the nodes that influence the particle and evaluate the shape
  // function at the particle position
  std::vector<double> S;
  findCellAndWeights(pPosition, numInfluenceNodes, S, pSize, pDefGrad);

  /*
  C=DECK HEXA8STIF
  C=PURPOSE Form elastic stiffness of 8-node hexahedron
  C=AUTHOR C. A. Felippa, June 1967 (Fortran IV, UCB 7094)
  C=VERSION May 1981 (Fortran 77, Lockheed Palo Alto VAX)
  C=EQUIPMENT Machine independent, except for long externals and
  lowercase
  C=KEYWORDS 3D linear elasticity eight node hexahedron
  C=KEYWORDS finite element stiffness matrix
  C=BLOCK ABSTRACT
  C
  C HEXA8STIF forms the element stiffness matrix of an 8-node
  C hexahedron ("brick") for linear elastic structural analysis.
  C
  C=END ABSTRACT
  C=BLOCK USAGE
  C
  C The calling sequence is
  C
  C call HEXA8STIF (OPT, X, Y, Z, C, P, SM, MD, STATUS)
  C
  C where the input arguments are
  C
  C OPT Option letter argument. Ignored in this
  implementation.
  C X (8 x 1) array of x coordinates of hexahedron nodes
  C Y (8 x 1) array of y coordinates of hexahedron nodes
  C Z (8 x 1) array of z coordinates of hexahedron nodes
  C C (6 x 6) constitutive material matrix 
  C P Gauss quadrature rule (1 to 4) assumed same in all
  iso-P
  C directions. Normally P=2 is used.
  C MD First dimension of SM in calling program.
  C
  C The outputs are:
  C
  C SM (24 x 24) computed element stiffness matrix. The
  C arrangement of rows and columns pertains to node
  C displacements arranged in the order
  C (vx1, vy1, vz1, vx2, vy2, vz2, ... vz8)
  C This particular ordering is obtained by setting
  array
  C LS to 1,4,7,10,13,16,19,22,2,5,8,11,14,
  C 17,20,23,3,6,9,12,15,18,21,24 (cf DATA statement for
  LS
  C below). To get other orderings reset that
  statement.
  C
  C STATUS Status character variable. Blank if no error
  C detected. Else error message.
  C
  C=END USAGE
  C=BLOCK FORTRAN
  subroutine HEXA8STIF
  & (opt, x, y, z, c, p, sm, md, status)
  C
  C A R G U M E N T S
  C
  character*(*) opt, status
  integer p, md
  double precision x(8), y(8), z(8), c(6,6)
  double precision sm(md,md)
  C
  C L O C A L V A R I A B L E S
  C
  double precision q(8), qx(8), qy(8), qz(8), cm(6,6)
  double precision xi, eta, mu, det, w
  double precision c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z
  double precision c4x, c4y, c4z, c5x, c5y, c5z, c6x, c6y, c6z
  integer i, ix, iy, iz, j, jx, jy, jz, k, l, m
  integer ls(24)
  C
  C D A T A
  C
  data ls /1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,
  & 3,6,9,12,15,18,21,24/
  C
  C L O G I C
  C
  status = ' '
  do 1200 j = 1,24
  do 1100 i = 1,24
  sm(i,j) = 0.0
  1100 continue
  1200 continue
  do 1400 j = 1,6
  do 1300 i = 1,6
  cm(i,j) = c(i,j)
  1300 continue
  1400 continue
  C
  do 3000 k = 1,p
  do 2500 l = 1,p
  do 2400 m = 1,p
  call HEXAGAUSSQ (p, k, p, l, p, m, xi, eta, mu, w)
  call HEXA8SHAPE (' ', xi, eta, mu, x, y, z, q, 
  & qx, qy, qz, det)
  if (det .le. 0.0) then
  status = 'Negative Jacobian determinant'
  if (det .eq. 0.0) then
  status = 'Zero Jacobian determinant'
  end if
  return
  end if
  w = w * det 
  C
  do 2000 j = 1,8
  jx = ls(j)
  jy = ls(j+8)
  jz = ls(j+16)
  c1x = (cm(1,1)*qx(j) + cm(1,4)*qy(j) + cm(1,6)*qz(j)) *w
  c1y = (cm(1,4)*qx(j) + cm(1,2)*qy(j) + cm(1,5)*qz(j)) *w
  c1z = (cm(1,6)*qx(j) + cm(1,5)*qy(j) + cm(1,3)*qz(j)) *w
  c2x = (cm(2,1)*qx(j) + cm(2,4)*qy(j) + cm(2,6)*qz(j)) *w
  c2y = (cm(2,4)*qx(j) + cm(2,2)*qy(j) + cm(2,5)*qz(j)) *w
  c2z = (cm(2,6)*qx(j) + cm(2,5)*qy(j) + cm(2,3)*qz(j)) *w
  c3x = (cm(3,1)*qx(j) + cm(3,4)*qy(j) + cm(3,6)*qz(j)) *w
  c3y = (cm(3,4)*qx(j) + cm(3,2)*qy(j) + cm(3,5)*qz(j)) *w
  c3z = (cm(3,6)*qx(j) + cm(3,5)*qy(j) + cm(3,3)*qz(j)) *w
  c4x = (cm(4,1)*qx(j) + cm(4,4)*qy(j) + cm(4,6)*qz(j)) *w
  c4y = (cm(4,4)*qx(j) + cm(4,2)*qy(j) + cm(4,5)*qz(j)) *w
  c4z = (cm(4,6)*qx(j) + cm(4,5)*qy(j) + cm(4,3)*qz(j)) *w
  c5x = (cm(5,1)*qx(j) + cm(5,4)*qy(j) + cm(5,6)*qz(j)) *w
  c5y = (cm(5,4)*qx(j) + cm(5,2)*qy(j) + cm(5,5)*qz(j)) *w
  c5z = (cm(5,6)*qx(j) + cm(5,5)*qy(j) + cm(5,3)*qz(j)) *w
  c6x = (cm(6,1)*qx(j) + cm(6,4)*qy(j) + cm(6,6)*qz(j)) *w
  c6y = (cm(6,4)*qx(j) + cm(6,2)*qy(j) + cm(6,5)*qz(j)) *w
  c6z = (cm(6,6)*qx(j) + cm(6,5)*qy(j) + cm(6,3)*qz(j)) *w
  do 1600 i = j,8
  ix = ls(i)
  iy = ls(i+8)
  iz = ls(i+16)
  sm(ix,jx) = sm(ix,jx)+qx(i)*c1x + qy(i)*c4x +
  qz(i)*c6x
  sm(jx,ix) = sm(ix,jx)
  sm(iy,jy) = sm(iy,jy)+qx(i)*c4y + qy(i)*c2y +
  qz(i)*c5y
  sm(jy,iy) = sm(iy,jy)
  sm(iz,jz) = sm(iz,jz)+qx(i)*c6z + qy(i)*c5z +
  qz(i)*c3z
  sm(jz,iz) = sm(iz,jz)
  sm(ix,jy) = sm(ix,jy)+qx(i)*c1y + qy(i)*c4y +
  qz(i)*c6y
  sm(iy,jx) = sm(iy,jx)+qx(i)*c4x + qy(i)*c2x +
  qz(i)*c5x
  sm(ix,jz) = sm(ix,jz)+qx(i)*c1z + qy(i)*c4z +
  qz(i)*c6z
  sm(iz,jx) = sm(iz,jx)+qx(i)*c6x + qy(i)*c5x +
  qz(i)*c3x
  sm(iy,jz) = sm(iy,jz)+qx(i)*c4z + qy(i)*c2z +
  qz(i)*c5z
  sm(iz,jy) = sm(iz,jy)+qx(i)*c6y + qy(i)*c5y +
  qz(i)*c3y
  sm(jx,iy) = sm(iy,jx)
  sm(jx,iz) = sm(iz,jx)
  sm(jy,ix) = sm(ix,jy)
  sm(jy,iz) = sm(iz,jy)
  sm(jz,ix) = sm(ix,jz)
  sm(jz,iy) = sm(iy,jz)
  1600 continue
  2000 continue
  2400 continue
  2500 continue
  3000 continue
  C
  return
  end
  */
}

/*!
 * PURPOSE : 
 *  Compute the value of the shape functions for a  eight-noded isoparametric 
 *  hexahedron ("brick element") and its C x-y-z derivatives, at a sample 
 *  point given by its hexahedron coordinates (xi,eta,zeta)
 * ORIGINAL AUTHOR :
 *  C. A. Felippa, June 1967 (Fortran IV, UCB 7094) 
 * INPUT :
 *  xi, eta,zeta : Isoparametric coordinates of given point
 *  x            : (8 x 1) array of x coordinates of hexahedron corners
 *  y            : (8 x 1) array of y coordinates of hexahedron corners
 *  z            : (8 x 1) array of z coordinates of hexahedron corners
 *  derivatives  : if False: return after computing S and detJ
 *                 if True: return after computing S, detJ, dSx, dSy, and dSz
 * OUTPUT :
 *  S            : (8 x 1) array of shape function values
 *  dSx          : (8 x 1) array of shape function x-derivatives (if DET>0)
 *  dSy          : (8 x 1) array of shape function y-derivatives (if DET>0)
 *  dSz          :  (8 x 1) array of shape function z-derivatives (if DET>0)
 *  detJ         : Value of Jacobian determinant. A nonpositive value should
 *                 be treated as an error in the calling program. In such
 *                 a case dSx,dSy,dSz are not computed.
 */
void shapeFunctionsHex8(double xi, double eta, double zeta, 
                        const Array8_D& x, const Array8_D& y, const Array8_D& z,
                        Array8_D& S, 
                        Array8_D& dSx, Array8_D& dSy, Array8_D& dSz,
                        double& detJ, 
                        bool derivatives = false)
{
  double d1 = 0.5 * (1.0+xi);
  double d2 = 0.5 * (1.0+eta);
  double d3 = 0.5 * (1.0+zeta);
  double d4 = 1.0 - d1;
  double d5 = 1.0 - d2;
  double d6 = 1.0 - d3;
  S[0] = d4 * d5 * d6;
  S[1] = d1 * d5 * d6;
  S[2] = d1 * d2 * d6;
  S[3] = d4 * d2 * d6;
  S[4] = d4 * d5 * d3;
  S[5] = d1 * d5 * d3;
  S[6] = d1 * d2 * d3;
  S[7] = d4 * d2 * d3;
  double d1h = 0.5 * d1;
  double d2h = 0.5 * d2;
  [[maybe_unused]] double d3h = 0.5 * d3;
  double d4h = 0.5 * d4;
  double d5h = 0.5 * d5;
  [[maybe_unused]] double d6h = 0.5 * d6;

  Array8_D s1, s2, s3;
  s1[0] = - d5h * d6;
  s1[1] = d5h * d6;
  s1[2] = d2h * d6;
  s1[3] = - d2h * d6;
  s1[4] = - d5h * d3;
  s1[5] = d5h * d3;
  s1[6] = d2h * d3;
  s1[7] = - d2h * d3;
  s2[0] = - d4h * d6;
  s2[1] = - d1h * d6;
  s2[2] = d1h * d6;
  s2[3] = d4h * d6;
  s2[4] = - d4h * d3;
  s2[5] = - d1h * d3;
  s2[6] = d1h * d3;
  s2[7] = d4h * d3;
  s3[0] = - d4h * d5;
  s3[1] = - d1h * d5;
  s3[2] = - d1h * d2;
  s3[3] = - d4h * d2;
  s3[4] = d4h * d5;
  s3[5] = d1h * d5;
  s3[6] = d1h * d2;
  s3[7] = d4h * d2;

  double xd1 = 0.0, xd2 = 0.0, xd3 = 0.0;
  double yd1 = 0.0, yd2 = 0.0, yd3 = 0.0;
  double zd1 = 0.0, zd2 = 0.0, zd3 = 0.0;
  for (int ii = 0; ii < 8; ++ii) {
    xd1 += x[ii] * s1[ii];
    xd2 += x[ii] * s2[ii];
    xd3 += x[ii] * s3[ii];
    yd1 += y[ii] * s1[ii];
    yd2 += y[ii] * s2[ii];
    yd3 += y[ii] * s3[ii];
    zd1 += z[ii] * s1[ii];
    zd2 += z[ii] * s2[ii];
    zd3 += z[ii] * s3[ii];
  }
  double a11 = yd2*zd3 - yd3*zd2;
  double a12 = yd3*zd1 - yd1*zd3;
  double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3;
  double a22 = xd1*zd3 - zd1*xd3;
  double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2;
  double a32 = yd1*xd3 - xd1*yd3;
  double a33 = xd1*yd2 - yd1*xd2;
  detJ = xd1*a11 + yd1*a21 + zd1*a31;

  if (!derivatives || !(detJ > 0.0)) return;

  double cdet = 1.0 / detJ;
  dSx[0] = cdet * ( a11 * s1[0] + a12 * s2[0] + a13 * s3[0]);
  dSy[0] = cdet * ( a21 * s1[0] + a22 * s2[0] + a23 * s3[0]);
  dSz[0] = cdet * ( a31 * s1[0] + a32 * s2[0] + a33 * s3[0]);
  dSx[1] = cdet * ( a11 * s1[1] + a12 * s2[1] + a13 * s3[1]);
  dSy[1] = cdet * ( a21 * s1[1] + a22 * s2[1] + a23 * s3[1]);
  dSz[1] = cdet * ( a31 * s1[1] + a32 * s2[1] + a33 * s3[1]);
  dSx[2] = cdet * ( a11 * s1[2] + a12 * s2[2] + a13 * s3[2]);
  dSy[2] = cdet * ( a21 * s1[2] + a22 * s2[2] + a23 * s3[2]);
  dSz[2] = cdet * ( a31 * s1[2] + a32 * s2[2] + a33 * s3[2]);
  dSx[3] = cdet * ( a11 * s1[3] + a12 * s2[3] + a13 * s3[3]);
  dSy[3] = cdet * ( a21 * s1[3] + a22 * s2[3] + a23 * s3[3]);
  dSz[3] = cdet * ( a31 * s1[3] + a32 * s2[3] + a33 * s3[3]);
  dSx[4] = cdet * ( a11 * s1[4] + a12 * s2[4] + a13 * s3[4]);
  dSy[4] = cdet * ( a21 * s1[4] + a22 * s2[4] + a23 * s3[4]);
  dSz[4] = cdet * ( a31 * s1[4] + a32 * s2[4] + a33 * s3[4]);
  dSx[5] = cdet * ( a11 * s1[5] + a12 * s2[5] + a13 * s3[5]);
  dSy[5] = cdet * ( a21 * s1[5] + a22 * s2[5] + a23 * s3[5]);
  dSz[5] = cdet * ( a31 * s1[5] + a32 * s2[5] + a33 * s3[5]);
  dSx[6] = cdet * ( a11 * s1[6] + a12 * s2[6] + a13 * s3[6]);
  dSy[6] = cdet * ( a21 * s1[6] + a22 * s2[6] + a23 * s3[6]);
  dSz[6] = cdet * ( a31 * s1[6] + a32 * s2[6] + a33 * s3[6]);
  dSx[7] = cdet * ( a11 * s1[7] + a12 * s2[7] + a13 * s3[7]);
  dSy[7] = cdet * ( a21 * s1[7] + a22 * s2[7] + a23 * s3[7]);
  dSz[7] = cdet * ( a31 * s1[7] + a32 * s2[7] + a33 * s3[7]);
}

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
ParticleInterpolator::getGaussPtsAndWeightsHex8(int numGaussPtsXi, int indexXi,
                                                int numGaussPtsEta, int indexEta,
                                                int numGaussPtsZeta, int indexZeta,
                                                double& xi, double& eta, double& zeta,
                                                double& weight) const
{
  double w1, w2, w3;
  lineGaussQuadratureWeights(numGaussPtsXi, indexXi, xi, w1);
  lineGaussQuadratureWeights(numGaussPtsEta, indexEta, eta, w2);
  lineGaussQuadratureWeights(numGaussPtsZeta, indexZeta, zeta, w3);
  weight = w1 * w2 * w3;
}

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
ParticleInterpolator::volumeHex8(const Array8_D& x, const Array8_D& y, const Array8_D& z,
                                 IntegrationRule integration_rule) const
{
  double vol = 0.0;
  if (integration_rule == IntegrationRule::ONE_POINT) {

    Array8_D dNxi, dNeta, dNzeta;
    naturalDerivativesHex8(0.0, 0.0, 0.0, dNxi, dNeta, dNzeta);

    double J11 = 0.0; double J12 = 0.0; double J13 = 0.0;
    double J21 = 0.0; double J22 = 0.0; double J23 = 0.0;
    double J31 = 0.0; double J32 = 0.0; double J33 = 0.0;
    for (int ii = 0; ii < 8; ++ii) {
      J11 += x[ii]*dNxi[ii];
      J12 += y[ii]*dNxi[ii];
      J13 += z[ii]*dNxi[ii];
      J21 += x[ii]*dNeta[ii];
      J22 += y[ii]*dNeta[ii];
      J23 += z[ii]*dNeta[ii];
      J31 += x[ii]*dNzeta[ii];
      J32 += y[ii]*dNzeta[ii];
      J33 += z[ii]*dNzeta[ii];
    }
    vol = (J11*J22*J33 + J21*J32*J13 + J31*J12*J23 
           - J31*J22*J13 - J11*J32*J23 - J21*J12*J33)*8.0;
  } else {

    if (integration_rule == IntegrationRule::EIGHT_POINT) {
      Array2_D gaussPoints;
      gaussPoints[1] = 1.0/std::sqrt(3.0);
      gaussPoints[0] = -gaussPoints[1];

      Array8_D dNxi, dNeta, dNzeta;
      for (auto xi : gaussPoints) {
        for (auto eta : gaussPoints) {
          for (auto zeta : gaussPoints) {

            naturalDerivativesHex8(xi, eta, zeta, dNxi, dNeta, dNzeta);

            double J11 = 0.0; double J12 = 0.0; double J13 = 0.0;
            double J21 = 0.0; double J22 = 0.0; double J23 = 0.0;
            double J31 = 0.0; double J32 = 0.0; double J33 = 0.0;
            for (int ii = 0; ii < 8; ++ii) {
              J11 += x[ii]*dNxi[ii];
              J12 += y[ii]*dNxi[ii];
              J13 += z[ii]*dNxi[ii];
              J21 += x[ii]*dNeta[ii];
              J22 += y[ii]*dNeta[ii];
              J23 += z[ii]*dNeta[ii];
              J31 += x[ii]*dNzeta[ii];
              J32 += y[ii]*dNzeta[ii];
              J33 += z[ii]*dNzeta[ii];
            }
            vol += (J11*J22*J33 + J21*J32*J13 + J31*J12*J23 
                    - J31*J22*J13 - J11*J32*J23 - J21*J12*J33);
          }
        }
      }
    }
  }
  return vol;
}

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
ParticleInterpolator::naturalDerivativesHex8(double xi, double eta, double zeta,
                                             Array8_D& dNxi, Array8_D& dNeta, 
                                             Array8_D& dNzeta) const
{
  dNxi[0] = -(1.0-eta)*(1.0-zeta)*0.125;
  dNxi[1] = (1.0-eta)*(1.0-zeta)*0.125;
  dNxi[2] = (1.0+eta)*(1.0-zeta)*0.125;
  dNxi[3] = -(1.0+eta)*(1.0-zeta)*0.125;
  dNxi[4] = -(1.0-eta)*(1.0+zeta)*0.125;
  dNxi[5] = (1.0-eta)*(1.0+zeta)*0.125;
  dNxi[6] = (1.0+eta)*(1.0+zeta)*0.125;
  dNxi[7] = -(1.0+eta)*(1.0+zeta)*0.125;
  dNeta[0] = -(1.0-xi)*(1.0-zeta) *0.125;
  dNeta[1] = -(1.0+xi)*(1.0-zeta) *0.125;
  dNeta[2] = (1.0+xi)*(1.0-zeta) *0.125;
  dNeta[3] = (1.0-xi)*(1.0-zeta) *0.125;
  dNeta[4] = -(1.0-xi)*(1.0+zeta) *0.125;
  dNeta[5] = -(1.0+xi)*(1.0+zeta) *0.125;
  dNeta[6] = (1.0+xi)*(1.0+zeta) *0.125;
  dNeta[7] = (1.0-xi)*(1.0+zeta) *0.125;
  dNzeta[0] = -(1.0-xi)*(1.0-eta) *0.125;
  dNzeta[1] = -(1.0+xi)*(1.0-eta) *0.125;
  dNzeta[2] = -(1.0+xi)*(1.0+eta) *0.125;
  dNzeta[3] = -(1.0-xi)*(1.0+eta) *0.125;
  dNzeta[4] = (1.0-xi)*(1.0-eta) *0.125;
  dNzeta[5] = (1.0+xi)*(1.0-eta) *0.125;
  dNzeta[6] = (1.0+xi)*(1.0+eta) *0.125;
  dNzeta[7] = (1.0-xi)*(1.0+eta) *0.125;
}

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
ParticleInterpolator::lineGaussQuadratureWeights(int numIntegrationPts, int index, 
                                                 double& xi, double& weight) const
{
  if (numIntegrationPts <= 1) {
    xi = 0.0;
    weight = 2.0;
  } else if (numIntegrationPts == 2) {
    if (index == 1) {
      xi = -1.0/std::sqrt(3.0);
    } else {
      xi = 1.0/std::sqrt(3.0);
    }
    weight = 1.0;
  } else if (numIntegrationPts == 3) {
    if (index == 1) {
      xi = -std::sqrt(0.6);
      weight = 5.0/9.0;
    } else if (index == 2) {
      xi = 0.0;
      weight = 8.0/9.0;
    } else {
      xi = std::sqrt(0.6);
      weight = 5.0/9.0;
    }
  } else {
    if (index == 1) {
      xi = -0.861136311594053;
      weight = 0.347854845137454;
    } else if (index == 2) {
      xi = -0.339981043584856;
      weight = 0.652145154862546;
    } else if (index == 3) {
      xi = 0.339981043584856;
      weight = 0.652145154862546;
    } else {
      xi = 0.861136311594053;
      weight = 0.347854845137454;
    }
  }
}


