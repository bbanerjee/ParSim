#ifndef ELLIP3D_DEM_TETRAHEDRON_H
#define ELLIP3D_DEM_TETRAHEDRON_H

#include <DiscreteElements/DEMContainers.h>

#include <Eigen/Dense>
#include <Core/Math/Matrix3.h>

#include <vector>

using Matrix4x3 = Eigen::Matrix<double, 4, 3>;
using Matrix3x4 = Eigen::Matrix<double, 3, 4>;
using NodeID = std::size_t;

namespace dem {

class DEMTetrahedron
{
private:

  std::vector<NodeID> d_nodeIDs;
  std::vector<DEMParticleP> d_particles;

  Matrix4x3 d_B;
  Matrix3 d_velGrad;
  Matrix3 d_defGrad;
  Matrix3 d_defGradRate;

  REAL volume(const Vec& p0, const Vec& p1, 
              const Vec& p2, const Vec& p3) const;
  REAL signedVolume() const;
  REAL signedInitialVolume() const;

  inline
  std::vector<double> referenceBasisFunctions(double r, double s, double t) const {
    auto N1 = 1 - r - s - t; auto N2 = r; auto N3 = s; auto N4 = t;
    return {{N1, N2, N3, N4}};
  }

  inline
  std::vector<Vec> referenceBasisDerivatives() const {
    Vec dN1_dr(-1, -1, -1);
    Vec dN2_dr( 1,  0,  0);
    Vec dN3_dr( 0,  1,  0);
    Vec dN4_dr( 0,  0,  1);
    return {{dN1_dr, dN2_dr, dN3_dr, dN4_dr}};
  }

  inline
  Vec transformToXYZ(double r, double s, double t, 
                     const Vec& p0, const Vec& p1, 
                     const Vec& p2, const Vec& p3) const {
    Matrix3 A = transformationMatrix(p0, p1, p2, p3);
    Vec rst(r, s, t);
    return transformToXYZ(rst, p0, A);
  }

  inline
  Matrix3 transformationMatrix(const Vec& p0, const Vec& p1, 
                               const Vec& p2, const Vec& p3) const {
    Vec a = p1 - p0; 
    Vec b = p2 - p0; 
    Vec c = p3 - p0; 
    Matrix3 A = columnMatrix3(a, b, c);
    return A;
  }

  inline
  Vec transformToXYZ(const Vec& rst, const Vec& p0, const Matrix3& A) const {
    return A * rst + p0;
  }

  inline
  Vec transformToRST(const Vec& xyz, 
                     const Vec& p0, const Vec& p1, 
                     const Vec& p2, const Vec& p3) const {
    Matrix3 A = transformationMatrix(p0, p1, p2, p3);
    Matrix3 Ainv = A.Inverse();
    return transformToRST(xyz, p0, Ainv);
  }

  inline
  Vec transformToRST(const Vec& xyz, const Vec& p0, const Matrix3& Ainv) const {
    return Ainv * (xyz - p0);
  }

public:

  DEMTetrahedron();
  DEMTetrahedron(NodeID n1, NodeID n2, NodeID n3, NodeID n4,
                 const DEMParticlePArray& particles);

  inline
  NodeID nodeID(int index) const { return d_nodeIDs[index]; }
  inline
  const DEMParticleP& particle(int index) const { return d_particles[index]; }

  inline
  Matrix4x3 getBMatrix() const { return d_B; }
  inline
  Matrix3 getVelGrad() const { return d_velGrad; }
  inline
  Matrix3 getDefGrad() const { return d_defGrad; }
  inline
  Matrix3 getDefGradRate() const { return d_defGradRate; }

  REAL volume() const;
  REAL initialVolume() const;

  inline
  std::vector<double> basisFunctions(const Vec& point) const;

  inline
  std::vector<Vec> basisDerivatives() const;

  inline std::pair<std::vector<double>, std::vector<Vec>> 
  basisFunctionsAndDerivatives(const Vec& point) const;

  void updateBMatrix();
  void updateVelGrad();
  void updateDefGrad();

};

} // namespace dem ends

#endif
