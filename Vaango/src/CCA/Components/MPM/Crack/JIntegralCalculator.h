/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#pragma once

#include <CCA/Components/MPM/Crack/JIntegralCalculatorDefs.h>

#include <expected> // C++23 std::expected

namespace Vaango {

class JIntegralCalculator
{
public:
  // Constructor would take references to static Crack class members or config
  // data
  JIntegralCalculator(const Uintah::MPMLabel* lb,
                      const Uintah::ParticleInterpolator* interpolator,
                      Uintah::MaterialManager* mat_manager,
                      double rJ_initial,
                      bool useVolumeIntegral,
                      double friction_coeff,
                      int cod_option,
                      double rdadx_val,
                      int stress_state,
                      int n8or27); // Pass necessary Crack member data

  // Main function to calculate J and K for a specific crack front node
  // Returns std::expected<std::pair<Vector, Vector>, JCalculationError>
  // Pair: {J_vector, K_vector}
  std::expected<std::pair<Uintah::Vector, Uintah::Vector>, JCalculationError>
  calculateJAndKForNode(
    int mat_id,
    int crack_node_idx,
    int crack_front_node_id,
    double crack_prop_velocity,
    double max_dx_cell,
    const NodalSolutionData& nodal_data,
    const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
    const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
    const Uintah::Point& crack_node_position, // d_cx[m][node]
    const Uintah::Vector& v1,
    const Uintah::Vector& v2,
    const Uintah::Vector& v3);

private:
  // Internal helper methods, potentially static or free functions if they don't
  // strictly depend on class members (but here, they often do)

  // Helper for Step 1
  CrackSegmentPosition
  getCrackSegmentPosition(int mat_id, int crack_node_id, int (&segs)[2]) const;
  Uintah::Point
  determineOrigin(CrackSegmentPosition pos,
                  const Uintah::Point& node_pos,
                  const Uintah::Point& neighbor_pos) const;

  // Helper for Step 2 & 3
  std::expected<Uintah::Point, JCalculationError>
  findJContourIntersection(int mat_id,
                           double& current_rJ,
                           const Uintah::Point& origin,
                           const Uintah::Vector& v1,
                           const Uintah::Vector& v2,
                           const Uintah::Vector& v3) const;

  // Helper for Step 4 & 5
  std::vector<JContourIntegrationPointData>
  evaluateSolutionsAtIntegrationPoints(
    const Uintah::Point& origin,
    const Uintah::Matrix3& T_global_to_local,
    const Uintah::Point& cross_pt,
    double rJ,
    const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
    const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
    const NodalSolutionData& nodal_data,
    int n_segments) const;

  // Helper for Step 6
  void
  transformTensorsToLocal(
    std::vector<JContourIntegrationPointData>& points_data,
    const Uintah::Matrix3& T_global_to_local) const;
  void
  transformVectorsToLocal(Uintah::Vector& vec_a,
                          Uintah::Vector& vec_b,
                          const Uintah::Matrix3& T_global_to_local) const;

  // Helper for Step 7 & 8
  std::pair<double, double>
  computeContourIntegral(
    double rJ,
    int n_segments,
    double xcprime,
    double ycprime,
    double scprime,
    const std::vector<JContourIntegrationPointData>& points_data) const;

  // Helper for Step 9
  std::pair<double, double>
  computeAreaIntegral(int mat_id,
                      double rJ,
                      double dx_max,
                      const Uintah::Point& origin,
                      const Uintah::Matrix3& TT_local_to_global,
                      const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
                      const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
                      const NodalSolutionData& nodal_data) const;

  // Helper for Step 10
  double
  computeFrictionalWork(const Uintah::Matrix3& sca,
                        const Uintah::Matrix3& scb,
                        const Uintah::Vector& uc) const;

  // Helper for Step 12a (COD)
  Uintah::Vector
  calculateCOD(int mat_id,
               CODOption option,
               double max_dx_cell,
               double rJ,
               const Uintah::Point& origin,
               const Uintah::Matrix3& T_global_to_local,
               const NodalSolutionData& nodal_data,
               const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
               const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad) const;

  // Assuming these are available globally or through member pointers
  const Uintah::MPMLabel* d_lb;
  const Uintah::Interpolator* d_interpolator;
  Uintah::MaterialManager* d_mat_manager;
  double d_rJ_initial;
  bool d_useVolumeIntegral;
  double d_cmu;
  CODOption d_CODOption;
  double d_rdadx;
  int d_stressState;
  int d_n8or27;

  // These would typically be accessed via member pointers to the Crack class,
  // or passed in as needed.
  // For this example, assuming they are accessible via the Crack instance
  // that creates this JIntegralCalculator.
  // std::vector<std::vector<int>> d_cfnset; // Need to be careful with these
  // std::vector<std::vector<int>> d_cfSegNodes;
  // std::vector<std::vector<int>> d_cfSegPreIdx;
  // std::vector<std::vector<Uintah::Point>> d_cx;
  // std::vector<std::vector<Uintah::Vector>> d_cfSegV1, d_cfSegV2, d_cfSegV3;
  // std::vector<std::vector<double>> d_cfSegVel;

  // Functions from original code, assumed to be member functions of Crack class
  // or passed as function pointers/lambdas
  // (These should ideally be moved into a separate helper or passed in)
  // void FindSegsFromNode(int m, int node, int (&segs)[2]) const;
  // void FindJIntegralPath(const Uintah::Point& origin, const Uintah::Vector&
  // v1, const Uintah::Vector& v2, const Uintah::Vector& v3, double (&A)[14])
  // const; bool FindIntersectionJPathAndCrackPlane(int m, double rJ, const
  // double (&A)[14], Uintah::Point& crossPt) const; void
  // GetPositionToComputeCOD(int m, const Uintah::Point& origin, const
  // Uintah::Matrix3& T, double& d) const;
};

} // namespace Vaango
