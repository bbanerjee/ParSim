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

#include <CCA/Components/MPM/Crack/JIntegralCalculator.h>

// Assuming these are globally available or passed in the constructor
extern Uintah::MPMLabel*
  lb; // Example, better to pass via JIntegralCalculator constructor
extern Uintah::Crack* d_crack_instance; // Access to Crack class members

namespace Vaango {

// NodalSolutionData Constructor Implementation
NodalSolutionData::NodalSolutionData(Uintah::DataWarehouse* new_dw,
                                     const Uintah::MPMLabel* lb,
                                     int dwi,
                                     const Uintah::Patch* patch,
                                     int NGC)
{
  // Need to handle potential errors if data is not found, or assert
  auto gac = Uintah::Ghost::AroundCells;
  new_dw->get(gMass, lb->gMassLabel, dwi, patch, gac, NGC);
  new_dw->get(Gmass, lb->GMassLabel, dwi, patch, gac, NGC);
  new_dw->get(GnumPatls, lb->GNumPatlsLabel, dwi, patch, gac, NGC);
  new_dw->get(gdisp, lb->gDisplacementLabel, dwi, patch, gac, NGC);
  new_dw->get(Gdisp, lb->GDisplacementLabel, dwi, patch, gac, NGC);
  new_dw->get(ggridStress, lb->gGridStressLabel, dwi, patch, gac, NGC);
  new_dw->get(GgridStress, lb->GGridStressLabel, dwi, patch, gac, NGC);
  new_dw->get(gdispGrads, lb->gDispGradsLabel, dwi, patch, gac, NGC);
  new_dw->get(GdispGrads, lb->GDispGradsLabel, dwi, patch, gac, NGC);
  new_dw->get(gW, lb->gStrainEnergyDensityLabel, dwi, patch, gac, NGC);
  new_dw->get(GW, lb->GStrainEnergyDensityLabel, dwi, patch, gac, NGC);
  new_dw->get(gK, lb->gKineticEnergyDensityLabel, dwi, patch, gac, NGC);
  new_dw->get(GK, lb->GKineticEnergyDensityLabel, dwi, patch, gac, NGC);
  new_dw->get(gacc, lb->gAccelerationLabel, dwi, patch, gac, NGC);
  new_dw->get(Gacc, lb->GAccelerationLabel, dwi, patch, gac, NGC);
  new_dw->get(gvel, lb->gVelocityLabel, dwi, patch, gac, NGC);
  new_dw->get(Gvel, lb->GVelocityLabel, dwi, patch, gac, NGC);
  new_dw->get(gvelGrads, lb->gVelGradsLabel, dwi, patch, gac, NGC);
  new_dw->get(GvelGrads, lb->GVelGradsLabel, dwi, patch, gac, NGC);
}

JIntegralCalculator::JIntegralCalculator(
  const Uintah::MPMLabel* lb_ptr,
  const Uintah::ParticleInterpolator* interpolator_ptr,
  Uintah::MaterialManager* mat_manager_ptr,
  double rJ_init,
  bool useVolIntegral,
  double friction_coeff,
  int cod_opt,
  double rdadx_val,
  int stress_state,
  int n8or27_val)
  : d_lb(lb_ptr)
  , d_interpolator(interpolator_ptr)
  , d_mat_manager(mat_manager_ptr)
  , d_rJ_initial(rJ_init)
  , d_useVolumeIntegral(useVolIntegral)
  , d_cmu(friction_coeff)
  , d_CODOption(static_cast<CODOption>(cod_opt))
  , d_rdadx(rdadx_val)
  , d_stressState(stress_state)
  , d_n8or27(n8or27_val)
{
}

std::expected<std::pair<Uintah::Vector, Uintah::Vector>, JCalculationError>
JIntegralCalculator::calculateJAndKForNode(
  int mat_id,
  int crack_node_idx,
  int crack_front_node_id,
  double crack_prop_velocity,
  double max_dx_cell,
  const NodalSolutionData& nodal_data,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
  const Uintah::Point& crack_node_position,
  const Uintah::Vector& v1,
  const Uintah::Vector& v2,
  const Uintah::Vector& v3)
{
  // Step 1: Define crack-front local coordinates
  int segs[2];
  d_crack_instance->FindSegsFromNode(
    mat_id,
    crack_front_node_id,
    segs); // Assuming FindSegsFromNode is accessible
  CrackSegmentPosition seg_pos =
    getCrackSegmentPosition(mat_id, crack_front_node_id, segs);

  Uintah::Point neighbor_pos;
  if (seg_pos == CrackSegmentPosition::SingleSegment) {
    int neighbor_node_id =
      (segs[0] < 0) ? d_crack_instance->d_cfSegNodes[mat_id][2 * segs[1]]
                    : d_crack_instance->d_cfSegNodes[mat_id][2 * segs[0] + 1];
    neighbor_pos = d_crack_instance->d_cx[mat_id][neighbor_node_id];
  } else if (seg_pos == CrackSegmentPosition::LeftEdge) {
    int neighbor_node_id =
      d_crack_instance
        ->d_cfSegNodes[mat_id][2 * segs[1]]; // Only right segment exists,
                                             // neighbor is its start node
    neighbor_pos = d_crack_instance->d_cx[mat_id][neighbor_node_id];
  } else if (seg_pos == CrackSegmentPosition::RightEdge) {
    int neighbor_node_id =
      d_crack_instance
        ->d_cfSegNodes[mat_id][2 * segs[0] + 1]; // Only left segment exists,
                                                 // neighbor is its end node
    neighbor_pos = d_crack_instance->d_cx[mat_id][neighbor_node_id];
  }
  // For middle, neighbor_pos is not directly used for origin,
  // crack_node_position is enough.

  Uintah::Point origin =
    determineOrigin(seg_pos, crack_node_position, neighbor_pos);
  Uintah::Matrix3 T_global_to_local(
    v1.x(), v1.y(), v1.z(), v2.x(), v2.y(), v2.z(), v3.x(), v3.y(), v3.z());
  Uintah::Matrix3 TT_local_to_global = T_global_to_local.Transpose();

  // Step 2 & 3: Find intersection (crossPt) between J-integral contour and
  // crack plane
  double rJ = d_rJ_initial;
  auto intersection_result =
    findJContourIntersection(mat_id, rJ, origin, v1, v2, v3);
  if (!intersection_result) {
    return std::unexpected(intersection_result.error());
  }
  Uintah::Point crossPt = intersection_result.value();

  double xc = crossPt.x(), yc = crossPt.y(), zc = crossPt.z();
  double xcprime = T_global_to_local(0, 0) * (xc - origin.x()) +
                   T_global_to_local(0, 1) * (yc - origin.y()) +
                   T_global_to_local(0, 2) * (zc - origin.z());
  double ycprime = T_global_to_local(1, 0) * (xc - origin.x()) +
                   T_global_to_local(1, 1) * (yc - origin.y()) +
                   T_global_to_local(1, 2) * (zc - origin.z());
  double scprime = std::sqrt(xcprime * xcprime + ycprime * ycprime);

  // Step 4 & 5: Set integral points and evaluate solutions
  int n_segments = 16; // Number of integration points on the contour
  std::vector<JContourIntegrationPointData> contour_points =
    evaluateSolutionsAtIntegrationPoints(origin,
                                         T_global_to_local,
                                         crossPt,
                                         rJ,
                                         pSize,
                                         pDefGrad,
                                         nodal_data,
                                         n_segments);

  // Get the (Sc, Uc) at crossPt (which is contour_points[0])
  Uintah::Vector Uca_at_cross, Ucb_at_cross;
  Uintah::Matrix3 Sca_at_cross, Scb_at_cross;

  // This part requires re-interpolation specific to the crossPt, as it uses the
  // original Sca/Scb/Uca/Ucb logic. The original code interpolates only once at
  // j=0 for (Sc, Uc) at crossPt. We need to re-create that logic or pass it
  // from evaluateSolutionsAtIntegrationPoints Let's assume for now that the
  // first element of contour_points has these properly filled (A redesign of
  // `evaluateSolutionsAtIntegrationPoints` could return these special values)
  Uintah::Vector temp_Uca(0., 0., 0.), temp_Ucb(0., 0., 0.);
  Uintah::Matrix3 temp_Sca(0.), temp_Scb(0.);
  std::vector<Uintah::IntVector> ni_cross(d_interpolator->size());
  std::vector<double> S_cross(d_interpolator->size());
  // Assuming pSize[0] and pDefGrad[0] are placeholder for particle size/defgrad
  // if needed at point
  d_interpolator->findCellAndWeights(
    crossPt, ni_cross, S_cross, pSize[0], pDefGrad[0]);

  double sumS = 0.;
  for (int k = 0; k < d_n8or27; ++k) {
    if (nodal_data.GnumPatls[ni_cross[k]] !=
        0) { // If it's a crack-affected node
      temp_Sca += nodal_data.ggridStress[ni_cross[k]] * S_cross[k];
      temp_Scb += nodal_data.GgridStress[ni_cross[k]] * S_cross[k];
      temp_Uca += nodal_data.gdisp[ni_cross[k]] * S_cross[k];
      temp_Ucb += nodal_data.Gdisp[ni_cross[k]] * S_cross[k];
      sumS += S_cross[k];
    }
  }
  if (sumS > 0) {
    Sca_at_cross = temp_Sca / sumS;
    Scb_at_cross = temp_Scb / sumS;
    Uca_at_cross = temp_Uca / sumS;
    Ucb_at_cross = temp_Ucb / sumS;
  }

  // Step 6: Transform the solutions to crack-front local coordinates (already
  // done in evaluateSolutionsAtIntegrationPoints) Now, transform Sca_at_cross,
  // Scb_at_cross, Uca_at_cross, Ucb_at_cross
  Uintah::Matrix3 sca_at_cross = Uintah::Matrix3(0.);
  Uintah::Matrix3 scb_at_cross = Uintah::Matrix3(0.);
  for (int i1 = 0; i1 < 3; ++i1) {
    for (int j1 = 0; j1 < 3; ++j1) {
      for (int i2 = 0; i2 < 3; ++i2) {
        for (int j2 = 0; j2 < 3; ++j2) {
          sca_at_cross(i1, j1) += T_global_to_local(i1, i2) *
                                  T_global_to_local(j1, j2) *
                                  Sca_at_cross(i2, j2);
          scb_at_cross(i1, j1) += T_global_to_local(i1, i2) *
                                  T_global_to_local(j1, j2) *
                                  Scb_at_cross(i2, j2);
        }
      }
    }
  }

  Uintah::Vector uca_at_cross, ucb_at_cross;
  transformVectorsToLocal(Uca_at_cross, Ucb_at_cross, T_global_to_local);
  uca_at_cross = Uca_at_cross; // Now holds the transformed value
  ucb_at_cross = Ucb_at_cross;

  // Step 7 & 8: Compute contour integral
  auto [Jx1, Jy1] = computeContourIntegral(
    rJ, n_segments, xcprime, ycprime, scprime, contour_points);

  // Step 9: Area integral (optional)
  double Jx2 = 0., Jy2 = 0.;
  if (d_useVolumeIntegral || Jx1 < 0.0)
    [[unlikely]] { // Unlikely hint for compiler
    auto [area_Jx2, area_Jy2] = computeAreaIntegral(mat_id,
                                                    rJ,
                                                    max_dx_cell,
                                                    origin,
                                                    TT_local_to_global,
                                                    pSize,
                                                    pDefGrad,
                                                    nodal_data);
    Jx2                       = area_Jx2;
    Jy2                       = area_Jy2;
  }

  // Step 10: Contribution of friction to energy release rate
  double fricWork = computeFrictionalWork(
    sca_at_cross, scb_at_cross, uca_at_cross - ucb_at_cross);

  // Step 11: J-integral vector
  Uintah::Vector cfJ_vector(Jx1 + Jx2 - fricWork, Jy1 + Jy2, 0.);

  // Step 12: Convert J-integral into stress intensity (K)
  Uintah::MPMMaterial* mpm_matl = static_cast<Uintah::MPMMaterial*>(
    d_mat_manager->getMaterial("MPM", mat_id));
  Uintah::ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

  Uintah::Vector D = calculateCOD(mat_id,
                                  d_CODOption,
                                  max_dx_cell,
                                  rJ,
                                  origin,
                                  T_global_to_local,
                                  nodal_data,
                                  pSize,
                                  pDefGrad);

  Uintah::Vector SIF;
  cm->convertJToK(
    mpm_matl, d_stressState, cfJ_vector, crack_prop_velocity, D, SIF);
  Uintah::Vector cfK_vector = SIF;

  return std::pair<Uintah::Vector, Uintah::Vector>{ cfJ_vector, cfK_vector };
}

// --- Private Helper Methods Implementations ---

CrackSegmentPosition
JIntegralCalculator::getCrackSegmentPosition(int mat_id,
                                             int crack_node_id,
                                             int (&segs)[2]) const
{
  d_crack_instance->FindSegsFromNode(
    mat_id,
    crack_node_id,
    segs); // Assume FindSegsFromNode is a member of Crack
  short singleSeg     = Uintah::NO;
  int neighbor_id     = -1;
  int segsNeighbor[2] = { -1, -1 };

  if (segs[1] < 0) { // Assuming R (right) is index 1, L (left) is index 0
    neighbor_id =
      d_crack_instance
        ->d_cfSegNodes[mat_id][2 * segs[0] + 1]; // End node of left segment
    d_crack_instance->FindSegsFromNode(mat_id, neighbor_id, segsNeighbor);
    if (segsNeighbor[0] < 0) { // If neighbor also only has one segment (left),
                               // it's a single crack
      singleSeg = Uintah::YES;
    }
  }
  if (segs[0] < 0) { // If L is invalid
    neighbor_id =
      d_crack_instance
        ->d_cfSegNodes[mat_id][2 * segs[1]]; // Start node of right segment
    d_crack_instance->FindSegsFromNode(mat_id, neighbor_id, segsNeighbor);
    if (segsNeighbor[1] < 0) { // If neighbor also only has one segment (right),
                               // it's a single crack
      singleSeg = Uintah::YES;
    }
  }

  if (singleSeg == Uintah::YES) {
    return CrackSegmentPosition::SingleSegment;
  }
  if (segs[1] < 0) {
    return CrackSegmentPosition::RightEdge; // Only left segment exists for this
                                            // node
  }
  if (segs[0] < 0) {
    return CrackSegmentPosition::LeftEdge; // Only right segment exists for this
                                           // node
  }
  return CrackSegmentPosition::Middle;
}

Uintah::Point
JIntegralCalculator::determineOrigin(CrackSegmentPosition pos,
                                     const Uintah::Point& node_pos,
                                     const Uintah::Point& neighbor_pos) const
{
  switch (pos) {
    case CrackSegmentPosition::SingleSegment:
      return node_pos + (neighbor_pos - node_pos) / 2.0;
    case CrackSegmentPosition::LeftEdge:  // Node is the end of the crack, shift
                                          // to neighbor
    case CrackSegmentPosition::RightEdge: // Node is the start of the crack,
                                          // shift to neighbor
      return neighbor_pos;
    case CrackSegmentPosition::Middle:
    default:
      return node_pos;
  }
}

std::expected<Uintah::Point, JCalculationError>
JIntegralCalculator::findJContourIntersection(int mat_id,
                                              double& current_rJ,
                                              const Uintah::Point& origin,
                                              const Uintah::Vector& v1,
                                              const Uintah::Vector& v2,
                                              const Uintah::Vector& v3) const
{
  double A[14];
  d_crack_instance->FindJIntegralPath(
    origin, v1, v2, v3, A); // Assuming this is a member of Crack

  Uintah::Point crossPt;
  double original_rJ = current_rJ;
  while (!d_crack_instance->FindIntersectionJPathAndCrackPlane(
    mat_id, current_rJ, A, crossPt)) {
    current_rJ *= 0.7;
    if (current_rJ / original_rJ < 0.01) [[unlikely]] {
      // Log error instead of exit
      std::cerr
        << "Error: J-integral radius (rJ) has been decreased 100 times "
        << "before finding the intersection between J-contour and crack plane. "
        << "J-integral calculation failed for this node.\n";
      return std::unexpected(JCalculationError::JIntegralRadiusTooSmall);
    }
  }
  return crossPt;
}

std::vector<JContourIntegrationPointData>
JIntegralCalculator::evaluateSolutionsAtIntegrationPoints(
  const Uintah::Point& origin,
  const Uintah::Matrix3& T_global_to_local,
  const Uintah::Point& cross_pt,
  double rJ,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
  const NodalSolutionData& nodal_data,
  int n_segments) const
{
  std::vector<JContourIntegrationPointData> points_data(n_segments + 1);
  double xc = cross_pt.x(), yc = cross_pt.y(), zc = cross_pt.z();
  double xcprime = T_global_to_local(0, 0) * (xc - origin.x()) +
                   T_global_to_local(0, 1) * (yc - origin.y()) +
                   T_global_to_local(0, 2) * (zc - origin.z());
  double ycprime = T_global_to_local(1, 0) * (xc - origin.x()) +
                   T_global_to_local(1, 1) * (yc - origin.y()) +
                   T_global_to_local(1, 2) * (zc - origin.z());
  double scprime  = std::sqrt(xcprime * xcprime + ycprime * ycprime);
  const double PI = M_PI; // Assuming M_PI from FastMath.h or <cmath>

  for (int j = 0; j <= n_segments; ++j) {
    double angle =
      2 * PI * static_cast<double>(j) / static_cast<double>(n_segments);
    double cosTheta =
      (xcprime * std::cos(angle) - ycprime * std::sin(angle)) / scprime;
    double sinTheta =
      (ycprime * std::cos(angle) + xcprime * std::sin(angle)) / scprime;

    double xprime = rJ * cosTheta;
    double yprime = rJ * sinTheta;

    // Convert back to global coordinates for interpolation
    Uintah::Point X_global = origin + (T_global_to_local.Transpose() *
                                       Uintah::Vector(xprime, yprime, 0.0))
                                        .asPoint();
    points_data[j].globalCoords = X_global;

    std::vector<Uintah::IntVector> ni(d_interpolator->size());
    std::vector<double> S(d_interpolator->size());
    // For pSize and pDefGrad, we might need a specific particle for each
    // integration point, but original code uses pSize[j] and pDefGrad[j],
    // implying a particle associated with each point. If there's only one pset,
    // then pSize[0] / pDefGrad[0] is used universally, as shown in the original
    // code's interpolation loop. Let's assume pSize[0] and pDefGrad[0] for
    // simplicity if one particle subset.
    d_interpolator->findCellAndWeights(X_global, ni, S, pSize[0], pDefGrad[0]);

    // Accumulate weighted sums
    for (int k = 0; k < d_n8or27; ++k) {
      bool is_below_crack = (nodal_data.GnumPatls[ni[k]] != 0 &&
                             j < n_segments / 2); // Original logic
      if (is_below_crack) {
        points_data[j].strainEnergyDensity += nodal_data.GW[ni[k]] * S[k];
        points_data[j].kineticEnergyDensity += nodal_data.GK[ni[k]] * S[k];
        points_data[j].globalStress += nodal_data.GgridStress[ni[k]] * S[k];
        points_data[j].globalDispGrads += nodal_data.GdispGrads[ni[k]] * S[k];
      } else {
        points_data[j].strainEnergyDensity += nodal_data.gW[ni[k]] * S[k];
        points_data[j].kineticEnergyDensity += nodal_data.gK[ni[k]] * S[k];
        points_data[j].globalStress += nodal_data.ggridStress[ni[k]] * S[k];
        points_data[j].globalDispGrads += nodal_data.gDispGrads[ni[k]] * S[k];
      }
    }
  }

  // Now transform all global tensors to local for all points
  transformTensorsToLocal(points_data, T_global_to_local);

  return points_data;
}

void
JIntegralCalculator::transformTensorsToLocal(
  std::vector<JContourIntegrationPointData>& points_data,
  const Uintah::Matrix3& T_global_to_local) const
{
  const Uintah::Matrix3& T = T_global_to_local;
  for (auto& point_data : points_data) {
    Uintah::Matrix3 st_val(0.), dg_val(0.);
    for (int i1 = 0; i1 < 3; ++i1) {
      for (int j1 = 0; j1 < 3; ++j1) {
        for (int i2 = 0; i2 < 3; ++i2) {
          for (int j2 = 0; j2 < 3; ++j2) {
            st_val(i1, j1) +=
              T(i1, i2) * T(j1, j2) * point_data.globalStress(i2, j2);
            dg_val(i1, j1) +=
              T(i1, i2) * T(j1, j2) * point_data.globalDispGrads(i2, j2);
          }
        }
      }
    }
    point_data.localStress    = st_val;
    point_data.localDispGrads = dg_val;
  }
}

void
JIntegralCalculator::transformVectorsToLocal(
  Uintah::Vector& vec_a,
  Uintah::Vector& vec_b,
  const Uintah::Matrix3& T_global_to_local) const
{
  const Uintah::Matrix3& T = T_global_to_local;
  vec_a                    = T * vec_a;
  vec_b                    = T * vec_b;
}

std::pair<double, double>
JIntegralCalculator::computeContourIntegral(
  double rJ,
  int n_segments,
  double xcprime,
  double ycprime,
  double scprime,
  const std::vector<JContourIntegrationPointData>& points_data) const
{
  double Jx1 = 0., Jy1 = 0.;
  const double PI = M_PI;

  for (int j = 0; j < n_segments; ++j) {
    double angle_j =
      2 * PI * static_cast<double>(j) / static_cast<double>(n_segments);
    double angle_j1 =
      2 * PI * static_cast<double>(j + 1) / static_cast<double>(n_segments);

    double cosTheta_j =
      (xcprime * std::cos(angle_j) - ycprime * std::sin(angle_j)) / scprime;
    double sinTheta_j =
      (ycprime * std::cos(angle_j) + xcprime * std::sin(angle_j)) / scprime;
    double t1_j = points_data[j].localStress(0, 0) * cosTheta_j +
                  points_data[j].localStress(0, 1) * sinTheta_j;
    double t2_j = points_data[j].localStress(1, 0) * cosTheta_j +
                  points_data[j].localStress(1, 1) * sinTheta_j;
    Uintah::Vector t123_j = Uintah::Vector(t1_j, t2_j, 0.); // Plane state
    Uintah::Vector dgx_j  = Uintah::Vector(points_data[j].localDispGrads(0, 0),
                                          points_data[j].localDispGrads(1, 0),
                                          points_data[j].localDispGrads(2, 0));
    Uintah::Vector dgy_j  = Uintah::Vector(points_data[j].localDispGrads(0, 1),
                                          points_data[j].localDispGrads(1, 1),
                                          points_data[j].localDispGrads(2, 1));
    double f1ForJx_j      = (points_data[j].strainEnergyDensity +
                        points_data[j].kineticEnergyDensity) *
                         cosTheta_j -
                       Dot(t123_j, dgx_j);
    double f1ForJy_j = (points_data[j].strainEnergyDensity +
                        points_data[j].kineticEnergyDensity) *
                         sinTheta_j -
                       Dot(t123_j, dgy_j);

    double cosTheta_j1 =
      (xcprime * std::cos(angle_j1) - ycprime * std::sin(angle_j1)) / scprime;
    double sinTheta_j1 =
      (ycprime * std::cos(angle_j1) + xcprime * std::sin(angle_j1)) / scprime;
    double t1_j1 = points_data[j + 1].localStress(0, 0) * cosTheta_j1 +
                   points_data[j + 1].localStress(0, 1) * sinTheta_j1;
    double t2_j1 = points_data[j + 1].localStress(1, 0) * cosTheta_j1 +
                   points_data[j + 1].localStress(1, 1) * sinTheta_j1;
    Uintah::Vector t123_j1 = Uintah::Vector(t1_j1, t2_j1, 0.); // Plane state
    Uintah::Vector dgx_j1 =
      Uintah::Vector(points_data[j + 1].localDispGrads(0, 0),
                     points_data[j + 1].localDispGrads(1, 0),
                     points_data[j + 1].localDispGrads(2, 0));
    Uintah::Vector dgy_j1 =
      Uintah::Vector(points_data[j + 1].localDispGrads(0, 1),
                     points_data[j + 1].localDispGrads(1, 1),
                     points_data[j + 1].localDispGrads(2, 1));
    double f1ForJx_j1 = (points_data[j + 1].strainEnergyDensity +
                         points_data[j + 1].kineticEnergyDensity) *
                          cosTheta_j1 -
                        Dot(t123_j1, dgx_j1);
    double f1ForJy_j1 = (points_data[j + 1].strainEnergyDensity +
                         points_data[j + 1].kineticEnergyDensity) *
                          sinTheta_j1 -
                        Dot(t123_j1, dgy_j1);

    Jx1 += f1ForJx_j + f1ForJx_j1;
    Jy1 += f1ForJy_j + f1ForJy_j1;
  }
  Jx1 *= rJ * PI / n_segments;
  Jy1 *= rJ * PI / n_segments;

  return { Jx1, Jy1 };
}

std::pair<double, double>
JIntegralCalculator::computeAreaIntegral(
  int mat_id,
  double rJ,
  double dx_max,
  const Uintah::Point& origin,
  const Uintah::Matrix3& TT_local_to_global,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad,
  const NodalSolutionData& nodal_data) const
{
  int nc = static_cast<int>(rJ / dx_max);
  if (rJ / dx_max - nc >= 0.5) {
    nc++;
  }
  if (nc < 2) {
    nc = 2; // Cell number J-path away from the origin
  }

  std::vector<double> c(4 * nc);
  for (int j = 0; j < 4 * nc; ++j) {
    c[j] = static_cast<double>(-4 * nc + 1 + 2 * j) / (4.0 * nc) * rJ;
  }

  std::vector<Uintah::Point> x_local_points;
  std::vector<Uintah::Point> X_global_points;

  for (double c_i : c) {
    for (double c_j : c) {
      Uintah::Point pq(c_i, c_j, 0.);
      if (pq.asVector().length() < rJ) {
        x_local_points.push_back(pq);
        X_global_points.push_back(origin +
                                  (TT_local_to_global * pq.asVector()));
      }
    }
  }

  if (X_global_points.empty()) {
    return { 0.0,
             0.0 }; // No integration points in area, return zero contribution
  }

  std::vector<Uintah::Vector> acceleration_local(X_global_points.size());
  std::vector<Uintah::Vector> velocity_local(X_global_points.size());
  std::vector<Uintah::Matrix3> disp_grad_local(X_global_points.size());
  std::vector<Uintah::Matrix3> vel_grad_local(X_global_points.size());

  for (size_t j = 0; j < X_global_points.size(); ++j) {
    Uintah::Vector acceleration_global(0.), velocity_global(0.);
    Uintah::Matrix3 disp_grad_global(0.), vel_grad_global(0.);

    std::vector<Uintah::IntVector> ni(d_interpolator->size());
    std::vector<double> S(d_interpolator->size());
    d_interpolator->findCellAndWeights(
      X_global_points[j], ni, S, pSize[0], pDefGrad[0]);

    for (int k = 0; k < d_n8or27; ++k) {
      // Original logic for below/above crack based on local y-coordinate
      // This assumes a straight crack aligned with the local x-axis
      bool is_below_crack =
        (nodal_data.GnumPatls[ni[k]] != 0 && x_local_points[j].y() < 0.);

      if (is_below_crack) {
        acceleration_global += nodal_data.Gacc[ni[k]] * S[k];
        velocity_global += nodal_data.Gvel[ni[k]] * S[k];
        disp_grad_global += nodal_data.GdispGrads[ni[k]] * S[k];
        vel_grad_global += nodal_data.GvelGrads[ni[k]] * S[k];
      } else {
        acceleration_global += nodal_data.gacc[ni[k]] * S[k];
        velocity_global += nodal_data.gvel[ni[k]] * S[k];
        disp_grad_global += nodal_data.gDispGrads[ni[k]] * S[k];
        vel_grad_global += nodal_data.gVelGrads[ni[k]] * S[k];
      }
    }

    acceleration_local[j] =
      TT_local_to_global.Transpose() *
      acceleration_global; // T_global_to_local * acceleration_global
    velocity_local[j] = TT_local_to_global.Transpose() * velocity_global;

    Uintah::Matrix3 disp_grad_temp(0.), vel_grad_temp(0.);
    const Uintah::Matrix3& T =
      TT_local_to_global.Transpose(); // T_global_to_local
    for (int i1 = 0; i1 < 3; ++i1) {
      for (int j1 = 0; j1 < 3; ++j1) {
        for (int i2 = 0; i2 < 3; ++i2) {
          for (int j2 = 0; j2 < 3; ++j2) {
            disp_grad_temp(i1, j1) +=
              T(i1, i2) * T(j1, j2) * disp_grad_global(i2, j2);
            vel_grad_temp(i1, j1) +=
              T(i1, i2) * T(j1, j2) * vel_grad_global(i2, j2);
          }
        }
      }
    }
    disp_grad_local[j] = disp_grad_temp;
    vel_grad_local[j]  = vel_grad_temp;
  }

  double f2ForJx = 0., f2ForJy = 0.;
    double rho = d_mat_manager->getMaterial("MPM", mat_id)->get         ConstitutiveModel()->get)InitialDensity(); // Corrected way to get density

    for (size_t j = 0; j < X_global_points.size(); ++j) {
      Uintah::Vector dgx = Uintah::Vector(disp_grad_local[j](0, 0),
                                          disp_grad_local[j](1, 0),
                                          0.); // Zero z for plane state
      Uintah::Vector dgy = Uintah::Vector(disp_grad_local[j](0, 1),
                                          disp_grad_local[j](1, 1),
                                          0.); // Zero z for plane state
      Uintah::Vector vgx = Uintah::Vector(vel_grad_local[j](0, 0),
                                          vel_grad_local[j](1, 0),
                                          0.); // Zero z for plane state
      Uintah::Vector vgy = Uintah::Vector(vel_grad_local[j](0, 1),
                                          vel_grad_local[j](1, 1),
                                          0.); // Zero z for plane state

      f2ForJx +=
        rho * (Dot(acceleration_local[j], dgx) - Dot(velocity_local[j], vgx));
      f2ForJy +=
        rho * (Dot(acceleration_local[j], dgy) - Dot(velocity_local[j], vgy));
    }

    double Jarea = M_PI * rJ * rJ;
    double Jx2_val =
      f2ForJx / static_cast<double>(X_global_points.size()) * Jarea;
    double Jy2_val =
      f2ForJy / static_cast<double>(X_global_points.size()) * Jarea;

    return { Jx2_val, Jy2_val };
}

double
JIntegralCalculator::computeFrictionalWork(const Uintah::Matrix3& sca,
                                           const Uintah::Matrix3& scb,
                                           const Uintah::Vector& uc) const
{
  if (d_cmu == 0.) {
    return 0.; // No friction
  }

  double ta          = std::sqrt(sca(1, 0) * sca(1, 0) + sca(1, 2) * sca(1, 2));
  double tb          = std::sqrt(scb(1, 0) * scb(1, 0) + scb(1, 2) * scb(1, 2));
  Uintah::Matrix3 tc = sca;
  if (tb > ta) {
    tc = scb;
  }

  double fricWork = std::fabs(tc(1, 0) * uc.x()) + std::fabs(tc(1, 2) * uc.z());
  if (tc(1, 1) < 0. &&
      uc.y() < 0.) { // Only if normal stress is compressive and crack closes
    fricWork += std::fabs(tc(1, 1) * uc.y());
  }
  return fricWork;
}

Uintah::Vector
JIntegralCalculator::calculateCOD(
  int mat_id,
  CODOption option,
  double max_dx_cell,
  double rJ,
  const Uintah::Point& origin,
  const Uintah::Matrix3& T_global_to_local,
  const NodalSolutionData& nodal_data,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pSize,
  const Uintah::ParticleVariable<Uintah::Matrix3>& pDefGrad) const
{
  double d_dist = 0.0;
  if (option == CODOption::FixedDistance ||
      option == CODOption::MaxDistanceOnCrackPlane) {
    if (d_crack_instance->d_doCrackPropagationStep) {
      d_dist = (d_rdadx < 1. ? 1. : d_rdadx) * max_dx_cell;
    } else {
      d_dist = rJ / 2.0;
    }

    if (option == CODOption::MaxDistanceOnCrackPlane) {
      d_crack_instance->GetPositionToComputeCOD(
        mat_id, origin, T_global_to_local, d_dist);
    }
  } else if (option == CODOption::IntersectionPoint) {
    // Here, we need the actual crossPt. Re-calculating or passing it in.
    // For simplicity, let's assume if IntersectionPoint is chosen, 'd_dist'
    // logic above is effectively bypassed for the point calculation. The
    // original code used crossPt directly, which is already `d_rJ` away in the
    // original. We'll use a placeholder `p_d` that we'd have to get from the
    // caller if `option == IntersectionPoint`. Or, assume if option ==
    // IntersectionPoint, d_dist is not directly used for point (-d,0,0) and a
    // different point (e.g., crossPt) is intended. Let's stick to the original
    // logic of (-d,0,0) with 'd' determined by option 0 or 1.
  }

  Uintah::Point p_d =
    origin + (T_global_to_local.Transpose() * Uintah::Vector(-d_dist, 0.0, 0.0))
               .asPoint();

  // Calculate displacements at point p_d
  Uintah::Vector disp_a = Uintah::Vector(0.);
  Uintah::Vector disp_b = Uintah::Vector(0.);
  std::vector<Uintah::IntVector> ni(d_interpolator->size());
  std::vector<double> S(d_interpolator->size());

  // Use pSize[0] and pDefGrad[0] as per the original code's single pset
  // assumption for interpolation
  d_interpolator->findCellAndWeights(p_d, ni, S, pSize[0], pDefGrad[0]);
  for (int k = 0; k < d_n8or27; ++k) {
    disp_a += nodal_data.gdisp[ni[k]] * S[k];
    disp_b += nodal_data.Gdisp[ni[k]] * S[k];
  }

  // Crack opening displacements in local coodinates
  return T_global_to_local * (disp_a - disp_b);
}

} // namespace Vaango

/*
// Crack.cpp (or relevant file where Crack class is defined)
#include "JIntegralCalculator.h" // Include the new header
#include <mpi.h> // For MPI_Bcast, MPI_Comm_size, MPI_Comm_rank
#include <iostream> // For std::cerr

void
Crack::CalculateFractureParameters(const ProcessorGroup* pg,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  // Make sure MPI is initialized and d_mpi_crack_comm is valid
  if (d_mpi_crack_comm == MPI_COMM_NULL) {
    // Handle error, e.g., throw exception or return
    std::cerr << "Error: d_mpi_crack_comm is not initialized.\n";
    return;
  }

  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, lb->simulationTimeLabel);
  double time = simTimeVar;

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  int pid, patch_size;
  MPI_Comm_size(d_mpi_crack_comm, &patch_size);
  MPI_Comm_rank(d_mpi_crack_comm, &pid);
  MPI_Datatype MPI_VECTOR = fun_getTypeDescription((Vector*)0)->getMPIType();

  for (int p = 0; p < patches->size(); ++p) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double dx_max      = Max(dx.x(), dx.y(), dx.z());

    // Use std::unique_ptr for interpolator for RAII
    std::unique_ptr<Interpolator> interpolator_ptr(
      d_flag->d_interpolator->clone(patch));
    // You might need to adjust Interpolator interface if clone returns raw
    // pointer

    int numMatls = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; ++m) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

      int dwi              = matls->get(m);
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Fetch all nodal data into a structured object
      const int NGC = d_NJ + d_NGN + 1;
      FractureAnalysis::NodalSolutionData nodal_data(
        new_dw, lb, dwi, patch, NGC);

      // Particle variables
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      old_dw->get(pSize, lb->pSizeLabel, pset);
      old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

      // Allocate memories for cfSegJ and cfSegK (resize existing vectors)
      // Assuming d_cfSegJ and d_cfSegK are std::vector<std::vector<Vector>>
      int cfNodeSize = static_cast<int>(d_cfSegNodes[m].size());
      d_cfSegJ[m].resize(cfNodeSize);
      d_cfSegK[m].resize(cfNodeSize);

      if (d_calFractParametersStep || d_doCrackPropagationStep) {
        // Initialize JIntegralCalculator (once per material or outer loop)
        // Pass all necessary Crack member variables
        FractureAnalysis::JIntegralCalculator j_calculator(
          lb,
          interpolator_ptr.get(),
          d_mat_manager,
          d_rJ,
          d_useVolumeIntegral,
          d_cmu[m],
          d_CODOption,
          d_rdadx,
          d_stressState[m],
          d_n8or27);

        for (int i = 0; i < patch_size;
             ++i) { // Loop over all patches for MPI distribution
          int num_crack_front_nodes_in_patch =
            static_cast<int>(d_cfnset[m][i].size());

          // Use std::vector instead of raw pointers for automatic memory
          // management
          std::vector<Vector> cfJ_local(num_crack_front_nodes_in_patch);
          std::vector<Vector> cfK_local(num_crack_front_nodes_in_patch);

          if (pid == i) [[likely]] { // Calculate J & K by processor i
            for (int l = 0; l < num_crack_front_nodes_in_patch; ++l) {
              int idx = d_cfnset[m][i][l]; // index of this node in d_cfSegNodes
              int node_id = d_cfSegNodes[m][idx]; // actual node number

              int preIdx_local = -1; // Index in the current local cf_nset array
              int preIdx_global = d_cfSegPreIdx[m][idx]; // Global index

              // Find preIdx_local in the current patch's d_cfnset[m][i]
              // This part of the logic needs to be careful: preIdx is meant to
              // check if a crack-front node is a duplicate in the global list
              // of crack-front nodes. If `preIdx_global` matches
              // `d_cfnset[m][i][ij]`, it's an already computed node within this
              // *local patch's* calculation.
              bool already_computed_locally = false;
              if (preIdx_global >=
                  0) { // If it's a "duplicate" (already processed globally)
                for (int ij = l - 1; ij >= 0; --ij) {
                  if (preIdx_global == d_cfnset[m][i][ij]) {
                    preIdx_local =
                      ij; // Found the local index of the already computed node
                    already_computed_locally = true;
                    break;
                  }
                }
              }

              if (!already_computed_locally) [[likely]] {
                // Original code had `if (preIdx < 0)` here. This means it only
                // calculates if it's not a duplicate *in the global list*. We
                // are replicating that logic: if d_cfSegPreIdx[m][idx] is -1,
                // it's a new node. If it's >=0, it's a duplicate of a node
                // potentially computed elsewhere. The inner `for` loop `for
                // (int ij = l - 1; ij >= 0; ij--)` was to check if a duplicate
                // *within the current local patch's processing batch* had been
                // done. To truly match original, this would be: if
                // (d_cfSegPreIdx[m][idx] < 0) { ... calculate ... } else { copy
                // from d_cfSegJ[m][d_cfSegPreIdx[m][idx]] } The original
                // `preIdx < 0` check seemed to imply local duplication within
                // `d_cfnset[m][i]`. Let's follow the original intent as `if
                // (preIdx < 0)` The new 'preIdx_local < 0' logic correctly
                // means 'this specific node isn't a duplicate in this
                // processing batch'

                Uintah::Point node_pos = d_cx[m][node_id];
                Uintah::Vector v1_val  = d_cfSegV1[m][idx];
                Uintah::Vector v2_val  = d_cfSegV2[m][idx];
                Uintah::Vector v3_val  = d_cfSegV3[m][idx];
                double crack_vel       = d_cfSegVel[m][idx];

                auto result = j_calculator.calculateJAndKForNode(m,
                                                                 idx,
                                                                 node_id,
                                                                 crack_vel,
                                                                 dx_max,
                                                                 nodal_data,
                                                                 pSize,
                                                                 pDefGrad,
                                                                 node_pos,
                                                                 v1_val,
                                                                 v2_val,
                                                                 v3_val);

                if (result) {
                  cfJ_local[l] = result->first;
                  cfK_local[l] = result->second;
                } else {
                  // Handle error, e.g., set to zero or a special error value
                  std::cerr
                    << "J-integral calculation failed for node " << node_id
                    << " with error: " << static_cast<int>(result.error())
                    << "\n";
                  cfJ_local[l] = Uintah::Vector(0., 0., 0.);
                  cfK_local[l] = Uintah::Vector(0., 0., 0.);
                }
              } else {
                // This node is a duplicate within this current patch's
                // processing batch Copy from already computed local result
                cfJ_local[l] = cfJ_local[preIdx_local];
                cfK_local[l] = cfK_local[preIdx_local];
              }
            } // End of loop over nodes (l)
          } // End if(pid==i)

          // Broadcast the results calculated by rank i to all the ranks
          // Use std::vector::data() for raw pointer to underlying array
          MPI_Bcast(cfJ_local.data(),
                    num_crack_front_nodes_in_patch,
                    MPI_VECTOR,
                    i,
                    d_mpi_crack_comm);
          MPI_Bcast(cfK_local.data(),
                    num_crack_front_nodes_in_patch,
                    MPI_VECTOR,
                    i,
                    d_mpi_crack_comm);

          // Save data in d_cfSegJ and d_cfSegK
          for (int l = 0; l < num_crack_front_nodes_in_patch; ++l) {
            int idx          = d_cfnset[m][i][l];
            d_cfSegJ[m][idx] = cfJ_local[l];
            d_cfSegK[m][idx] = cfK_local[l];
          }
        } // End of loop over ranks (i)

        if (pid == 0) {
          OutputCrackFrontResults(
            m, time, delT); // Assuming this is a member function of Crack
        }
      } // End if(calFractParametersStep || doCrackPropagationStep)
    } // End of loop over matls
  }
}
*/