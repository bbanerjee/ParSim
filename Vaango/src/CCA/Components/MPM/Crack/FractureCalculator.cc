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

#include <CCA/Components/MPM/Crack/FractureCalculator.h>

#include <algorithm>
#include <array>
#include <format>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

namespace Vaango {

// Implementation using C++23 features
std::expected<CrackFractureCalculator::LocalCoordinates, std::string>
CrackFractureCalculator::setupLocalCoordinates(int material_id,
                                               int node_id) const
{
  LocalCoordinates coords;

  // Find segments connected to this node
  std::array<int, 2> segments;
  FindSegsFromNode(material_id, node_id, segments.data());

  // Detect single segment crack using modern C++23 pattern matching
  bool is_single_segment = false;
  int neighbor           = -1;

  if (segments[1] < 0) {
    neighbor = d_cfSegNodes[material_id][2 * segments[0] + 1];
    std::array<int, 2> neighbor_segs = { -1, -1 };
    FindSegsFromNode(material_id, neighbor, neighbor_segs.data());
    is_single_segment = (neighbor_segs[0] < 0);
  } else if (segments[0] < 0) {
    neighbor = d_cfSegNodes[material_id][2 * segments[1]];
    std::array<int, 2> neighbor_segs = { -1, -1 };
    FindSegsFromNode(material_id, neighbor, neighbor_segs.data());
    is_single_segment = (neighbor_segs[1] < 0);
  }

  // Calculate origin using structured binding (C++23)
  if (is_single_segment) {
    auto [pt1, pt2] =
      std::make_pair(d_cx[material_id][node_id], d_cx[material_id][neighbor]);
    coords.origin = pt1 + (pt2 - pt1) / 2.0;
  } else {
    if (segments[1] < 0 || segments[0] < 0) {
      coords.origin = d_cx[material_id][neighbor];
    } else {
      coords.origin = d_cx[material_id][node_id];
    }
  }

  // Set up local coordinate system
  int idx   = /* calculate index */;
  coords.v1 = d_cfSegV1[material_id][idx];
  coords.v2 = d_cfSegV2[material_id][idx];
  coords.v3 = d_cfSegV3[material_id][idx];

  // Create transformation matrices using modern initialization
  coords.T = Matrix3{ coords.v1.x(), coords.v1.y(), coords.v1.z(),
                      coords.v2.x(), coords.v2.y(), coords.v2.z(),
                      coords.v3.x(), coords.v3.y(), coords.v3.z() };

  coords.TT = Matrix3{ coords.v1.x(), coords.v2.x(), coords.v3.x(),
                       coords.v1.y(), coords.v2.y(), coords.v3.y(),
                       coords.v1.z(), coords.v2.z(), coords.v3.z() };

  return coords;
}

std::expected<CrackFractureCalculator::JIntegralData, std::string>
CrackFractureCalculator::setupJIntegralContour(
  const LocalCoordinates& coords) const
{
  JIntegralData j_data;

  // Find J-integral path parameters
  std::array<double, 14> path_params;
  FindJIntegralPath(
    coords.origin, coords.v1, coords.v2, coords.v3, path_params.data());

  // Find intersection with crack plane
  double radius = d_rJ;
  auto intersection_result =
    findJIntegralCrackPlaneIntersection(path_params, radius);
  if (!intersection_result) {
    return std::unexpected(
      std::format("Failed to find J-integral intersection: {}",
                  intersection_result.error()));
  }

  Point cross_point = intersection_result.value();
  j_data.radius     = radius;

  // Set up integration points using ranges (C++23)
  constexpr int num_segments = 16;
  constexpr double PI        = std::numbers::pi; // C++23 math constants

  j_data.integration_points.reserve(num_segments + 1);
  j_data.strain_energy.reserve(num_segments + 1);
  j_data.kinetic_energy.reserve(num_segments + 1);
  j_data.stresses.reserve(num_segments + 1);
  j_data.displacement_grads.reserve(num_segments + 1);

  // Calculate local coordinates of intersection
  auto [x0, y0, z0] =
    std::make_tuple(coords.origin.x(), coords.origin.y(), coords.origin.z());
  double xc_prime = coords.v1.x() * (cross_point.x() - x0) +
                    coords.v1.y() * (cross_point.y() - y0) +
                    coords.v1.z() * (cross_point.z() - z0);
  double yc_prime = coords.v2.x() * (cross_point.x() - x0) +
                    coords.v2.y() * (cross_point.y() - y0) +
                    coords.v2.z() * (cross_point.z() - z0);
  double sc_prime = std::sqrt(xc_prime * xc_prime + yc_prime * yc_prime);

  // Generate integration points using ranges
  for (int j : std::views::iota(0, num_segments + 1)) {
    double angle = 2.0 * PI * static_cast<double>(j) / num_segments;
    double cos_theta =
      (xc_prime * std::cos(angle) - yc_prime * std::sin(angle)) / sc_prime;
    double sin_theta =
      (yc_prime * std::cos(angle) + xc_prime * std::sin(angle)) / sc_prime;

    // Local coordinates
    double x_prime = radius * cos_theta;
    double y_prime = radius * sin_theta;

    // Global coordinates
    double x = coords.v1.x() * x_prime + coords.v2.x() * y_prime + x0;
    double y = coords.v1.y() * x_prime + coords.v2.y() * y_prime + y0;
    double z = coords.v1.z() * x_prime + coords.v2.z() * y_prime + z0;

    j_data.integration_points.emplace_back(x, y, z);
    j_data.strain_energy.emplace_back(0.0);
    j_data.kinetic_energy.emplace_back(0.0);
    j_data.stresses.emplace_back(Matrix3::zero());
    j_data.displacement_grads.emplace_back(Matrix3::zero());
  }

  return j_data;
}

std::expected<CrackFractureCalculator::CrackSurfaceData, std::string>
CrackFractureCalculator::calculateCrackSurfaceData(
  const LocalCoordinates& coords) const
{
  CrackSurfaceData surface_data;

  // Determine the intersection point based on COD option
  Point calculation_point;
  if (d_CODOption == 0 || d_CODOption == 1) {
    // Calculate distance from crack tip
    double distance =
      d_doCrackPropagationStep ? std::max(1.0, d_rdadx) * d_dx_max : d_rJ / 2.0;

    // Adjust distance if point is not on crack plane (COD option 1)
    if (d_CODOption == 1) {
      auto adjusted_distance = adjustDistanceForCrackPlane(coords, distance);
      if (!adjusted_distance) {
        return std::unexpected(
          std::format("Failed to adjust distance for crack plane: {}",
                      adjusted_distance.error()));
      }
      distance = adjusted_distance.value();
    }

    // Calculate global coordinates of the point (-distance, 0, 0) in local
    // system
    auto [x0, y0, z0] =
      std::make_tuple(coords.origin.x(), coords.origin.y(), coords.origin.z());
    calculation_point = Point(-distance * coords.v1.x() + x0,
                              -distance * coords.v1.y() + y0,
                              -distance * coords.v1.z() + z0);
  } else if (d_CODOption == 2) {
    // Use the intersection point between J-integral contour and crack plane
    std::array<double, 14> path_params;
    FindJIntegralPath(
      coords.origin, coords.v1, coords.v2, coords.v3, path_params.data());

    double temp_radius = d_rJ;
    auto intersection_result =
      findJIntegralCrackPlaneIntersection(path_params, temp_radius);
    if (!intersection_result) {
      return std::unexpected(
        std::format("Failed to find intersection for COD calculation: {}",
                    intersection_result.error()));
    }
    calculation_point = intersection_result.value();
  }

  surface_data.intersection_point = calculation_point;

  // Set up interpolation for the calculation point
  auto interpolator = d_flag->d_interpolator->clone(d_current_patch);
  if (!interpolator) {
    return std::unexpected(
      "Failed to create interpolator for crack surface data");
  }

  // Modern C++23 RAII wrapper for interpolator cleanup
  std::unique_ptr<Interpolator> interp_guard(interpolator);

  std::vector<IntVector> node_indices(interpolator->size());
  std::vector<double> shape_functions(interpolator->size());

  // Find interpolation weights at the calculation point
  try {
    // Note: Using placeholder values for particle size and deformation gradient
    // In practice, these would come from the nearest particle or be
    // appropriately chosen
    Matrix3 dummy_size   = Matrix3::identity();
    Matrix3 dummy_deform = Matrix3::identity();

    interpolator->findCellAndWeights(calculation_point,
                                     node_indices,
                                     shape_functions,
                                     dummy_size,
                                     dummy_deform);
  } catch (const std::exception& e) {
    return std::unexpected(
      std::format("Interpolation failed at calculation point: {}", e.what()));
  }

  // Calculate displacements and stresses at the calculation point
  Vector displacement_above = Vector::zero();
  Vector displacement_below = Vector::zero();
  Matrix3 stress_above      = Matrix3::zero();
  Matrix3 stress_below      = Matrix3::zero();

  double total_weight_above = 0.0;
  double total_weight_below = 0.0;

  // Interpolate values from grid nodes using modern range-based approach
  for (auto&& [node_idx, weight] :
       std::views::zip(node_indices, shape_functions)) {
    if (weight < std::numeric_limits<double>::epsilon()) {
      continue; // Skip negligible contributions
    }

    // Check if this node is in the crack zone
    bool is_crack_zone = (d_GnumPatls[node_idx] != 0);

    if (is_crack_zone) {
      // Above crack surface (original grid)
      displacement_above += d_gDisp[node_idx] * weight;
      stress_above += d_gGridStress[node_idx] * weight;
      total_weight_above += weight;

      // Below crack surface (ghost grid)
      displacement_below += d_GDisp[node_idx] * weight;
      stress_below += d_GGridStress[node_idx] * weight;
      total_weight_below += weight;
    } else {
      // Non-crack zone - use regular grid values for both surfaces
      displacement_above += d_gDisp[node_idx] * weight;
      displacement_below += d_gDisp[node_idx] * weight;
      stress_above += d_gGridStress[node_idx] * weight;
      stress_below += d_gGridStress[node_idx] * weight;
      total_weight_above += weight;
      total_weight_below += weight;
    }
  }

  // Normalize by total weights to ensure proper interpolation
  if (total_weight_above > std::numeric_limits<double>::epsilon()) {
    displacement_above /= total_weight_above;
    stress_above /= total_weight_above;
  }

  if (total_weight_below > std::numeric_limits<double>::epsilon()) {
    displacement_below /= total_weight_below;
    stress_below /= total_weight_below;
  }

  // Store the results
  surface_data.relative_displacement_a = displacement_above;
  surface_data.relative_displacement_b = displacement_below;
  surface_data.stress_traction_a       = stress_above;
  surface_data.stress_traction_b       = stress_below;

  // Transform to local coordinates for consistency with the original algorithm
  Vector displacement_diff       = displacement_above - displacement_below;
  Vector local_displacement_diff = coords.T * displacement_diff;

  // Transform stress tensors to local coordinates
  Matrix3 local_stress_a = transformStressToLocal(stress_above, coords.T);
  Matrix3 local_stress_b = transformStressToLocal(stress_below, coords.T);

  // Store transformed values (these are what the original code used)
  surface_data.relative_displacement_a = coords.T * displacement_above;
  surface_data.relative_displacement_b = coords.T * displacement_below;
  surface_data.stress_traction_a       = local_stress_a;
  surface_data.stress_traction_b       = local_stress_b;

  return surface_data;
}

std::expected<CrackFractureCalculator::FractureResult, std::string>
CrackFractureCalculator::computeFractureParameters(
  const LocalCoordinates& coords,
  const JIntegralData& j_data,
  const CrackSurfaceData& surface_data,
  int material_id,
  int node_idx) const
{
  FractureResult result;

  // Calculate contour integral
  auto [jx_contour, jy_contour] = calculateContourIntegral(j_data, coords);

  // Calculate area integral if needed
  double jx_area = 0.0, jy_area = 0.0;
  if (d_useVolumeIntegral || jx_contour < 0.0) {
    auto area_result = calculateAreaIntegral(coords, jx_contour);
    if (!area_result) {
      return std::unexpected(std::format("Area integral calculation failed: {}",
                                         area_result.error()));
    }
    std::tie(jx_area, jy_area) = area_result.value();
  }

  // Calculate friction work
  result.friction_work = calculateFrictionWork(surface_data, material_id);

  // Compute J-integral vector
  result.energy_release_rate = jx_contour + jx_area - result.friction_work;
  result.J_integral =
    Vector(result.energy_release_rate, jy_contour + jy_area, 0.0);

  // Convert to stress intensity factors
  // This would call the constitutive model conversion method
  MPMMaterial* mpm_matl =
    static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", material_id));
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

  // Get crack opening displacement and velocity
  Vector cod            = /* calculate COD */;
  double crack_velocity = d_cfSegVel[material_id][node_idx];

  cm->convertJToK(mpm_matl,
                  d_stressState[material_id],
                  result.J_integral,
                  crack_velocity,
                  cod,
                  result.stress_intensity_factors);

  return result;
}

std::pair<double, double>
CrackFractureCalculator::calculateContourIntegral(
  const JIntegralData& j_data,
  const LocalCoordinates& coords) const
{
  constexpr double PI = std::numbers::pi;
  const int num_segments =
    static_cast<int>(j_data.integration_points.size()) - 1;

  std::vector<double> f1_for_jx(num_segments + 1);
  std::vector<double> f1_for_jy(num_segments + 1);

  // Calculate integrand values using ranges
  for (auto&& [j, point] : j_data.integration_points | std::views::enumerate) {
    double angle = 2.0 * PI * static_cast<double>(j) / num_segments;

    cosTheta  = (xcprime * cos(angle) - ycprime * sin(angle)) / scprime;
    sinTheta  = (ycprime * cos(angle) + xcprime * sin(angle)) / scprime;
    double t1 = st[j](0, 0) * cosTheta + st[j](0, 1) * sinTheta;
    double t2 = st[j](1, 0) * cosTheta + st[j](1, 1) * sinTheta;

    Vector t123 = Vector(t1, t2, 0. /*t3*/); // plane state
    Vector dgx  = Vector(dg[j](0, 0), dg[j](1, 0), dg[j](2, 0));
    Vector dgy  = Vector(dg[j](0, 1), dg[j](1, 1), dg[j](2, 1));

    f1ForJx[j] = (W[j] + K[j]) * cosTheta - Dot(t123, dgx);
    f1ForJy[j] = (W[j] + K[j]) * sinTheta - Dot(t123, dgy);
  }

  // Integrate using modern algorithms
  double jx = std::transform_reduce(f1_for_jx.begin(),
                                    f1_for_jx.end() - 1,
                                    f1_for_jx.begin() + 1,
                                    0.0,
                                    std::plus<>{},
                                    [](double a, double b) { return a + b; });

  double jy = std::transform_reduce(f1_for_jy.begin(),
                                    f1_for_jy.end() - 1,
                                    f1_for_jy.begin() + 1,
                                    0.0,
                                    std::plus<>{},
                                    [](double a, double b) { return a + b; });

  jx *= j_data.radius * PI / num_segments;
  jy *= j_data.radius * PI / num_segments;

  return { jx, jy };
}

void
CrackFractureCalculator::CalculateFractureParameters(
  const ProcessorGroup* pg,
  const PatchSubset* patches,
  const MaterialSubset* matls,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  // Get simulation time using modern variable handling
  simTime_vartype sim_time_var;
  old_dw->get(sim_time_var, lb->simulationTimeLabel);
  double time = sim_time_var;

  delt_vartype del_t;
  old_dw->get(del_t, lb->delTLabel, getLevel(patches));

  // Process patches using ranges
  for (auto&& [patch_idx, patch] :
       patches->getVector() | std::views::enumerate) {
    Vector dx     = patch->dCell();
    double dx_max = std::max({ dx.x(), dx.y(), dx.z() });

    // MPI setup with modern initialization
    int patch_size, processor_id;
    MPI_Comm_size(d_mpi_crack_comm, &patch_size);
    MPI_Comm_rank(d_mpi_crack_comm, &processor_id);

    // Process materials using ranges
    int num_materials = d_mat_manager->getNumMaterials("MPM");
    for (int material_id : std::views::iota(0, num_materials)) {
      if (!d_calFractParametersStep && !d_doCrackPropagationStep) {
        continue;
      }

      // Resize data structures
      int crack_front_size = static_cast<int>(d_cfSegNodes[material_id].size());
      d_cfSegJ[material_id].resize(crack_front_size);
      d_cfSegK[material_id].resize(crack_front_size);

      // Process crack front nodes across all processors
      for (int proc_id : std::views::iota(0, patch_size)) {
        int num_nodes = static_cast<int>(d_cfnset[material_id][proc_id].size());
        if (num_nodes == 0) {
          continue;
        }

        std::vector<Vector> crack_front_j(num_nodes);
        std::vector<Vector> crack_front_k(num_nodes);

        if (processor_id == proc_id) {
          // Process nodes in this processor using ranges and modern error
          // handling
          for (auto&& [local_idx, global_idx] :
               d_cfnset[material_id][proc_id] | std::views::enumerate) {

            auto result =
              processCrackFrontNode(material_id, global_idx, old_dw, new_dw);
            if (!result) {
              std::print(stderr,
                         "Error processing crack front node {}: {}\n",
                         global_idx,
                         result.error());
              continue;
            }

            crack_front_j[local_idx] = result.value().J_integral;
            crack_front_k[local_idx] = result.value().stress_intensity_factors;
          }
        }

        // Broadcast results using modern MPI handling
        MPI_Datatype mpi_vector =
          fun_getTypeDescription(static_cast<Vector*>(nullptr))->getMPIType();
        MPI_Bcast(crack_front_j.data(),
                  num_nodes,
                  mpi_vector,
                  proc_id,
                  d_mpi_crack_comm);
        MPI_Bcast(crack_front_k.data(),
                  num_nodes,
                  mpi_vector,
                  proc_id,
                  d_mpi_crack_comm);

        // Store results using ranges
        for (auto&& [local_idx, global_idx] :
             d_cfnset[material_id][proc_id] | std::views::enumerate) {
          d_cfSegJ[material_id][global_idx] = crack_front_j[local_idx];
          d_cfSegK[material_id][global_idx] = crack_front_k[local_idx];
        }
      }

      // Output results from processor 0
      if (processor_id == 0) {
        OutputCrackFrontResults(material_id, time, del_t);
      }
    }
  }
}

std::expected<CrackFractureCalculator::FractureResult, std::string>
CrackFractureCalculator::processCrackFrontNode(
  int material_id,
  int node_idx,
  const DataWarehouse* old_dw,
  const DataWarehouse* new_dw) const
{
  // Setup local coordinates
  auto coords_result = setupLocalCoordinates(material_id, node_idx);
  if (!coords_result) {
    return std::unexpected(std::format("Failed to setup local coordinates: {}",
                                       coords_result.error()));
  }

  // Setup J-integral contour
  auto j_data_result = setupJIntegralContour(coords_result.value());
  if (!j_data_result) {
    return std::unexpected(std::format("Failed to setup J-integral contour: {}",
                                       j_data_result.error()));
  }

  // Calculate crack surface data
  auto surface_data_result = calculateCrackSurfaceData(coords_result.value());
  if (!surface_data_result) {
    return std::unexpected(
      std::format("Failed to calculate crack surface data: {}",
                  surface_data_result.error()));
  }

  // Compute fracture parameters
  return computeFractureParameters(coords_result.value(),
                                   j_data_result.value(),
                                   surface_data_result.value(),
                                   material_id,
                                   node_idx);
}

} // end namespace Vaango