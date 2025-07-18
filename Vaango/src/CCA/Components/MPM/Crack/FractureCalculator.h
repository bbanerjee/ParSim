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

#ifndef __MPM_CRACK_FRACTURE_CALCULATOR__
#define __MPM_CRACK_FRACTURE_CALCULATOR__

#include <Core/Geometry/Vector.h>
#include <Core/Math/Matrix3.h>

#include <expected>

// Forward declarations to avoid circular includes
namespace Uintah {
  class ProcessorGroup;
  class PatchSubset;
  class MaterialSubset;
  class DataWarehouse;
}

namespace Vaango {

// C++23 features used:
// - std::expected for error handling
// - std::span for safer array handling
// - ranges and views for cleaner iteration
// - std::format for better string formatting
// - structured bindings and auto improvements

class FractureCalculator
{
public:
  struct LocalCoordinates
  {
    Uintah::Vector v1, v2, v3;
    Uintah::Matrix3 T, TT; // Transformation matrices
    Uintah::Point origin;
  };

  struct JIntegralData
  {
    std::vector<Uintah::Point> integration_points;
    std::vector<double> strain_energy;
    std::vector<double> kinetic_energy;
    std::vector<Uintah::Matrix3> stresses;
    std::vector<Uintah::Matrix3> displacement_grads;
    double radius;
  };

  struct CrackSurfaceData
  {
    Uintah::Vector relative_displacement_a, relative_displacement_b;
    Uintah::Matrix3 stress_traction_a, stress_traction_b;
    Uintah::Point intersection_point;
  };

  struct FractureResult
  {
    Uintah::Vector J_integral;
    Uintah::Vector stress_intensity_factors;
    double energy_release_rate;
    double friction_work;
  };

private:
  // Modern C++23 error handling with expected
  std::expected<LocalCoordinates, std::string>
  setupLocalCoordinates(int material_id, int node_id) const;

  std::expected<JIntegralData, std::string>
  setupJIntegralContour(const LocalCoordinates& coords) const;

  std::expected<CrackSurfaceData, std::string>
  calculateCrackSurfaceData(const LocalCoordinates& coords) const;

  std::expected<FractureResult, std::string>
  computeFractureParameters(const LocalCoordinates& coords,
                            const JIntegralData& j_data,
                            const CrackSurfaceData& surface_data,
                            int material_id,
                            int node_idx) const;

  // Helper methods using modern C++23 features
  std::expected<Point, std::string>
  findJIntegralCrackPlaneIntersection(std::span<const double, 14> path_params,
                                      double& radius) const;

  void
  evaluateFieldsAtIntegrationPoints(JIntegralData& j_data) const;

  std::pair<double, double>
  calculateContourIntegral(const JIntegralData& j_data,
                           const LocalCoordinates& coords) const;

  std::expected<std::pair<double, double>, std::string>
  calculateAreaIntegral(const LocalCoordinates& coords,
                        double contour_jx) const;

  double
  calculateFrictionWork(const CrackSurfaceData& surface_data,
                        int material_id) const;

public:
  // Main refactored function using C++23 features
  void
  CalculateFractureParameters(const Uintah::ProcessorGroup* pg,
                              const Uintah::PatchSubset* patches,
                              const Uintah::MaterialSubset* matls,
                              Uintah::DataWarehouse* old_dw,
                              Uintah::DataWarehouse* new_dw);

private:
  // Process a single crack front node using modern C++23
  std::expected<FractureResult, std::string>
  processCrackFrontNode(int material_id,
                        int node_idx,
                        const Uintah::DataWarehouse* old_dw,
                        const Uintah::DataWarehouse* new_dw) const;
};

} // end namespace Vaango

#endif // __MPM_CRACK_FRACTURE_CALCULATOR__