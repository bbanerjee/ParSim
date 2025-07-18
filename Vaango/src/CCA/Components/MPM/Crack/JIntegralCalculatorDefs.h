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

#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Math/Matrix3.h>

#include <vector>

// Forward declarations to avoid circular includes
namespace Uintah {
class ProcessorGroup;
} // namespace Uintah

namespace Vaango {

// Enums for better readability
enum class CrackSegmentPosition
{
  LeftEdge,
  RightEdge,
  Middle,
  SingleSegment
};

enum class CODOption
{
  FixedDistance           = 0,
  MaxDistanceOnCrackPlane = 1,
  IntersectionPoint       = 2
};

// Represents all necessary nodal data
struct NodalSolutionData
{
  // Using constNCVariable directly might be tricky for copying;
  // ideally, these would be views or references to avoid deep copies.
  // For simplicity, let's assume direct access or pass a wrapper.
  Uintah::constNCVariable<double> gMass, Gmass;
  Uintah::constNCVariable<int> GnumPatls;
  Uintah::constNCVariable<Uintah::Vector> gdisp, Gdisp;
  Uintah::constNCVariable<Uintah::Matrix3> ggridStress, GgridStress;
  Uintah::constNCVariable<Uintah::Matrix3> gdispGrads, GdispGrads;
  Uintah::constNCVariable<double> gW, GW;
  Uintah::constNCVariable<double> gK, GK;
  Uintah::constNCVariable<Uintah::Vector> gacc, Gacc;
  Uintah::constNCVariable<Uintah::Vector> gvel, Gvel;
  Uintah::constNCVariable<Uintah::Matrix3> gvelGrads, GvelGrads;

  // A way to initialize this struct
  NodalSolutionData(Uintah::DataWarehouse* new_dw,
                    const Uintah::MPMLabel* lb,
                    int dwi,
                    const Uintah::Patch* patch,
                    int NGC);
};

// Represents data at a single integration point on the J-contour
struct JContourIntegrationPointData
{
  Uintah::Point X;     // Integration points
  double W;            // Strain energy density
  double K;            // Kinetic energy density
  Uintah::Matrix3 ST;  // Stress in global coordinates
  Uintah::Matrix3 DG;  // DispGrad in global coordinates
  Uintah::Matrix3 st;  // Stress in local coordinates
  Uintah::Matrix3 dg;  // DispGrad in local coordinates
};

// Error codes for std::expected
enum class JCalculationError
{
  JIntegralRadiusTooSmall,
  MPICommunicationFailed, // Though MPI is external, could be an internal error
                          // representation
  InvalidCrackGeometry    // For FindSegsFromNode or GetPositionToComputeCOD
};

} // namespace Vaango