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

#include <Core/GeometryPiece/DynamicSDFGeometry.h>

#include <Core/GeometryPiece/GeometryPieceFactory.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Grid/Box.h>

#include <Core/Exceptions/ProblemSetupException.h>

#include <iostream>



namespace Uintah {



DynamicSDFGeometry::DynamicSDFGeometry(ProblemSpecP& ps,

                                       const std::string& input_ups_dir)

{


  ps->get("res", d_res); // Grid resolution for SDF
  
  // Find the tri geometry piece inside
  ProblemSpecP tri_ps = ps->findBlock("tri");
  if (!tri_ps) {
    throw ProblemSetupException("DynamicSDFGeometry requires a <tri> block", __FILE__, __LINE__);
  }
  d_mesh = std::make_shared<TriGeometryPiece>(tri_ps, input_ups_dir);

  Box box = d_mesh->getBoundingBox();
  // Pad the box slightly to ensure SDF is defined outside
  Vector pad = (box.upper() - box.lower()) * 0.1;
  // Ensure pad has some minimum thickness to avoid degenerate boxes in 2D
  if (pad.x() < 1e-6) pad.x(1e-6);
  if (pad.y() < 1e-6) pad.y(1e-6);
  if (pad.z() < 1e-6) pad.z(1e-6);
  Box paddedBox(box.lower() - pad, box.upper() + pad);

  d_sdf = std::make_unique<LocalSDF>(paddedBox, d_res);

  // Generate SDF
  Vector spacing =
    (paddedBox.upper() - paddedBox.lower()) / (d_res - IntVector(1, 1, 1));
  auto triangles = d_mesh->getTriangles();
  auto points    = d_mesh->getPoints();

  for (int i = 0; i < d_res.x(); ++i) {
    for (int j = 0; j < d_res.y(); ++j) {
      for (int k = 0; k < d_res.z(); ++k) {
        Point p      = paddedBox.lower() + spacing * Vector(i, j, k);
        bool isInside = d_mesh->inside(p);

        // Brute force distance to triangles
        double minDistSq = 1.0e30;
        for (const auto& triIdx : triangles) {
          Point v0 = points[triIdx.x()];
          Point v1 = points[triIdx.y()];
          Point v2 = points[triIdx.z()];

          // Simplified distance to centroid for now
          // TODO: Implement real point-to-triangle distance
          Point centroid =
            Point(0, 0, 0) + (v0.asVector() + v1.asVector() + v2.asVector()) / 3.0;
          double dSq = (p - centroid).length2();
          if (dSq < minDistSq) {
            minDistSq = dSq;
          }
        }

        double dist = std::sqrt(minDistSq);
        d_sdf->setDistance(IntVector(i, j, k), isInside ? -dist : dist);
      }
    }
    }
  }
  
  DynamicSDFGeometry::DynamicSDFGeometry(const DynamicSDFGeometry& other)
    : GeometryPiece(other)
  {
    d_res  = other.d_res;
    d_mesh = other.d_mesh;
    if (other.d_sdf) {
      d_sdf = other.d_sdf->clone();
    }
  }
  
  GeometryPieceP
  DynamicSDFGeometry::clone() const
  {
  
  return std::make_shared<DynamicSDFGeometry>(*this);
}

bool DynamicSDFGeometry::inside(const Point& p) const {
  return d_sdf->getDistance(p) < 0;
}

double DynamicSDFGeometry::getSDF(const Point& p) const {
  return d_sdf->getDistance(p);
}

Vector DynamicSDFGeometry::getSDFGradient(const Point& p) const {
  return d_sdf->getGradient(p);
}

Box DynamicSDFGeometry::getBoundingBox() const {
  return d_sdf->getBoundingBox();
}

void DynamicSDFGeometry::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("res", d_res);
  // Mesh output would go here
}

} // namespace Uintah
