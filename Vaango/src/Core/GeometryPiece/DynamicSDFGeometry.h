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

#ifndef __CORE_GEOMETRYPIECE_DYNAMIC_SDF_GEOMETRY_H__
#define __CORE_GEOMETRYPIECE_DYNAMIC_SDF_GEOMETRY_H__

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/LocalSDF.h>
#include <Core/GeometryPiece/TriGeometryPiece.h>
#include <memory>

namespace Uintah {

class DynamicSDFGeometry : public GeometryPiece {
public:
  inline static const std::string TYPE_NAME{ "dynamic_sdf" };

  DynamicSDFGeometry(ProblemSpecP& ps, const std::string& input_ups_dir = "");
  DynamicSDFGeometry(const DynamicSDFGeometry& other);
  virtual ~DynamicSDFGeometry() = default;

  virtual std::string getType() const override { return TYPE_NAME; }
  virtual GeometryPieceP clone() const override;

  virtual bool inside(const Point& p) const override;
  virtual double getSDF(const Point& p) const override;
  virtual Vector getSDFGradient(const Point& p) const override;
  virtual Box getBoundingBox() const override;

  const LocalSDF& getSDF() const { return *d_sdf; }

private:
  virtual void outputHelper(ProblemSpecP& ps) const override;

  std::unique_ptr<LocalSDF> d_sdf;
  std::shared_ptr<TriGeometryPiece> d_mesh;
  IntVector d_res;
};

} // namespace Uintah

#endif
