/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __CORE_GRID_MATERIAL_H__
#define __CORE_GRID_MATERIAL_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <string_view>
#include <memory>

namespace Uintah {

class Material
{
public:
  Material(ProblemSpecP& ps);
  Material() = default;

  virtual ~Material(); 

  Material(const Material& mat) = delete;
  Material& operator=(const Material& mat) = delete;

  virtual ProblemSpecP outputProblemSpec(ProblemSpecP& ps);

  // Return index associated with this material's
  // location in the data warehouse
  int getDWIndex() const;

  void setDWIndex(int);

  const MaterialSubset* thisMaterial() const { return d_mat_subset.get(); }

  bool hasName() const { return d_have_name; }
  std::string_view getName() const { return d_name; }

protected:
  // Index associated with this material's spot in the DW
  int d_dwindex{ -1 };
  std::unique_ptr<MaterialSubset> d_mat_subset{ nullptr };

private:
  bool d_have_name{ false };
  std::string d_name{ "" };

};
} // End namespace Uintah

#endif // __CORE_GRID_MATERIAL_H__
