/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef SGP_DATA_H
#define SGP_DATA_H

#include <CCA/Components/MPM/ParticleCreator/ParticleCreatorStructs.h>
#include <Core/GeometryPiece/SpecialGeomPiece.h>

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace Vaango {
namespace Helpers {

// Helper for getting and iterating over SpecialGeometryPiece data
template<typename T>
struct SGPData
{
  const std::vector<T>* data = nullptr;
  typename std::vector<T>::const_iterator iterator;

  SGPData() = default;

  void
  initialize(Uintah::SpecialGeomPiece* sgp,
             const Uintah::ObjectVars& obj_vars,
             const std::string& name,
             Uintah::GeometryObject* obj_ptr)
  {
    if (sgp) {
      if constexpr (std::is_same_v<T, double>) {
        data     = sgp->getScalar(name);
        auto key = std::make_pair(name, obj_ptr);
        if (auto map_iter = obj_vars.scalars.find(key);
            map_iter != obj_vars.scalars.end()) {
          iterator = map_iter->second.begin();
        }
      } else if constexpr (std::is_same_v<T, Uintah::Vector>) {
        data     = sgp->getVector(name);
        auto key = std::make_pair(name, obj_ptr);
        if (auto map_iter = obj_vars.vectors.find(key);
            map_iter != obj_vars.vectors.end()) {
          iterator = map_iter->second.begin();
        }
      } else if constexpr (std::is_same_v<T, Uintah::Matrix3>) {
        data     = sgp->getTensor(name);
        auto key = std::make_pair(name, obj_ptr);
        if (auto map_iter = obj_vars.tensors.find(key);
            map_iter != obj_vars.tensors.end()) {
          iterator = map_iter->second.begin();
        }
      }
    }
  }

  std::optional<T>
  get_and_advance()
  {
    if (data && iterator != data->end()) {
      return *iterator++;
    }
    return std::nullopt;
  }
};

} // namespace Helpers
} // end namespace Vaango

#endif // SGP_DATA_H
