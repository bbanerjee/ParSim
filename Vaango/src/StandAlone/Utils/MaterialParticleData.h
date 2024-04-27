/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#ifndef __VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_DATA_H__
#define __VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_DATA_H__

#include <StandAlone/Utils/MaterialParticleVarData.h>

#include <Core/Grid/Variables/ParticleSubset.h>

#include <map>
#include <string>

namespace Vaango {
namespace Utils {
namespace CompareUda {

class MaterialParticleData
{
public:
  MaterialParticleData()
    : matl_(-999)
    , particleIDs_(0)
  {
  }
  MaterialParticleData(int matl)
    : matl_(matl)
    , particleIDs_(0)
  {
  }
  MaterialParticleData(const MaterialParticleData& copy)
    : matl_(copy.matl_)
    , vars_(copy.vars_)
    , particleIDs_(copy.particleIDs_)
  {
  }

  ~MaterialParticleData() {}

  MaterialParticleVarData&
  operator[](const std::string& varName)
  {
    MaterialParticleVarData& result = vars_[varName];
    vars_[varName].setVarName(varName);

    if (varName == "p.particleID") {
      particleIDs_ = &result;
    }
    return result;
  }

  void
  compare(MaterialParticleData& data2,
          double time1,
          double time2,
          double abs_tolerance,
          double rel_tolerance);

  void
  setMatl(int matl)
  {
    matl_ = matl;
  }

private:
  void
  createPatchMap();

  void
  gather(Uintah::ParticleSubset* gatherSubset);

  void
  sort();

  int matl_;
  std::map<std::string, MaterialParticleVarData> vars_;
  MaterialParticleVarData* particleIDs_; // will point to one of vars_
};

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_DATA_H__
