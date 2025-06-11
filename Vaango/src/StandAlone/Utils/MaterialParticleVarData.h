/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef __VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_VAR_DATA_H__
#define __VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_VAR_DATA_H__

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <map>
#include <string>

namespace Vaango {
namespace Utils {
namespace CompareUda {

/**********************************************************************
 * MaterialParticleVarData and MaterialParticleData are for comparing
 * ParticleVariables when the patch distributions are different in the
 * different uda's -- p.particleID must be a supplied variable for this
 * to work.
 *********************************************************************/

class MaterialParticleVarData
{
public:
  MaterialParticleVarData()
    : m_name("")
    , m_particleIDData(0)
    , m_patchMap(0)
  {
  }
  MaterialParticleVarData(const std::string& name)
    : m_name(name)
    , m_particleIDData(0)
    , m_patchMap(0)
  {
  }

  ~MaterialParticleVarData(); // needs to delete the m_particleVars
                              // and patchMap if "p.particleID"

  void
  setVarName(const std::string& varName)
  {
    m_name = varName;
  }

  const std::string&
  getName()
  {
    return m_name;
  }

  // add for each patch
  void
  add(Uintah::ParticleVariableBase* pvb, const Uintah::Patch* patch);

  // gather all patches (vectors) into one
  void
  gather(Uintah::ParticleSubset* gatherSubset);

  bool
  compare(MaterialParticleVarData& data2,
          int matl,
          double time1,
          double time2,
          double abs_tolerance,
          double rel_tolerance);

  void
  createPatchMap();

  void
  setParticleIDData(MaterialParticleVarData* particleIDData)
  {
    m_particleIDData = particleIDData;
    m_patchMap       = particleIDData->m_patchMap;
  }

  const std::vector<Uintah::ParticleVariableBase*>&
  getParticleVars()
  {
    return m_particleVars;
  }

  Uintah::long64
  getParticleID(Uintah::particleIndex index);
  const Uintah::Patch*
  getPatch(Uintah::particleIndex index);

private:
  template<class T>
  bool
  compare(MaterialParticleVarData& data2,
          Uintah::ParticleVariable<T>* value1,
          Uintah::ParticleVariable<T>* value2,
          int matl,
          double time1,
          double time2,
          double abs_tolerance,
          double rel_tolerance);

  std::string m_name;
  // vector elements each represent a patch -- doesn't matter which
  std::vector<Uintah::ParticleVariableBase*> m_particleVars;
  std::vector<Uintah::ParticleSubset*> subsets_;
  std::vector<const Uintah::Patch*> m_patches;

  MaterialParticleVarData* m_particleIDData;
  std::map<Uintah::long64, const Uintah::Patch*>* m_patchMap;
};

} // namespace CompareUda
} // namespace Utils
} // namespace Vaango

#endif //__VAANGO_STANDALONE_UTILS_MATERIAL_PARTICLE_VAR_DATA_H__
