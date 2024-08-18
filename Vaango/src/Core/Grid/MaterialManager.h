/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __CORE_GRID_MATERIAL_MANAGER_H__
#define __CORE_GRID_MATERIAL_MANAGER_H__

#include <Core/Util/RefCounted.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Util/DOUT.hpp>

#include <sci_defs/uintah_defs.h>

#include <map>
#include <memory>
#include <vector>

#define OVERHEAD_WINDOW 40

namespace Uintah {

class VarLabel;
class Material;
class EmptyMaterial;

class MaterialManager : public RefCounted
{
public:
  MaterialManager();
  ~MaterialManager() override;

  MaterialManager(const MaterialManager&) = delete;

  auto
  operator=(const MaterialManager&) -> MaterialManager& = delete;

  void
  clearMaterials();

  void
  finalizeMaterials();

  void
  registerMaterial(const std::string& name, std::shared_ptr<Material> material);

  void
  registerMaterial(const std::string& name,
                   std::shared_ptr<Material> material,
                   std::uint32_t index);
  void
  registerEmptyMaterial(std::shared_ptr<EmptyMaterial> material);

  auto
  allMaterials(const std::string& name) const -> const MaterialSet*;

  auto
  allMaterials() const -> const MaterialSet*;

  auto
  getNumMaterials(const std::string& name) const -> size_t;

  auto
  getNumMaterials() const -> size_t
  {
    return d_all_materials.size();
  }

  auto
  getMaterial(const std::string& name, std::uint32_t index) const -> Material*;

  auto
  getMaterial(std::uint32_t index) const -> Material*
  {
    return index < d_all_materials.size() ? d_all_materials[index].get()
                                          : nullptr;
  }

  auto
  getAllInOneMaterial() const -> const MaterialSubset*
  {
    return d_all_in_one_material;
  }

  auto
  originalAllMaterials() const -> const MaterialSet*;

  void
  setOriginalMatlsFromRestart(MaterialSet* matls);

  auto
  parseAndLookupMaterial(ProblemSpecP& params, const std::string& name) const
    -> Material*;

private:
  void
  registerMaterial(std::shared_ptr<Material> material);

  void
  registerMaterial(std::shared_ptr<Material> material, std::uint32_t index);

  auto
  getMaterialByName(const std::string& name) const -> Material*;

  // Named materials
  std::map<std::string, std::shared_ptr<Material>> d_named_materials;

  // All empty materials
  std::vector<std::shared_ptr<EmptyMaterial>> d_empty_materials;

  // All materials from simulator components
  std::vector<std::shared_ptr<Material>> d_all_materials;
  MaterialSet* d_all_materialsets{ nullptr };

  // Materials from each simulator component
  std::map<std::string, std::vector<std::shared_ptr<Material>>> d_materials;
  std::map<std::string, MaterialSet*> d_materialsets;

  // The switcher needs to clear the materials, but don't
  // delete them or there might be VarLabel problems when
  // CMs are destroyed.  Store them here until the end.
  std::vector<std::shared_ptr<Material>> d_materials_old;

  // Keep track of all the original materials if switching
  MaterialSet* d_all_materialsets_old{ nullptr };
  MaterialSubset* d_all_in_one_material{ nullptr };

  inline static int s_count = 0;

}; // end class MaterialManager

} // End namespace Uintah

#endif //__CORE_GRID_MATERIAL_MANAGER_H__
