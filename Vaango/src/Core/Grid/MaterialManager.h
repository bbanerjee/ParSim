/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __CORE_GRID_MATERIAL_MANAGER_H__
#define __CORE_GRID_MATERIAL_MANAGER_H__

#include <Core/Geometry/Vector.h>
#include <Core/Grid/Ghost.h>
#include <Core/Grid/SimulationTime.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Math/MinMax.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/RefCounted.h>

#include <iostream>
#include <map>
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
  ~MaterialManager();

  void clearMaterials();
  void finalizeMaterials();

  void registerMaterial(const std::string& name, Material* material);
  void registerMaterial(const std::string& name, Material* material,
                        unsigned int index);
  void registerEmptyMaterial(EmptyMaterial* material);

  const MaterialSet* allMaterials(const std::string& name) const;
  const MaterialSet* allMaterials() const;

  size_t getNumMaterials(const std::string& name) const;
  size_t getNumMaterials() const
  {
    return d_all_materials.size();
  }

  Material* getMaterial(const std::string& name, int index) const;
  Material* getMaterial(int index) const { return d_all_materials[index]; }

  Material* getMaterialByName(const std::string& name) const;

  MaterialSubset* getAllInOneMaterial() { return d_all_in_one_material; }
  const MaterialSet* originalAllMaterials() const;

  void setOriginalMaterialsFromRestart(MaterialSet* materials);

  Material* parseAndLookupMaterial(ProblemSpecP& params,
                                   const std::string& name) const;

private:
  MaterialManager(const MaterialManager&) = delete;
  MaterialManager& operator=(const MaterialManager&) = delete;

  void registerMaterial(Material* material);
  void registerMaterial(Material* material, unsigned int index);

  std::vector<Material*> d_all_materials;
  MaterialSet* d_all_materialsets{ nullptr };
  std::vector<EmptyMaterial*> d_empty_materials;
  std::map<std::string, Material*> d_named_materials;

  std::map<std::string, std::vector<Material*>> d_materials;
  std::map<std::string, MaterialSet*> d_materialsets;

  // The switcher needs to clear the materials, but don't
  // delete them or there mihgt be VarLabel problems when
  // CMs are destroyed.  Store them here until the end.
  std::vector<Material*> d_materials_old;

  // Keep track of all the original materials if switching
  MaterialSet* d_all_materials_old{ nullptr };
  MaterialSubset* d_all_in_one_material{ nullptr };

  static int s_count;

}; // end class MaterialManager

} // End namespace Uintah

#endif //__CORE_GRID_MATERIAL_MANAGER_H__
