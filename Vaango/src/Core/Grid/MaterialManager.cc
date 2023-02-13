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

#include <Core/Grid/MaterialManager.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/EmptyMaterial.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <sstream>
#include <iostream>

namespace Uintah {

MaterialManager::MaterialManager()
{
  if (s_count > 0) {
    std::ostringstream out;
    out << "**ERROR** Multiple material managers created";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  s_count++;
}

MaterialManager::~MaterialManager()
{
  clearMaterials();

  if (d_all_materialsets_old && d_all_materialsets_old->removeReference()) {
    delete d_all_materialsets_old;
  }

  s_count--;
}

void
MaterialManager::clearMaterials()
{
  for (auto& material : d_all_materials) {
    d_materials_old.push_back(material);
  }

  if (d_all_materialsets && d_all_materialsets->removeReference()) {
    delete d_all_materialsets;
  }
  d_all_materialsets = nullptr;

  if (d_all_in_one_material && d_all_in_one_material->removeReference()) {
    delete d_all_in_one_material;
  }
  d_all_in_one_material = nullptr;

  for (auto& [material_name, material] : d_materials) {
    if (d_materialsets[material_name] &&
        d_materialsets[material_name]->removeReference()) {
      delete d_materialsets[material_name];
    }
    d_materials[material_name].clear();
    d_materialsets[material_name] = nullptr;
  }

  d_all_materials.clear();
  d_empty_materials.clear();
  d_named_materials.clear();
}

void
MaterialManager::finalizeMaterials()
{
  // All materials
  if (d_all_materialsets && d_all_materialsets->removeReference()) {
    delete d_all_materialsets;
  }
  d_all_materialsets = scinew MaterialSet();
  d_all_materialsets->addReference();

  std::vector<int> all_materialsvec(d_all_materials.size());
  for (size_t i = 0; i < d_all_materials.size(); i++) {
    all_materialsvec[i] = d_all_materials[i]->getDWIndex();
  }
  d_all_materialsets->addAll(all_materialsvec);

  // Original all_materials (before switch using Switcher)
  if (d_all_materialsets_old == nullptr) {
    d_all_materialsets_old = scinew MaterialSet();
    d_all_materialsets_old->addReference();
    d_all_materialsets_old->addAll(all_materialsvec);
  }

  // All materials in one:  a material that represents all materials
  // (i.e. summed over all materials -- the whole enchilada)
  if (d_all_in_one_material && d_all_in_one_material->removeReference()) {
    delete d_all_in_one_material;
  }
  d_all_in_one_material = scinew MaterialSubset();
  d_all_in_one_material->addReference();
  d_all_in_one_material->add(static_cast<int>(d_all_materials.size()));

  for (auto& [material_name, materialset] : d_materialsets) {
    if (materialset && materialset->removeReference()) {
      delete materialset;
    }
  }

  for (auto& [material_name, material] : d_materials) {
    d_materialsets[material_name] = scinew MaterialSet();
    d_materialsets[material_name]->addReference();

    std::vector<int> materialsvec(material.size());
    for (size_t i = 0; i < material.size(); i++) {
      materialsvec[i] = material[i]->getDWIndex();
    }
    d_materialsets[material_name]->addAll(materialsvec);
  }
}

void
MaterialManager::registerMaterial(const std::string& name,
                                  std::shared_ptr<Material> material)
{
  d_materials[name].push_back(material);
  registerMaterial(material);
}

void
MaterialManager::registerMaterial(const std::string& name,
                                  std::shared_ptr<Material> material,
                                  std::uint32_t index)
{
  d_materials[name].push_back(material);
  registerMaterial(material, index);
}

void
MaterialManager::registerMaterial(std::shared_ptr<Material> material)
{
  material->setDWIndex((int)d_all_materials.size());
  d_all_materials.push_back(material);
  if (material->hasName()) {
    d_named_materials[static_cast<std::string>(material->getName())] =
      material;
  }
}

void
MaterialManager::registerMaterial(std::shared_ptr<Material> material,
                                  std::uint32_t index)
{
  if (d_all_materials.size() <= index) {
    index = d_all_materials.size();
    d_all_materials.resize(index + 1);
  }
  material->setDWIndex(index);
  d_all_materials[index] = material;

  if (material->hasName()) {
    d_named_materials[static_cast<std::string>(material->getName())] =
      material;
  }
}

void
MaterialManager::registerEmptyMaterial(std::shared_ptr<EmptyMaterial> material)
{
  d_empty_materials.push_back(material);
  registerMaterial(material);
}

auto
MaterialManager::allMaterials(const std::string& name) const
  -> const MaterialSet*
{
  if (d_materialsets.find(name) != d_materialsets.end()) {
    return d_materialsets.at(name);
  }
  return nullptr;
}

auto
MaterialManager::allMaterials() const -> const MaterialSet*
{
  ASSERT(d_all_materialsets != nullptr);
  return d_all_materialsets;
}

auto
MaterialManager::getNumMaterials(const std::string& name) const -> size_t
{
  if (d_materials.find(name) != d_materials.end()) {
    return d_materials.at(name).size();
  }
  return 0;
}

auto
MaterialManager::getMaterial(const std::string& name, std::uint32_t index) const
  -> Material*
{
  if (d_materials.find(name) != d_materials.end()) {
    if (index < d_materials.at(name).size()) {
      return d_materials.at(name)[index].get();
    }
  }
  return nullptr;
}

auto
MaterialManager::originalAllMaterials() const -> const MaterialSet*
{
  ASSERT(d_all_materialsets_old != nullptr);
  return d_all_materialsets_old;
}

void
MaterialManager::setOriginalMatlsFromRestart(MaterialSet* matls)
{
  d_all_materialsets_old->removeReference();
  d_all_materialsets_old = matls;
  d_all_materialsets_old->addReference();
}

auto
MaterialManager::parseAndLookupMaterial(ProblemSpecP& params,
                                        const std::string& name) const
  -> Material*
{
  // for single material problems return matl 0
  Material* result = getMaterial(0);
  if (getNumMaterials() > 1) {
    std::string material_name;
    if (!params->get(name, material_name)) {
      std::ostringstream out;
      out << "**ERROR** Cannot find <material> section";
      throw ProblemSetupException(out.str(), __FILE__, __LINE__);
    }

    // remove quotation, a common user input error
    material_name.erase(
      std::remove(material_name.begin(), material_name.end(), '\"'),
      material_name.end());

    result = getMaterialByName(material_name);
    if (!result) {
      std::ostringstream out;
      out << "**ERROR** Cannot find a material named:" << material_name;
      throw ProblemSetupException(out.str(), __FILE__, __LINE__);
    }
  }
  return result;
}

auto
MaterialManager::getMaterialByName(const std::string& name) const -> Material*
{
  auto iter = d_named_materials.find(name);
  if (iter != d_named_materials.end()) {
    return iter->second.get();
  }
  return nullptr;
}

} // end namespace Uintah