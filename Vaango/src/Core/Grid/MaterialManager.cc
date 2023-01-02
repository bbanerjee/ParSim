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
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <sstream>
#include <string>

namespace Uintah {

// static vars
int MaterialManager::s_count = 0;

MaterialManager::MaterialManager()
{
  if (s_count++ >= 1) {
    std::ostringstream out;
    out << "**ERROR** Multiple material managers created";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
}

MaterialManager::~MaterialManager()
{
  clearMaterials();

  for (size_t i = 0; i < d_materials_old.size(); i++) {
    delete d_materials_old[i];
  }

  if (d_all_materials_old && d_all_materials_old->removeReference()) {
    delete d_all_materials_old;
  }
}

void
MaterialManager::clearMaterials()
{
  for (size_t i = 0; i < d_all_materials.size(); i++) {
    d_materials_old.push_back(d_all_materials[i]);
  }

  if (d_all_materialsets && d_all_materialsets->removeReference()) {
    delete d_all_materialsets;
  }
  d_all_materialsets = nullptr;

  if (d_all_in_one_material && d_all_in_one_material->removeReference()) {
    delete d_all_in_one_material;
  }
  d_all_in_one_material = nullptr;

  d_all_materials.clear();
  d_empty_materials.clear();
  d_named_materials.clear();

  for (auto& material : d_materials) {
    std::string name = material.first;
    if (d_materialsets[name] && d_materialsets[name]->removeReference()) {
      delete d_materialsets[name];
    }
    d_materials[name].clear();
    d_materialsets[name] = nullptr;
  }
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
  std::vector<int> all_materials(d_all_materials.size());
  for (size_t i = 0; i < d_all_materials.size(); i++) {
    all_materials[i] = d_all_materials[i]->getDWIndex();
  }
  d_all_materialsets->addAll(all_materials);

  // Original all_materials (before switch using Switcher)
  if (d_all_materials_old == nullptr) {
    d_all_materials_old = scinew MaterialSet();
    d_all_materials_old->addReference();
    d_all_materials_old->addAll(all_materials);
  }

  // All materials in one:  a material that represents all materials
  // (i.e. summed over all materials -- the whole enchilada)
  if (d_all_in_one_material && d_all_in_one_material->removeReference()) {
    delete d_all_in_one_material;
  }
  d_all_in_one_material = scinew MaterialSubset();
  d_all_in_one_material->addReference();
  d_all_in_one_material->add((int)d_all_materials.size());

  for (auto& material : d_materialsets) {
    if (material.second && material.second->removeReference()) {
      delete material.second;
    }
  }

  for (auto& material : d_materials) {
    std::string name = material.first;
    d_materialsets[name] = scinew MaterialSet();
    d_materialsets[name]->addReference();

    std::vector<int> materials(d_materials[name].size());
    for (size_t i = 0; i < d_materials[name].size(); i++) {
      materials[i] = d_materials[name][i]->getDWIndex();
    }
    d_materialsets[name]->addAll(materials);
  }
}

void
MaterialManager::registerMaterial(const std::string& name, Material* material)
{
  d_materials[name].push_back(material);
  registerMaterial(material);
}

void
MaterialManager::registerMaterial(const std::string& name, Material* material,
                                  unsigned int index)
{
  d_materials[name].push_back(material);
  registerMaterial(material, index);
}

void
MaterialManager::registerMaterial(Material* material)
{
  material->setDWIndex((int)d_all_materials.size());
  d_all_materials.push_back(material);
  if (material->hasName()) {
    d_named_materials[material->getName()] = material;
  }
}

void
MaterialManager::registerMaterial(Material* material, unsigned int index)
{
  material->setDWIndex(index);

  if (d_all_materials.size() <= index) {
    d_all_materials.resize(index + 1);
  }
  d_all_materials[index] = material;

  if (material->hasName()) {
    d_named_materials[material->getName()] = material;
  }
}

void
MaterialManager::registerEmptyMaterial(EmptyMaterial* material)
{
  d_empty_materials.push_back(material);
  registerMaterial(material);
}

const MaterialSet*
MaterialManager::allMaterials(const std::string& name) const
{
  if (d_materialsets.find(name) != d_materialsets.end()) {
    return d_materialsets.at(name);
  }
  return nullptr;
}

const MaterialSet*
MaterialManager::allMaterials() const
{
  ASSERT(d_all_materialsets != nullptr);
  return d_all_materialsets;
}

size_t
MaterialManager::getNumMaterials(const std::string& name) const
{
  if (d_materials.find(name) != d_materials.end()) {
    return d_materials.at(name).size();
  }
  return 0;
}

Material*
MaterialManager::getMaterial(const std::string& name, int index) const
{
  if (d_materials.find(name) != d_materials.end()) {
    return d_materials.at(name)[index];
  }
  return nullptr;
}

Material*
MaterialManager::getMaterialByName(const std::string& name) const
{
  auto iter = d_named_materials.find(name);
  if (iter != d_named_materials.end()) {
    return iter->second;
  }
  return nullptr;
}

const MaterialSet*
MaterialManager::originalAllMaterials() const
{
  ASSERT(d_all_materials_old != nullptr);
  return d_all_materials_old;
}

void
MaterialManager::setOriginalMaterialsFromRestart(MaterialSet* materials)
{
  if (d_all_materials_old && d_all_materials_old->removeReference())
    delete d_all_materials_old;
  d_all_materials_old = materials;
}

Material*
MaterialManager::parseAndLookupMaterial(ProblemSpecP& params,
                                        const std::string& name) const
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


} // end namespace Uintah