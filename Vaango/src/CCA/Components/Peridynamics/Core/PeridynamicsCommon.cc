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

#include <CCA/Components/Peridynamics/Core/PeridynamicsCommon.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/ProblemSpec/ProblemSpec.h>

namespace Vaango {

template<class T>
std::map<int, T>
PeridynamicsCommon::initializeMap(const T& val)
{
  std::map<int, T> myMap;
  int numPeridynamicsMatls = s_mat_manager->getNumMaterials("Peridynamics");

  for (int m = 0; m < numPeridynamicsMatls; m++) {
    PeridynamicsMaterial* mpm_matl = static_cast<PeridynamicsMaterial*>(
      s_mat_manager->getMaterial("Peridynamics", m));

    int dwi    = mpm_matl->getDWIndex();
    myMap[dwi] = val;
  }
  return myMap;
}

template<>
std::map<int, double>
PeridynamicsCommon::initializeMap(const double& val);

template<>
std::map<int, Uintah::Vector>
PeridynamicsCommon::initializeMap(const Uintah::Vector& val);

PeridynamicsCommon::PeridynamicsCommon(
  const Uintah::MaterialManagerP matManager)
{
  d_pd_labels   = std::make_unique<PeridynamicsLabel>();
  s_mat_manager = matManager;
}

void
PeridynamicsCommon::materialProblemSetup(const Uintah::ProblemSpecP& prob_spec,
                                         PeridynamicsFlags* flags,
                                         bool is_restart)
{
  d_pd_flags = flags;

  // so all components can know how many particle ghost cells to ask for
  d_pd_flags->particle_ghost_type  = d_particle_ghost_type;
  d_pd_flags->particle_ghost_layer = d_particle_ghost_layer;

  // Search for the MaterialProperties block and then get the Peridynamics
  // section
  Uintah::ProblemSpecP mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");
  Uintah::ProblemSpecP mpm_mat_ps = mat_ps->findBlock("Peridynamics");
  for (Uintah::ProblemSpecP ps = mpm_mat_ps->findBlock("material");
       ps != nullptr;
       ps = ps->findNextBlock("material")) {
    std::string index("");
    ps->getAttribute("index", index);

    constexpr int DEFAULT_VALUE = -1;
    std::stringstream id(index);
    int index_val = DEFAULT_VALUE;
    id >> index_val;

    if (!id) {
      //  std::stringstream parsing failed... on many (most) systems, the
      // original value assigned to index_val would be left
      // intact... but on some systems (redstorm) it inserts garbage,
      // so we have to manually restore the value.
      index_val = DEFAULT_VALUE;
    }
    // std::cout << "Material attribute = " << index_val << ", " << index << ",
    // " << id << "\n";

    // Create and register as an Peridynamics material
    std::shared_ptr<PeridynamicsMaterial> mat =
      std::make_shared<PeridynamicsMaterial>(
        ps, s_mat_manager, d_pd_flags, is_restart);

    // When doing restart, we need to make sure that we load the materials
    // in the same order that they were initially created.  Restarts will
    // ALWAYS have an index number as in <material index = "0">.
    // Index_val = -1 means that we don't register the material by its
    // index number.
    if (index_val > -1) {
      s_mat_manager->registerMaterial("Peridynamics", mat, index_val);
    } else {
      s_mat_manager->registerMaterial("Peridynamics", mat);
    }

    // If new particles are to be created, create a copy of each material
    // without the associated geometry
    if (d_pd_flags->d_createNewParticles) {
      std::shared_ptr<PeridynamicsMaterial> mat_copy =
        std::make_shared<PeridynamicsMaterial>();
      mat_copy->copyWithoutGeom(ps, mat.get(), d_pd_flags);
      s_mat_manager->registerMaterial("Peridynamics", mat_copy);
    }
  }
}

} // namespace Vaango