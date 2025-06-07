/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/ProblemSpec/ProblemSpec.h>

namespace Uintah {

template<class T>
std::map<int, T>
MPMCommon::initializeMap(const T& val)
{
  std::map<int, T> myMap;
  int numMPMMatls = s_materialManager->getNumMaterials("MPM");

  for (int m = 0; m < numMPMMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(s_materialManager->getMaterial("MPM", m));

    int dwi    = mpm_matl->getDWIndex();
    myMap[dwi] = val;
  }
  return myMap;
}

template std::map<int, double>
MPMCommon::initializeMap(const double& val);

template std::map<int, Vector>
MPMCommon::initializeMap(const Vector& val);

MPMCommon::MPMCommon(const MaterialManagerP matManager)
{
  d_lb              = std::make_shared<MPMLabel>();
  s_materialManager = matManager;
}

void
MPMCommon::materialProblemSetup(const ProblemSpecP& prob_spec,
                                MPMFlags* flags,
                                bool is_restart,
                                const std::string& input_ups_dir)
{
  d_flags = flags;

  // so all components can know how many particle ghost cells to ask for
  d_flags->particle_ghost_type  = d_particle_ghost_type;
  d_flags->particle_ghost_layer = d_particle_ghost_layer;

  // Search for the MaterialProperties block and then get the MPM section
  ProblemSpecP mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");
  ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");
  for (ProblemSpecP ps = mpm_mat_ps->findBlock("material"); ps != nullptr;
       ps              = ps->findNextBlock("material")) {
    std::string index("");
    ps->getAttribute("index", index);

    const int DEFAULT_VALUE = -1;
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

    // Create and register as an MPM material
    std::shared_ptr<MPMMaterial> mat = std::make_shared<MPMMaterial>(
      ps, s_materialManager, d_flags, is_restart, input_ups_dir);

    // Add particle state
    mat->registerParticleState(d_particleState, d_particleState_preReloc);

    // When doing restart, we need to make sure that we load the materials
    // in the same order that they were initially created.  Restarts will
    // ALWAYS have an index number as in <material index = "0">.
    // Index_val = -1 means that we don't register the material by its
    // index number.
    if (index_val > -1) {
      s_materialManager->registerMaterial("MPM", mat, index_val);
    } else {
      s_materialManager->registerMaterial("MPM", mat);
    }

    // If new particles are to be created, create a copy of each material
    // without the associated geometry
    if (flags->d_createNewParticles) {
      std::shared_ptr<MPMMaterial> mat_copy = std::make_shared<MPMMaterial>();
      mat_copy->copyWithoutGeom(ps, mat.get(), d_flags);
      s_materialManager->registerMaterial("MPM", mat_copy);
    }
  }
}

} // end namespace Uintah