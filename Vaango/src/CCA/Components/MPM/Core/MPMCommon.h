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

#ifndef __CCA_COMPONENTS_MPM_MPMCOMMON_H__
#define __CCA_COMPONENTS_MPM_MPMCOMMON_H__

#include <Core/Grid/Ghost.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <map>
#include <memory>
#include <vector>

namespace Uintah {

class ProcessorGroup;

class MPMFlags;
class MPMLabel;

class MPMCommon
{
public:
  template<class T>
  static std::map<int, T>
  initializeMap(const T& val);

public:
  MPMCommon(const MaterialManagerP materialManager);
  ~MPMCommon() {}

  virtual void
  materialProblemSetup(const ProblemSpecP& prob_spec,
                       MPMFlags* flags,
                       bool is_restart,
                       const std::string& input_ups_dir = "");

  // Used by the switcher
  virtual void
  setupForSwitching()
  {
    d_particleState.clear();
    d_particleState_preReloc.clear();
  }

  inline void
  setParticleGhostLayer(Ghost::GhostType type, int num_ghost_cells)
  {
    d_particle_ghost_type  = type;
    d_particle_ghost_layer = num_ghost_cells;
  }

  inline void
  getParticleGhostLayer(Ghost::GhostType& type, int& num_ghost_cells)
  {
    type            = d_particle_ghost_type;
    num_ghost_cells = d_particle_ghost_layer;
  }

public:
  // Particle state
  std::vector<std::vector<const VarLabel*>> d_particleState;
  std::vector<std::vector<const VarLabel*>> d_particleState_preReloc;

  std::shared_ptr<MPMLabel> d_lb{ nullptr };

protected:
  //! so all components can know how many particle ghost cells to ask for
  Ghost::GhostType d_particle_ghost_type{ Ghost::None };
  int d_particle_ghost_layer{ 0 };

private:
  inline static MaterialManagerP s_materialManager{ nullptr };

  MPMFlags* d_flags{ nullptr };
};
} // namespace Uintah

#endif //__CCA_COMPONENTS_MPM_MPMCOMMON_H__
