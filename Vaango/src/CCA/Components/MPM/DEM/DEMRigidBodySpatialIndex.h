/*
 * The MIT License
 *
 * Copyright (c) 2013-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __VAANGO_CCA_COMPONENTS_MPM_DEM_SPATIAL_INDEX_H__
#define __VAANGO_CCA_COMPONENTS_MPM_DEM_SPATIAL_INDEX_H__

#include <CCA/Components/MPM/DEM/DEMCommon.h>
#include <unordered_map>

namespace Uintah {
class Patch;
}

namespace Vaango {

class DEMRigidBodySpatialIndex
{
public:
  DEMRigidBodySpatialIndex(int mat_j,
                           const DEMParticleSets& psets,
                           const DEMParticleInputData& inputs,
                           const Uintah::Patch* patch,
                           bool is_dem_material);

  ParticleIDToCurrentIdxMap
  query(const Uintah::Point& pos_i,
        int mat_i,
        Uintah::particleIndex pidx_i,
        int cell_radius = 1) const;

private:
  int d_mat_j;
  bool d_is_dem;
  const DEMParticleInputData& d_inputs;
  const Uintah::Patch* d_patch;

  // cell -> particles whose pX falls in that cell
  std::unordered_map<Uintah::IntVector,
                     std::vector<Uintah::particleIndex>,
                     IntVectorHash> d_cell_to_particles;

  // Non-DEM cross-material only: identity map built once, returned in O(1)
  ParticleIDToCurrentIdxMap d_non_dem_full_map;
};

} // namespace Vaango

#endif // __VAANGO_CCA_COMPONENTS_MPM_DEM_SPATIAL_INDEX_H__