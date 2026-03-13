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

#include <CCA/Components/MPM/DEM/DEMRigidBodySpatialIndex.h>
#include <CCA/Components/MPM/DEM/DEMCommon.h>
#include <Core/Util/DOUT.hpp>

using namespace Uintah;

// Debug streams:  To turn on use  export SCI_DEBUG="DEM:+,SerialMPM:+" 
extern Dout dem_doing;
extern Dout dem_doing_fn;
extern Dout dem_dbg;

namespace Vaango {

// -----------------------------------------------------------------------
// Build: O(n_j) — bucket each mat_j particle into its patch cell,
// exactly as the MPM numInCell pattern does.
// -----------------------------------------------------------------------
DEMRigidBodySpatialIndex::DEMRigidBodySpatialIndex(
  int mat_j,
  const DEMParticleSets& psets,
  const DEMParticleInputData& inputs,
  const Uintah::Patch* patch,
  bool is_dem_material)
  : d_mat_j(mat_j)
  , d_is_dem(is_dem_material)
  , d_inputs(inputs)
  , d_patch(patch)
{
  for (auto pidx_j : *psets.psets_all[mat_j]) {
    IntVector c;
    d_patch->findCell(inputs.pX_old[mat_j][pidx_j], c);
    d_cell_to_particles[c].push_back(pidx_j);

    if (!is_dem_material) {
      d_non_dem_full_map[(Uintah::long64)pidx_j] = pidx_j;
    }
  }
}

// -----------------------------------------------------------------------
// Query: walk the (2*cell_radius+1)^3 neighbourhood around pos_i's cell,
// applying identical skip logic to the original per-particle loops.
// -----------------------------------------------------------------------
ParticleIDToCurrentIdxMap
DEMRigidBodySpatialIndex::query(const Uintah::Point& pos_i,
                                int mat_i,
                                Uintah::particleIndex pidx_i,
                                int cell_radius) const
{
  // Cross-material non-DEM: no per-particle filtering needed — O(1)
  if (!d_is_dem && mat_i != d_mat_j) {
    return d_non_dem_full_map;
  }

  // Find the cell that pos_i lives in, matching patch->findCell convention
  IntVector center;
  d_patch->findCell(pos_i, center);

  const Uintah::long64 rbID_i =
    d_is_dem ? d_inputs.pRigidBodyID_old[mat_i][pidx_i] : -1;

  ParticleIDToCurrentIdxMap unique_bodies;

  for (int dx = -cell_radius; dx <= cell_radius; ++dx) {
    for (int dy = -cell_radius; dy <= cell_radius; ++dy) {
      for (int dz = -cell_radius; dz <= cell_radius; ++dz) {

        IntVector cell = center + IntVector(dx, dy, dz);

        auto cell_it = d_cell_to_particles.find(cell);
        if (cell_it == d_cell_to_particles.end()) {
          continue;
        }

        for (auto pidx_j : cell_it->second) {

          // ---- same skip logic as the original buildUniqueBodies ----
          if (d_is_dem) {
            if (mat_i == d_mat_j &&
                (pidx_j <= pidx_i ||
                 rbID_i == d_inputs.pRigidBodyID_old[d_mat_j][pidx_j])) {
              // DOUT(dem_dbg, "[DEMRigidBodySpatialIndex::query] Skipping mat_i=" 
              //      << mat_i << " pidx_i=" << pidx_i 
              //      << " mat_j=" << d_mat_j << " pidx_j=" << pidx_j);
              continue;
            }
          } else {
            // Same-material non-DEM: only process upper triangle
            if (mat_i == d_mat_j && pidx_i <= pidx_j) {
              continue;
            }
          }
          // -----------------------------------------------------------

          if (!d_is_dem) {
            unique_bodies[(Uintah::long64)pidx_j] = pidx_j;
            continue;
          }

          // DEM: keep the closest representative per rigid body
          const Uintah::long64 rbID_j = d_inputs.pRigidBodyID_old[d_mat_j][pidx_j];
          DOUT(dem_dbg, "[DEMRigidBodySpatialIndex::query] " 
                        << "mat_i=" << mat_i << " mat_j=" << d_mat_j 
                        << " pidx_i=" << pidx_i << " pidx_j=" << pidx_j 
                        << " rbID[j]=" << rbID_j);
          auto rb_it          = unique_bodies.find(rbID_j);
          if (rb_it == unique_bodies.end()) {
            unique_bodies[rbID_j] = pidx_j;
          } else {
            double d2_new =
              (pos_i - d_inputs.pX_old[d_mat_j][pidx_j]).length2();
            double d2_old =
              (pos_i - d_inputs.pX_old[d_mat_j][rb_it->second]).length2();
            if (d2_new < d2_old) {
              rb_it->second = pidx_j;
            }
          }
        }
      }
    }
  }
  return unique_bodies;
}

void 
DEMRigidBodySpatialIndex::print(std::ostream& out) const
{
  out << "DEMRigidBodySpatialIndex:"
      << " mat_j=" << d_mat_j
      << " is_dem=" << std::boolalpha << d_is_dem
      << " num_cells=" << d_cell_to_particles.size()
      << "\n";

  // Sort cells for deterministic output
  std::vector<IntVector> cells;
  cells.reserve(d_cell_to_particles.size());
  for (const auto& [cell, _] : d_cell_to_particles) {
    cells.push_back(cell);
  }
  std::sort(cells.begin(), cells.end(), [](const IntVector& a, const IntVector& b) {
    if (a.x() != b.x()) return a.x() < b.x();
    if (a.y() != b.y()) return a.y() < b.y();
    return a.z() < b.z();
  });

  for (const auto& cell : cells) {
    const auto& particles = d_cell_to_particles.at(cell);
    out << "  cell=(" << cell.x() << "," << cell.y() << "," << cell.z() << ")"
        << " num_particles=" << particles.size() << "\n";
    for (auto pidx_j : particles) {
      out << "    pidx_j=" << pidx_j;
      if (d_is_dem) {
        out << " rbID=" << d_inputs.pRigidBodyID_old[d_mat_j][pidx_j]
            << " pos="  << d_inputs.pX_old[d_mat_j][pidx_j];
      }
      out << "\n";
    }
  }

  if (!d_is_dem && !d_non_dem_full_map.empty()) {
    out << "  non_dem_full_map size=" << d_non_dem_full_map.size() << "\n";
    for (const auto& [rbID, pidx] : d_non_dem_full_map) {
      out << "    rbID=" << rbID << " pidx=" << pidx << "\n";
    }
  }
}

std::ostream& 
operator<<(std::ostream& out, const DEMRigidBodySpatialIndex& idx)
{
  idx.print(out);
  return out;
}

} // end namespace Vaango
