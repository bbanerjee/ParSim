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

// Debug streams:  To turn on use  export SCI_DEBUG="DEMdbg:+,SerialMPM:+" 
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
  bool j_is_dem_material)
  : d_mat_j(mat_j)
  , d_j_is_dem(j_is_dem_material)
  , d_inputs(inputs)
  , d_patch(patch)
{
  for (auto pidx_j : *psets.psets_all[mat_j]) {
    IntVector c;
    d_patch->findCell(inputs.pX_old[mat_j][pidx_j], c);
    d_cell_to_particles[c].push_back(pidx_j);

    if (!d_j_is_dem) {
      d_non_dem_full_map[(Uintah::long64)pidx_j] = pidx_j;
    }
  }
}

// -----------------------------------------------------------------------
// Query: walk the (2*cell_radius+1)^3 neighbourhood around pos_i's cell,
// Returns a map from DEM ID -> representative particle index for
// all j-material bodies that particle pidx_i should interact with. 
// -----------------------------------------------------------------------
DEMBodyIDToCurrentParticleIdxMap
DEMRigidBodySpatialIndex::query(const Uintah::Point& pos_i,
                                int mat_i,
                                Uintah::particleIndex pidx_i,
                                int cell_radius) const
{
  // DEM ID for particle i
  const Uintah::long64 pDEMID_i = d_inputs.pDEMBodyID_old[mat_i][pidx_i];

  // Find the cell that pos_i lives in
  IntVector center;
  d_patch->findCell(pos_i, center);

  DEMBodyIDToCurrentParticleIdxMap close_particles;

  for (int dx = -cell_radius; dx <= cell_radius; ++dx) {
    for (int dy = -cell_radius; dy <= cell_radius; ++dy) {
      for (int dz = -cell_radius; dz <= cell_radius; ++dz) {

        IntVector cell = center + IntVector(dx, dy, dz);
        auto cell_it = d_cell_to_particles.find(cell);
        if (cell_it == d_cell_to_particles.end()) {
          continue;
        }

        for (auto pidx_j : cell_it->second) {

          // DEM ID for particle j
          const auto pDEMID_j = d_inputs.pDEMBodyID_old[d_mat_j][pidx_j];
          bool same_dem_body_id = (pDEMID_i == pDEMID_j);
          // if (mat_i != d_mat_j) {
          //   DOUT(dem_dbg, "[DEMRigidBodySpatialIndex::query] DEM IDs: "
          //                << "i: [mat, pidx, demID]=" << mat_i << " " << pidx_i << " " << pDEMID_i
          //                << " j: [mat, pidx, demID]=" << d_mat_j << " " << pidx_j << " " << pDEMID_j);
          // }

          if (same_dem_body_id) {
            continue;
          }

          if (!d_j_is_dem) {
            close_particles[(Uintah::long64)pidx_j] = pidx_j;
            continue;
          }

          // DEM: keep the closest representative per rigid body
          // DOUT(dem_dbg, "[DEMRigidBodySpatialIndex::query] " 
          //               << "mat_i=" << mat_i << " mat_j=" << d_mat_j 
          //               << " pidx_i=" << pidx_i << " pidx_j=" << pidx_j 
          //               << " pDEMID[j]=" << pDEMID_j);
          auto rb_it          = close_particles.find(pDEMID_j);
          if (rb_it == close_particles.end()) {
            close_particles[pDEMID_j] = pidx_j;
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
  return close_particles;
}

void 
DEMRigidBodySpatialIndex::print(std::ostream& out) const
{
  out << "DEMRigidBodySpatialIndex:"
      << " mat_j=" << d_mat_j
      << " j_is_dem=" << std::boolalpha << d_j_is_dem
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
      if (d_j_is_dem) {
        out << " pDEMBodyID=" << d_inputs.pDEMBodyID_old[d_mat_j][pidx_j]
            << " pos="  << d_inputs.pX_old[d_mat_j][pidx_j];
      }
      out << "\n";
    }
  }

  if (!d_j_is_dem && !d_non_dem_full_map.empty()) {
    out << "  non_dem_full_map size=" << d_non_dem_full_map.size() << "\n";
    for (const auto& [pDEMID, pidx] : d_non_dem_full_map) {
      out << "    pDEMID=" << pDEMID << " pidx=" << pidx << "\n";
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
