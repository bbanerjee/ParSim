/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/ParticleCreator/SGPDataManager.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreatorStructs.h>
#include <format>
#include <iostream>

namespace Vaango::Helpers {

void SGPDataManager::initialize(Uintah::SpecialGeomPiece* sgp,
                                const Uintah::ObjectVars& obj_vars,
                                Uintah::GeometryObject* obj,
                                bool with_color,
                                bool do_scalar_diffusion,
                                bool with_gauss_solver) {
  if (!sgp) {
    d_has_data = false;
    return;
  }
  
  d_has_data = true;
  
  // Initialize core particle properties
  d_volume.initialize(sgp, obj_vars, "p.volume", obj);
  d_temperature.initialize(sgp, obj_vars, "p.temperature", obj);
  d_external_force.initialize(sgp, obj_vars, "p.externalforce", obj);
  d_fiber_dir.initialize(sgp, obj_vars, "p.fiberdirs", obj);
  d_velocity.initialize(sgp, obj_vars, "p.velocity", obj);
  d_size.initialize(sgp, obj_vars, "p.size", obj);
  
  // Initialize optional properties based on simulation flags
  if (with_color) {
    d_color.initialize(sgp, obj_vars, "p.color", obj);
  }
  
  if (do_scalar_diffusion) {
    d_concentration.initialize(sgp, obj_vars, "p.concentration", obj);
    d_area.initialize(sgp, obj_vars, "p.area", obj);
  }
  
  if (with_gauss_solver) {
    d_pos_charge.initialize(sgp, obj_vars, "p.poscharge", obj);
    d_neg_charge.initialize(sgp, obj_vars, "p.negcharge", obj);
    d_permittivity.initialize(sgp, obj_vars, "p.permittivity", obj);
  }
}

} // end namespace Vaango::Helpers