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

#ifndef SGP_DATA_MANAGER_H
#define SGP_DATA_MANAGER_H

#include <CCA/Components/MPM/ParticleCreator/SGPData.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/SpecialGeomPiece.h>
#include <memory>

namespace Vaango::Helpers {

/**
 * Manages all special geometry piece particle data in one place.
 * Replaces the need for 12+ separate SGPData declarations.
 */
class SGPDataManager {
public:
  SGPDataManager() = default;
  
  /**
   * Initialize all SGP data structures for a special geometry piece.
   * Only initializes data that's actually present in the geometry piece.
   * 
   * @param sgp The special geometry piece
   * @param obj_vars The object variables containing particle data
   * @param obj The geometry object
   * @param with_color Whether to initialize color data
   * @param do_scalar_diffusion Whether to initialize diffusion data
   * @param with_gauss_solver Whether to initialize electrostatics data
   */
  void initialize(Uintah::SpecialGeomPiece* sgp,
                  const Uintah::ObjectVars& obj_vars,
                  Uintah::GeometryObject* obj,
                  bool with_color,
                  bool do_scalar_diffusion,
                  bool with_gauss_solver);
  
  /**
   * Check if this geometry piece has any special data.
   */
  bool hasData() const { return d_has_data; }
  
  // Accessor methods that return data and advance iterators
  auto getVolume() -> std::optional<double> { return d_volume.get_and_advance(); }
  auto getTemperature() -> std::optional<double> { return d_temperature.get_and_advance(); }
  auto getExternalForce() -> std::optional<Vector> { return d_external_force.get_and_advance(); }
  auto getFiberDir() -> std::optional<Uintah::Vector> { return d_fiber_dir.get_and_advance(); }
  auto getVelocity() -> std::optional<Uintah::Vector> { return d_velocity.get_and_advance(); }
  auto getSize() -> std::optional<Uintah::Matrix3> { return d_size.get_and_advance(); }
  auto getColor() -> std::optional<double> { return d_color.get_and_advance(); }
  auto getArea() -> std::optional<Uintah::Vector> { return d_area.get_and_advance(); }
  auto getConcentration() -> std::optional<double> { return d_concentration.get_and_advance(); }
  auto getPosCharge() -> std::optional<double> { return d_pos_charge.get_and_advance(); }
  auto getNegCharge() -> std::optional<double> { return d_neg_charge.get_and_advance(); }
  auto getPermittivity() -> std::optional<double> { return d_permittivity.get_and_advance(); }

private:
  bool d_has_data = false;
  
  SGPData<double> d_volume;
  SGPData<double> d_temperature;
  SGPData<Uintah::Vector> d_external_force;
  SGPData<Uintah::Vector> d_fiber_dir;
  SGPData<Uintah::Vector> d_velocity;
  SGPData<Uintah::Matrix3> d_size;
  SGPData<double> d_color;
  SGPData<Uintah::Vector> d_area;
  SGPData<double> d_concentration;
  SGPData<double> d_pos_charge;
  SGPData<double> d_neg_charge;
  SGPData<double> d_permittivity;
};

} // end namespace Vaango::Helpers

#endif // SGP_DATA_MANAGER_H
