/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __VAANGO_PERIDYNAMICS_MATERIAL_H__
#define __VAANGO_PERIDYNAMICS_MATERIAL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>

#include <Core/Grid/Material.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <vector>

namespace Uintah {
  class Patch;
  class DataWarehouse;
  class GeometryObject;
}

namespace Vaango 
{
  class PeridynamicsMaterial : public Uintah::Material 
  {

  public:

    // Default Constructor
    PeridynamicsMaterial();

    // Standard Peridynamics Material Constructor
    PeridynamicsMaterial(Uintah::ProblemSpecP& ps, 
                         const Uintah::GridP grid,
                         Uintah::SimulationStateP& ss, 
                         PeridynamicsFlags* flags);
         
   ~PeridynamicsMaterial();

   virtual void registerParticleState(Uintah::SimulationState* ss);

   virtual Uintah::ProblemSpecP outputProblemSpec(Uintah::ProblemSpecP& ps);

   // Return correct constitutive model pointer for this material
   PeridynamicsMaterialModel* getMaterialModel() const;

   // Return correct basic damage model pointer for this material
   PeridynamicsDamageModel* getDamageModel() const;

   double getInitialDensity() const {return d_density;}
   
   Uintah::particleIndex countParticles(const Uintah::Patch* patch);
   void createParticles(Uintah::particleIndex numParticles,
                        Uintah::CCVariable<short int>& cellNAPID,
                        const Uintah::Patch* patch,
                        Uintah::DataWarehouse* new_dw);

   ParticleCreator* getParticleCreator();

 private:

   PeridynamicsLabel* d_varLabel;
   PeridynamicsMaterialModel* d_materialModel;
   PeridynamicsDamageModel* d_damageModel;
   ParticleCreator* d_particle_creator;

   double d_density;

   std::vector<Uintah::GeometryObject*> d_geom_objs;

   // Prevent copying of this class
   // copy constructor
   PeridynamicsMaterial(const PeridynamicsMaterial &mpmm);
   PeridynamicsMaterial& operator=(const PeridynamicsMaterial &mpmm);

   ///////////////////////////////////////////////////////////////////////////
   //
   // The standard set of initialization actions except particlecreator
   //
   void standardInitialization(Uintah::ProblemSpecP& ps, 
                               Uintah::GridP grid,
                               PeridynamicsFlags* flags);
 };

} // End namespace Vaango

#endif // __VAANGO_PERIDYNAMICS_MATERIAL_H__
