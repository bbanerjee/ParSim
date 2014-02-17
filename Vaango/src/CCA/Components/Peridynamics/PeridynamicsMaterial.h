#ifndef __VAANGO_PERIDYNAMICS_MATERIAL_H__
#define __VAANGO_PERIDYNAMICS_MATERIAL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>

#include <Core/Grid/Material.h>
#include <Core/Grid/SimulationState.h>
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
  class ParticleCreator;

  class PeridynamicsMaterial : public Uintah::Material 
  {

  public:

    // Default Constructor
    PeridynamicsMaterial();

    // Standard Peridynamics Material Constructor
    PeridynamicsMaterial(Uintah::ProblemSpecP& ps, 
                         Uintah::SimulationStateP& ss, 
                         PeridynamicsFlags* flags);
         
   ~PeridynamicsMaterial();

   virtual void registerParticleState(Uintah::SimulationState* ss);

   virtual Uintah::ProblemSpecP outputProblemSpec(Uintah::ProblemSpecP& ps);

   // Return correct constitutive model pointer for this material
   PeridynamicsMaterialModel* getMaterialModel() const;

   // Return correct basic damage model pointer for this material
   PeridynamicsFailureModel* getFailureModel() const;

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
   PeridynamicsFailureModel* d_failureModel;
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
   void standardInitialization(Uintah::ProblemSpecP& ps, PeridynamicsFlags* flags);
 };

} // End namespace Vaango

#endif // __VAANGO_PERIDYNAMICS_MATERIAL_H__
