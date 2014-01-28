#ifndef __VAANGO_PERIDYNAMICS_MATERIAL_H__
#define __VAANGO_PERIDYNAMICS_MATERIAL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsSimulationState.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>

#include <Core/Grid/Material.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <vector>

namespace Vaango 
{

  class Uintah::Patch;
  class Uintah::DataWarehouse;
  class Uintah::GeometryObject;
  class Uintah::ParticleCreator;

  class PeridynamicsMaterial : public Uintah::Material 
  {

  public:

    // Default Constructor
    PeridynamicsMaterial();

    // Standard Peridynamics Material Constructor
    PeridynamicsMaterial(Uintah::ProblemSpecP& ps, 
                         PeridynamicsSimulationStateP& ss, 
                         PeridynamicsFlags* flags);
         
   ~PeridynamicsMaterial();

   virtual void registerParticleState(PeridynamicsSimulationState* ss);

   virtual Uintah::ProblemSpecP outputProblemSpec(Uintah::ProblemSpecP& ps);

   // Return correct constitutive model pointer for this material
   PeridynamicsMaterialModel* getMaterialModel() const;

   // Return correct basic damage model pointer for this material
   PeridynamicsFailureModel* getFailureModel() const;

   Uintah::particleIndex countParticles(const Patch* patch);
   void createParticles(Uintah::particleIndex numParticles,
                        Uintah::CCVariable<short int>& cellNAPID,
                        const Uintah::Patch* patch,
                        Uintah::DataWarehouse* new_dw);

   Uintah::ParticleCreator* getParticleCreator();

 private:

   PeridynamicsLabel* d_varLabel;
   PeridynamicsMaterialModel* d_materialModel;
   PeridynamicsFailureModel* d_failureModel;
   Uintah::ParticleCreator* d_particle_creator;
   std::vector<GeometryObject*> d_geom_objs;

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
