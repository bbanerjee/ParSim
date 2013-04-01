#ifndef MATITI_PERIMATERIAL_H
#define MATITI_PERIMATERIAL_H

namespace Uintah {

using namespace SCIRun;

 class MPMMaterial : public Material {
 public:

   // Default Constructor
   MPMMaterial();

   // Standard MPM Material Constructor
   MPMMaterial(ProblemSpecP&, SimulationStateP& ss, MPMFlags* flags);
         
   ~MPMMaterial();

   virtual void registerParticleState(SimulationState* ss);

   virtual ProblemSpecP outputProblemSpec(ProblemSpecP& ps);

   /*!  Create a copy of the material without the associated geometry */
   void copyWithoutGeom(ProblemSpecP& ps,const MPMMaterial* mat,
                        MPMFlags* flags);
         
   //////////
   // Return correct constitutive model pointer for this material
   ConstitutiveModel* getConstitutiveModel() const;

   //////////
   // Return correct basic damage model pointer for this material
   Vaango::BasicDamageModel* getBasicDamageModel() const;

   // Return correct burn model pointer for this material
   particleIndex countParticles(const Patch* patch);

   void createParticles(particleIndex numParticles,
                        CCVariable<short int>& cellNAPID,
                        const Patch*,
                        DataWarehouse* new_dw);


   ParticleCreator* getParticleCreator();

   double getInitialDensity() const;

   // Get the specific heats at room temperature
   double getInitialCp() const;
   double getInitialCv() const;

   // for temperature dependent plasticity models
   double getRoomTemperature() const;
   double getMeltTemperature() const;

   bool getIsRigid() const;

   bool getIncludeFlowWork() const;
   double getSpecificHeat() const;
   double getThermalConductivity() const;

   int nullGeomObject() const;


   // For MPMICE
   double getGamma() const;
   void initializeCCVariables(CCVariable<double>& rhom,
                              CCVariable<double>& rhC,
                              CCVariable<double>& temp,   
                              CCVariable<Vector>& vCC,
                              CCVariable<double>& vfCC,
                              const Patch* patch);

   void initializeDummyCCVariables(CCVariable<double>& rhom,
                                   CCVariable<double>& rhC,
                                   CCVariable<double>& temp,   
                                   CCVariable<Vector>& vCC,
                                   CCVariable<double>& vfCC,
                                   const Patch* patch);

   bool d_doBasicDamage;

 private:

   MPMLabel* d_lb;
   ConstitutiveModel* d_cm;
   ParticleCreator* d_particle_creator;

   double d_density;
   bool d_includeFlowWork;
   double d_specificHeat;
   double d_thermalConductivity;

   // Specific heats at constant pressure and constant volume
   // (values at room temperature - [273.15 + 20] K)
   double d_Cp, d_Cv;

   // for temperature dependent plasticity models
   double d_troom;
   double d_tmelt;

   // for implicit rigid body contact
   bool d_is_rigid;

   // For basic damage computations
   Vaango::BasicDamageModel* d_basicDamageModel;

   std::vector<GeometryObject*> d_geom_objs;

   // Prevent copying of this class
   // copy constructor
   MPMMaterial(const MPMMaterial &mpmm);
   MPMMaterial& operator=(const MPMMaterial &mpmm);

   ///////////////////////////////////////////////////////////////////////////
   //
   // The standard set of initialization actions except particlecreator
   //
   void standardInitialization(ProblemSpecP& ps, MPMFlags* flags);
 };

} // End namespace Uintah

#endif // MATITI_PERIMATERIAL_H
