#ifndef EMU2DC_MATERIAL_H
#define EMU2DC_MATERIAL_H

#include <DamageModel.h>
#include <DamageModelUP.h>

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <memory>
#include <iostream>

namespace Emu2DC {

  class Material 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Material& dam);

  public:

    enum class MicroModulusModel {
      Constant=0,
      Conical=1
    };
  
    /**
     * Constructor/destructor of a bond material
     * (Use the clone method rather than a copy contructor or operator=)
     */
    Material();
    virtual ~Material();
    void clone(const Material* mat);

    /**
     *  Initialize material properties of the bonds from the input file
     */
    void initialize(Uintah::ProblemSpecP& ps);

    /**
     *  Compute the force in a bond given the displacement
     *  The bond material stores the history if needed.
     */
    void computeForce(const Point3D& nodePos,
                      const Point3D& familyPos,
                      const Vector3D& nodeDisp,
                      const Vector3D& familyDisp,
                      const double& horizonSize,
                      Vector3D& force);

    /**
     *  Purpose : Compute micromodulus
     *  Options :
     *  2D constant micromodulus   
     *     dmicroF =  6.0d0*young/(pi*thickness*(horizon**3)*(1.d0/3.d0)*(4.d0/3.d0))  
     *        ==> thickness effect will be vanished in volume integration
     *     dmicroF =  13.5d0*young/(pi*(horizon**3)) ==> thickness = 1
     *  2D canonical micromodulus
     *     dmicroF = 24.0d0*young*(1.0d0-bondlength/horizon)/
     *                   (pi*thickness*(horizon**3)*(1.d0/3.d0)*(4.d0/3.d0))  
     *        ==> thickness effect will be vanished in volume integration
     *     dmicroF = 54.0d0*young*(1.0d0-bondlength/horizon)/(pi*(horizon**3))
     */
    double computeMicroModulus(const double& bondLength, const double& horizonSize);

    /**
     *  Compute the critical strain in a bond 
     */
    double criticalStrain(const double& horizonSize) const;

    // Access methods
    inline void id(const int& ID) {d_id = ID;}
    inline int id() const {return d_id;}
    inline bool hasName() const {return d_have_name;}
    inline std::string name() const {return d_name;}
    void microModulusModel(const MicroModulusModel& microModulus) {d_micro_modulus_model = microModulus;}
    MicroModulusModel microModulusModel() const {return d_micro_modulus_model;}


    inline double density() const {return d_density;}
    inline double youngModulus() const {return d_young_modulus;}
    inline double fractureEnergy() const {return d_fracture_energy;}

    double strain() const {return d_strain;}
    double strainEnergy() const {return d_strain_energy;}
    double microModulus() const {return d_micro_modulus;}

    DamageModel* damageModel() const
    {
      return d_damage_model.get();
    }

  protected:

    int d_id;
    bool d_have_name;
    std::string d_name;
    MicroModulusModel d_micro_modulus_model;

    double d_density;
    double d_young_modulus;
    double d_fracture_energy;

    double d_micro_modulus;
    double d_strain;
    double d_strain_energy;

    DamageModelUP d_damage_model;

  private:
    Material(const Material& mat);
    Material& operator=(const Material& mat);

  }; // end class

} // end namespace

#endif
