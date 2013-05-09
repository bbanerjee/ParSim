#ifndef EMU2DC_MATERIAL_H
#define EMU2DC_MATERIAL_H

#include <DamageModel.h>
#include <DamageModelUP.h>

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

    enum class MicroModulusType {
      Constant=0,
      Conical=1
    };
  
    Material();
    Material(const Material& mat);
    virtual ~Material();

    Material& operator=(const Material& mat);

    void initialize(Uintah::ProblemSpecP& ps);

    double criticalStrain(const double horizonSiaze) const;

    inline void id(const int& ID) {d_id = ID;}
    inline int id() const {return d_id;}
    inline bool hasName() const {return d_have_name;}
    inline std::string name() const {return d_name;}
    void microModulus(const MicroModulusType& microModulus) {d_micro_modulus = microModulus;}
    MicroModulusType microModulus() const {return d_micro_modulus;}


    inline double density() const {return d_density;}
    inline double youngModulus() const {return d_young_modulus;}
    inline double fractureEnergy() const {return d_fracture_energy;}

    DamageModel* damageModel() const
    {
      return d_damage_model.get();
    }

  protected:

    int d_id;
    bool d_have_name;
    std::string d_name;
    MicroModulusType d_micro_modulus;

    double d_density;
    double d_young_modulus;
    double d_fracture_energy;

    DamageModelUP d_damage_model;

  }; // end class

} // end namespace

#endif
