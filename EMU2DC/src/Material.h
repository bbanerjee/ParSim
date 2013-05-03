#ifndef EMU2DC_MATERIAL_H
#define EMU2DC_MATERIAL_H

#include <DamageModel.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <memory>
#include <iostream>

namespace Emu2DC {

  typedef std::unique_ptr<DamageModel> DamageModelUP;

  class Material 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Material& dam);

  public:
  
    Material();
    virtual ~Material();

    void initialize(Uintah::ProblemSpecP& ps);

    inline void id(const int& ID) {d_id = ID;}
    inline int id() const {return d_id;}
    inline bool hasName() const {return d_have_name;}
    inline std::string name() const {return d_name;}
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

    double d_density;
    double d_young_modulus;
    double d_fracture_energy;

    DamageModelUP d_damage_model;

  }; // end class

} // end namespace

#endif
