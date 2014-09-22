#ifndef MATITI_MATERIAL_H
#define MATITI_MATERIAL_H

#include <MaterialModels/DamageModelBase.h>
#include <Pointers/DamageModelSP.h>

//#include <MaterialModels/DamageModelSimple.h>
//#include <Pointers/DamageModelSimpleSP.h>

#include <MaterialModels/Density.h>
#include <Pointers/DensitySP.h>

#include <Woods/Wood.h>
#include <Pointers/WoodSP.h>
 
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <memory>
#include <iostream>
#include <vector>

namespace Matiti {

  class Material 
  {
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Material& dam);

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
    void clone(const Material* mat, double randomNum, double coeffOfVar);
    void cloneAverage(const Material* mat1, const Material* mat2);

    /**
     *  Initialize material properties of the bonds from the input file
     */
    void initialize(Uintah::ProblemSpecP& ps);
    void initialize(int id, MicroModulusModel model, 
                    double density, double modulus, double fractureEnergy);

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

    void computeForce(const Point3D& nodePos,
                      const Point3D& familyPos,
                      const Vector3D& nodeDisp,
                      const Vector3D& familyDisp,
                      const double& horizonSize,
                      const DensitySP density,
                      const WoodSP wood,
                      const Vector3D& gridSize,
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
    double computeCriticalStrain(const double& horizonSize);
//    double computeCriticalStrain(const NodeP node1, const NodeP node2) const;
    /**
     * Compute damage factor
     */
    double computeDamageFactor(const double& damage_index) const;

    bool earlywoodPoint(const Point3D& xi, const DensitySP density, const WoodSP wood);

    // Access methods
    inline void id(const int& ID) {d_id = ID;}
    inline int id() const {return d_id;}
    inline bool hasName() const {return d_have_name;}
    inline std::string name() const {return d_name;}
    inline std::string densityType() const {return d_density_type;}
    void microModulusModel(const MicroModulusModel& microModulus) {d_micro_modulus_model = microModulus;}
    MicroModulusModel microModulusModel() const {return d_micro_modulus_model;}


    inline double density() const {return d_density;}
    inline double youngModulus() const {return d_young_modulus;}
    inline double fractureEnergy() const {return d_fracture_energy;}

    inline const double& ringWidth() const { return d_ring; }
    inline void ringWidth(const double& width) { d_ring = width; }

    const DensitySP getDensity() const {return d_node_density;}
    void setDensity(const DensitySP den) {d_node_density = den; } 

    const WoodSP getWood() const { return d_wood; }
    void setWood(const WoodSP wood) {d_wood = wood; }

    double strain() const {return d_strain;}
    double criticalStrain() const {return d_critical_strain;}
    void setCriticalStrain(const double critStr) {d_critical_strain = critStr;}
    double strainEnergy() const {return d_strain_energy;}
    double microModulus() const {return d_micro_modulus;}

    DamageModelSP damageModel() const
    {
      return d_damage_model;
    }

    std::vector<double> densityPolyCoeff() const
    {
      return d_coeffs;
    }

  protected:

    int d_id;
    bool d_have_name;
    std::string d_name;
    std::string d_density_type;
    MicroModulusModel d_micro_modulus_model;

    double d_density;
    double d_young_modulus;
    double d_fracture_energy;
    double d_example;

    double d_micro_modulus;
    double d_strain;
    double d_critical_strain;
    double d_strain_energy;
    double d_ring;
    double d_earlywood_fraction;
 
    DamageModelSP d_damage_model;

    DensitySP d_node_density;

    WoodSP d_wood;


  private:
    Material(const Material& mat);
    Material& operator=(const Material& mat);

    std::vector<double> d_coeffs;

  }; // end class

} // end namespace

#endif
