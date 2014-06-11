#ifndef __MATITI_WOOD_H
#define __MATITI_WOOD_H

#include <MaterialModels/Material.h>
#include <GeometryPiece/BoxGeometryPiece.h>

#include <Pointers/NodeP.h>
#include <Pointers/WoodSP.h>

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>

#include <iostream>
#include <memory>

#include <Core/ProblemSpec/ProblemSpecP.h>


namespace Matiti {

  class Wood
  {
  public:
      Wood();
      Wood(const Wood& wood);
      Wood(const WoodSP wood); 
      ~Wood();

      void clone(const WoodSP wood);    
       
      void initialize(Uintah::ProblemSpecP& ps);
    
      double computeMicroModulus(const double& bondLength,
                                 const double& horizonSize,
                                 const bool& fiberBond,
                                 const bool& earlywoodNode,
                                 const bool& earlywoodFamily,
                                 const Vector3D& gridSize);

       double computeCriticalStrain(const NodeP node1, const NodeP node2,
                                    const bool& earlywoodNode1,
                                    const bool& earlywoodNode2);

       inline const double& earlywoodFraction() const { return d_earlywood_fraction;}
       inline void earlywoodFraction(const double& fraction) { d_earlywood_fraction = fraction;}
    




  protected:

//****************** earlywood properties ***********************************
      double d_earlywood_young_modulus_radial;
      double d_earlywood_young_modulus_tangential;
      double d_earlywood_young_modulus_longitudinal;

//      double d_earlywood_young_modulus_fiber;
//      double d_earlywood_young_modulus_matrix;

      double d_earlywood_poission_ratio_tangential_radial;
      double d_earlywood_poission_ratio_tangential_longitudinal;
      double d_earlywood_poission_ratio_radial_longitudinal;
      double d_earlywood_poission_ratio_longitudinal_radial;
     
//      double d_earlywood_longitudinal_poission_ratio;
//      double d_earlywood_transverse_poission_ratio;

//********************  latewood properties  *********************************
      double d_latewood_young_modulus_radial;
      double d_latewood_young_modulus_tangential;
      double d_latewood_young_modulus_longitudinal;

//      double d_latewood_young_modulus_fiber;
//      double d_latewood_young_modulus_matrix;

      double d_latewood_poission_ratio_tangential_radial;
      double d_latewood_poission_ratio_tangential_longitudinal;
      double d_latewood_poission_ratio_radial_longitudinal;
      double d_latewood_poission_ratio_longitudinal_radial;
     
//      double d_latewood_longitudinal_poission_ratio;
//      double d_latewood_transverse_poission_ratio;

//*****************************************************************************

      double d_fracture_energy_fiber;  // fracture energy in fiber direction
      double d_fracture_energy_matrix;  // frcture energy in perpendicular to the fiber direction

//      double d_fiber_fraction;
      double d_earlywood_fraction;

      BoxGeometryPiece* d_box;

     


















 };   //end of the class
};   //end of the namespace












#endif
