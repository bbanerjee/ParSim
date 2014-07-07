#include <Woods/Wood.h>
#include <Core/Node.h>
//#include <GeometryPiece/BoxGeometryPiece.h>

//#include <Core/ProblemSpec/ProblemSpec.h>

//#include <iostream>
#include <cmath>

using namespace Matiti;

Wood::Wood() 
     :d_earlywood_young_modulus_radial(0.0)
//    , d_earlywood_young_modulus_tangential(0.0)
    , d_earlywood_young_modulus_longitudinal(0.0)
//    , d_earlywood_poission_ratio_tangential_radial(0.0)
//    , d_earlywood_poission_ratio_tangential_longitudinal(0.0)
    , d_earlywood_poission_ratio_radial_longitudinal(0.0)
    , d_earlywood_poission_ratio_longitudinal_radial(0.0)
//    , d_earlywood_longitudinal_poission_ratio(0.0)
//    , d_earlywood_transverse_poission_ratio(0.0)

    , d_latewood_young_modulus_radial(0.0)
//    , d_latewood_young_modulus_tangential(0.0)
    , d_latewood_young_modulus_longitudinal(0.0)
//    , d_latewood_poission_ratio_tangential_radial(0.0)
//    , d_latewood_poission_ratio_tangential_longitudinal(0.0)
    , d_latewood_poission_ratio_radial_longitudinal(0.0)
    , d_latewood_poission_ratio_longitudinal_radial(0.0)
//    , d_latewood_longitudinal_poission_ratio(0.0)
//    , d_latewood_transverse_poission_ratio(0.0)

    , d_fracture_energy_fiber(0.0), d_fracture_energy_matrix(0.0)
//    , d_young_modulus_fiber(0.0), d_young_modulus_matrix(0.0)
//    , d_fiber_fraction(0.5)
    , d_earlywood_fraction(0.5) {}
     

Wood::Wood(const Wood& wood)
    : d_earlywood_young_modulus_radial(wood.d_earlywood_young_modulus_radial),
//      d_earlywood_young_modulus_tangential(wood.d_earlywood_young_modulus_tangential), 
      d_earlywood_young_modulus_longitudinal(wood.d_earlywood_young_modulus_longitudinal),
//      d_earlywood_poission_ratio_tangential_radial(wood.d_earlywood_poission_ratio_tangential_radial),
//      d_earlywood_poission_ratio_tangential_longitudinal
//                 (wood.d_earlywood_poission_ratio_tangential_longitudinal),
      d_earlywood_poission_ratio_radial_longitudinal
                 (wood.d_earlywood_poission_ratio_radial_longitudinal),
//      d_earlywood_poission_ratio_longitudinal_radial
//                (wood.d_earlywood_poission_ratio_longitudinal_radial),
    
      d_latewood_young_modulus_radial(wood.d_latewood_young_modulus_radial),
//      d_latewood_young_modulus_tangential(wood.d_latewood_young_modulus_tangential), 
      d_latewood_young_modulus_longitudinal(wood.d_latewood_young_modulus_longitudinal),
//      d_latewood_poission_ratio_tangential_radial(wood.d_latewood_poission_ratio_tangential_radial),
//      d_latewood_poission_ratio_tangential_longitudinal
//                 (wood.d_latewood_poission_ratio_tangential_longitudinal),
      d_latewood_poission_ratio_radial_longitudinal
                 (wood.d_latewood_poission_ratio_radial_longitudinal),
//      d_latewood_poission_ratio_longitudinal_radial
//                (wood.d_latewood_poission_ratio_longitudinal_radial),      

      d_fracture_energy_fiber(wood.d_fracture_energy_fiber),
      d_fracture_energy_matrix(wood.d_fracture_energy_matrix),
//      d_fiber_fraction(wood.d_fiber_fraction),
      d_earlywood_fraction(wood.d_earlywood_fraction)
      {
       d_earlywood_poission_ratio_longitudinal_radial =
       d_earlywood_young_modulus_radial*d_earlywood_poission_ratio_radial_longitudinal
                                              /d_earlywood_young_modulus_longitudinal;

       d_latewood_poission_ratio_longitudinal_radial =
       d_latewood_young_modulus_radial*d_latewood_poission_ratio_radial_longitudinal
                                              /d_latewood_young_modulus_longitudinal;


       /* Haplin-Tsai relationships  
       d_young_modulus_fiber = d_fiber_fraction*d_young_modulus_longitudinal+
                               (1-d_fiber_fraction)*d_young_modulus_radial;
       d_young_modulus_matrix = d_young_modulus_longitudinal*d_young_modulus_radial/
                               (d_fiber_fraction*d_young_modulus_radial+
                               (1-d_fiber_fraction)*d_young_modulus_longitudinal);
       d_longitudinal_poission_ratio = d_fiber_fraction*d_poission_ratio_radial_longitudinal+
                                      (1-d_fiber_fraction)*d_poission_ratio_longitudinal_radial;
       d_transverse_poission_ratio =
                     d_young_modulus_matrix*d_longitudinal_poission_ratio/d_young_modulus_fiber;*/
      }




Wood::Wood(const WoodSP wood)
   :  d_earlywood_young_modulus_radial(wood->d_earlywood_young_modulus_radial),
//      d_earlywood_young_modulus_tangential(wood->d_earlywood_young_modulus_tangential), 
      d_earlywood_young_modulus_longitudinal(wood->d_earlywood_young_modulus_longitudinal),
//      d_earlywood_poission_ratio_tangential_radial
//                 (wood->d_earlywood_poission_ratio_tangential_radial),
//      d_earlywood_poission_ratio_tangential_longitudinal
//                 (wood->d_earlywood_poission_ratio_tangential_longitudinal),
      d_earlywood_poission_ratio_radial_longitudinal
                 (wood->d_earlywood_poission_ratio_radial_longitudinal),
//      d_earlywood_poission_ratio_longitudinal_radial
//                (wood->d_earlywood_poission_ratio_longitudinal_radial),
    
      d_latewood_young_modulus_radial(wood->d_latewood_young_modulus_radial),
//      d_latewood_young_modulus_tangential(wood->d_latewood_young_modulus_tangential), 
      d_latewood_young_modulus_longitudinal(wood->d_latewood_young_modulus_longitudinal),
//      d_latewood_poission_ratio_tangential_radial
//                 (wood->d_latewood_poission_ratio_tangential_radial),
//      d_latewood_poission_ratio_tangential_longitudinal
//                 (wood->d_latewood_poission_ratio_tangential_longitudinal),
      d_latewood_poission_ratio_radial_longitudinal
                 (wood->d_latewood_poission_ratio_radial_longitudinal),
//      d_latewood_poission_ratio_longitudinal_radial
//                (wood->d_latewood_poission_ratio_longitudinal_radial),      

      d_fracture_energy_fiber(wood->d_fracture_energy_fiber),
      d_fracture_energy_matrix(wood->d_fracture_energy_matrix),
//      d_fiber_fraction(wood->d_fiber_fraction),
      d_earlywood_fraction(wood->d_earlywood_fraction)
      {
       d_earlywood_poission_ratio_longitudinal_radial =
       d_earlywood_young_modulus_radial*d_earlywood_poission_ratio_radial_longitudinal
                                              /d_earlywood_young_modulus_longitudinal;

       d_latewood_poission_ratio_longitudinal_radial =
       d_latewood_young_modulus_radial*d_latewood_poission_ratio_radial_longitudinal
                                              /d_latewood_young_modulus_longitudinal;


       /* Haplin-Tsai relationships  
       d_young_modulus_fiber = d_fiber_fraction*d_young_modulus_longitudinal+
                               (1-d_fiber_fraction)*d_young_modulus_radial;
       d_young_modulus_matrix = d_young_modulus_longitudinal*d_young_modulus_radial/
                               (d_fiber_fraction*d_young_modulus_radial+
                               (1-d_fiber_fraction)*d_young_modulus_longitudinal);
       d_longitudinal_poission_ratio = d_fiber_fraction*d_poission_ratio_radial_longitudinal+
                                      (1-d_fiber_fraction)*d_poission_ratio_longitudinal_radial;
       d_transverse_poission_ratio =
                     d_young_modulus_matrix*d_longitudinal_poission_ratio/d_young_modulus_fiber;*/
      }
 

Wood::~Wood()
{
}


void
Wood::clone(const WoodSP wood)
   {
      d_earlywood_young_modulus_radial = wood->d_earlywood_young_modulus_radial;
//      d_earlywood_young_modulus_tangential = wood->d_earlywood_young_modulus_tangential; 
      d_earlywood_young_modulus_longitudinal = wood->d_earlywood_young_modulus_longitudinal;
      d_earlywood_poission_ratio_tangential_radial = 
                           wood->d_earlywood_poission_ratio_tangential_radial;
      d_earlywood_poission_ratio_tangential_longitudinal =
                           wood->d_earlywood_poission_ratio_tangential_longitudinal;
      d_earlywood_poission_ratio_radial_longitudinal = 
                           wood->d_earlywood_poission_ratio_radial_longitudinal;
//      d_earlywood_poission_ratio_longitudinal_radial = 
                         //  wood->d_earlywood_poission_ratio_longitudinal_radial;

      d_latewood_young_modulus_radial = wood->d_latewood_young_modulus_radial;
//      d_latewood_young_modulus_tangential = wood->d_latewood_young_modulus_tangential; 
      d_latewood_young_modulus_longitudinal = wood->d_latewood_young_modulus_longitudinal;
      d_latewood_poission_ratio_tangential_radial = 
                           wood->d_latewood_poission_ratio_tangential_radial;
      d_latewood_poission_ratio_tangential_longitudinal =
                           wood->d_latewood_poission_ratio_tangential_longitudinal;
      d_latewood_poission_ratio_radial_longitudinal = 
                           wood->d_latewood_poission_ratio_radial_longitudinal;
//      d_latewood_poission_ratio_longitudinal_radial = 
                         //  wood->d_latewood_poission_ratio_longitudinal_radial;


      d_fracture_energy_fiber = wood->d_fracture_energy_fiber;
      d_fracture_energy_matrix = wood->d_fracture_energy_matrix;
//      d_fiber_fraction = wood->d_fiber_fraction;
      d_earlywood_fraction = wood->d_earlywood_fraction;

      d_earlywood_poission_ratio_longitudinal_radial =
      d_earlywood_young_modulus_radial*d_earlywood_poission_ratio_radial_longitudinal
                                   /d_earlywood_young_modulus_longitudinal;

      d_latewood_poission_ratio_longitudinal_radial =
      d_latewood_young_modulus_radial*d_latewood_poission_ratio_radial_longitudinal
                                   /d_latewood_young_modulus_longitudinal;


       /* Haplin-Tsai relationships  
       d_young_modulus_fiber = d_fiber_fraction*d_young_modulus_longitudinal+
                               (1-d_fiber_fraction)*d_young_modulus_radial;
       d_young_modulus_matrix = d_young_modulus_longitudinal*d_young_modulus_radial/
                               (d_fiber_fraction*d_young_modulus_radial+
                               (1-d_fiber_fraction)*d_young_modulus_longitudinal);
       d_longitudinal_poission_ratio = d_fiber_fraction*d_poission_ratio_radial_longitudinal+
                                      (1-d_fiber_fraction)*d_poission_ratio_longitudinal_radial;
       d_transverse_poission_ratio =
                     d_young_modulus_matrix*d_longitudinal_poission_ratio/d_young_modulus_fiber; */
      }





void 
Wood::initialize(Uintah::ProblemSpecP& ps)
{
//  std::cout << "first of initial wood" << std::endl;
  ps->require("fracture_energy_fiber", d_fracture_energy_fiber);
  ps->require("fracture_energy_matrix", d_fracture_energy_matrix);
//  ps->require("fiber_fraction", d_fiber_fraction);
  ps->require("earlywood_fraction", d_earlywood_fraction);


  Uintah::ProblemSpecP earlywood_ps = ps->findBlock("Earlywood_Properties");
  if (!earlywood_ps) return;
  earlywood_ps->require("young_modulus_radial", d_earlywood_young_modulus_radial);
//  earlywood_ps->require("young_modulus_tangential", d_earlywood_young_modulus_tangential);
  earlywood_ps->require("young_modulus_longitudinal", d_earlywood_young_modulus_longitudinal);
  earlywood_ps->require("poission_ratio_tangential_radial",
                         d_earlywood_poission_ratio_tangential_radial);
  earlywood_ps->require("poission_ratio_tangential_longitudinal",
                         d_earlywood_poission_ratio_tangential_longitudinal);
  earlywood_ps->require("poission_ratio_radial_longitudinal", 
                         d_earlywood_poission_ratio_radial_longitudinal);
//  earlywood_ps->require("poission_ratio_longitudinal_radial",
//                         d_earlywood_poission_ratio_longitudinal_radial);


  Uintah::ProblemSpecP latewood_ps = ps->findBlock("Latewood_Properties");
  if (!latewood_ps) return;
  latewood_ps->require("young_modulus_radial", d_latewood_young_modulus_radial);
//  latewood_ps->require("young_modulus_tangential", d_latewood_young_modulus_tangential);
  latewood_ps->require("young_modulus_longitudinal", d_latewood_young_modulus_longitudinal);
  latewood_ps->require("poission_ratio_tangential_radial",
                         d_latewood_poission_ratio_tangential_radial);
  latewood_ps->require("poission_ratio_tangential_longitudinal",
                         d_latewood_poission_ratio_tangential_longitudinal);
  latewood_ps->require("poission_ratio_radial_longitudinal", 
                         d_latewood_poission_ratio_radial_longitudinal);
//  latewood_ps->require("poission_ratio_longitudinal_radial",
//                         d_latewood_poission_ratio_longitudinal_radial);

  
  d_earlywood_poission_ratio_longitudinal_radial =
  d_earlywood_young_modulus_radial*d_earlywood_poission_ratio_radial_longitudinal
                               /d_earlywood_young_modulus_longitudinal;

  d_latewood_poission_ratio_longitudinal_radial =
  d_latewood_young_modulus_radial*d_latewood_poission_ratio_radial_longitudinal
                               /d_latewood_young_modulus_longitudinal;


/*   Haplin-Tsai relationships  
  d_young_modulus_fiber = d_fiber_fraction*d_young_modulus_longitudinal+
                          (1-d_fiber_fraction)*d_young_modulus_radial;
  d_young_modulus_matrix = d_young_modulus_longitudinal*d_young_modulus_radial/
                          (d_fiber_fraction*d_young_modulus_radial+
                          (1-d_fiber_fraction)*d_young_modulus_longitudinal);
  d_longitudinal_poission_ratio = d_fiber_fraction*d_poission_ratio_radial_longitudinal+
                          (1-d_fiber_fraction)*d_poission_ratio_longitudinal_radial;
  d_transverse_poission_ratio =
           d_young_modulus_matrix*d_longitudinal_poission_ratio/d_young_modulus_fiber; */
// std::cout << "end of initial wood" << std::endl; 
}




double
Wood::computeMicroModulus(const double& bondLength,
                    const double& horizonSize,
                    const bool& fiberBond,
                    const bool& earlywoodNode,
                    const bool& earlywoodFamily,
                    const Vector3D& gridSize)
{

    double grid_size = gridSize.y();
    double m = horizonSize/grid_size;

//    std::cout << "grid size= " << "(" << gridSize.x() << ", " << gridSize.y() << ", " << gridSize.z() << ")" << std::endl;
//    std::cout << "grid_size= " << grid_size << " m= " << m << std::endl;

//  double m = 5;
  double fiber_coeff = M_PI*m/2;
  double matrix_coeff = M_PI*m/(M_PI*m-2);
 

  double conical_coeff = 6*(1-bondLength/horizonSize)/(horizonSize*horizonSize);


  double earlywood_denom = (1- d_earlywood_poission_ratio_radial_longitudinal*
                                   d_earlywood_poission_ratio_longitudinal_radial);

//*********************** compute earlywood fiber micromodulus *******************************
  double earlywood_c_fiber = (d_earlywood_young_modulus_longitudinal+
                 d_earlywood_poission_ratio_radial_longitudinal*d_earlywood_young_modulus_radial)*
                 conical_coeff/earlywood_denom;
         earlywood_c_fiber *= fiber_coeff;
//********************************************************************************************


  double latewood_denom = (1- d_latewood_poission_ratio_radial_longitudinal*
                                   d_latewood_poission_ratio_longitudinal_radial);


//*********************** compute latewood fiber micromodulus ********************************
  double latewood_c_fiber = (d_latewood_young_modulus_longitudinal+
                 d_latewood_poission_ratio_radial_longitudinal*d_latewood_young_modulus_radial)*
                 conical_coeff/latewood_denom;
         latewood_c_fiber *= fiber_coeff;
//********************************************************************************************

        
   
  conical_coeff *= 2/(M_PI*horizonSize);

//*********************** compute earlywood matrix micromodulus ******************************* 
  double earlywood_c_matrix = (d_earlywood_young_modulus_radial+
                 d_earlywood_poission_ratio_radial_longitudinal*d_earlywood_young_modulus_radial)*
                 conical_coeff/earlywood_denom;
         earlywood_c_matrix *= matrix_coeff;
//*********************************************************************************************


//*********************** compute latewood matrix micromodulus ********************************

  double latewood_c_matrix = (d_latewood_young_modulus_radial+
                 d_latewood_poission_ratio_radial_longitudinal*d_latewood_young_modulus_radial)*
                 conical_coeff/latewood_denom;
         latewood_c_fiber *= matrix_coeff;
//*********************************************************************************************



 if (fiberBond)
  {
    if (earlywoodNode) return earlywood_c_fiber;
    else return latewood_c_fiber;
  }
 else 
  {
    if ((earlywoodNode) && (earlywoodFamily)) return earlywood_c_matrix;
    else if (!(earlywoodNode) && !(earlywoodFamily)) return latewood_c_matrix;
    else return std::min(earlywood_c_matrix, earlywood_c_fiber);
  }

 
       
}


double
Wood::computeCriticalStrain(const NodeP node1, const NodeP node2,
                      const bool& earlywoodNode1,
                      const bool& earlywoodNode2)
{ 
  Vector3D xi = node2->position() - node1->position();
  double horizonSize = node1->horizonSize();
  double initial_length = xi.length();
  double conical_coeff = 12*(1-initial_length/horizonSize)/(M_PI*horizonSize*horizonSize*horizonSize);
  double earlywood_denom = (1- d_earlywood_poission_ratio_radial_longitudinal*
                                   d_earlywood_poission_ratio_longitudinal_radial);
  double latewood_denom = (1- d_latewood_poission_ratio_radial_longitudinal*
                                   d_latewood_poission_ratio_longitudinal_radial);
  double earlywood_c_fiber_iso = (d_earlywood_young_modulus_longitudinal+
                 d_earlywood_poission_ratio_radial_longitudinal*d_earlywood_young_modulus_radial)*
                 conical_coeff/earlywood_denom;
  double earlywood_c_matrix_iso = (d_earlywood_young_modulus_radial+
                 d_earlywood_poission_ratio_radial_longitudinal*d_earlywood_young_modulus_radial)*
                 conical_coeff/earlywood_denom;
  double latewood_c_fiber_iso = (d_latewood_young_modulus_longitudinal+
                 d_latewood_poission_ratio_radial_longitudinal*d_latewood_young_modulus_radial)*
                 conical_coeff/latewood_denom;
  double latewood_c_matrix_iso = (d_latewood_young_modulus_radial+
                 d_latewood_poission_ratio_radial_longitudinal*d_latewood_young_modulus_radial)*
                 conical_coeff/latewood_denom;


//************************** calculate critical relative elongations *************************

  double h = horizonSize;
  double fiber_coefficient = 20*d_fracture_energy_fiber/(h*h*h*h);
  double matrix_coefficient = 20*d_fracture_energy_matrix/(h*h*h*h);
  

  double earlywood_fiber_critical_elongation = std::sqrt(fiber_coefficient/earlywood_c_fiber_iso);
  double earlywood_matrix_critical_elongation = std::sqrt(matrix_coefficient/earlywood_c_matrix_iso);
  double latewood_fiber_critical_elongation = std::sqrt(fiber_coefficient/latewood_c_fiber_iso);
  double lateywood_matrix_critical_elongation = std::sqrt(matrix_coefficient/latewood_c_fiber_iso);

  bool fiber_bond = ((xi.y() == 0) && (xi.z() == 0));

  if (fiber_bond)
  {
    if (earlywoodNode1) return earlywood_fiber_critical_elongation;
    else return latewood_fiber_critical_elongation;
  }
 else 
  {
    if ((earlywoodNode1) && (earlywoodNode2)) return earlywood_matrix_critical_elongation;
    else if (!(earlywoodNode1) && !(earlywoodNode2)) return lateywood_matrix_critical_elongation;
    else return std::max
               (earlywood_fiber_critical_elongation, earlywood_matrix_critical_elongation);
  }

}

  
  










   
