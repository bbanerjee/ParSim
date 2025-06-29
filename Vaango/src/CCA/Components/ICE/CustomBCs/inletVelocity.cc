/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#include <CCA/Components/ICE/CustomBCs/inletVelocity.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <Core/Math/MersenneTwister.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Variables/CellIterator.h>

using namespace std;

static Uintah::DebugStream cout_BC_CC("ICE_BC_CC", false);
namespace Uintah {
/* ______________________________________________________________________
 Purpose~   -returns (true) if the inletVel BC is specified on any face,
            -reads input parameters needed setBC routines
 ______________________________________________________________________  */
bool read_inletVel_BC_inputs(const ProblemSpecP& prob_spec,
                             MaterialManagerP& materialManager,
                             inletVel_globalVars* global,
                             GridP& grid)
{
  //__________________________________
  // search the BoundaryConditions problem spec
  // determine if inletVelocity bcs are specified
  ProblemSpecP grid_ps= prob_spec->findBlock("Grid");
  ProblemSpecP bc_ps  = grid_ps->findBlock("BoundaryConditions");
 
  bool usingBC = false;
  string whichProfile; 
  
  for( ProblemSpecP face_ps = bc_ps->findBlock( "Face" ); face_ps != nullptr; face_ps=face_ps->findNextBlock( "Face" ) ) {
    map<string,string> face;
    face_ps->getAttributes(face);
    bool setThisFace = false;
    
    for( ProblemSpecP bc_iter = face_ps->findBlock( "BCType" ); bc_iter != nullptr; bc_iter = bc_iter->findNextBlock( "BCType" ) ) {
      map<string,string> bc_type;
      bc_iter->getAttributes(bc_type);
      
      whichProfile = bc_type["var"];
      if ( (whichProfile == "powerLawProfile" || whichProfile == "logWindProfile") && !setThisFace ) {
        usingBC = true;
        setThisFace = true;
        
        //__________________________________
        // bulletproofing
        if (bc_type["id"] == "all"){
          string warn="ERROR:\n Inputs:inletVelocity Boundary Conditions: You've specified the 'id' = all \n The 'id' must be the ice material.";
          throw ProblemSetupException(warn, __FILE__, __LINE__);  
        }
      }
    }
  }
  //__________________________________
  //  read in variables required by the boundary
  //  conditions and put them in the variable basket
  if(usingBC ){
  
    // set default values
    global->vonKarman = 0.4;
   
    ProblemSpecP inlet_ps = bc_ps->findBlock("inletVelocity");
    if (!inlet_ps) {
      string warn="ERROR:\n Inputs:Boundary Conditions: Cannot find inletVelocity_BC block";
      throw ProblemSetupException(warn, __FILE__, __LINE__);
    }
    
    //__________________________________
    //  Find the material associated with this BC
    Material* matl = materialManager->parseAndLookupMaterial(inlet_ps, "material");
    global->iceMatl_indx = matl->getDWIndex();
    
    int numICEMatls = materialManager->getNumMaterials( "ICE" );
    bool foundMatl = false;
    ostringstream indicies;
    
    for(int m = 0; m < numICEMatls; m++){
      ICEMaterial* matl = (ICEMaterial*) materialManager->getMaterial( "ICE",  m );
      int indx = matl->getDWIndex();
      indicies << " " << indx ;
      if(indx == global->iceMatl_indx){
        foundMatl = true;
      }
    }
    
    if(foundMatl==false){
      ostringstream warn;                                                                                              
      warn << "ERROR:\n Inputs: inletVelocity Boundary Conditions: The ice_material_index: "<< global->iceMatl_indx<< " is not "    
           << "\n an ICE material: " << indicies.str() << endl;
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
    
    //__________________________________
    // computational domain
    BBox b;
    grid->getInteriorSpatialRange(b);
    global->gridMin = b.min();
    global->gridMax = b.max();
    
    inlet_ps -> get( "roughness",             global->roughness   );
    inlet_ps -> get( "vonKarmanConstant",     global->vonKarman   );
    inlet_ps -> get( "exponent",              global->exponent    ); 
    inlet_ps -> require( "verticalDirection", global->verticalDir );
    
    Vector tmp = b.max() - b.min();
    double maxHeight = tmp[ global->verticalDir ];   // default value
    
    inlet_ps -> get( "maxHeight",             maxHeight   );
    Vector lo = b.min().asVector();
    global->maxHeight = maxHeight - lo[ global->verticalDir ];

  }
  return usingBC;
}

/* ______________________________________________________________________ 
 Function~  addRequires--   
 ______________________________________________________________________  */
void addRequires_inletVel(Task* t, 
                          const string& where,
                          ICELabel* lb,
                          const MaterialSubset* ice_matls,
                          const bool recursive)
{
  cout_BC_CC<< "Doing addRequires_inletVel: \t\t" <<t->getName()
            << " " << where << endl;
  
  Ghost::GhostType  gn  = Ghost::None;

  //std::cout << " addRequires_inletVel: " << recursive <<  " where: " << where << endl;
  
  if(where == "implicitPressureSolve"){
    t->needs(Task::OldDW, lb->velocity_CCLabel, ice_matls, gn);
  }
  else if(where == "velFC_Exchange"){
    
    //__________________________________
    // define parent data warehouse
    Task::WhichDW pOldDW = Task::OldDW;
    if(recursive) {
      pOldDW  = Task::ParentOldDW;
    }
  
    t->needs(pOldDW, lb->velocity_CCLabel, ice_matls, gn);
  }
}

//______________________________________________________________________ 
void  preprocess_inletVelocity_BCs(DataWarehouse* old_dw,
                                   ICELabel* lb,
                                   const int indx,
                                   const Patch* patch,
                                   const string& where,
                                   bool& set_BCs,
                                   const bool recursive,
                                   inletVel_globalVars* global,
                                   inletVel_localVars* local )
{
  if( indx != global->iceMatl_indx ){
    return;
  }
  
  set_BCs = false;
  
  //std::cout << " preprocess_inletVelocity_BCs: " << where << " recursive: " << recursive << endl;
  
  if(where == "velFC_Exchange"){
    set_BCs = true;
    
    // change the definition of parent(old)DW
    // when implicit
    DataWarehouse* pOldDW = old_dw;
    if(recursive) {
      pOldDW  = old_dw->getOtherDataWarehouse(Task::ParentOldDW); 
    }
    
    pOldDW->get(local->vel_CC, lb->velocity_CCLabel, indx, patch,Ghost::None,0);
  }
  
  if( where == "CC_Exchange" ) {
    set_BCs = true;
  }
  if( where == "Advection" ){
    set_BCs = true;
  }
}
/*_________________________________________________________________
 Purpose~  Set inlet velocity boundary conditions
___________________________________________________________________*/
int  set_inletVelocity_BC(const Patch* patch,
                          const Patch::FaceType face,
                          CCVariable<Vector>& vel_CC,
                          const string& var_desc,
                          Iterator& bound_ptr,
                          const string& bc_kind,
                          const Vector& bc_value,
                          inletVel_globalVars* global,
                          [[maybe_unused]] inletVel_localVars* local )
{
  int nCells = 0;
  
  if (var_desc == "Velocity" && (bc_kind == "powerLawProfile" || bc_kind == "logWindProfile") ) {
    cout_BC_CC << "    Vel_CC (" << bc_kind << ") \t\t" <<patch->getFaceName(face)<< endl;

    // bulletproofing
    if (!global ){
      throw InternalError("set_inletVelocity_BC", __FILE__, __LINE__);
    }
    const Level* level = patch->getLevel();
    
    int nDir = patch->getFaceAxes(face)[0];  //normal velocity direction
    int vDir = global->verticalDir;              // vertical direction
    

    //__________________________________
    // compute the velocity in the normal direction
    // u = U_infinity * pow( h/height )^n
    if( bc_kind == "powerLawProfile" ){
      double d          =  global->gridMin(vDir);
      double gridHeight =  global->gridMax(vDir);
      double height     =  global->maxHeight;
      Vector U_infinity =  bc_value;
      double n          =  global->exponent;
      
      //std::cout << "     height: " << height << " exponent: " << n << " U_infinity: " << U_infinity 
      //     << " nDir: " << nDir << " vDir: " << vDir << endl;
           
      for (bound_ptr.reset(); !bound_ptr.done(); bound_ptr++)   {
        IntVector c = *bound_ptr; 
        
        Point here   = level->getCellPosition(c);
        double h     = here.asVector()[vDir] ;
        
        vel_CC[c]    = U_infinity;             // set the components that are not normal to the face           
        double ratio = (h - d)/height;
        ratio = Clamp(ratio,0.0,1.0);
        
        if( h > d && h < height){
          vel_CC[c][nDir] = U_infinity[nDir] * pow(ratio, n);
        }else{                                // if height < h < gridHeight
          vel_CC[c][nDir] = U_infinity[nDir];
        }
        
        // Clamp edge/corner values 
        if( h < d || h > gridHeight ){
          vel_CC[c] = Vector(0,0,0);
        }
        
        // std::cout << "        " << c <<  " h " << h  << " h/height  " << ratio << " vel_CC: " << vel_CC[c] <<endl;                              
      }
      nCells += bound_ptr.size();
    }
    
    //__________________________________
    //   u = U_star * (1/vonKarman) * ln( (z-d)/roughness)
    else if( bc_kind == "logWindProfile" ){
    
      double inv_K       = 1.0/global->vonKarman;
      double d           = global->gridMin(vDir);  // origin
      double gridMax     = global->gridMax(vDir);
      Vector frictionVel = bc_value;
      double roughness   = global->roughness;
      
//      std::cout << "     d: " << d << " frictionVel: " << frictionVel << " roughness: " << roughness 
//                << " nDir: " << nDir << " vDir: " << vDir << endl;
    
      for (bound_ptr.reset(); !bound_ptr.done(); bound_ptr++)   {
        IntVector c = *bound_ptr;
        
        Point here = level->getCellPosition(c);
        double z   = here.asVector()[vDir] ;
        
        vel_CC[c]    = frictionVel;            // set the components that are not normal to the face
        double ratio = (z - d)/roughness;
        
        vel_CC[c][nDir] = frictionVel[nDir] * inv_K * log(ratio);
        
        // Clamp edge/corner values 
        if(z < d || z > gridMax){
          vel_CC[c] = Vector(0,0,0);
        }
//        std::cout << "        " << c <<  " z " << z  << " z-d " << z-d << " ratio " << ratio << " vel_CC: " << vel_CC[c] <<endl;
      }
      nCells += bound_ptr.size();
    }else{
      ostringstream warn;
      warn << "ERROR ICE::set_inletVelocity_BC  This type of boundary condition has not been implemented ("
           << bc_kind << ")\n" << endl; 
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  }

  return nCells; 
}
 
}  // using namespace Uintah
