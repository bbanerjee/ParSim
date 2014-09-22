#include <Peridynamics.h>
#include <Core/FamilyComputer.h>
#include <MaterialModels/Material.h>
#include <Woods/Wood.h>
#include <Pointers/WoodSP.h>
#include <Core/Body.h>
#include <Core/Node.h>
#include <Core/Bond.h>
#include <Core/Exception.h>

#include <Pointers/MaterialUP.h>
#include <Pointers/MaterialSP.h>

#include <Geometry/Vector3D.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>
#include <chrono>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>

#include <errno.h>
#include <fenv.h>

using namespace Matiti;

Peridynamics::Peridynamics() 
{
//myfile.open("checkMatiti.txt", ios::out);
}

Peridynamics::~Peridynamics() 
{
}

 // ofstream myfile;


 void
 Peridynamics::problemSetup(Uintah::ProblemSpecP& ps)
 {
// myfile.open("checkMatiti.txt", ios::out);
 //if (myfile.is_open()) {
     //myfile.precision(20);
  // Set up the simulation state data
  d_state.initialize(ps);
  // std::cout << d_state ;

  // Set up the time information
  d_time.initialize(ps);
  // std::cout << d_time ;

  // Set up the output information
  d_output.initialize(ps);
  // std::cout << d_output ;

  // Set up the domain
  d_domain.initialize(ps);
  // std::cout << d_domain ;

  // Set up the initial material list
  int count = 0;
  for (Uintah::ProblemSpecP mat_ps = ps->findBlock("Material"); mat_ps != 0;
       mat_ps = mat_ps->findNextBlock("Material")) {
    MaterialSP mat = std::make_shared<Material>();
    mat->initialize(mat_ps);
    d_wood_bool = (mat->hasName() == true) && (mat->name() == "wood");
    if (d_wood_bool)  {
        WoodSP wood = std::make_shared<Wood>();
        wood->initialize(mat_ps);
        d_wood_list.emplace_back(wood);
    } 
    mat->id(count);
    d_mat_list.emplace_back(mat);
    ++count;
    //std::cout << count << std::endl;
   // std::cout << *mat ;
  }


/*   for (auto iter = d_mat_list.begin(); iter != d_mat_list.end(); ++iter) {
        MaterialSP mat = *iter;
        if ((mat->name() == "wood")) */

  // Set up the body information
  count = 0;
  for (Uintah::ProblemSpecP body_ps = ps->findBlock("Body"); body_ps != 0;
       body_ps = body_ps->findNextBlock("Body")) {
//    std::cout << count << endl;
    // Initialize the body (nodes, elements, cracks)
    BodySP body = std::make_shared<Body>();
    body->initialize(body_ps, d_domain, d_state, d_mat_list); 
    body->id(count);
    d_body_list.emplace_back(body);
    
    ++count;
//     std::cout << *body;
 
    // Compute the horizons of nodes in the body
    HorizonComputer compute_horizon;
    compute_horizon(body, d_state);
  }

  // Create the initial family of each node and remove any bonds
  // intersected by initial cracks
  //feclearexcept(FE_ALL_EXCEPT);
  for (auto iter = d_body_list.begin(); iter != d_body_list.end(); ++iter) {

    // Create family
    (*iter)->createInitialFamily(d_domain);

    const NodePArray& node_list = (*iter)->nodes();

//    std::cout << *(*iter);
//    myfile << "Before Removing the bonds intersected by cracks: ";
/*    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
        NodeP cur_node = *node_iter;
        if (cur_node->getID() < 13) {
//           myfile << "\n\nNode(" << cur_node->getID() << ")= [" << cur_node->x() << ", " << cur_node->y() 
//                  << ", " << cur_node->z() << "]\n";
//           myfile << "Bonds: ";
           BondPArray& bonds = cur_node->getBonds();
           for (auto family_iter = bonds.begin(); family_iter != bonds.end(); family_iter++) {
               BondP bond = *family_iter;                           
//               myfile << "[" << bond->first()->getID() << ", " << bond->second()->getID() << "],   ";
           }
         }
    } */
            

    // Remove bonds
    (*iter)->removeBondsIntersectedByCracks();

    // After all bond deletions have been completed update nodal damage indices
    (*iter)->updateDamageIndex();

    // Print body information
//    std::cout << "After Removing the bonds: " << std::endl;
//    myfile << "\n\nAfter Removing the bonds: ";
  }

//  std::cout << "After Removing the bonds: " << std::endl;
/*   for (auto iter = d_body_list.begin(); iter != d_body_list.end(); ++iter) { 
     const NodePArray& node_list = (*iter)->nodes();
     for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
        NodeP cur_node = *node_iter;
        if (cur_node->getID() < 13) {
           myfile << "\n\nNode(" << cur_node->getID() << ")= [" << cur_node->x() << ", " << cur_node->y() 
                  << ", " << cur_node->z() << "]\n";
           myfile << "Bonds: ";
           BondPArray& bonds = cur_node->getBonds();
           for (auto family_iter = bonds.begin(); family_iter != bonds.end(); family_iter++) {
               BondP bond = *family_iter;                           
               myfile << "[" << bond->first()->getID() << ", " << bond->second()->getID() << "],   ";
           }
         }
     }
   } */

  
  

  //nodeNominate = ((cur_node->getID() == minNodeId) || (cur_node->getID() == maxNodeId));
     
  //int fe = fetestexcept(FE_ALL_EXCEPT);
  //std::cout << "Floating point exception somewhere in the code ?" << fe
  //          << std::endl;
  //std::cout << "   The following exceptions are set" << std::endl;
  //if (fe & FE_DIVBYZERO) puts ("FE_DIVBYZERO");
  //if (fe & FE_INEXACT)   puts ("FE_INEXACT");
  //if (fe & FE_INVALID)   puts ("FE_INVALID");
  //if (fe & FE_OVERFLOW)  puts ("FE_OVERFLOW");
  //if (fe & FE_UNDERFLOW) puts ("FE_UNDERFLOW");
 

  d_num_broken_bonds = 0;

// }else {
//   cout << "Unable to open the file";  
// } 

}
//  

void 
Peridynamics::problemSetup(Time& time,
                           OutputVTK& output,
                           SimulationState& state,
                           Domain& domain,
                           MaterialSPArray& matList,
                           BodySPArray& bodyList)
{
  d_time.clone(time);
  d_output.clone(output);
  d_state.clone(state);
  d_domain.clone(domain);
  d_mat_list = matList;
  d_body_list = bodyList;
  d_num_broken_bonds = 0;
}



void
Peridynamics::run()
{
  std::cerr << "....Begin Solver...." << std::endl;

  // Constants
  Vector3D Zero(0.0, 0.0, 0.0);

  // Write the output at the beginning of the simulation
  d_output.write(d_time, d_domain, d_body_list);

  // Displacement driven computation. Body alreay has initial velocity.  
  // Compute initial displacement
  applyInitialConditions();

  // Do time incrementation
  double cur_time = 0.0;
  int cur_iter = 1;
  bool nodeNominates = false;
  int minNodeId = 100000;
  int maxNodeId = -1;

  for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {    
         const NodePArray& node_list = (*body_iter)->nodes(); 
         for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
              NodeP cur_node = *node_iter;
              nodeNominates = (cur_node->initialFamilySize() != cur_node->currentFamilySize());
              if (nodeNominates) 
              {
               if (minNodeId > cur_node->getID())
                   minNodeId = cur_node->getID();
               if (maxNodeId < cur_node->getID())
                   maxNodeId = cur_node->getID();
              }
          }         
  }


  /*   
  ofstream myfile;
  myfile.open("nodalInformation.txt", ios::out);
  if (myfile.is_open()) {
    myfile << "Number of the node    " << "Volume of the node    " << "Mass of the node    " << 
          "Family size of the node \n";
  } else {
    cout << "Unable to open the file";  
  }

  // Find the total initial momentom
  Vector3D initialTotalMomentum(0.0, 0.0, 0.0);
  for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {
    const NodePArray& node_list = (*body_iter)->nodes();
    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

      // Get the node
      NodeP cur_node = *node_iter;
        
      initialTotalMomentum +=cur_node->velocity()*(cur_node->volume()*cur_node->densityNode());
   
      if (myfile.is_open()) {
        myfile << "    " << cur_node->getID() << "        " << cur_node->volume() << "        " << 
        cur_node->volume()*cur_node->densityNode() << "        " << cur_node->numAdjacentElements() << "\n";
      } else {
        cout << "Unable to open the file";
      }
    }
  }
  */
  
  /*
  fstream myfile;
  myfile.open("finalMomentum.txt", ios::out);
  if (myfile.is_open()) {
    myfile <<"coeff_of_restitution= " << d_velocitybc.getRestitution() << "   horizon_factor= "
           << d_state.horizonFactor() << "\n" ;
  } else {
    cout << "Unable to open the file"; 
  }
  */
//  myfile.open("checkMatiti.txt", ios::out);
  myfile.open("checkMatiti.txt", ios::out);
  if (myfile.is_open()) {
     myfile.precision(20);

     int nodeNumber = 740;
    
     while (cur_time < d_time.maxTime() && cur_iter < d_time.maxIter()) {
   
       auto t1 = std::chrono::high_resolution_clock::now();

       // Get the current delT
       double delT = d_time.delT();
       //bool condition = ((cur_iter == 1) || (cur_iter == 5) || (cur_iter == 500 ));
       //bool condition = ((cur_iter == 1) || (cur_iter%20 == 0));
       
/*       if (condition) 
          myfile << "\n Iteration Number= " << cur_iter << "\n";*/
   
       //std::cout << std::endl << " Iteration Number " << cur_iter << std::endl << std::endl;

       // Do the computations separately for each body
       for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {    
          const NodePArray& node_list = (*body_iter)->nodes();
          // Get acceleration due to body force
          Vector3D body_force = (*body_iter)->bodyForce();
          // Get the grid's size
          Vector3D grid_size = (*body_iter)->gridSize();
  
//          std::cout << " STAGE ONE VV" << std::endl;
          // Step 1:
          // Update nodal velocity and nodal displacement
          // 1. v(n+1/2) = v(n) + dt/2m * f(u(n))
          // 2. u(n+1) = u(n) + dt * v(n+1/2)

          // Compute internal force (Step 1)
          //    and add body force to internal force
          for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

              // Get the node
              NodeP cur_node = *node_iter;
               // Compute the internal force at the node
               if (cur_iter == 1) { 
                  Vector3D internal_force_density(0.0, 0.0, 0.0);
                  computeInternalForceDensity(cur_node, internal_force_density, grid_size);
                  cur_node->internalForce(internal_force_density);
               }
           }

           for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

               // Get the node
               NodeP cur_node = *node_iter;
                if (cur_node->omit()) {
//                  cur_node->velocity(Zero);
//                  cur_node->displacement(Zero);
//                  cur_node->newVelocity(Zero);
//                  cur_node->newDisplacement(Zero);
                  cur_node->newVelocity(Zero);
                  cur_node->newDisplacement(cur_node->displacement());
                  continue;  // skip this node
                } else {

                Vector3D internal_force_density = cur_node->internalForce();

               // Compute body force density
                Vector3D body_force_density = body_force*cur_node->densityNode();

               // Get external force
                Vector3D external_force_density = cur_node->externalForce();
               
               // Compute acceleration (F_ext - F_int = m a)
               // **TODO** Make sure mass is conserved
                Vector3D acceleration = (external_force_density + internal_force_density + body_force_density)/cur_node->densityNode();
        
//              int id = cur_node->getID();
               //if (id == 1 || id == 2 || id == 104 || id == 105) {
//                  std::cout << "Node = " << id << " F_ext = " << external_force_density << " F_int = " 
//                  << internal_force_density << std::endl;
//                  std::cout << " body_force=  " << body_force_density << " density = " << cur_node->densityNode() 
//                  << " volume = " << cur_node->volume() 
//                  << " area = " << cur_node->area() 
//                  << " acceleration= " << acceleration << std::endl;
//                  std::cout << *cur_node ;
            // }
        

            // Integrate acceleration with velocity Verlet algorithm
            // and Update nodal velocity
            // 1. v(n+1/2) = v(n) + dt/2m * f(u(n))
                Vector3D velocity_mid(0.0, 0.0, 0.0);
                integrateNodalAcceleration(cur_node->velocity(), acceleration, 0.5*delT, velocity_mid);
                cur_node->newVelocity(velocity_mid);
                cur_node->velocity(cur_node->newVelocity());

            // Integrate the mid step velocity
            // and Update nodal displacement
            // 2. u(n+1) = u(n) + dt * v(n+1/2)
                Vector3D displacement_new(0.0, 0.0, 0.0);
                integrateNodalVelocity(cur_node->displacement(), cur_node->velocity(), delT, displacement_new);
//              cur_node->displacement(displacement_new);
                cur_node->newDisplacement(displacement_new);

//              std::cout << "horizon_size= " << cur_node->horizonSize() << std::endl;

//              std::cout << "After Stage 1: Node = " << cur_node->getID() 
//                        << " vel_old = " << cur_node->velocity() 
//                        << " vel_new = " << cur_node->newVelocity() 
//                        << " disp_old = " << cur_node->displacement() 
//                        << " disp_new = " << cur_node->newDisplacement() 
//                        << std::endl;

              } //end of else

            } //end of for (node)

           // **TODO** Apply any external forces due to contact  
           // (*body_iter)->applyContactForces();

       } // end of for (body)
    
    
       //int fe = fetestexcept(FE_ALL_EXCEPT);
       //std::cout << "Floating point exception somewhere in the code ?" << fe
       //          << std::endl;
       //std::cout << "   The following exceptions are set" << std::endl;
       //if (fe & FE_DIVBYZERO) puts ("FE_DIVBYZERO");
       //if (fe & FE_INEXACT)   puts ("FE_INEXACT");
       //if (fe & FE_INVALID)   puts ("FE_INVALID");
       //if (fe & FE_OVERFLOW)  puts ("FE_OVERFLOW");
       //if (fe & FE_UNDERFLOW) puts ("FE_UNDERFLOW");

//       std::cout << " STAGE 2 VV " << std::endl;
       // Stage 2 of Velocity Verlet
       //   Update nodal velocity
       //     3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
       
       for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {

           // Get the nodes in the body
           const NodePArray& node_list = (*body_iter)->nodes();

           // Get body force per unit volume
           Vector3D body_force = (*body_iter)->bodyForce();

           // Get the grid's size
           Vector3D grid_size = (*body_iter)->gridSize();
   //      std::cout << "grid size= " << "(" << grid_size.x() << ", " << grid_size.y() << ", " << grid_size.z() << ")" << std::endl;

        // Compute internal force  (Step 2)
           for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

               // Get the node
               NodeP cur_node = *node_iter;
               cur_node->displacement(cur_node->newDisplacement());
           }

           for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

               // Get the node
               NodeP cur_node = *node_iter;

     //        std::cout << "Before stage 2 : Node = " << cur_node->getID() 
//             << " vel = " << cur_node->velocity() << " disp = " << cur_node->displacement() << std::endl;

//             Compute updated internal force from updated nodal displacements
               std::vector <int> posBrokenBond;
               Vector3D internal_force_density1(0.0, 0.0, 0.0);
               computeInternalForceDensity(cur_node, internal_force_density1, grid_size, cur_iter, posBrokenBond);
               cur_node->internalForce(internal_force_density1);
        
            // Compute body force density
               Vector3D body_force_density = body_force*cur_node->densityNode();

            // Get external force
               Vector3D external_force_density = cur_node->externalForce();

               
               

            // Compute acceleration (F_ext - F_int = m a) 
               Vector3D acceleration = (external_force_density + internal_force_density1 + body_force_density)/cur_node->densityNode();

            // Integrate acceleration with velocity Verlet algorithm
            // and Update nodal velocity
            // 3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
               Vector3D velocity_new(0.0, 0.0, 0.0);
               integrateNodalAcceleration(cur_node->velocity(), acceleration, 0.5*delT, velocity_new);
               cur_node->newVelocity(velocity_new);
         
              // bool condition = (cur_node->getID() == 1) && ((cur_iter == 1) || (cur_iter == 5) || (cur_iter == 500 ));    
/*               if ((cur_node->getID() == nodeNumber) && (condition)) {
                 myfile << "\nNode= " << cur_node->getID() << " , (" << cur_node->x() << ", " << cur_node->y()
                                      << ", " << cur_node->z() << ") , volume= " << cur_node->volume()
                                      << ", horizon_size= " << cur_node->horizonSize() << "\n";

                 myfile << "external_force= " << "  (" << cur_node->externalForce().x() << ", " << cur_node->externalForce().y()
                                              << ", " << cur_node->externalForce().z() << ") \n";
                 myfile << "internal_force= " << "  (" << cur_node->internalForce().x() << ", " << cur_node->internalForce().y()
                                              << ", " << cur_node->internalForce().z() << ") \n";
                 myfile << "displacement= " << "  (" << cur_node->displacement().x() << ", " << cur_node->displacement().y()
                                            << ", " << cur_node->displacement().z() << ") \n";               
                 myfile << "velocity= " << "  (" << cur_node->newVelocity().x() << ", " << cur_node->newVelocity().y()
                                        << ", " << cur_node->newVelocity().z() << ") \n";
               }  */
              
//               if (cur_node->getID() == nodeNumber)
//               std::cout << "Iteration: " << cur_iter << std::endl;
//               nodeNominates = ((cur_node->getID() == minNodeId) || (cur_node->getID() == maxNodeId));
//               if (nodeNominates) 
//               {
//               myfile << "Iteration: " << cur_iter << std::endl;

//               std::cout << "Node " << cur_node->getID() << std::endl;
//               myfile << "Node " << cur_node->getID() << std::endl;
//               }

               BondPArray& bonds = cur_node->getBonds();
               for (auto family_iter = bonds.begin(); family_iter != bonds.end(); family_iter++) {
                    BondP bond = *family_iter;
                    const NodeP node1 = bond->first();
                    const NodeP node2 = bond->second();
//                    const Material* mat = bond->material();
                    const MaterialSP mat = bond->materialS();
//                    if((cur_node->getID() == nodeNumber) && (condition))
//                    if (cur_node->getID() == nodeNumber)
//                    if (nodeNominates) {
              
                    const double dmgij = std::max(node1->damageIndex(), node2->damageIndex());
                    double damage_fac = mat->computeDamageFactor(dmgij);
                    const double horizon = node1->horizonSize();
                    double critical_strain_cur;
                    critical_strain_cur = mat->computeCriticalStrain(horizon);
                    double ecr2 = critical_strain_cur*damage_fac;
                    double str = mat->strain();
//                    myfile << "str= " << str << "\n";
/*                    if ((cur_node->getID() == nodeNumber) && (condition)) {
                       if ((str > ecr2) && (node1->failureAllowed() == true) && (node2->failureAllowed() == true)) {
                           myfile << ",****Breaking bond between nodes " << node1->getID() << " and " << node2->getID() 
                                  << ", critical strain = " << critical_strain_cur << "\n, ecr2 = "
                                  << ecr2 << ", strain = " << str << std::endl;
                       }
                    } */
               }

//               if ((cur_iter > 490) && (cur_node->getID() == 1513)) {
//                   myfile << "Iteration= " << cur_iter << std::endl;
//                   myfile << *cur_node << std::endl;
//                   myfile << "posBrokenBond= ";
                  int modifyPos = 0;                 
                  for (auto pos_iter = posBrokenBond.begin(); pos_iter != posBrokenBond.end(); pos_iter++) {
                      int position = *pos_iter;
//                     if (cur_node->getID() == nodeNumber)
//                      if (nodeNominates)
                        if ((cur_iter == 600 ) && (cur_node->getID() == 500))
                           std::cout << "posBrokenBond= " << position << ", ";
                      myfile << position << ", ";  
                      position-=modifyPos;
                      bonds.erase(bonds.begin()+position);
                      modifyPos++;
                  }
//               }
//               std::cout << std::endl;
               if ((cur_iter > 490) && (cur_node->getID() == 1513))
               myfile << std::endl;               

               if (cur_node->omit()) continue;  // skip this node
               cur_node->updateDamageIndex();
  
//               if ((cur_node->getID() == nodeNumber) && (condition)) {
                 if ((cur_iter > 490) && (cur_node->getID() == 1513))
                      myfile << "+++++Damage_index= " << cur_node->damageIndex() << "\n\n";
//               }
               

//             std::cout << "After stage 2 : Node = " << cur_node->getID() << " internal force= " << cur_node->internalForce() 
//             << " vel = " << cur_node->newVelocity() << " disp = " << cur_node->newDisplacement() << std::endl;
           } // end of for (node)

        // Apply displacement (and actually velocity too) boundary conditions
        // before the first stage of time integration
       //  Adds a reaction force to the external force in the direction opposite the
       //  internal force and with the same magnitude
           (*body_iter)->applyDisplacementBC();

        // Print body information
//         std::cout << *(*body_iter);
       } // end of for (body)
 
       // Apply domain boundary conditions to the body
       for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {
           d_domain.applyVelocityBC(*body_iter);
       }

       // Update kinematic quantities and delT
       //feclearexcept(FE_ALL_EXCEPT);
       Vector3D finalTotalMomentum(0.0, 0.0, 0.0);
       double kineticEnergy(0.0);
       double delT_new = delT;
       for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {

         // Loop through nodes in the body
            const NodePArray& node_list = (*body_iter)->nodes();
            for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {

             // Get the node
                NodeP cur_node = *node_iter;

             // Update kinematic quantities
                cur_node->displacement(cur_node->newDisplacement());
                cur_node->velocity(cur_node->newVelocity());

             // Compute momentum and KE
                finalTotalMomentum +=cur_node->velocity()*(cur_node->volume()*cur_node->densityNode());
                kineticEnergy +=cur_node->velocity().lengthSq()*1/2*cur_node->volume()
                *cur_node->densityNode();

             // Compute stable delT
                delT_new = std::min(cur_node->computeStableTimestep(d_time.timeStepFactor()), delT_new);
             }
        }

//       if (cur_iter > 1) {
//          for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {
//          myfile << "\n Body= " << *(*body_iter) << "\n";
//          }
//       } 

    //  Break bonds that are stretched (only after first iteration is complete)
//      if (cur_iter > 1) {
//          for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {

//               const NodePArray& node_list = (*body_iter)->nodes();
//               breakBonds(node_list);
//          }
//      }

    //  Write the finalTotalMomentum in an output file named finalMomentum.txt
    /*
        if (myfile.is_open()) {
            myfile << "   " << cur_iter<< "    " << finalTotalMomentum[0]
            << "   " << finalTotalMomentum[1] << "   " 
            << finalTotalMomentum[2] << "    " << finalTotalMomentum.length() << "    " << kineticEnergy << "\n";
        } else {
            std::cout << "Unable to open the file"; 
        }
    */
    
      //check the momentum conservation

     // const double err=0.0001;
     // const vector3D errVec=(err, err, err);
     // Vector3D diffFinalMomentum=finalTotalMomentum-initialTotalMomentum;
     // if (abs(diffFinalMomentum[0])>err || abs(diffFinalMomentum[1])>err || abs(diffFinalMomentum[2])>err) {
    //      std::ostringstream out; 
    //      out << "The momentum conservation is not true. Time step=" << cur_iter 
    //      << "Initial Momentum=" << initialTotalMomentum << "Final   Momentum=" << finalTotalMomentum;
    //      throw Exception(out.str(), __FILE__, __LINE__); 
    //  }   


    //  Get memory usage
        double res_mem = 0.0, shar_mem = 0.0;
        checkMemoryUsage(res_mem, shar_mem);

    //  Print out current time and iteration
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "Time = " << cur_time << " Iteration = " << cur_iter 
                  << " Compute time = " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
                  << " ms"
                  << " Memory = " << (res_mem-shar_mem)/1024 << " MB" << std::endl;

     // Update the current time and current iteration
        cur_time = d_time.incrementTime(delT);
        d_time.setDelT(delT_new);
        ++cur_iter;
        
//        if(condition)
//        std::cout << std::endl << "ITERATION= " << cur_iter << std::endl << std::endl;

     // Output nodal information every snapshots_frequency iteration   
        int output_freq = d_output.outputIteratonInterval();
        if (cur_iter%output_freq == 0) {
           d_output.write(d_time, d_domain, d_body_list);
           //std::cout << "Wrote out data at time " << d_time << std::endl;
        }
   }

 }else {
   cout << "Unable to open the file";  
 }
 myfile.close();
}
 
  
  
  
//} //end of while

//Is the material a wood?
//*TODO it should be completed 

bool
Peridynamics::isMaterialWood()
{
 MaterialSP mat = d_mat_list.front();
 return ((mat->hasName() == true) && (mat->name() == "wood"));
}

// Apply the initial velocity 
void
Peridynamics::applyInitialConditions()
{
  // Compute initial displacement
  for (auto body_iter = d_body_list.begin(); body_iter != d_body_list.end(); ++body_iter) {

    // Get the initial velocity
    const Vector3D& init_vel = (*body_iter)->initialVelocity();
    
    // Loop through nodes in the body
    const NodePArray& node_list = (*body_iter)->nodes();
    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
      if ((*node_iter)->omit()) continue;
      (*node_iter)->velocity(init_vel);
    }
  }
}

// Compute the internal force per unit volume
void
Peridynamics::computeInternalForceDensity(const NodeP cur_node,
                                          Vector3D& internalForce, const Vector3D& gridSize)
{

  // Initialize strain energy and spsum
  double strain_energy = 0.0;
  double spsum = 0.0;

  // **WARNING** For now Family is fixed at start and is not recomputed each time step
  // Get the family of node mi (all the nodes within its horizon, delta).
  const BondPArray& bonds = cur_node->getBonds();
//  std::cout << " Node = " << cur_node->getID() << "  Position = (" << cur_node->x() << ", " << cur_node->y() << ", " << cur_node->z()
//  << ")  " << "Bonds = " << std::endl;
//  for (auto iter = bonds.begin(); iter != bonds.end(); ++iter) {
//      std::cout << *(*iter);
//  }
  
  //const NodePArray& family_nodes = cur_node->getFamily();
  //const MaterialSPArray& bond_materials = cur_node->getBondMaterials();

  // Loop over the family of current node mi.
  for (auto family_iter = bonds.begin(); family_iter != bonds.end(); family_iter++) {

    BondP bond = *family_iter;

    const NodeP family_node = bond->second();
//    std::cout << "Bond= " << family_node << std::endl;
    if (family_node->omit()) continue;  // skip this node

//    Find the peridynamic interparticle force.
    // = force density per unit volume due to peridynamic interaction between nodes
    
    
    if (isMaterialWood() == true) {
 
    bond->computeInternalForce(d_mat_list, gridSize);
    }
    else {
    bond->computeInternalForce();
    }

//    std::cout << "bond internal force= " << bond->internalForce() << std::endl;
    // Sum up the force on node mi due to all the attached bonds.
    // force at the current configuration (n+1)
    internalForce += bond->internalForce();

    strain_energy += bond->computeStrainEnergy();
    spsum += (bond->computeMicroModulus()/(cur_node->densityNode()));
  }
  cur_node->internalForce(internalForce);
  cur_node->strainEnergy(strain_energy);
  cur_node->spSum(spsum);  
//  std::cout << *cur_node;
}


// Compute the internal force per unit volume
void
Peridynamics::computeInternalForceDensity(const NodeP cur_node,
                                          Vector3D& internalForce, const Vector3D& gridSize, const int& iter,
                                          std::vector <int>& posBrokenBond)
{

  // Initialize strain energy and spsum
  double strain_energy = 0.0;
  double spsum = 0.0;

  // **WARNING** For now Family is fixed at start and is not recomputed each time step
  // Get the family of node mi (all the nodes within its horizon, delta).
//   const BondPArray& bonds = cur_node->getBonds();
     const BondPArray& bonds = cur_node->getBonds();
//  std::cout << " Node = " << cur_node->getID() << "  Position = (" << cur_node->x() << ", " << cur_node->y() << ", " << cur_node->z()
//  << ")  " << std::endl; // << "Bonds = " << std::endl;
//  for (auto iter = bonds.begin(); iter != bonds.end(); ++iter) {
//      std::cout << *(*iter);
//  }
  
  //const NodePArray& family_nodes = cur_node->getFamily();
  //const MaterialSPArray& bond_materials = cur_node->getBondMaterials();

  // Loop over the family of current node mi.
  
  int pos = 0;
//  std::vector <int> posBrokenBond;

  for (auto family_iter = bonds.begin(); family_iter != bonds.end(); family_iter++) {

     BondP bond = *family_iter;

     const NodeP family_node = bond->second();
 //    std::cout << "Bond= " << family_node << std::endl;
     if (family_node->omit()) continue;  // skip this node

     // Find the peridynamic interparticle force.
     // = force density per unit volume due to peridynamic interaction between nodes
    
    
     if (isMaterialWood() == true) {
         bond->computeInternalForce(d_mat_list, gridSize);
     }
     else {
         bond->computeInternalForce();
     }

 //    std::cout << "bond internal force= " << bond->internalForce() << std::endl;
    // Sum up the force on node mi due to all the attached bonds.
    // force at the current configuration (n+1)
    
//     if (iter > 1) {
     bool brokenBond = false;
     if (isMaterialWood() == true) {
       brokenBond = bond->checkAndFlagBrokenBond(d_mat_list);
     }
     else {
       brokenBond = bond->checkAndFlagBrokenBond();
     }
     if (brokenBond == true) {
       posBrokenBond.emplace_back(pos);
     }
     pos++;
//    if (brokenBond == true) {
//       Vector3D zeroForce(0.0, 0.0, 0.0);
//       bond->internalForce(zeroForce);       
//    }

//    std::cout << "bond internal force= " << bond->internalForce() << std::endl;
//    std::cout << "Broken Bond= " << bond->isBroken() << "= " << brokenBond << std::endl;

     internalForce += bond->internalForce();

     strain_energy += bond->computeStrainEnergy();
     spsum += (bond->computeMicroModulus()/(cur_node->densityNode()));
  
     cur_node->internalForce(internalForce);
     cur_node->strainEnergy(strain_energy);
     cur_node->spSum(spsum);  
  }      
//  std::cout << *cur_node;
}

// Integrates the node accelerations due to peridynamic ("structured") interaction.
// Update the node velocities
// [one-step Velocity-Verlet formulation]
// 1. v(n+1/2) = v(n) + dt/2m * f(q(n))
// 2. q(n+1) = q(n) + dt * v(n+1/2)
// 3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
void
Peridynamics::integrateNodalAcceleration(const Vector3D& velOld,
                                         const Vector3D& accOld,
                                         double delT,
                                         Vector3D& velNew)
{
  velNew = velOld + accOld*delT;
  if (accOld.isnan() || velOld.isnan()) {
    std::ostringstream out;
    out << "Acceleration/old velocity is nan.  Vel = " << velOld << " acc = " << accOld ;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}

void
Peridynamics::integrateNodalVelocity(const Vector3D& dispOld,
                                     const Vector3D& velOld,
                                     double delT,
                                     Vector3D& dispNew)
{
  dispNew = dispOld + velOld*delT;
}

void 
Peridynamics::breakBonds(const NodePArray& nodes)
{
  // First flag and delete broken bonds attached to each node
  for (auto iter = nodes.begin(); iter != nodes.end(); iter++) {
    NodeP cur_node = *iter;
    if (cur_node->omit()) continue;  // skip this node

    // Break bonds and update the damage index
    if (isMaterialWood() == true) {
    cur_node->findAndDeleteBrokenBonds(d_mat_list);
    }
    else {
    cur_node->findAndDeleteBrokenBonds();
//    std::cout << std::endl << "*********************check if the Bond is Broken**************************" << std::endl << std::endl;
    }
  }

  // After all bond deletions have been completed update nodal damage indices
  for (auto iter = nodes.begin(); iter != nodes.end(); iter++) {
    NodeP cur_node = *iter;
    if (cur_node->omit()) continue;  // skip this node

    cur_node->updateDamageIndex();
  }
}

void
Peridynamics::checkMemoryUsage(double& resident_mem, double& shared_mem)
{
  int tSize = 0, resident = 0, share = 0;
  std::ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  resident_mem = resident * page_size_kb;
  shared_mem = share * page_size_kb;
}

