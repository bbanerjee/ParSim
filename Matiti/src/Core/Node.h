/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MATITI_NODE_H__
#define __MATITI_NODE_H__

#include <Containers/ElementPArray.h>
#include <Containers/BondPArray.h>
#include <Containers/MaterialSPArray.h>
#include <Pointers/MaterialUP.h>
#include <Pointers/ElementP.h>
#include <MaterialModels/Material.h>

//#include <NodePArray.h>
//#include <MaterialUPArray.h>

#include <Types/Types.h>
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <iostream>
#define _USE_MATH_DEFINE
#include <cmath>

namespace Matiti {
    
  // This structure defines the node type
  class Node {

    public:

      friend std::ostream& operator<<(std::ostream& out, const Matiti::Node& node);

    public:

      Node();
      Node(const int id, const double xx, const double yy, const double zz, const int surfaceNode);
      Node(const Node& node);
      ~Node();

      bool operator<(const Node& node) const { return (d_id < node.d_id); }

      /**
       * Compute stable timestep
       * Inputs: timestep reduction factor
       */
      double computeStableTimestep(const double& fac) const;

      /**
       * Compute initial displacement
       * Inputs: initial velocity
       *         time increment
       */
      void computeInitialDisplacement(const Vector3D& initVel, double delT);
      
//      void getInterval(Array3& interval);

      inline void omit(const bool& omit) { d_omit = omit; }
      inline bool omit() const { return d_omit; }

      inline void onSurface(const bool& flag) { d_surfaceNode = flag; }
      inline bool onSurface() const { return d_surfaceNode; }
      
      inline const long64& getID() const { return d_id; }
      inline void setID(const long64& id) { d_id = id; }

      inline const int& matType() const { return d_mat_type; }
      inline void matType(const int& mat_type) { d_mat_type = mat_type; }

      inline const double& horizonSize() const { return d_horizon_size; }
      inline void horizonSize(const double& horizon_size) { d_horizon_size = horizon_size; }

      inline const double& area() const { return d_area; }
      inline void area(const double& area) { d_area = area; }

      inline const double& densityNode() const { return d_density; }
      inline void densityNode(const double& density) { d_density = density; }

      inline const double& volume() const { return d_volume; }
      inline const double& radius() const { return d_radius; }
      inline void volume(const double& volume) { 
         d_volume = volume; 
         d_radius = std::pow(0.75*d_volume/M_PI, (1.0/3.0));
      }


      /**
       * Assign node material
       * ** WARNING** This is only for the initial setup before bonds are computed.
       */
      void assignMaterial(const Material* mat) { d_material->clone(mat);}
      void assignMaterial(const Material* mat,
                          double randomNum,
                          double coeffOfVar) { 
         d_material->clone(mat, randomNum, coeffOfVar);
      }
      const Material* material() const {return d_material.get();}
      double density() const {return d_material->density();}

      inline double x() const {return d_pos.x();}
      inline double y() const {return d_pos.y();}
      inline double z() const {return d_pos.z();}

      inline const Point3D& position() const { return d_pos; }
      inline void position(const Point3D& pos)  { d_pos = pos; }

      inline const Vector3D& displacement() const { return d_disp; }
      inline void displacement(const Vector3D& disp)  { d_disp = disp; }
      inline void xDisplacement(double disp)  { d_disp[0] = disp; }
      inline void yDisplacement(double disp)  { d_disp[1] = disp; }
      inline void zDisplacement(double disp)  { d_disp[2] = disp; }

      inline const Vector3D& newDisplacement() const { return d_disp_new; }
      inline void newDisplacement(const Vector3D& disp)  { d_disp_new = disp; }
      inline void xNewDisplacement(double disp)  { d_disp_new[0] = disp; }
      inline void yNewDisplacement(double disp)  { d_disp_new[1] = disp; }
      inline void zNewDisplacement(double disp)  { d_disp_new[2] = disp; }

      inline const Vector3D& velocity() const { return d_vel; }
      inline void velocity(const Vector3D& vel)  { d_vel = vel; }
      inline void xVelocity(double vel)  { d_vel[0] = vel; }
      inline void yVelocity(double vel)  { d_vel[1] = vel; }
      inline void zVelocity(double vel)  { d_vel[2] = vel; }

      inline const Vector3D& newVelocity() const { return d_vel_new; }
      inline void newVelocity(const Vector3D& vel)  { d_vel_new = vel; }
      inline void xNewVelocity(double vel)  { d_vel_new[0] = vel; }
      inline void yNewVelocity(double vel)  { d_vel_new[1] = vel; }
      inline void zNewVelocity(double vel)  { d_vel_new[2] = vel; }

      inline const Vector3D& acceleration() const { return d_accel; }
      inline void acceleration(const Vector3D& accel)  { d_accel = accel; }

      inline const Vector3D& internalForce() const { return d_int_force; } 
      inline void internalForce(const Vector3D& internalForce) {d_int_force = internalForce;}
      inline void xInternalForce(double force)  { d_int_force[0] = force; }
      inline void yInternalForce(double force)  { d_int_force[1] = force; }
      inline void zInternalForce(double force)  { d_int_force[2] = force; }

      inline const Vector3D& externalForce() const { return d_ext_force; }
      inline void externalForce(const Vector3D& extForce)  { d_ext_force = extForce; }
      inline void xExternalForce(double force)  { d_ext_force[0] = force; }
      inline void yExternalForce(double force)  { d_ext_force[1] = force; }
      inline void zExternalForce(double force)  { d_ext_force[2] = force; }

      inline const Array3& getInterval() const { return d_interval;}
      inline void setInterval(const Array3& interval) {d_interval[0] = interval[0];
                                                       d_interval[1] = interval[1];
                                                       d_interval[2] = interval[2];} 

      inline int numAdjacentElements() const { return d_adjacent_elements.size(); }

      inline double distance(const Node& node) const
      {
        return d_pos.distance(node.d_pos);
      }

      const ElementPArray& getAdjacentElements() const
      {
        return d_adjacent_elements;
      }

      const ElementP& getAdjacentElement(const int& index) const
      {
        return d_adjacent_elements[index];
      }

      inline void addAdjacentElement(const ElementP& elem) 
      {
        d_adjacent_elements.push_back(elem);
      }

      // Node "family" = neighbor list access methods
      void setBonds(const BondPArray& fam) { d_bonds.clear(); d_bonds = fam; }
      BondPArray& getBonds() {return d_bonds;}

      // void setFamily(const NodePArray& fam);
      // const NodePArray& getFamily() const {return d_neighbor_list;}
      // const MaterialUPArray& getBondMaterials() const {return d_bond_materials;}

      void initialFamilySize(const int size) {d_initial_family_size = size;}
      int initialFamilySize() const {return d_initial_family_size;}
      int currentFamilySize() const {return (int) d_bonds.size();}

      /**
       *  Find and delete broken bonds
       */
      void findAndDeleteBrokenBonds();
      void findAndDeleteBrokenBonds(const MaterialSPArray& matList);

      /**
       *  Damage index access methods
       */
      void updateDamageIndex();
      double damageIndex() const {return d_damage_index;}

      /**
       * Store some data in case it's needed later
       */
      void strainEnergy(double energy) {d_strain_energy = energy;}
      void spSum(double spsum) {d_sp_sum = spsum;}

      /**
       * Hack to stop nodes from flying apart when subjected to an external load
       */
      inline void allowFailure(bool flag) {d_allow_failure = flag;}
      inline bool failureAllowed() const {return d_allow_failure;}


    private:

      long64 d_id;
      int d_mat_type;
      double d_horizon_size;
      bool d_omit;                  // Omit this node from the computation if true
      bool d_surfaceNode;           // This node is on the surface of the body if true
                                    // TODO: The surface can be a crack surface.
      double d_area;                // Zero if inside, non-zero on the surface 
      double d_volume;              // Volume of the node
      double d_radius;              // Radius of the ball containing the volume of the node

      double d_density;

      MaterialUP d_material;  // For initial setup  **WARNING** Potential problems.

      ElementPArray d_adjacent_elements; // The elements adjacent to this node, 
      
      // Using array of structs instead of struct of arrays
      BondPArray d_bonds;                // The bonds attached to this node
      int d_initial_family_size;

      // NodePArray d_neighbor_list;        // The nodes inside the horizon of this node
      // MaterialUPArray d_bond_materials;  // One material per bond to store history

      // Ideally all these variables should be arrays indexed by node ID for faster access
      // TODO: At some point in the future refactor the entire code to make these into arrays.
      Point3D  d_pos; 
      Vector3D d_disp;  
      Vector3D d_disp_new;  
      Vector3D d_vel;  
      Vector3D d_vel_new;  
      Vector3D d_accel;  
      Vector3D d_int_force;
      Vector3D d_ext_force;
      Array3 d_interval;  

      // Not really necessary but storing for now
      double d_strain_energy;
      double d_sp_sum;

      double d_damage_index;
      bool d_allow_failure;  // Hack to stop boundary nodes at which external forces are applied
                             // from flying apart.
  };

} // end namespace

#endif
