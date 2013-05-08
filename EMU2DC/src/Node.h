#ifndef __EMU2DC_NODE_H__
#define __EMU2DC_NODE_H__

#include <ElementPArray.h>
#include <NodePArray.h>
#include <Types.h>
#include <iostream>
#include <cmath>

namespace Emu2DC {
    
  // This structure defines the node type
  class Node {

    public:

      friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Node& node);

    public:

      Node();
      Node(const int id, const double xx, const double yy, const double zz, const int surfaceNode);
      ~Node();

      bool operator<(const Node& node) const;

      inline void dimension(const int dim) {d_dimension = dim;}
      inline int dimension() const {return d_dimension;}

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

      inline const double& volume() const { return d_volume; }
      inline void volume(const double& volume) { d_volume = volume; }

      inline const double& density() const { return d_density; }
      inline void density(const double& density) { d_density = density; }

      inline const double& youngModulus() const { return d_young; }
      inline void youngModulus(const double& young) { d_young = young; }

      inline const double& strainEnergy() const { return d_strain_energy; }
      inline void strainEnergy(const double& strain_energy) { d_strain_energy = strain_energy; }

      inline const double& damageIndex() const { return d_damage_index; }
      inline void damageIndex(const double& damage_index) { d_damage_index = damage_index; }

      inline const Array3& position() const { return d_pos; }
      inline void position(const Array3& pos)  { d_pos = pos; }

      inline const Array3& displacement() const { return d_disp; }
      inline void displacement(const Array3& disp)  { d_disp = disp; }

      inline const Array3& oldDisplacement() const { return d_old_disp; }
      inline void oldDisplacement(const Array3& disp)  { d_old_disp = disp; }

      inline const Array3& newDisplacement() const { return d_new_disp; }
      inline void newDisplacement(const Array3& disp)  { d_new_disp = disp; }

      inline const Array3& velocity() const { return d_veloc; }
      inline void velocity(const Array3& veloc)  { d_veloc = veloc; }

      inline const Array3& acceleration() const { return d_accel; }
      inline void acceleration(const Array3& accel)  { d_accel = accel; }

      inline const Array3& force() const { return d_force; }
      inline void force(const Array3& force)  { d_force = force; }

      inline const Array3& externalForce() const { return d_force; }
      inline void externalForce(const Array3& extForce)  { d_force = extForce; }

      inline int numAdjacentElements() const { return d_adjacent_elements.size(); }

      inline double distance(const Node& node) const
      {
        double dx = d_pos[0] - node.d_pos[0];
        double dy = d_pos[1] - node.d_pos[1];
        double dz = d_pos[2] - node.d_pos[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
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
      void setFamily(const NodePArray& fam) {d_neighbor_list = fam;}
      const NodePArray& getFamily() const {return d_neighbor_list;}

    private:

      int d_dimension;
      long64 d_id;
      int d_mat_type;
      double d_horizon_size;
      bool d_omit;         // Omit this node from the computation if true
      bool d_surfaceNode;  // This node is on the surface of the body if true

      double d_volume;
      double d_density;
      double d_young;

      double d_strain_energy;
      double d_damage_index;

      ElementPArray d_adjacent_elements; // The elements adjacent to this node, 
      NodePArray d_neighbor_list;        // The nodes inside the horizon of this node

      Array3 d_pos;  // array 
      Array3 d_disp;  // array
      Array3 d_veloc;  // array
      Array3 d_accel;  // array
      Array3 d_new_veloc;  // array
      Array3 d_new_disp;  // array
      Array3 d_old_disp;  // array
      Array3 d_force;  // array

      int* d_bc;
  };

} // end namespace

#endif
