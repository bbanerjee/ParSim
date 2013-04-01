#ifndef EMU2DC_NODE_H
#define EMU2DC_NODE_H

#include <array>
#include <vector>
#include <iostream>


namespace Emu2DC {
  typedef std::array<double, 3> Array3;

  class Vector3 : public Array3 {
    public:
    
    Vector3 operator+(const Vector3 vec) {
      Vector3 out_vec;
      out_vec[0] = (*this)[0] + vec[0]; 
      out_vec[1] = (*this)[1] + vec[1]; 
      out_vec[2] = (*this)[2] + vec[2]; 
      return out_vec;
    }
  };
}

namespace Emu2DC {
    
  // This structure defines the node type
  class Node {

    public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Node& node);

      Node();
      ~Node();

      bool operator<(const Node& node) const;

      inline void omit(const bool& omit)
      {
        d_omit = omit;
      }

      inline bool omit() const
      {
        return d_omit;
      }
      
      inline int getID() const
      {
        return d_id;
      }

      inline void setID(const int& id)
      {
        d_id = id;
      }

      inline int getMatType() const
      {
        return d_mat_type;
      }

      inline void setMatType(const int& mat_type)
      {
        d_mat_type = mat_type;
      }

      inline double getHorizonSize() const
      {
        return d_horizon_size;
      }

      inline void setHorizonSize(const double& horizon_size)
      {
        d_horizon_size = horizon_size;
      }

      inline bool isHangingNode() const
      {
        return d_iflag;
      }

      inline void makeHangingNode(const bool& flag)
      {
        d_iflag = flag;
      }

      inline double getVolume() const
      {
        return d_volume;
      }

      inline void setVolume(const double& volume)
      {
        d_volume = volume;
      }

      inline double getDensity() const
      {
        return d_density;
      }

      inline void setDensity(const double& density)
      {
        d_density = density;
      }

      inline double getYoung() const
      {
        return d_young;
      }

      inline void setYoung(const double& young)
      {
        d_young = young;
      }

      inline double getStrainEnergy() const
      {
        return d_strain_energy;
      }

      inline void setStrainEnergy(const double& strain_energy)
      {
        d_strain_energy = strain_energy;
      }

      inline double getDamageIndex() const
      {
        return d_damage_index;
      }

      inline void setDamageIndex(const double& damage_index)
      {
        d_damage_index = damage_index;
      }

      inline void getPosition(Array3& pos) const 
      {
        pos = d_pos;
      }

      inline void setPosition(const Array3& pos)  
      {
        d_pos = pos;
      }

      inline void getDisplacement(Array3& disp) const 
      {
        disp = d_disp;
      }

      inline void setDisplacement(const Array3& disp)  
      {
        d_disp = disp;
      }

      inline void getOldDisplacement(Array3& disp) const 
      {
        disp = d_old_disp;
      }

      inline void setOldDisplacement(const Array3& disp)  
      {
        d_old_disp = disp;
      }

      inline void getNewDisplacement(Array3& disp) const 
      {
        disp = d_new_disp;
      }

      inline void setNewDisplacement(const Array3& disp)  
      {
        d_new_disp = disp;
      }

      inline void getVelocity(Array3& veloc) const 
      {
        veloc = d_veloc;
      }

      inline void setVelocity(const Array3& veloc)  
      {
        d_veloc = veloc;
      }

      inline void getAcceleration(Array3& accel) const 
      {
        accel = d_accel;
      }

      inline void setAcceleration(const Array3& accel)  
      {
        d_accel = accel;
      }

      inline void getForce(Array3& force) const 
      {
        force = d_force;
      }

      inline void setForce(const Array3& force)  
      {
        d_force = force;
      }

      inline int getNumNodeElements() const
      {
        return d_nnodeelements;
      }

      inline void getAdjacentElements(std::vector<int>& nodeelements) const
      {
        nodeelements = d_nodeelements;
      }

      inline int getAdjacentElement(const int& index) const
      {
        return d_nodeelements[index];
      }

      inline void addAdjacentElement(const int& elem_id) 
      {
        d_nodeelements.push_back(elem_id);
        d_nodeelements_size.push_back(0.0);
        d_nodeelements_depth.push_back(0);
        d_nodeelements_nhanging_nodes.push_back(0);
        d_nnodeelements = d_nodeelements.size();
      }

      inline void getAdjacentElementsSize(std::vector<double>& nodeelements_size) const
      {
        nodeelements_size = d_nodeelements_size;
      }

      inline double getAdjacentElementSize(const int& index) const
      {
        return d_nodeelements_size[index];
      }

      inline void updateAdjacentElementSize(const int& index, const double& elem_size) 
      {
        d_nodeelements_size[index] = elem_size;
      }

      inline void getAdjacentElementsDepth(std::vector<int>& nodeelements_depth) const
      {
        nodeelements_depth = d_nodeelements_depth;
      }

      inline int getAdjacentElementDepth(const int& index) const
      {
        return d_nodeelements_depth[index];
      }

      inline void updateAdjacentElementDepth(const int& index, const int& elem_depth) 
      {
        d_nodeelements_depth[index] = elem_depth;
      }

      inline void getAdjacentElementsNumHangingNodes(std::vector<int>& nodeelements_nhanging_nodes) const
      {
        nodeelements_nhanging_nodes = d_nodeelements_nhanging_nodes;
      }

      inline int getAdjacentElementNumHangingNodes(const int& index) const
      {
        return d_nodeelements_nhanging_nodes[index];
      }

      inline void updateAdjacentElementNumHangingNodes(const int& index) 
      {
        d_nodeelements_nhanging_nodes[index] = d_nodeelements_nhanging_nodes[index] + 1; 
      }

    private:

      bool d_omit;  // Omit this node from the computation

      int d_id;
      int d_mat_type;
      double d_horizon_size;
      bool d_iflag;  // iflag = 1: node is a hanging node; 
                   // iflag = 0: node is not a hanging node

      double d_volume;
      double d_density;
      double d_young;

      double d_strain_energy;
      double d_damage_index;

      Array3 d_pos;  // array 
      Array3 d_disp;  // array
      Array3 d_veloc;  // array
      Array3 d_accel;  // array
      Array3 d_new_veloc;  // array
      Array3 d_new_disp;  // array
      Array3 d_old_disp;  // array
      Array3 d_force;  // array

      int* d_bc;
    
      int d_nnodeelements;
      std::vector<int> d_nodeelements; // The id of the elements adjacent to this node, 
                                       // can be 1,2,3 or 4 elements
      std::vector<double> d_nodeelements_size; // The size(area) of the elements adjacent to this node, 
                                               // can be 1,2,3 or 4 elements
      std::vector<int> d_nodeelements_depth; // The depth of each element adjacent to this node
      std::vector<int> d_nodeelements_nhanging_nodes; // The number of hanging nodes for each element 
                                                      // connected to this node
  };

} // end namespace

#endif
