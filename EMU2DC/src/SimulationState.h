#ifndef EMU2DC_SIMULATION_STATE_H
#define EMU2DC_SIMULATION_STATE_H

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>
#include <string>

namespace Emu2DC {

  class SimulationState {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::SimulationState& state);

  public:
 
    SimulationState();
    ~SimulationState();

    void initialize(const Uintah::ProblemSpecP& ps);

    inline double maxTime() const { return d_max_time;}
    inline double delT() const { return d_delT;}
    inline int maxIter() const { return d_max_iter;}

    inline void outputFolder(const std::string& folder) {d_output_folder_name = folder;}
    inline std::string outputFolder() const {return d_output_folder_name;}
    inline std::string outputFile() const {return d_output_file_name;}
    inline int outputIteratonInterval() const {return d_output_iter_interval;}

    enum class ModulusType {
      Constant=0,
      Conical=1
    };
    inline bool isDynamic() const {return d_is_dynamic;}
    inline ModulusType modulusType() const {return d_modulus_type;}
    inline double horizonFactor() const {return d_horizon_factor;}
    inline int dimensions() const {return d_dimensions;}

  private:

    // Simulation time 
    double d_max_time;  // Maximum time for which simulation is to be run
    double d_delT;      // Timestep size
    int d_max_iter;     // Maximum number of iterations

    //  Output file folder and name 
    std::string d_output_folder_name;
    std::string d_output_file_name;
    int d_output_iter_interval;

    // Peridynamics flags
    int d_dimensions;           // 2D or 3D
    bool d_is_dynamic;          // static or dyanamic simulation
    ModulusType d_modulus_type; // constant = 0 or conical = 1
    double d_horizon_factor;    // Factor with which to multiply the nodal volume to get
                                // horizon size 
  }; // end class
   
} // end namespace

#endif
