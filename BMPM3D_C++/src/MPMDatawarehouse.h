#ifndef __BRMPM_MPMDATAWAREHOUSE__
#define __BRMPM_MPMDATAWAREHOUSE__

#include <MPMDataTypes.h>
#include <MPMPatchP.h>
//#include <Output.h>
//#include <OutputVTK.h>
//#include <MPMTime.h>
// #include <MPMShapeFunction.h>

#include <map>

namespace BrMPM {

  class MPMDatawarehouse {

  public:
    MPMDatawarehouse();
//     MPMDatawarehouse(Uintah::ProblemSpecP& ps);
    ~MPMDatawarehouse();

    /*
    void saveData(double dt, MaterialSPArray& matlist);

    void dumpData(double dt, MaterialSPArray& matlist);

    bool checkSave(double dt);
    */

    // Add a variable to the data warehouse
    template<typename T>
    void add(const std::string& label, int dwi, T& val);

    // Zero the variable data in the data warehouse
    void zero(const std::string& label, int dwi);

    // Get the particle data
    template<typename T>
    void get(const std::string& label, int dwi, T& val);

    // Update the datawarehouse with new the particle data
    template<typename T>
    void put(const std::string& label, int dwi, T& val);

    // Add particles to the datawarehouse
    void addParticles(const int& dwi, Point3DParticleData& pX, DoubleParticleData& pVol,
                      Vector3DParticleData& pN, DoubleParticleData& density,
                      const int& shSize);

    void createGrid(int dwi, MPMPatchP& patch);

    void zeroGrid(int dwi);

  private:

    // Datawarehouse data map
    // The string contains the join of the variable name and the 
    // material index (dwi)
    int d_id;
    std::map<std::string, MPMVar> d_var;

//     OutputVTK d_out;
//     MPMTime d_time;
//     MPMShapeFunction d_shapefunction;

  }; // end class

}

 // end namespace

#endif
