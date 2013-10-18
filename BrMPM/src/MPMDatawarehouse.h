#ifndef __BRMPM_MPMDATAWAREHOUSE__
#define __BRMPM_MPMDATAWAREHOUSE__

//#include <Output.h>
//#include <OutputVTK.h>
//#include <MPMTime.h>
// #include <MPMShapeFunction.h>

#include <vector>
#include <string>
#include <map>

#include <boost/variant.hpp>

namespace BrMPM {

  // Define types
  typedef std::vector<int>          IntegerParticleData;
  typedef std::vector<double>       DoubleParticleData;
  typedef std::vector<std::string>  StringParticleData;
  typedef std::vector<Point3D>      Point3DParticleData;
  typedef std::vector<Vector3D>     Vector3DParticleData;
  typedef std::vector<Matrix3D>     Matrix3DParticleData;

  typedef std::vector<CellIndexVector>          VectorIntParticleData;
  typedef std::vector<CellInterpolationVector>  VectorDoubleParticleData;

  typedef std::vector<double>    DoubleNodeData;
  typedef std::vector<Vector3D>  Vector3DNodeData;

  // Define the variants
  typedef boost::variant<IntegerParticleData, DoubleParticleData, Point3DParticleData,
                         Vector3DParticleData, Matrix3DParticleData> MPMParticleVar;
  typedef boost::variant<DoubleNodeData, Vector3DNodeData> MPMNodeVar;
  typedef boost::variant<VectorIntParticleData, VectorDoubleParticleData> MPMInterpolationVar;

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
    void addParticleVar(const std::string& label, int dwi, MPMParticleVar& val);
    void addNodeVar(const std::string& label, int dwi, MPMNodeVar& val);
    void addInterpolationVar(const std::string& label, int dwi, MPMInterpolationVar& val);

    // Zero the variable data in the data warehouse
    void zeroParticleVar(const std::string& label, int dwi);
    void zeroNodeVar(const std::string& label, int dwi);
    void zeroInterpolationVar(const std::string& label, int dwi);

    // Get the particle data
    void getParticleVar(const std::string& label, int dwi, MPMParticleVar& val);
    void getNodeVar(const std::string& label, int dwi, MPMNodeVar& val);
    void getInterpolationVar(const std::string& label, int dwi, MPMInterpolationVar& val);

    // Add particles to the datawarehouse
    void addParticles(const int& dwi, Point3DParticleData& pX, DoubleParticleData& pVol,
                      Vector3DParticleData& pN, DoubleParticleData& density,
                      const int& shSize);

    /*
    void createGrid(int dwi, MPMPatchP& patch);

    void zeroGrid(int dwi);
    */

  private:

    // Datawarehouse data map
    // The string contains the join of the variable name and the 
    // material index (dwi)
    int d_id;
    std::map<std::string, MPMParticleVar> d_particles;
    std::map<std::string, MPMNodeVar> d_nodes;
    std::map<std::string, MPMInterpolationVar> d_interp;

//     OutputVTK d_out;
//     MPMTime d_time;
//     MPMShapeFunction d_shapefunction;

  }; // end class

}

 // end namespace

#endif
