#include <cmath>
#include <math.h>
#include <stdlib.h>

#include <MPMDatawarehouse.h>

using namespace BrMPM;

MPMDatawarehouse::MPMDatawarehouse()
  : d_id(0)
//             : d_id(0), d_time(new MPMTime())
{
}

//MPMDatawarehouse::MPMDatawarehouse(Uintah::ProblemSpecP& ps)
//{
//   d_shapefunction.initialise(ps);
//}


MPMDatawarehouse::~MPMDatawarehouse()
{
}

/*
void
MPMDatawarehouse::saveData(double dt, MaterialSPArray& matlist)
{
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.saveData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id += 1;
}

 
void
MPMDatawarehouse::dumpData(double dt, MaterialSPArray& matlist)
{
   if (checkSave(dt)) {
           d_out.outputFileCount(d_save.dumpData(d_out.outputFileCount(), matlist));
        }
   incrementTime(dt);
   d_id+=1;
}

bool
MPMDatawarehouse::checkSave(double dt)
{
  double dr=d_time.currentTime()/dt;
  double dt0=dt*std::min(dr-std::floor(dr), std::ceil(dr)-dr);
  return dt0<dt/2;
}
*/

void MPMDatawarehouse::addParticleVar(const std::string& label, int dwi,
                                      MPMParticleVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  d_particles[label_dwi] = val;
}

void MPMDatawarehouse::addNodeVar(const std::string& label, int dwi,
                                  MPMNodeVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  d_nodes[label_dwi] = val;
}

void MPMDatawarehouse::addInterpolationVar(const std::string& label, int dwi,
                                           MPMInterpolationVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  d_interp[label_dwi] = val;
}

template<typename T>
void MPMDatawarehouse::zeroParticleVar(const std::string& label, int dwi)
{
  std::string label_dwi = label + std::to_string(dwi);
  MPMData<T> zero(0.0);
  boost::get<T>(d_particles[label_dwi]) = zero;
}

template<typename T>
void MPMDatawarehouse::zeroNodeVar(const std::string& label, int dwi)
{
  std::string label_dwi = label + std::to_string(dwi);
  MPMData<T> zero(0.0);
  boost::get<T>(d_nodes[label_dwi]) = zero;
}

template<typename T>
void MPMDatawarehouse::zeroInterpolationVar(const std::string& label, int dwi)
{
  std::string label_dwi = label + std::to_string(dwi);
  MPMData<T> zero(0.0);
  boost::get<T>(d_interp[label_dwi]) = zero;
}

template<typename T>
void MPMDatawarehouse::getParticleVar(const std::string& label, int dwi,
                                      MPMParticleVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = d_particles[label_dwi];
}

template<typename T>
void MPMDatawarehouse::getNodeVar(const std::string& label, int dwi,
                                  MPMNodeVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = d_nodes[label_dwi];
}

template<typename T>
void MPMDatawarehouse::getInterpolationVar(const std::string& label, int dwi,
                                           MPMInterpolationVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = d_interp[label_dwi];
}

void MPMDatawarehouse::addParticles(const int& dwi,
                                    Point3DParticleData& pX,
                                    DoubleParticleData& pVol,
                                    Vector3DParticleData& pN,
                                    DoubleParticleData& density,
                                    const int& numNearNodes)
{
  // get the number of particles
  int numPart = pX.size();

  // Add initial position, position, volume, mass from inputs
  addParticleVar("pX", dwi, pX);
  Point3DParticleData px = pX.clone();
  addParticleVar("px", dwi, px);
  addParticleVar("pN", dwi, pN);
  Vector3DParticleData pn = pN.clone();
  addParticleVar("pn", dwi, pn);
  addParticleVar("pVol", dwi, pVol);
  DoubleParticleData pm = pVol*density;
  addParticleVar("pm", dwi, pm);

  // Create default values of other particle variables
  Vector3DParticleData pw(numPart, Vector3D(0.0));
  addParticleVar("pw", dwi, pw);
  Vector3DParticleData pvI(numPart, Vector3D(0.0));
  addParticleVar("pvI", dwi, pvI);
  Vector3DParticleData pxI(numPart, Vector3D(0.0));
  addParticleVar("pxI", dwi, pxI);
  Vector3DParticleData pfe(numPart, Vector3D(0.0));
  addParticleVar("pfe", dwi, pfe);
  Vector3DParticleData pfi(numPart, Vector3D(0.0));
  addParticleVar("pfi", dwi, pfi);
  Vector3DParticleData pfc(numPart, Vector3D(0.0));
  addParticleVar("pfc", dwi, pfc);
  Vector3DParticleData pwc(numPart, Vector3D(0.0));
  addParticleVar("pwc", dwi, pwc);
  Matrix3DParticleData pGv(numPart, Matrix3D(0.0));
  addParticleVar("pGv", dwi, pGv);
  Matrix3DParticleData pVS(numPart, Matrix3D(0.0));
  addParticleVar("pVS", dwi, pVS);
  Matrix3D one; one.Identity();
  Matrix3DParticleData pF(numPart, one);
  addParticleVar("pF", dwi, pF);

  // Create the interpolation information
  std::vector<int> zeroVecInt(numNearNodes, 0);
  VectorIntParticleData cIdx(numPart, zeroVecInt);
  addParticleVar("cIdx", dwi, cIdx);
  std::vector<double> zeroVecDouble(numNearNodes, 0.0);
  VectorDoubleParticleData cW(numPart, zeroVecDouble);
  addParticleVar("cW", dwi, cW);
  VectorDoubleParticleData cGradx(numPart, zeroVecDouble);
  addParticleVar("cGradx", dwi, cGradx);
  VectorDoubleParticleData cGrady(numPart, zeroVecDouble);
  addParticleVar("cGrady", dwi, cGrady);
  VectorDoubleParticleData cGradz(numPart, zeroVecDouble);
  addParticleVar("cGradz", dwi, cGradz);
}

/*
void MPMDatawarehouse::createGrid(int dwi, MPMPatchP& patch)
{
  DoubleNodeData gx = patch.initGrid();
  addNodeVar("gx", dwi, gx);
  zeroGrid(dwi);
}

void MPMDatawarehouse::zeroGrid(int dwi)
{
  DoubleNodeData gx;
  getNodeVar("gx", dwi, gx);
  // TODO: Add others after patch is done
}
*/

