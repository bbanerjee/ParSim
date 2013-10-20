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

template<typename T>
void MPMDatawarehouse::add(const std::string& label, int dwi,
                           T& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  MPMVar var(val);
  d_var[label_dwi] = var;
}

void MPMDatawarehouse::zero(const std::string& label, int dwi)
{
  std::string label_dwi = label + std::to_string(dwi);
  boost::apply_visitor(ZeroVisitor(), d_var[label_dwi]);
}

template<typename T>
void MPMDatawarehouse::get(const std::string& label, int dwi, T& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = boost::get<T>(d_var[label_dwi]);
}

template<typename T>
void MPMDatawarehouse::put(const std::string& label, int dwi, T& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  MPMVar var(val);
  boost::get<T>(d_var[label_dwi]) = var;
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
  add("pX", dwi, pX);
  Point3DParticleData px(pX);
  add("px", dwi, px);
  add("pN", dwi, pN);
  Vector3DParticleData pn(pN);
  add("pn", dwi, pn);
  add("pVol", dwi, pVol);
  DoubleParticleData pm = pVol*density;
  add("pm", dwi, pm);

  // Create default values of other particle variables
  Vector3DParticleData pw(numPart, Vector3D(0.0));
  add("pw", dwi, pw);
  Vector3DParticleData pvI(numPart, Vector3D(0.0));
  add("pvI", dwi, pvI);
  Vector3DParticleData pxI(numPart, Vector3D(0.0));
  add("pxI", dwi, pxI);
  Vector3DParticleData pfe(numPart, Vector3D(0.0));
  add("pfe", dwi, pfe);
  Vector3DParticleData pfi(numPart, Vector3D(0.0));
  add("pfi", dwi, pfi);
  Vector3DParticleData pfc(numPart, Vector3D(0.0));
  add("pfc", dwi, pfc);
  Vector3DParticleData pwc(numPart, Vector3D(0.0));
  add("pwc", dwi, pwc);
  Matrix3DParticleData pGv(numPart, Matrix3D(0.0));
  add("pGv", dwi, pGv);
  Matrix3DParticleData pVS(numPart, Matrix3D(0.0));
  add("pVS", dwi, pVS);
  Matrix3D one; one.Identity();
  Matrix3DParticleData pF(numPart, one);
  add("pF", dwi, pF);

  // Create the interpolation information
  std::vector<int> zeroVecInt(numNearNodes, 0);
  VectorIntParticleData cIdx(numPart, zeroVecInt);
  add("cIdx", dwi, cIdx);
  std::vector<double> zeroVecDouble(numNearNodes, 0.0);
  VectorDoubleParticleData cW(numPart, zeroVecDouble);
  add("cW", dwi, cW);
  VectorDoubleParticleData cGradx(numPart, zeroVecDouble);
  add("cGradx", dwi, cGradx);
  VectorDoubleParticleData cGrady(numPart, zeroVecDouble);
  add("cGrady", dwi, cGrady);
  VectorDoubleParticleData cGradz(numPart, zeroVecDouble);
  add("cGradz", dwi, cGradz);
}

/*
void MPMDatawarehouse::createGrid(int dwi, MPMPatchP& patch)
{
  DoubleNodeData gx = patch.initGrid();
  add("gx", dwi, gx);
  zeroGrid(dwi);
}

void MPMDatawarehouse::zeroGrid(int dwi)
{
  DoubleNodeData gx;
  get("gx", dwi, gx);
  // TODO: Add others after patch is done
}
*/

