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

void MPMDatawarehouse::add(const std::string& label, int dwi,
                           MPMVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  d_var[label_dwi] = val;
}

void MPMDatawarehouse::zero(const std::string& label, int dwi)
{
  std::string label_dwi = label + std::to_string(dwi);
  boost::apply_visitor(ZeroVisitor(), d_var[label_dwi]);
}

void MPMDatawarehouse::get(const std::string& label, int dwi,
                           MPMVar& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = d_var[label_dwi];
}

template<typename T>
void MPMDatawarehouse::get(const std::string& label, int dwi, T& val)
{
  std::string label_dwi = label + std::to_string(dwi);
  val = boost::get<T>(d_var[label_dwi]);
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
  MPMVar pXVar(pX);
  add("pX", dwi, pXVar);
  Point3DParticleData px(pX);
  MPMVar pxVar(px);
  add("px", dwi, pxVar);
  MPMVar pNVar(pN);
  add("pN", dwi, pNVar);
  Vector3DParticleData pn(pN);
  MPMVar pnVar(pn);
  add("pn", dwi, pnVar);
  MPMVar pVolVar(pVol);
  add("pVol", dwi, pVolVar);
  DoubleParticleData pm = pVol*density;
  MPMVar pmVar(pm);
  add("pm", dwi, pmVar);

  // Create default values of other particle variables
  Vector3DParticleData pw(numPart, Vector3D(0.0));
  MPMVar pwVar(pw);
  add("pw", dwi, pwVar);
  Vector3DParticleData pvI(numPart, Vector3D(0.0));
  MPMVar pvIVar(pvI);
  add("pvI", dwi, pvIVar);
  Vector3DParticleData pxI(numPart, Vector3D(0.0));
  MPMVar pxIVar(pxI);
  add("pxI", dwi, pxIVar);
  Vector3DParticleData pfe(numPart, Vector3D(0.0));
  MPMVar pfeVar(pfe);
  add("pfe", dwi, pfeVar);
  Vector3DParticleData pfi(numPart, Vector3D(0.0));
  MPMVar pfiVar(pfi);
  add("pfi", dwi, pfiVar);
  Vector3DParticleData pfc(numPart, Vector3D(0.0));
  MPMVar pfcVar(pfc);
  add("pfc", dwi, pfcVar);
  Vector3DParticleData pwc(numPart, Vector3D(0.0));
  MPMVar pwcVar(pwc);
  add("pwc", dwi, pwcVar);
  Matrix3DParticleData pGv(numPart, Matrix3D(0.0));
  MPMVar pGvVar(pGv);
  add("pGv", dwi, pGvVar);
  Matrix3DParticleData pVS(numPart, Matrix3D(0.0));
  MPMVar pVSVar(pVS);
  add("pVS", dwi, pVSVar);
  Matrix3D one; one.Identity();
  Matrix3DParticleData pF(numPart, one);
  MPMVar pFVar(pF);
  add("pF", dwi, pFVar);

  // Create the interpolation information
  std::vector<int> zeroVecInt(numNearNodes, 0);
  VectorIntParticleData cIdx(numPart, zeroVecInt);
  MPMVar cIdxVar(cIdx);
  add("cIdx", dwi, cIdxVar);
  std::vector<double> zeroVecDouble(numNearNodes, 0.0);
  VectorDoubleParticleData cW(numPart, zeroVecDouble);
  MPMVar cWVar(cW);
  add("cW", dwi, cWVar);
  VectorDoubleParticleData cGradx(numPart, zeroVecDouble);
  MPMVar cGradxVar(cGradx);
  add("cGradx", dwi, cGradxVar);
  VectorDoubleParticleData cGrady(numPart, zeroVecDouble);
  MPMVar cGradyVar(cGrady);
  add("cGrady", dwi, cGradyVar);
  VectorDoubleParticleData cGradz(numPart, zeroVecDouble);
  MPMVar cGradzVar(cGradz);
  add("cGradz", dwi, cGradzVar);
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

