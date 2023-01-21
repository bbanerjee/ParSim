/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef ____CCA_COMPONENTS_MPM_PARTICLE_CREATOR_PARTICLECREATOR_H____
#define ____CCA_COMPONENTS_MPM_PARTICLE_CREATOR_PARTICLECREATOR_H____

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <map>
#include <vector>

namespace Uintah {

using particleIndex = int;
using particleId    = int;

class GeometryObject;
class Patch;
class DataWarehouse;
class MPMFlags;
class MPMMaterial;
class MPMLabel;
class AMRMPMLabel;
class HydroMPMLabel;
class ParticleSubset;
class VarLabel;

class ParticleCreator
{
public:
  ParticleCreator(MPMMaterial* matl, MPMFlags* flags);

  virtual ~ParticleCreator() = default;

  using VecGeometryObjectSP = std::vector<std::shared_ptr<GeometryObject>>;

  virtual particleIndex
  createParticles(MPMMaterial* matl,
                  CCVariable<short int>& cellNAPID,
                  const Patch*,
                  DataWarehouse* new_dw,
                  const VecGeometryObjectSP& objects);

  virtual void
  registerPermanentParticleState(MPMMaterial* matl);

  std::vector<const VarLabel*>
  returnParticleState();

  std::vector<const VarLabel*>
  returnParticleStatePreReloc();

  using GeomName   = std::pair<std::string, GeometryObject*>;
  using GeomPoint  = std::map<GeometryObject*, std::vector<Point>>;
  using GeomScalar = std::map<GeomName, std::vector<double>>;
  using GeomVector = std::map<GeomName, std::vector<Vector>>;
  using GeomTensor = std::map<GeomName, std::vector<Matrix3>>;

  struct ObjectVars
  {
    GeomPoint points;
    GeomScalar scalars;
    GeomVector vectors;
    GeomTensor tensors;
  };

  struct ParticleVars
  {
    ParticleVariable<Point> position;
    ParticleVariable<Vector> pDisp, pVelocity, pAcc, pExternalForce;
    ParticleVariable<Matrix3> pSize;
    ParticleVariable<double> pMass, pVolume, pTemperature, pSpecificVolume,
      pErosion;
    ParticleVariable<double> pColor, pTempPrevious, p_q;
    ParticleVariable<long64> pParticleID;
    ParticleVariable<Vector> pFiberDir;
    ParticleVariable<int> pLoadCurveID;
    ParticleVariable<IntVector> pLoadCurveIDVector;

    // Body forces
    ParticleVariable<Vector> pBodyForceAcc;
    ParticleVariable<double> pCoriolisImportance;

    // ImplicitParticleCreator
    ParticleVariable<double> pVolumeold;

    // MembraneParticleCreator
    ParticleVariable<Vector> pTang1, pTang2, pNorm;

    // AMR
    ParticleVariable<int> pRefined;
    ParticleVariable<int> pLastLevel;

    // Switch between explicit and implicit MPM
    ParticleVariable<double> pExternalHeatFlux;

    // For friction contact
    ParticleVariable<double> pSurface;
    ParticleVariable<double> pSurfaceGrad;

    // Scalar Diffusion
    ParticleVariable<double> pConcentration;
    ParticleVariable<double> pConcPrevious;
    ParticleVariable<Vector> pConcGrad;
    ParticleVariable<double> pExternalScalarFlux;
    ParticleVariable<double> pPosCharge;
    ParticleVariable<double> pNegCharge;
    ParticleVariable<Vector> pPosChargeGrad;
    ParticleVariable<Vector> pNegChargeGrad;
    ParticleVariable<double> pPermittivity;
    ParticleVariable<Vector> pArea;

    // Hydro-mechanical coupling MPM
    ParticleVariable<double> pFluidMass, pSolidMass, pPorePressure, pPorosity;
    ParticleVariable<Vector> pFluidVelocity, pFluidAcceleration;
    ParticleVariable<Vector> pPrescribedPorePressure;
  };

protected:
  virtual ParticleSubset*
  allocateVariables(particleIndex numParticles,
                    int dwi,
                    const Patch* patch,
                    DataWarehouse* new_dw,
                    ParticleVars& pvars);

  virtual particleIndex
  countAndCreateParticles(const Patch*, GeometryObject* obj, ObjectVars& vars);

  void
  createPoints(const Patch* patch, GeometryObject* obj, ObjectVars& vars);

  virtual void
  initializeParticle(const Patch* patch,
                     GeometryObject* obj,
                     MPMMaterial* matl,
                     Point p,
                     IntVector cell_idx,
                     particleIndex i,
                     CCVariable<short int>& cellNAPI,
                     ParticleVars& pvars);

  //////////////////////////////////////////////////////////////////////////
  /*! Get the LoadCurveID applicable for this material point */
  //////////////////////////////////////////////////////////////////////////
  int
  getLoadCurveID(const Point& pp, const Vector& dxpp);

  IntVector
  getLoadCurveID(const Point& pp,
                 const Vector& dxpp,
                 Vector& areacomps,
                 int mat_id);

  //////////////////////////////////////////////////////////////////////////
  /*! Print MPM physical boundary condition information */
  //////////////////////////////////////////////////////////////////////////
  void
  printPhysicalBCs();

  //////////////////////////////////////////////////////////////////////////
  /*! Calculate the external force to be applied to a particle */
  //////////////////////////////////////////////////////////////////////////
  virtual void
  applyForceBC(const Vector& dxpp,
               const Point& pp,
               const double& pMass,
               Vector& pExtForce);

  int
  checkForSurface(const GeometryPieceP piece, const Point p, const Vector dxpp);

  double
  checkForSurface2(const GeometryPieceP piece,
                   const Point p,
                   const Vector dxpp);

  std::unique_ptr<MPMLabel> d_mpm_labels;
  std::unique_ptr<AMRMPMLabel> d_amrmpm_labels;
  std::unique_ptr<HydroMPMLabel> d_hydrompm_labels;
  MPMFlags* d_flags;

  bool d_useLoadCurves;
  bool d_useLoadCurvesVector;
  bool d_withColor;
  bool d_doScalarDiffusion;
  bool d_artificialViscosity;
  bool d_computeScaleFactor;
  bool d_useCPTI;
  bool d_withGaussSolver;
  bool d_coupledFlow;

  std::vector<const VarLabel*> particle_state, particle_state_preReloc;
};

} // End of namespace Uintah

#endif // ____CCA_COMPONENTS_MPM_PARTICLE_CREATOR_PARTICLECREATOR_H____
