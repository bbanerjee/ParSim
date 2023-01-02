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

#ifndef __PARTICLE_CREATOR_H__
#define __PARTICLE_CREATOR_H__

#include <Core/Thread/CrowdMonitor.h>

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/SimulationStateP.h>
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
class ParticleSubset;
class VarLabel;

class ParticleCreator
{
public:
  ParticleCreator(MPMMaterial* matl, MPMFlags* flags);

  virtual ~ParticleCreator() = default;

  virtual particleIndex
  createParticles(MPMMaterial* matl,
                  CCVariable<short int>& cellNAPID,
                  const Patch*,
                  DataWarehouse* new_dw,
                  const std::vector<std::shared_ptr<GeometryObject>>& objects);

  virtual void
  registerPermanentParticleState(MPMMaterial* matl);

  std::vector<const VarLabel*>
  returnParticleState();
  std::vector<const VarLabel*>
  returnParticleStatePreReloc();

  using GeomPoint   = std::map<GeometryObject*, std::vector<Point>>;
  using GeomReal    = std::map<GeometryObject*, std::vector<double>>;
  using GeomVector  = std::map<GeometryObject*, std::vector<Vector>>;
  using GeomMatrix3 = std::map<GeometryObject*, std::vector<Matrix3>>;

  struct ObjectVars
  {
    GeomPoint d_object_points;
    GeomReal d_object_vols;
    GeomReal d_object_temps;
    GeomReal d_object_colors;
    GeomVector d_object_forces;
    GeomVector d_object_fibers;
    GeomVector d_object_velocity; // gcd add
    GeomMatrix3 d_object_size;
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

  std::unique_ptr<MPMLabel> d_lb;
  MPMFlags* d_flags;

  bool d_useLoadCurves;
  bool d_withColor;
  bool d_doScalarDiffusion;
  bool d_artificialViscosity;
  bool d_computeScaleFactor;
  bool d_useCPTI;

  std::vector<const VarLabel*> particle_state, particle_state_preReloc;

  mutable CrowdMonitor d_lock;

public:
  /*! For material addition capability */
  virtual void
  allocateVariablesAddRequires(Task* task,
                               const MPMMaterial* matl,
                               const PatchSet* patch) const;

  /*! For material addition capability */
  virtual void
  allocateVariablesAdd(DataWarehouse* new_dw,
                       ParticleSubset* addset,
                       map<const VarLabel*, ParticleVariableBase*>* newState,
                       ParticleSubset* delset,
                       DataWarehouse* old_dw);
};

} // End of namespace Uintah

#endif // __PARTICLE_CREATOR_H__
