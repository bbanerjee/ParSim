/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/SimulationState.h>
#include <vector>
#include <map>

namespace Uintah {
  typedef int particleIndex;
  typedef int particleId;

  class GeometryObject;
  class Patch;
  class DataWarehouse;
  class MPMFlags;
  class MPMMaterial;
  class MPMLabel;
  class ParticleSubset;
  class VarLabel;

  class ParticleCreator {

  public:
    
    ParticleCreator(MPMMaterial* matl, MPMFlags* flags);


    virtual ~ParticleCreator();


    virtual particleIndex createParticles(MPMMaterial* matl,
                                          CCVariable<short int>& cellNAPID,
                                          const Patch*,DataWarehouse* new_dw,
                                          std::vector<GeometryObject*>&);


    virtual void registerPermanentParticleState(MPMMaterial* matl);

    std::vector<const VarLabel* > returnParticleState();
    std::vector<const VarLabel* > returnParticleStatePreReloc();

    typedef std::map<GeometryObject*, std::vector<Point> > geomPoints;
    typedef std::map<GeometryObject*, std::vector<double> > geomVols;
    typedef std::map<GeometryObject*, std::vector<Vector> > geomVecs;
    typedef std::map<GeometryObject*, std::vector<Matrix3> > geomMat3s;

    typedef struct {
      geomPoints d_object_points;
      geomVols d_object_vols;
      geomVols d_object_temps;
      geomVols d_object_colors;
      geomVecs d_object_forces;
      geomVecs d_object_fibers;  
      geomVecs d_object_velocity; // gcd add
      geomMat3s d_object_size;
    } ObjectVars;

    typedef struct {
      ParticleVariable<Point> position;
      ParticleVariable<Vector> pDisp, pVelocity, pAcc, pExternalForce;
      ParticleVariable<Matrix3> pSize;
      ParticleVariable<double> pMass, pVolume, pTemperature, 
                               pSpecificVolume, pErosion;
      ParticleVariable<double> pColor,pTempPrevious,p_q;
      ParticleVariable<long64> pParticleID;
      ParticleVariable<Vector> pFiberDir; 
      ParticleVariable<int> pLoadCurveID;
      // Body forces
      ParticleVariable<Vector> pBodyForceAcc;
      ParticleVariable<double> pCoriolisImportance;
      // ImplicitParticleCreator
      ParticleVariable<double> pVolumeold;
      //MembraneParticleCreator
      ParticleVariable<Vector> pTang1, pTang2, pNorm;
      // AMR
      ParticleVariable<int> pRefined;
      ParticleVariable<int> pLastLevel;

      // Switch between explicit and implicit MPM
      ParticleVariable<double> pExternalHeatFlux;

      // For friction contact
      ParticleVariable<double> pSurface;

    } ParticleVars;

  protected:

    virtual ParticleSubset* allocateVariables(particleIndex numParticles,
                                              int dwi, const Patch* patch,
                                              DataWarehouse* new_dw,
                                              ParticleVars& pvars);

    virtual particleIndex countAndCreateParticles(const Patch*,
                                                  GeometryObject* obj,
                                                  ObjectVars& vars);


    void createPoints(const Patch* patch, GeometryObject* obj, ObjectVars& vars);



    virtual void initializeParticle(const Patch* patch,
                                    GeometryObject* obj,
                                    MPMMaterial* matl,
                                    Point p, IntVector cell_idx,
                                    particleIndex i,
                                    CCVariable<short int>& cellNAPI,
                                    ParticleVars& pvars);
    
    //////////////////////////////////////////////////////////////////////////
    /*! Get the LoadCurveID applicable for this material point */
    //////////////////////////////////////////////////////////////////////////
    int getLoadCurveID(const Point& pp, const Vector& dxpp);

    //////////////////////////////////////////////////////////////////////////
    /*! Print MPM physical boundary condition information */
    //////////////////////////////////////////////////////////////////////////
    void printPhysicalBCs();

    //////////////////////////////////////////////////////////////////////////
    /*! Calculate the external force to be applied to a particle */
    //////////////////////////////////////////////////////////////////////////
    virtual void applyForceBC(const Vector& dxpp,  const Point& pp,
                              const double& pMass,  Vector& pExtForce);
    
    int checkForSurface(const GeometryPieceP piece, const Point p,
                        const Vector dxpp);
    double checkForSurface2(const GeometryPieceP piece, const Point p,
                            const Vector dxpp);
    


    MPMLabel* d_lb;
    MPMFlags* d_flags;

    bool d_useLoadCurves;
    bool d_withColor;
    bool d_doScalarDiffusion;
    bool d_artificialViscosity;
    bool d_computeScaleFactor;
    bool d_useCPTI;

    std::vector<const VarLabel* > particle_state, particle_state_preReloc;
    
    mutable CrowdMonitor   d_lock;

  public:

    /*! For material addition capability */
    virtual void allocateVariablesAddRequires(Task* task, 
                                              const MPMMaterial* matl,
                                              const PatchSet* patch) const;

    /*! For material addition capability */
    virtual void allocateVariablesAdd(DataWarehouse* new_dw,
                                      ParticleSubset* addset,
                                      map<const VarLabel*,ParticleVariableBase*>* newState,
                                      ParticleSubset* delset,
                                      DataWarehouse* old_dw);

  };



} // End of namespace Uintah

#endif // __PARTICLE_CREATOR_H__
