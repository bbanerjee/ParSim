/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/ParticleCreator/ParticleCreatorStructs.h>
#include <CCA/Components/MPM/ParticleCreator/SGPDataManager.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/GeometryPiece/SpecialGeomPiece.h>

#include <map>
#include <optional>
#include <string_view>
#include <vector>


namespace Uintah {

class ParticleCreator
{
public:
  ParticleCreator(MPMMaterial* matl, MPMFlags* flags);

  virtual ~ParticleCreator() = default;


  virtual auto
  createParticles(MPMMaterial* matl,
                  CCVariable<short int>& cellNAPID,
                  const Patch*,
                  DataWarehouse* new_dw,
                  const VecGeometryObjectSP& objects) -> particleIndex;

  virtual void
  registerPermanentParticleState(MPMMaterial* matl);

  auto
  returnParticleState() -> std::vector<const VarLabel*>;

  auto
  returnParticleStatePreReloc() -> std::vector<const VarLabel*>;


protected:
  virtual auto
  allocateVariables(particleIndex numParticles,
                    int dwi,
                    const Patch* patch,
                    DataWarehouse* new_dw,
                    ParticleVars& pvars) -> ParticleSubset*;

  virtual auto
  countAndCreateParticles(const Patch*, GeometryObject* obj, ObjectVars& vars)
    -> particleIndex;

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
  auto
  getLoadCurveID(const Point& pp, const Vector& dxpp) -> int;

  auto
  getLoadCurveID(const Point& pp,
                 const Vector& dxpp,
                 Vector& areacomps,
                 int mat_id) -> IntVector;

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

  auto
  checkForSurface(const GeometryPieceP piece, const Point p, const Vector dxpp)
    -> int;

  auto
  checkForSurface2(const GeometryPieceP piece, const Point p, const Vector dxpp)
    -> double;

  MPMMaterial* d_matl{ nullptr };
  std::unique_ptr<MPMLabel> d_mpm_labels{ nullptr };
  std::unique_ptr<AMRMPMLabel> d_amrmpm_labels{ nullptr };
  std::unique_ptr<HydroMPMLabel> d_hydrompm_labels{ nullptr };
  MPMFlags* d_flags{ nullptr };

  bool d_useLoadCurves{ false };
  bool d_useLoadCurvesVector{ false };
  bool d_withColor{ false };
  bool d_doScalarDiffusion{ false };
  bool d_artificialViscosity{ false };
  bool d_computeScaleFactor{ false };
  bool d_useCPTI{ false };
  bool d_withGaussSolver{ false };
  bool d_coupledFlow{ false };

  std::vector<const VarLabel*> particle_state, particle_state_preReloc;

private:

  /**
   * Handles particle creation for discrete materials.
   * Creates one master particle at center plus geometric proxy particles.
   */
  particleIndex handleDEMMaterial(const Patch* patch,
                                      GeometryObject* obj,
                                      ObjectVars& obj_vars);

  /**
   * Handles particle creation for rigid materials.
   * Uses normal discretization plus adds a master particle at center.
   */
  particleIndex handleRigidMaterial(const Patch* patch,
                                   GeometryObject* obj,
                                   ObjectVars& obj_vars);

  /**
   * Handles particle creation for special geometry pieces.
   * Uses the particle creators from the SpecialGeomPiece class.
   */
  particleIndex handleSpecialGeometryPiece(const Patch* patch,
                                          GeometryObject* obj,
                                          ObjectVars& obj_vars);

  /**
   * Creates a master particle at the center of the geometry's bounding box.
   */
  void createMasterParticle(const Patch* patch,
                           GeometryPieceP piece,
                           GeometryObject* obj,
                           ObjectVars& obj_vars);

  /**
   * Creates surface visualization particles based on geometry type.
   * Handles cylinders, spheres, and provides fallback for other geometries.
   */
  void createSurfaceParticles(const Patch* patch,
                             GeometryPieceP piece,
                             GeometryObject* obj,
                             ObjectVars& obj_vars);

  /**
   * Creates corner particles at bounding box corners.
   * These ensure the object is visible to neighbor patches.
   */
  void createCornerParticles(const Patch* patch,
                            GeometryPieceP piece,
                            GeometryObject* obj,
                            ObjectVars& obj_vars);

  /**
   * Populates particle variables from special geometry piece data.
   */
  void populateParticleVariables(SpecialGeomPiece* sgp,
                                const Patch* patch,
                                GeometryObject* obj,
                                ObjectVars& obj_vars,
                                const ParticleVarPairs& pairs,
                                int numPts);
  /**
   * Initializes a single particle with position, cell, and basic properties.
   */
  void initializeParticleBasics(const Patch* patch,
                               GeometryObject* obj,
                               MPMMaterial* matl,
                               const Point& point,
                               const IntVector& cell_idx,
                               particleIndex pidx,
                               CCVariable<short int>& cellNAPID,
                               ParticleVars& pvars);
  
  /**
   * Applies special geometry piece data to a particle.
   */
  void applySpecialGeometryData(Vaango::Helpers::SGPDataManager& sgp_data,
                               MPMMaterial* matl,
                               const Patch* patch,
                               particleIndex pidx,
                               ParticleVars& pvars);
  
  /**
   * Applies load curve if particle is on surface.
   */
  void applyLoadCurve(GeometryPieceP piece,
                     const Point& point,
                     const Vector& dxpp,
                     particleIndex pidx,
                     ParticleVars& pvars);
  
  /**
   * Processes and initializes all particles for a single geometry object.
   */
  void processGeometryObjectParticles(GeometryObject* obj,
                                     GeometryPieceP piece,
                                     const Patch* patch,
                                     MPMMaterial* matl,
                                     const ObjectVars& obj_vars,
                                     const Vector& dxpp,
                                     particleIndex& current_particle_index,
                                     CCVariable<short int>& cellNAPID,
                                     ParticleVars& pvars,
                                     Vaango::Helpers::SGPDataManager& sgp_data);
};

} // End of namespace Uintah

#endif // ____CCA_COMPONENTS_MPM_PARTICLE_CREATOR_PARTICLECREATOR_H____
