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

#ifndef __CONSTITUTIVE_MODEL_H__
#define __CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <vector>

namespace Uintah {

class Task;
class Patch;
class VarLabel;
class MPMLabel;
class MPMFlags;
class MPMMaterial;
class DataWarehouse;
class ParticleSubset;
class ParticleVariableBase;

//////////////////////////////////////////////////////////////////////////
/*!
  \class ConstitutiveModel

  \brief Base class for contitutive models.

  \author Steven G. Parker \n
  Department of Computer Science \n
  University of Utah \n
  Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n

  Long description...
*/
//////////////////////////////////////////////////////////////////////////

// Needed for adding and deleting materials during particle conversion
// Defined in Ports/DataWarehouse.h
// typedef ParticleLabelVariableMap ParticleLabelVariableMap;

class ConstitutiveModel
{
public:
  enum class ModelType
  {
    TOTAL_FORM,
    RATE_FORM,
    INCREMENTAL
  };

  ConstitutiveModel(MPMFlags* MFlag);
  ConstitutiveModel(const ConstitutiveModel* cm);
  virtual ~ConstitutiveModel();

  void
  setMaterialManager(MaterialManagerP& manager)
  {
    d_materialManager = manager;
  }

  virtual ModelType
  modelType() const = 0;
  virtual void
  outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) = 0;

  // Basic constitutive model calculations
  virtual void
  computeStressTensor(const PatchSubset* patches,
                      const MPMMaterial* matl,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  virtual void
  computeStressTensorImplicit(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  ///////////////////////////////////////////////////////////////////////
  /*! Initial computes and requires for the constitutive model */
  ///////////////////////////////////////////////////////////////////////
  virtual void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches) const;

  ///////////////////////////////////////////////////////////////////////
  /*! Initialize the variables used in the CM */
  ///////////////////////////////////////////////////////////////////////
  virtual void
  initializeCMData(const Patch* patch,
                   const MPMMaterial* matl,
                   DataWarehouse* new_dw) = 0;

  ///////////////////////////////////////////////////////////////////////
  /*!
   * Actually initialize the stress and deformation gradient assuming linear
   * elastic behavior after computing the body force acceleration
   *
   * **WARNING** 1) Assumes zero shear stresses and that body forces are aligned
   *                with coordinate directions
   *             2) Needs the model to have a "initializeWithBodyForce" flag
   *                set as true.  A more general implementation is not worth
   *                the significant extra effort.
   */
  ///////////////////////////////////////////////////////////////////////
  virtual void
  initializeStressAndDefGradFromBodyForce(const Patch* patch,
                                          const MPMMaterial* matl,
                                          DataWarehouse* new_dw) const;

  ///////////////////////////////////////////////////////////////////////
  /*! Set up the computes and requires for the task that computes the
      stress tensor and associated kinematic and thermal quantities */
  ///////////////////////////////////////////////////////////////////////
  virtual void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches) const;

  virtual void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches,
                         const bool recursion,
                         const bool SchedParent) const;

  virtual void
  scheduleCheckNeedAddMPMMaterial(Task* task,
                                  const MPMMaterial* matl,
                                  const PatchSet* patches) const;

  // Determine if addition of an acceptor material is needed
  virtual void
  checkNeedAddMPMMaterial(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  /////////////////////////////////////////////////////////////////
  /*! Add particle conversion related requires to the task graph */
  /////////////////////////////////////////////////////////////////
  virtual void
  allocateCMDataAddRequires(Task* task,
                            const MPMMaterial* matl,
                            const PatchSet* patch,
                            MPMLabel* lb) const;

  /////////////////////////////////////////////////////////////////
  /*! Copy the data from the particle to be deleted to the particle
      to be added */
  /////////////////////////////////////////////////////////////////
  virtual void
  allocateCMDataAdd(DataWarehouse* new_dw,
                    ParticleSubset* addset,
                    ParticleLabelVariableMap* newState,
                    ParticleSubset* delset,
                    DataWarehouse* old_dw) = 0;

  virtual void
  addParticleState(std::vector<const VarLabel*>& from,
                   std::vector<const VarLabel*>& to) = 0;

  ////////////////////////////////////////////////////////////////////////
  /*! Carry forward CM variables for RigidMPM */
  ////////////////////////////////////////////////////////////////////////
  virtual void
  carryForward(const PatchSubset* patches,
               const MPMMaterial* matl,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  virtual double
  computeRhoMicroCM(double pressure,
                    const double p_ref,
                    const MPMMaterial* matl,
                    double temperature,
                    double rho_guess) = 0;

  virtual void
  computePressEOSCM(double rho_m,
                    double& press_eos,
                    double p_ref,
                    double& dp_drho,
                    double& ss_new,
                    const MPMMaterial* matl,
                    double temperature) = 0;

  virtual double
  getCompressibility() = 0;

  double
  computeRateOfWork(const Matrix3& stress,
                    const Matrix3& rateOfDeformation) const;

  virtual Vector
  getInitialFiberDir();

  double
  computeRhoMicro(double press,
                  double gamma,
                  double cv,
                  double Temp,
                  double rho_guess);

  void
  computePressEOS(double rhoM,
                  double gamma,
                  double cv,
                  double Temp,
                  double& press,
                  double& dp_drho,
                  double& dp_de);

  //////////
  // Convert J-integral into stress intensity for hypoelastic materials
  // (for FRACTURE)
  virtual void
  convertJToK(const MPMMaterial* matl,
              const string& stressState,
              const Vector& J,
              const double& C,
              const Vector& V,
              Vector& SIF);

  //////////
  // Detect if crack propagates and the direction (for FRACTURE)
  virtual short
  crackPropagates(const double& Vc,
                  const double& KI,
                  const double& KII,
                  double& theta);

  virtual void
  addRequiresDamageParameter(Task* task,
                             const MPMMaterial* matl,
                             const PatchSet* patches) const;

  virtual void
  getDamageParameter(const Patch* patch,
                     ParticleVariable<int>& damage,
                     int dwi,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  inline void
  setWorld(const ProcessorGroup* myworld)
  {
    d_world = myworld;
  }

  inline void
  setSharedState(MaterialManager* sharedState)
  {
    d_mat_manager = sharedState;
  }

  // Make a clone of the constitutive model

  virtual std::unique_ptr<ConstitutiveModel>
  clone() = 0;

protected:
  ///////////////////////////////////////////////////////////////////////
  /*! Initialize the common quantities that all the explicit constituive
   *  models compute : called by initializeCMData */
  ///////////////////////////////////////////////////////////////////////
  void
  initSharedDataForExplicit(const Patch* patch,
                            const MPMMaterial* matl,
                            DataWarehouse* new_dw);

  /////////////////////////////////////////////////////////////////
  /*! Computes and Requires common to all hypo-elastic constitutive models
   *  that do explicit time stepping : called by addComputesAndRequires */
  /////////////////////////////////////////////////////////////////
  void
  addSharedCRForHypoExplicit(Task* task,
                             const MaterialSubset* matlset,
                             const PatchSet* patches) const;

  /*!
   *  Computes and requires used by all models where stress has been
   *  pre-rotated.
   */
  void
  addComputesAndRequiresForRotatedExplicit(Task* task,
                                           const MaterialSubset* matlset,
                                           const PatchSet*) const;

  /////////////////////////////////////////////////////////////////
  /*! Computes and Requires common to all constitutive models that
   *  do explicit time stepping : called by addComputesAndRequires */
  /////////////////////////////////////////////////////////////////
  void
  addSharedCRForExplicit(Task* task,
                         const MaterialSubset* matlset,
                         const PatchSet* patches) const;

  /////////////////////////////////////////////////////////////////
  /*! Particle conversion related requires common to all constitutive
      models that do explicit time stepping : called by
      allocateCMDataAddRequires */
  /////////////////////////////////////////////////////////////////
  void
  addSharedRForConvertExplicit(Task* task,
                               const MaterialSubset* matlset,
                               const PatchSet*) const;

  /////////////////////////////////////////////////////////////////
  /*! Copy the data common to all constitutive models from the
      particle to be deleted to the particle to be added.
      Called by allocateCMDataAdd */
  /////////////////////////////////////////////////////////////////
  void
  copyDelToAddSetForConvertExplicit(DataWarehouse* new_dw,
                                    ParticleSubset* delset,
                                    ParticleSubset* addset,
                                    ParticleLabelVariableMap* newState);

  /////////////////////////////////////////////////////////////////
  /*! Carry forward the data common to all constitutive models
      when using RigidMPM.
      Called by carryForward */
  /////////////////////////////////////////////////////////////////
  void
  carryForwardSharedData(ParticleSubset* pset,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw,
                         const MPMMaterial* matl);

  /*!
    \brief Calculate the artificial bulk viscosity (q)

    \f[
    q = \rho (A_1 | c D_{kk} dx | + A_2 D_{kk}^2 dx^2)
    ~~\text{if}~~ D_{kk} < 0
    \f]
    \f[
    q = 0 ~~\text{if}~~ D_{kk} >= 0
    \f]

    where \f$ \rho \f$ = current density \n
    \f$ dx \f$ = characteristic length = (dx+dy+dz)/3 \n
    \f$ A_1 \f$ = Coeff1 (default = 0.2) \n
    \f$ A_2 \f$ = Coeff2 (default = 2.0) \n
    \f$ c \f$ = Local bulk sound speed = \f$ \sqrt{K/\rho} \f$ \n
    \f$ D_{kk} \f$ = Trace of rate of deformation tensor \n
  */
  double
  artificialBulkViscosity(double Dkk, double c, double rho, double dx) const;

  /*! Calculate gradient of vector field for 8 noded interpolation, B matrix
      for Kmat and B matrix for Kgeo */
  void
  computeGradAndBmats(Matrix3& grad,
                      std::vector<IntVector>& ni,
                      std::vector<Vector>& d_S,
                      const double* oodx,
                      constNCVariable<Vector>& gVec,
                      const Array3<int>& l2g,
                      double B[6][24],
                      double Bnl[3][24],
                      int* dof);

  MPMLabel* lb;
  MPMFlags* flag;
  ModelType d_modelType;

  int NGP;
  int NGN;
  const ProcessorGroup* d_world;

  // don't store MaterialManagerP or it will add a reference
  // that will never be removed
  MaterialManager* d_mat_manager;

  MaterialManagerP d_materialManager;
};

std::ostream&
operator<<(std::ostream& out, const ConstitutiveModel::ModelType& mt);

} // End namespace Uintah

#endif // __CONSTITUTIVE_MODEL_H__
