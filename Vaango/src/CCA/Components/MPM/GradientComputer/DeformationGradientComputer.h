#ifndef __DEFORMATION_GRADIENT_COMPUTER_H__
#define __DEFORMATION_GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/ParticleVariableBase.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <vector>

namespace Uintah {

  class Task;
  class Patch;

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class DeformationGradientComputer
    \brief Class for computing deformation gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class DeformationGradientComputer {


  public:
         
    DeformationGradientComputer(MPMFlags* MFlag, SimulationStateP& ss);
    DeformationGradientComputer(const DeformationGradientComputer* gc);
    virtual ~DeformationGradientComputer();

    // Make a clone of the gradient computer
    DeformationGradientComputer* clone();

    // Computes and requires     
    void addInitialComputesAndRequires(Task* task,
                                       const MPMMaterial* mpm_matl,
                                       const PatchSet*);

    void addComputesAndRequires(Task* task,
                                const MPMMaterial* mpm_matl,
                                const PatchSet*);

    void addComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches,
                                const bool /*recurse*/,
                                const bool SchedParent) const;

    void addRequiresForConvert(Task* task,
                               const MPMMaterial* mpm_matl);

    void copyAndDeleteForConvert(DataWarehouse* new_dw,
                                 ParticleSubset* addset,
                                 map<const VarLabel*,
                                 ParticleVariableBase*>* newState,
                                 ParticleSubset* delset,
                                 DataWarehouse* old_dw );

    void initializeGradient(const Patch* patch,
                            const MPMMaterial* mpm_matl,
                            DataWarehouse* new_dw);

    void computeDeformationGradient(const PatchSubset* patches,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw);

    void computeDeformationGradient(const PatchSubset* patches,
                                    const MPMMaterial* mpm_matl,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw,
                                    const bool recurse);

  protected:

    void addComputesAndRequiresExplicit(Task* task,
                                        const MPMMaterial* mpm_matl);
   
    void addComputesAndRequiresImplicit(Task* task,
                                        const MPMMaterial* mpm_matl);

    void initializeGradientExplicit(const Patch* patch,
                                    const MPMMaterial* mpm_matl,
                                    DataWarehouse* new_dw);

    void initializeGradientImplicit(const Patch* patch,
                                    const MPMMaterial* mpm_matl,
                                    DataWarehouse* new_dw);

    void computeDeformationGradientExplicit(const Patch* patch,
                                            const MPMMaterial* mpm_matl,
                                            const double& delT,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw);

    void computeDeformationGradientImplicit(const Patch* patch,
                                            const MPMMaterial* mpm_matl,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw);

    void computeDeformationGradientImplicit(const Patch* patch,
                                            const MPMMaterial* mpm_matl,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* parent_old_dw,
                                            DataWarehouse* new_dw);

    void computeDeformationGradientFromVelocity(const Matrix3& velGrad_old,
                                                const Matrix3& velGrad_new,
                                                const Matrix3& defGrad_old,
                                                const double& delT,
                                                Matrix3& defGrad_new,
                                                Matrix3& defGrad_inc);

    void computeDeformationGradientFromTotalDisplacement(const Matrix3& dispGrad_new,
                                                         const Matrix3& defGrad_old,
                                                         Matrix3& defGrad_new,
                                                         Matrix3& defGrad_inc);

    void seriesUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                     const Matrix3& defGrad_old,
                                     const double& delT,
                                     Matrix3& defGrad_new,
                                     Matrix3& defGrad_inc);

    void seriesUpdateLinearVelGrad(const Matrix3& velGrad_old,
                                   const Matrix3& velGrad_new,
                                   const Matrix3& defGrad_old,
                                   const double& delT,
                                   Matrix3& defGrad_new,
                                   Matrix3& defGrad_inc);

    void subcycleUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                       const Matrix3& defGrad_old,
                                       const double& delT,
                                       Matrix3& defGrad_new,
                                       Matrix3& defGrad_inc);

    void computeDeformationGradientFromIncrementalDisplacement(const Matrix3& dispGrad_new,
                                                               const Matrix3& defGrad_old,
                                                               Matrix3& defGrad_new,
                                                               Matrix3& defGrad_inc);

  protected:

    MPMLabel* lb;
    MPMFlags* flag;
    int NGP;
    int NGN;
    SimulationStateP d_sharedState;

    static const Matrix3 Identity;
    static const Matrix3 Zero;
  };

} // End namespace Uintah
      


#endif  // __DEFORMATION_GRADIENT_COMPUTER_H__

