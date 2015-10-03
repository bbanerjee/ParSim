/*
 * The MIT License
 *
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __VISCOELASTIC_FORTRAN_CONSTITUTIVE_MODEL_H__
#define __VISCOELASTIC_FORTRAN_CONSTITUTIVE_MODEL_H__

#include <cmath>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <vector>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class ViscoElasticFortran
    \brief Linear viscoelastic relaxation model that calls Tim Fuller's FORTRAN90 code
    \author Biswajit Banerjee \n

    \warning Only isotropic materials
  */
  /////////////////////////////////////////////////////////////////////////////

  class ViscoElasticFortran : public Uintah::ConstitutiveModel {

  public:

    struct ModelParameters {
      double K;     // Bulk modulus
      double G;     // Shear modulus = G_\infty + sum_{i=1}^10 G_i
      double G00;   // Shear modulus G_\infty/G
      double G01;   // Normalized Prony series coefficients (G_i/G)
      double G02;
      double G03;
      double G04;
      double G05;
      double G06;
      double G07;
      double G08;
      double G09;
      double G10;
      double Tau01;   // Prony series shear relaxation times
      double Tau02;
      double Tau03;
      double Tau04;
      double Tau05;
      double Tau06;
      double Tau07;
      double Tau08;
      double Tau09;
      double Tau10;
      double C1_WLF;   // WLF time-temperature superposition parameter C1
      double C2_WLF;   // WLF time-temperature superposition parameter C1
      double Tref_WLF; // WLF time-temperature superposition reference temperature
    };

  private:

    ModelParameters d_param;
    int d_nProp;
    std::vector<double> d_props;

    int d_nStateV;
    std::vector<const Uintah::VarLabel*> stateVLabels;
    std::vector<const Uintah::VarLabel*> stateVLabels_preReloc;

    // Prevent copying of this class
    ViscoElasticFortran& operator=(const ViscoElasticFortran &cm);

  public:
    // constructors
    ViscoElasticFortran(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
    ViscoElasticFortran(const ViscoElasticFortran* cm);
       
    // destructor
    virtual ~ViscoElasticFortran();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps,bool output_cm_tag = true);

    // clone
    ViscoElasticFortran* clone();

    // Initialize local labels
    void initializeLocalMPMLabels();

    // compute stable timestep for this patch
    virtual void computeStableTimestep(const Uintah::Patch* patch,
                                       const Uintah::MPMMaterial* matl,
                                       Uintah::DataWarehouse* new_dw);

    // compute stress at each particle in the patch
    virtual void computeStressTensor(const Uintah::PatchSubset* patches,
                                     const Uintah::MPMMaterial* matl,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw);


    // carry forward CM data for RigidMPM
    virtual void carryForward(const Uintah::PatchSubset* patches,
                              const Uintah::MPMMaterial* matl,
                              Uintah::DataWarehouse* old_dw,
                              Uintah::DataWarehouse* new_dw);

    // initialize  each particle's constitutive model data
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const Uintah::MPMMaterial* matl,
                                       const Uintah::PatchSet* ) const;
    
    virtual void initializeCMData(const Uintah::Patch* patch,
                                  const Uintah::MPMMaterial* matl,
                                  Uintah::DataWarehouse* new_dw);

    virtual void allocateCMDataAddRequires(Uintah::Task* task, const Uintah::MPMMaterial* matl,
                                           const Uintah::PatchSet* patch, 
                                           Uintah::MPMLabel* lb) const;


    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* subset,
                                   SCIRun::ParticleLabelVariableMap* newState,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw);


    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches,
                                        const bool recursion,
                                        const bool schedParent=true) const;

    virtual double computeRhoMicroCM(double pressure,
                                     const double p_ref,
                                     const Uintah::MPMMaterial* matl, 
                                     double temperature,
                                     double rho_guess);

    virtual void computePressEOSCM(double rho_m, double& press_eos,
                                   double p_ref,
                                   double& dp_drho, double& ss_new,
                                   const Uintah::MPMMaterial* matl, 
                                   double temperature);

    virtual double getCompressibility();


    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to);
  };

} // End namespace Uintah

#endif  // __VISCOELASTIC_FORTRAN_CONSTITUTIVE_MODEL_H__ 
