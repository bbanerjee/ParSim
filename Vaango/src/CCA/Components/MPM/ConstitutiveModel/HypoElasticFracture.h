/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Callaghan Innovation, New Zealand
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

#ifndef __HYPOELASTIC_FRACTURE_CONSTITUTIVE_MODEL_H__
#define __HYPOELASTIC_FRACTURE_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/HypoElastic.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {

class HypoElasticFracture : public HypoElastic
{
public:

  struct Toughness
  {
    double Vc, KIc, KIIc;

    Toughness(double Vc_in, double KIc_in, double KIIc_in) {
      Vc = Vc_in; KIc = KIc_in; KIIc = KIIc_in;
    }
  };

  struct FractureParams
  {
    // Parameters in the empirical criterion (KI/KIc)^p+(KII/KIIc)^q=1 
    // for crack initialization (KIIc=r*KIc)
    double p, q, r;

    // Fracture toughness at various velocities in the format Vector(Vc,KIc,KIIc)
    std::vector<Toughness> Kc;
  };

  // constructors
  HypoElasticFracture(ProblemSpecP& ps, MPMFlags* flag);
  HypoElasticFracture(const HypoElasticFracture* cm);
  HypoElasticFracture& operator=(const HypoElasticFracture& cm) = delete;
  virtual ~HypoElasticFracture() override = default;
  HypoElasticFracture* clone() override;

  ModelType modelType() const override
  {
    return ModelType::INCREMENTAL;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  /*  initialize */
  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  /* compute stress */
  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, const bool recursion,
                              const bool schedParent = true) const override;

  /* for material conversion */
  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;
  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  // Convert J-integral into stress intensity factors
  // for hypoelastic materials (for FRACTURE)
  void convertJToK(const MPMMaterial* matl, const std::string& stressState,
                   const Vector& J, const double& C, const Vector& V,
                   Vector& SIF) override;

  // Detect if crack propagates and the propagation direction (for FRACTURE)
  short crackPropagates(const double& Vc, const double& KI, const double& KII,
                        double& theta) override;

private:

  FractureParams d_fracParam;
  std::string d_crackPropagationCriterion;

  double crackPropagationAngleFromStrainEnergyDensityCriterion(const double&,
                                                               const double&,
                                                               const double&);

};

} // End namespace Uintah

#endif // __HYPOELASTIC_FRACTURE_CONSTITUTIVE_MODEL_H__
