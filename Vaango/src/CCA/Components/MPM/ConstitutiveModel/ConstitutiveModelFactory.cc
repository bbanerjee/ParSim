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

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModelFactory.h>

#include <sci_defs/uintah_defs.h> // For NO_FORTRAN

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/CompMooneyRivlin.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFracture.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticImplicit.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/PolarOrthotropicHypoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyper.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyperImplicit.h>

#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/JWLppMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/MurnaghanMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/ProgramBurn.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/ViscoSCRAMHotSpot.h>

#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/ElasticPlasticHP.h>
#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/UCNH.h>

#include <CCA/Components/MPM/ConstitutiveModel/PorousModels/P_Alpha.h>
#include <CCA/Components/MPM/ConstitutiveModel/PorousModels/SoilFoam.h>

#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/ArenaMixture.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Arenisca.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Arenisca3.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Arenisca4.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/CamClay.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/MohrCoulomb.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/NonLocalDruckerPrager.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/SoilModelBrannon.h>

#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/IdealGasMP.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/Membrane.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/RigidMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/ShellMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/Water.h>

#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularPlasticity.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularPlasticityCap.h>

#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/MWViscoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoScram.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoScramImplicit.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoTransIsoHyper.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoTransIsoHyperImplicit.h>

#include <CCA/Components/MPM/ConstitutiveModel/ViscoPlasticModels/ViscoPlastic.h>

#include <CCA/Components/MPM/ConstitutiveModel/ManufacturedSolutions/CNH_MMS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ManufacturedSolutions/HypoElastic_MMS.h>

#if !defined(NO_FORTRAN)
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFortran.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Diamm.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Kayenta.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoElasticFortran.h>
#endif

#include <CCA/Components/MPM/Core/MPMFlags.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <fstream>
#include <iostream>
#include <string>

namespace Uintah {

std::unique_ptr<ConstitutiveModel>
ConstitutiveModelFactory::create(ProblemSpecP& ps, MPMFlags* flags)
{
  ProblemSpecP child = ps->findBlock("constitutive_model");
  if (!child)
    throw ProblemSetupException(
      "Cannot find constitutive_model tag", __FILE__, __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "No type for constitutive_model", __FILE__, __LINE__);

  if (flags->d_integratorType != "implicit" &&
      flags->d_integratorType != "explicit" &&
      flags->d_integratorType != "fracture") {
    string txt = "MPM: time integrator [explicit or implicit] hasn't been set.";
    throw ProblemSetupException(txt, __FILE__, __LINE__);
  }

  if (flags->d_integratorType == "implicit" &&
      (mat_type == "comp_neo_hook_plastic")) {
    string txt = "MPM:  You cannot use implicit MPM and comp_neo_hook_plastic";
    throw ProblemSetupException(txt, __FILE__, __LINE__);
  }

  if (mat_type == "rigid")
    return std::make_unique<RigidMaterial>(child, flags);

  else if (mat_type == "comp_mooney_rivlin")
    return std::make_unique<CompMooneyRivlin>(child, flags);

  else if (mat_type == "nonlocal_drucker_prager")
    return std::make_unique<NonLocalDruckerPrager>(child, flags);

  else if (mat_type == "Arenisca")
    return std::make_unique<Arenisca>(child, flags);

  else if (mat_type == "Arenisca3")
    return std::make_unique<Arenisca3>(child, flags);

  else if (mat_type == "arena")
    return std::make_unique<Vaango::Arena>(child, flags);

  else if (mat_type == "arena_mixture")
    return std::make_unique<Vaango::ArenaMixture>(child, flags);

  else if (mat_type == "Arenisca4")
    return std::make_unique<Arenisca4>(child, flags);

  else if (mat_type == "soil_model_brannon")
    return std::make_unique<SoilModelBrannon>(child, flags);

  else if (mat_type == "tabular_eos")
    return std::make_unique<Vaango::TabularEquationOfState>(child, flags);

  else if (mat_type == "tabular_plasticity")
    return std::make_unique<Vaango::TabularPlasticity>(child, flags);

  else if (mat_type == "tabular_plasticity_cap")
    return std::make_unique<Vaango::TabularPlasticityCap>(child, flags);

  else if (mat_type == "mohr_coulomb")
    return std::make_unique<MohrCoulomb>(child, flags);

  else if (mat_type == "comp_neo_hook") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return std::make_unique<UCNH>(child, flags, false, false);
    else if (flags->d_integratorType == "implicit")
      return std::make_unique<UCNH>(child, flags);
  } else if (mat_type == "cnh_damage")
    return std::make_unique<UCNH>(child, flags, false, true);

  else if (mat_type == "UCNH")
    return std::make_unique<UCNH>(child, flags);

  else if (mat_type == "cnh_mms")
    return std::make_unique<CNH_MMS>(child, flags);

  else if (mat_type == "cnhp_damage")
    return std::make_unique<UCNH>(child, flags, true, true);

  else if (mat_type == "trans_iso_hyper") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return std::make_unique<TransIsoHyper>(child, flags);
    else if (flags->d_integratorType == "implicit")
      return std::make_unique<TransIsoHyperImplicit>(child, flags);
  }

  else if (mat_type == "visco_trans_iso_hyper") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return std::make_unique<ViscoTransIsoHyper>(child, flags);
    else if (flags->d_integratorType == "implicit")
      return std::make_unique<ViscoTransIsoHyperImplicit>(child, flags);
  }

  else if (mat_type == "ideal_gas")
    return std::make_unique<IdealGasMP>(child, flags);

  else if (mat_type == "p_alpha")
    return std::make_unique<P_Alpha>(child, flags);

  else if (mat_type == "water")
    return std::make_unique<Water>(child, flags);

  else if (mat_type == "comp_neo_hook_plastic")
    return std::make_unique<UCNH>(child, flags, true, false);

  else if (mat_type == "visco_scram") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return std::make_unique<ViscoScram>(child, flags);
    else if (flags->d_integratorType == "implicit")
      return std::make_unique<ViscoScramImplicit>(child, flags);
  }

  else if (mat_type == "viscoSCRAM_hs")
    return std::make_unique<ViscoSCRAMHotSpot>(child, flags);

  else if (mat_type == "hypo_elastic") {
    if (flags->d_integratorType == "explicit")
      return std::make_unique<HypoElastic>(child, flags);

    else if (flags->d_integratorType == "implicit") {
      if (!flags->d_doGridReset) {
         std::ostringstream msg;
        msg << "\n ERROR: One may not use HypoElastic along with \n"
            << " <do_grid_reset>false</do_grid_reset> \n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
      return std::make_unique<HypoElasticImplicit>(child, flags);
    }
  }

  else if (mat_type == "hypo_elastic_fracture")
    return std::make_unique<HypoElasticFracture>(child, flags);

#if !defined(NO_FORTRAN)
  else if (mat_type == "hypo_elastic_fortran")
    return std::make_unique<HypoElasticFortran>(child, flags);

  else if (mat_type == "kayenta")
    return std::make_unique<Kayenta>(child, flags);

  else if (mat_type == "diamm")
    return std::make_unique<Diamm>(child, flags);

  else if (mat_type == "viscoelastic_fortran")
    return std::make_unique<Vaango::ViscoElasticFortran>(child, flags);
#endif

  else if (mat_type == "mw_visco_elastic")
    return std::make_unique<MWViscoElastic>(child, flags);

  else if (mat_type == "membrane")
    return std::make_unique<Membrane>(child, flags);

  else if (mat_type == "murnaghanMPM")
    return std::make_unique<MurnaghanMPM>(child, flags);

  else if (mat_type == "program_burn")
    return std::make_unique<ProgramBurn>(child, flags);

  else if (mat_type == "shell_CNH")
    return std::make_unique<ShellMaterial>(child, flags);

  else if (mat_type == "elastic_plastic")
    return std::make_unique<ElasticPlasticHP>(child, flags);

  else if (mat_type == "elastic_plastic_hp")
    return std::make_unique<ElasticPlasticHP>(child, flags);

  else if (mat_type == "soil_foam")
    return std::make_unique<SoilFoam>(child, flags);

  else if (mat_type == "visco_plastic")
    return std::make_unique<ViscoPlastic>(child, flags);

  else if (mat_type == "murnaghanMPM")
    return std::make_unique<MurnaghanMPM>(child, flags);

  else if (mat_type == "jwlpp_mpm")
    return std::make_unique<JWLppMPM>(child, flags);

  else if (mat_type == "camclay")
    return std::make_unique<CamClay>(child, flags);

  else if (mat_type == "polar_orthotropic_hypoelastic")
    return std::make_unique<PolarOrthotropicHypoElastic>(child, flags);

  else if (mat_type == "hypo_elastic_mms") {
    return std::make_unique<HypoElastic_MMS>(child, flags);
  }

  else
    throw ProblemSetupException(
      "Unknown Material Type R (" + mat_type + ")", __FILE__, __LINE__);

  return nullptr;
}

} // end namespace Uintah