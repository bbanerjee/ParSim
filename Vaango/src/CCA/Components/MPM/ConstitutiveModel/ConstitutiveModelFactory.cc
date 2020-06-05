/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/PolarOrthotropicHypoElastic.h>

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/CompMooneyRivlin.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFracture.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticImplicit.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyper.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyperImplicit.h>

#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/MurnaghanMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/ProgramBurn.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/JWLppMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/ViscoSCRAMHotSpot.h>

#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/ElasticPlasticHP.h>
#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/UCNH.h>
#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/ViscoPlastic.h>

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

#include <CCA/Components/MPM/ConstitutiveModel/ManufacturedSolutions/CNH_MMS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ManufacturedSolutions/HypoElastic_MMS.h>

#if !defined(NO_FORTRAN)
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFortran.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Diamm.h>
#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/Kayenta.h>
#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticModels/ViscoElasticFortran.h>
#endif

#include <CCA/Components/MPM/MPMFlags.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <fstream>
#include <iostream>
#include <string>

using std::cerr;
using std::ifstream;
using std::ofstream;
using std::ostringstream;

using namespace Uintah;

ConstitutiveModel*
ConstitutiveModelFactory::create(ProblemSpecP& ps, MPMFlags* flags)
{
  ProblemSpecP child = ps->findBlock("constitutive_model");
  if (!child)
    throw ProblemSetupException("Cannot find constitutive_model tag", __FILE__,
                                __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for constitutive_model", __FILE__,
                                __LINE__);

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
    return (scinew RigidMaterial(child, flags));

  else if (mat_type == "comp_mooney_rivlin")
    return (scinew CompMooneyRivlin(child, flags));

  else if (mat_type == "nonlocal_drucker_prager")
    return (scinew NonLocalDruckerPrager(child, flags));

  else if (mat_type == "Arenisca")
    return (scinew Arenisca(child, flags));

  else if (mat_type == "Arenisca3")
    return (scinew Arenisca3(child, flags));

  else if (mat_type == "arena")
    return (scinew Vaango::Arena(child, flags));

  else if (mat_type == "arena_mixture")
    return (scinew Vaango::ArenaMixture(child, flags));

  else if (mat_type == "Arenisca4")
    return (scinew Arenisca4(child, flags));

  else if (mat_type == "soil_model_brannon")
    return (scinew SoilModelBrannon(child, flags));

  else if (mat_type == "tabular_eos")
    return (scinew Vaango::TabularEquationOfState(child, flags));

  else if (mat_type == "tabular_plasticity")
    return (scinew Vaango::TabularPlasticity(child, flags));

  else if (mat_type == "tabular_plasticity_cap")
    return (scinew Vaango::TabularPlasticityCap(child, flags));

  else if (mat_type == "mohr_coulomb") 
    return (scinew MohrCoulomb(child, flags));

  else if (mat_type == "comp_neo_hook") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return (scinew UCNH(child, flags, false, false));
    else if (flags->d_integratorType == "implicit")
      return (scinew UCNH(child, flags));
  } else if (mat_type == "cnh_damage")
    return (scinew UCNH(child, flags, false, true));

  else if (mat_type == "UCNH")
    return (scinew UCNH(child, flags));

  else if (mat_type == "cnh_mms")
    return (scinew CNH_MMS(child, flags));

  else if (mat_type == "cnhp_damage")
    return (scinew UCNH(child, flags, true, true));

  else if (mat_type == "trans_iso_hyper") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return (scinew TransIsoHyper(child, flags));
    else if (flags->d_integratorType == "implicit")
      return (scinew TransIsoHyperImplicit(child, flags));
  }

  else if (mat_type == "visco_trans_iso_hyper") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return (scinew ViscoTransIsoHyper(child, flags));
    else if (flags->d_integratorType == "implicit")
      return (scinew ViscoTransIsoHyperImplicit(child, flags));
  }

  else if (mat_type == "ideal_gas")
    return (scinew IdealGasMP(child, flags));

  else if (mat_type == "p_alpha")
    return (scinew P_Alpha(child, flags));

  else if (mat_type == "water")
    return (scinew Water(child, flags));

  else if (mat_type == "comp_neo_hook_plastic")
    return (scinew UCNH(child, flags, true, false));

  else if (mat_type == "visco_scram") {
    if (flags->d_integratorType == "explicit" ||
        flags->d_integratorType == "fracture")
      return (scinew ViscoScram(child, flags));
    else if (flags->d_integratorType == "implicit")
      return (scinew ViscoScramImplicit(child, flags));
  }

  else if (mat_type == "viscoSCRAM_hs")
    return (scinew ViscoSCRAMHotSpot(child, flags));

  else if (mat_type == "hypo_elastic") {
    if (flags->d_integratorType == "explicit")
      return (scinew HypoElastic(child, flags));

    else if (flags->d_integratorType == "implicit") {
      if (!flags->d_doGridReset) {
        ostringstream msg;
        msg << "\n ERROR: One may not use HypoElastic along with \n"
            << " <do_grid_reset>false</do_grid_reset> \n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
      return (scinew HypoElasticImplicit(child, flags));
    }
  }

  else if (mat_type == "hypo_elastic_fracture") 
    return (scinew HypoElasticFracture(child, flags));

#if !defined(NO_FORTRAN)
  else if (mat_type == "hypo_elastic_fortran")
    return (scinew HypoElasticFortran(child, flags));

  else if (mat_type == "kayenta")
    return (scinew Kayenta(child, flags));

  else if (mat_type == "diamm")
    return (scinew Diamm(child, flags));

  else if (mat_type == "viscoelastic_fortran")
    return (scinew Vaango::ViscoElasticFortran(child, flags));
#endif

  else if (mat_type == "mw_visco_elastic")
    return (scinew MWViscoElastic(child, flags));

  else if (mat_type == "membrane")
    return (scinew Membrane(child, flags));

  else if (mat_type == "murnaghanMPM")
    return (scinew MurnaghanMPM(child, flags));

  else if (mat_type == "program_burn")
    return (scinew ProgramBurn(child, flags));

  else if (mat_type == "shell_CNH")
    return (scinew ShellMaterial(child, flags));

  else if (mat_type == "elastic_plastic")
    return (scinew ElasticPlasticHP(child, flags));

  else if (mat_type == "elastic_plastic_hp")
    return (scinew ElasticPlasticHP(child, flags));

  else if (mat_type == "soil_foam")
    return (scinew SoilFoam(child, flags));

  else if (mat_type == "visco_plastic")
    return (scinew ViscoPlastic(child, flags));

  else if (mat_type == "murnaghanMPM")
    return (scinew MurnaghanMPM(child, flags));

  else if (mat_type == "jwlpp_mpm")
    return (scinew JWLppMPM(child, flags));

  else if (mat_type == "camclay")
    return (scinew CamClay(child, flags));

  else if (mat_type == "polar_orthotropic_hypoelastic")
    return (scinew PolarOrthotropicHypoElastic(child, flags));

  else if (mat_type == "hypo_elastic_mms") {
    return (scinew HypoElastic_MMS(child, flags));
  }

  else
    throw ProblemSetupException("Unknown Material Type R (" + mat_type + ")",
                                __FILE__, __LINE__);

  return nullptr;
}
