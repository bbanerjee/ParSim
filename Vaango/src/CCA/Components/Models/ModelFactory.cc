/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#include <CCA/Components/Models/ModelFactory.h>

#include <sci_defs/uintah_defs.h>

#if !defined(NO_ICE)
#include <CCA/Components/Models/FluidsBased/AdiabaticTable.h>
#include <CCA/Components/Models/FluidsBased/BinaryProperties.h>
#include <CCA/Components/Models/FluidsBased/MassMomEng_src.h>
#include <CCA/Components/Models/FluidsBased/Mixing.h>
#include <CCA/Components/Models/FluidsBased/PassiveScalar.h>
#include <CCA/Components/Models/FluidsBased/TestModel.h>
#include <CCA/Components/Models/FluidsBased/flameSheet_rxn.h>
#include <CCA/Components/Models/ParticleBased/TracerParticles.h>
#endif

#if !defined(NO_ICE) && !defined(NO_MPM)
#include <CCA/Components/Models/HEChem/DDT0.h>
#include <CCA/Components/Models/HEChem/DDT1.h>
#include <CCA/Components/Models/HEChem/IandG.h>
#include <CCA/Components/Models/HEChem/JWLpp.h>
#include <CCA/Components/Models/HEChem/LightTime.h>
#include <CCA/Components/Models/HEChem/MesoBurn.h>
#include <CCA/Components/Models/HEChem/Simple_Burn.h>
#include <CCA/Components/Models/HEChem/Steady_Burn.h>
#include <CCA/Components/Models/HEChem/Unsteady_Burn.h>
#include <CCA/Components/Models/HEChem/ZeroOrder.h>

#include <CCA/Components/Models/SolidReactionModel/SolidReactionModel.h>
#endif

#include <CCA/Ports/ModelInterface.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>
#include <vector>

using namespace Uintah;

std::vector<std::unique_ptr<ModelInterface>>
ModelFactory::makeModels(const ProcessorGroup* myworld,
                         const MaterialManagerP materialManager,
                         const ProblemSpecP& restart_prob_spec,
                         const ProblemSpecP& prob_spec)
{
  std::vector<std::unique_ptr<ModelInterface>> d_models;

  ProblemSpecP model_spec = restart_prob_spec->findBlock("Models");

  if (!model_spec) {
    return d_models;
  }

  for (ProblemSpecP model_ps = model_spec->findBlock("Model");
       model_ps != nullptr;
       model_ps = model_ps->findNextBlock("Model")) {

    std::string type;

    if (!model_ps->getAttribute("type", type)) {
      throw ProblemSetupException(
        "\nERROR<Model>: Could not determine the type of the model.\n",
        __FILE__,
        __LINE__);
    }

    // A no-op is need to make the conditionals work correctly if
    // there are no models.
    if (0) {
    }

#if !defined(NO_ICE)
    // ICE turned on
    else if (type == "BinaryProperties") {
      d_models.emplace_back(
        std::make_unique<BinaryProperties>(myworld, materialManager, model_ps));
    } else if (type == "AdiabaticTable") {
      d_models.emplace_back(
        std::make_unique<AdiabaticTable>(myworld, materialManager, model_ps));
    } else if (type == "Mixing") {
      d_models.emplace_back(
        std::make_unique<Mixing>(myworld, materialManager, model_ps));
    } else if (type == "flameSheet_rxn") {
      d_models.emplace_back(
        std::make_unique<flameSheet_rxn>(myworld, materialManager, model_ps));
    } else if (type == "mass_momentum_energy_src") {
      d_models.emplace_back(
        std::make_unique<MassMomEng_src>(myworld, materialManager, model_ps));
    } else if (type == "PassiveScalar") {
      d_models.emplace_back(
        std::make_unique<PassiveScalar>(myworld, materialManager, model_ps));
    } else if (type == "TracerParticles") {
      d_models.emplace_back(
        std::make_unique<TracerParticles>(myworld, materialManager, model_ps));
    }
#endif

#if !defined(NO_ICE) && !defined(NO_MPM)
    else if (type == "Test") {
      d_models.emplace_back(
        std::make_unique<TestModel>(myworld, materialManager, model_ps));
    }

    else if (type == "Simple_Burn") {
      d_models.emplace_back(std::make_unique<Simple_Burn>(
        myworld, materialManager, model_ps, prob_spec));
    } else if (type == "Steady_Burn") {
      d_models.emplace_back(std::make_unique<Steady_Burn>(
        myworld, materialManager, model_ps, prob_spec));
    } else if (type == "Unsteady_Burn") {
      d_models.emplace_back(std::make_unique<Unsteady_Burn>(
        myworld, materialManager, model_ps, prob_spec));
    } else if (type == "MesoBurn") {
      d_models.emplace_back(std::make_unique<MesoBurn>(
        myworld, materialManager, model_ps, prob_spec));
    } else if (type == "IandG") {
      d_models.emplace_back(
        std::make_unique<IandG>(myworld, materialManager, model_ps));
    } else if (type == "JWLpp") {
      d_models.emplace_back(
        std::make_unique<JWLpp>(myworld, materialManager, model_ps, prob_spec));
    } else if (type == "ZeroOrder") {
      d_models.emplace_back(std::make_unique<ZeroOrder>(
        myworld, materialManager, model_ps, prob_spec));
    } else if (type == "LightTime") {
      d_models.emplace_back(
        std::make_unique<LightTime>(myworld, materialManager, model_ps));
    } else if (type == "DDT0") {
      d_models.emplace_back(
        std::make_unique<DDT0>(myworld, materialManager, model_ps, prob_spec));
    } else if (type == "DDT1") {
      d_models.emplace_back(
        std::make_unique<DDT1>(myworld, materialManager, model_ps, prob_spec));
    } else if (type == "SolidReactionModel") {
      d_models.emplace_back(std::make_unique<SolidReactionModel>(
        myworld, materialManager, model_ps, prob_spec));
    }
#endif

    else {
      std::ostringstream msg;
      msg << "\nERROR<Model>: Unknown model : " << type << ".\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  }

  return d_models;
}
