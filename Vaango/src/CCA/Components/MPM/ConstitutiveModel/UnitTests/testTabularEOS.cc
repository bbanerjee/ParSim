#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>

#include <Core/Parallel/Parallel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;
using Uintah::ProblemSetupException;
using Uintah::InvalidValue;
using nlohmann::json;

class VaangoEnv : public ::testing::Environment {
public:

  int d_argc;
  char** d_argv;

  explicit VaangoEnv(int argc, char**argv) {
    d_argc = argc;
    d_argv = argv;
  }

  virtual ~VaangoEnv() {}

  virtual void SetUp() {
    Uintah::Parallel::determineIfRunningUnderMPI(d_argc, d_argv);
    Uintah::Parallel::initializeManager(d_argc, d_argv);
    auto ps = createInput();


  }

  virtual void TearDown() {
    Uintah::Parallel::finalizeManager();

  }

  ProblemSpecP createInput() {

    // Create a new document
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

    // Create root node
    xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
    xmlDocSetRootElement(doc, rootNode);

    // Meta
    auto meta = xmlNewChild(rootNode, nullptr, BAD_CAST "Meta", BAD_CAST "");
    xmlNewChild(meta, nullptr, BAD_CAST "title", BAD_CAST "Unit test Tabular EOS");

    // Simulation component 
    auto simComp = xmlNewChild(rootNode, nullptr, BAD_CAST "SimulationComponent",
                               BAD_CAST "");
    xmlNewProp(simComp, BAD_CAST "type", BAD_CAST "mpm");

    // Time
    auto time = xmlNewChild(rootNode, nullptr, BAD_CAST "Time", BAD_CAST "");
    xmlNewChild(time, nullptr, BAD_CAST "maxTime", BAD_CAST "1.0");
    xmlNewChild(time, nullptr, BAD_CAST "initTime", BAD_CAST "0.0");
    xmlNewChild(time, nullptr, BAD_CAST "delt_min", BAD_CAST "1.0e-6");
    xmlNewChild(time, nullptr, BAD_CAST "delt_max", BAD_CAST "0.01");
    xmlNewChild(time, nullptr, BAD_CAST "timestep_multiplier", BAD_CAST "0.3");

    // DataArchiver
    auto da = xmlNewChild(rootNode, nullptr, BAD_CAST "DataArchiver", BAD_CAST "");
    xmlNewChild(da, nullptr, BAD_CAST "outputTimestepInterval", BAD_CAST "1");
    auto save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "g.mass");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.x");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.color");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.temperature");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.velocity");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.particleID");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.stress");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.deformationGradient");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST ""); 
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "g.acceleration");
    save = xmlNewChild(da, nullptr, BAD_CAST "checkpoint", BAD_CAST "");
    xmlNewProp(save, BAD_CAST "cycle", BAD_CAST "2");
    xmlNewProp(save, BAD_CAST "timestepInterval", BAD_CAST "4000");

    // MPM
    auto mpm = xmlNewChild(rootNode, nullptr, BAD_CAST "MPM", BAD_CAST "");
    xmlNewChild(mpm, nullptr, BAD_CAST "time_integrator", BAD_CAST "explicit");
    xmlNewChild(mpm, nullptr, BAD_CAST "interpolator", BAD_CAST "linear");
    xmlNewChild(mpm, nullptr, BAD_CAST "use_load_curves", BAD_CAST "false");
    xmlNewChild(mpm, nullptr, BAD_CAST "minimum_particle_mass", BAD_CAST "1.0e-15");
    xmlNewChild(mpm, nullptr, BAD_CAST "minimum_mass_for_acc", BAD_CAST "1.0e-15");
    xmlNewChild(mpm, nullptr, BAD_CAST "maximum_particle_velocity", BAD_CAST "1.0e5");
    xmlNewChild(mpm, nullptr, BAD_CAST "artificial_damping_coeff", BAD_CAST "0.0");
    xmlNewChild(mpm, nullptr, BAD_CAST "artificial_viscosity", BAD_CAST "true");
    xmlNewChild(mpm, nullptr, BAD_CAST "artificial_viscosity_heating", BAD_CAST "false");
    xmlNewChild(mpm, nullptr, BAD_CAST "do_contact_friction_heating", BAD_CAST "false");
    xmlNewChild(mpm, nullptr, BAD_CAST "create_new_particles", BAD_CAST "false");
    xmlNewChild(mpm, nullptr, BAD_CAST "use_momentum_form", BAD_CAST "false");
    xmlNewChild(mpm, nullptr, BAD_CAST "with_color", BAD_CAST "true");
    xmlNewChild(mpm, nullptr, BAD_CAST "use_prescribed_deformation", BAD_CAST "true");
    xmlNewChild(mpm, nullptr, BAD_CAST "prescribed_deformation_file", BAD_CAST "HydrostaticCompression_PrescribedDeformation.inp");
    xmlNewChild(mpm, nullptr, BAD_CAST "minimum_subcycles_for_F", BAD_CAST "-2");
    auto ero = xmlNewChild(mpm, nullptr, BAD_CAST "erosion", BAD_CAST "");
    xmlNewProp(ero, BAD_CAST "algorithm", BAD_CAST "none");

    // Physical constants
    auto pc = xmlNewChild(rootNode, nullptr, BAD_CAST "PhysicalConstants", BAD_CAST "");
    xmlNewChild(pc, nullptr, BAD_CAST "gravity", BAD_CAST "[0,0,0]");

    // Material properties
    auto matProp = xmlNewChild(rootNode, nullptr, BAD_CAST "MaterialProperties", BAD_CAST "");
    mpm = xmlNewChild(matProp, nullptr, BAD_CAST "MPM", BAD_CAST "");
    auto mat = xmlNewChild(mpm, nullptr, BAD_CAST "material", BAD_CAST "");
    xmlNewProp(mat, BAD_CAST "name", BAD_CAST "TabularEOS");
    xmlNewChild(mat, nullptr, BAD_CAST "density", BAD_CAST "1050");
    xmlNewChild(mat, nullptr, BAD_CAST "melt_temp", BAD_CAST "3695.0");
    xmlNewChild(mat, nullptr, BAD_CAST "room_temp", BAD_CAST "294.0");
    xmlNewChild(mat, nullptr, BAD_CAST "thermal_conductivity", BAD_CAST "174.0e-7");
    xmlNewChild(mat, nullptr, BAD_CAST "specific_heat", BAD_CAST "134.0e-8");
    auto cm = xmlNewChild(mat, nullptr, BAD_CAST "constitutive_model", BAD_CAST "");
    xmlNewProp(cm, BAD_CAST "type", BAD_CAST "tabular_eos");
    xmlNewChild(cm, nullptr, BAD_CAST "filename", BAD_CAST "tabular_eos.json");
    xmlNewChild(cm, nullptr, BAD_CAST "independent_variables", BAD_CAST "DensityRatio");
    xmlNewChild(cm, nullptr, BAD_CAST "dependent_variables", BAD_CAST "Pressure");
    auto interp = xmlNewChild(cm, nullptr, BAD_CAST "interpolation", BAD_CAST "");
    xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");
    auto geom = xmlNewChild(mat, nullptr, BAD_CAST "geom_object", BAD_CAST "");
    auto box = xmlNewChild(geom, nullptr, BAD_CAST "box", BAD_CAST "");
    xmlNewProp(box, BAD_CAST "label", BAD_CAST "Plate1");
    xmlNewChild(box, nullptr, BAD_CAST "min", BAD_CAST "[0.0,0.0,0.0]");
    xmlNewChild(box, nullptr, BAD_CAST "max", BAD_CAST "[1.0,1.0,1.0]");
    xmlNewChild(geom, nullptr, BAD_CAST "res", BAD_CAST "[1,1,1]");
    xmlNewChild(geom, nullptr, BAD_CAST "velocity", BAD_CAST "[0,0,0]");
    xmlNewChild(geom, nullptr, BAD_CAST "temperature", BAD_CAST "294");
    xmlNewChild(geom, nullptr, BAD_CAST "color", BAD_CAST "0");
    auto contact = xmlNewChild(mpm, nullptr, BAD_CAST "contact", BAD_CAST "");
    xmlNewChild(contact, nullptr, BAD_CAST "type", BAD_CAST "null");
    xmlNewChild(contact, nullptr, BAD_CAST "materials", BAD_CAST "[0]");
    xmlNewChild(contact, nullptr, BAD_CAST "mu", BAD_CAST "0.1");

    // Grid
    auto grid = xmlNewChild(rootNode, nullptr, BAD_CAST "Grid", BAD_CAST "");
    xmlNewChild(grid, nullptr, BAD_CAST "BoundaryConditions", BAD_CAST "");
    auto level = xmlNewChild(grid, nullptr, BAD_CAST "Level", BAD_CAST "");
    box = xmlNewChild(level, nullptr, BAD_CAST "box", BAD_CAST "");
    xmlNewProp(box, BAD_CAST "label", BAD_CAST "1");
    xmlNewChild(box, nullptr, BAD_CAST "lower", BAD_CAST "[-2,-2,-2]");
    xmlNewChild(box, nullptr, BAD_CAST "upper", BAD_CAST "[3,3,3]");
    xmlNewChild(box, nullptr, BAD_CAST "resolution", BAD_CAST "[5,5,5]");
    xmlNewChild(box, nullptr, BAD_CAST "extraCells", BAD_CAST "[0,0,0]");
    xmlNewChild(box, nullptr, BAD_CAST "patches", BAD_CAST "[1,1,1]");

    // Print the document to stdout
    xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    if (!ps) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }

    return ps;
  }
};

int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new VaangoEnv(argc, argv));
  return RUN_ALL_TESTS();
}

TEST(TabularDataTest, parseVariableNames)
{
}
