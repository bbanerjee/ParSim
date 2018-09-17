#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/SimulationController/AMRSimulationController.h>
#include <CCA/Components/Regridder/RegridderCommon.h>
#include <CCA/Components/Solvers/SolverFactory.h>
#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>
#include <CCA/Components/LoadBalancers/LoadBalancerFactory.h>
#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <CCA/Components/Schedulers/SchedulerCommon.h>
#include <CCA/Components/Schedulers/SchedulerFactory.h>

#include <CCA/Ports/SolverInterface.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/Output.h>

#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/Exception.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/Environment.h>
#include <Core/OS/Dir.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using nlohmann::json;
using Uintah::Dir;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;
using Uintah::ProblemSpecReader;
using Uintah::Exception;
using Uintah::ProblemSetupException;
using Uintah::InvalidValue;
using Uintah::ProcessorGroup;
using Uintah::SimulationController;
using Uintah::AMRSimulationController;
using Uintah::RegridderCommon;
using Uintah::SolverInterface;
using Uintah::SolverFactory;
using Uintah::UintahParallelComponent;
using Uintah::ComponentFactory;
using Uintah::SimulationInterface;
using Uintah::LoadBalancerCommon;
using Uintah::LoadBalancerFactory;
using Uintah::DataArchiver;
using Uintah::Output;
using Uintah::SchedulerCommon;
using Uintah::SchedulerFactory;

class VaangoEnv : public ::testing::Environment {
public:

  int d_argc;
  char** d_argv;
  char** d_env;

  explicit VaangoEnv(int argc, char** argv, char* env[]) {
    d_argc = argc;
    d_argv = argv;
    d_env = env;
  }

  virtual ~VaangoEnv() {}

  virtual void SetUp() {
    Uintah::Parallel::determineIfRunningUnderMPI(d_argc, d_argv);
    Uintah::Parallel::initializeManager(d_argc, d_argv);
    Uintah::create_sci_environment(d_env, 0, true );
  }

  virtual void TearDown() {
    Uintah::Parallel::finalizeManager();
  }

  static ProblemSpecP createInput() {

    char currPath[2000];
    if (!getcwd(currPath, sizeof(currPath))) {
      std::cout << "Current path not found\n";
    }
    //std::cout << "Dir = " << currPath << std::endl;

    // Create a new document
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

    // Create root node
    xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
    xmlDocSetRootElement(doc, rootNode);

    // Meta
    auto meta = xmlNewChild(rootNode, nullptr, BAD_CAST "Meta", BAD_CAST "");
    xmlNewChild(meta, nullptr, BAD_CAST "title", BAD_CAST "Unit test Tabular Plasticity");

    // Simulation component 
    auto simComp = xmlNewChild(rootNode, nullptr, BAD_CAST "SimulationComponent",
                               BAD_CAST "");
    xmlNewProp(simComp, BAD_CAST "type", BAD_CAST "mpm");

    // Time
    auto time = xmlNewChild(rootNode, nullptr, BAD_CAST "Time", BAD_CAST "");
    xmlNewChild(time, nullptr, BAD_CAST "maxTime", BAD_CAST "1.0");
    xmlNewChild(time, nullptr, BAD_CAST "initTime", BAD_CAST "0.0");
    xmlNewChild(time, nullptr, BAD_CAST "delt_min", BAD_CAST "1.0e-6");
    xmlNewChild(time, nullptr, BAD_CAST "delt_max", BAD_CAST "0.04");
    xmlNewChild(time, nullptr, BAD_CAST "timestep_multiplier", BAD_CAST "0.3");
    xmlNewChild(time, nullptr, BAD_CAST "max_Timesteps", BAD_CAST "5");

    // DataArchiver
    auto da = xmlNewChild(rootNode, nullptr, BAD_CAST "DataArchiver", BAD_CAST "");
    xmlNewChild(da, nullptr, BAD_CAST "filebase", BAD_CAST "UniaxialStrainRotateTabularPlasticityCap.uda");
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
    std::string prescribed = 
      std::string(currPath) + "/" + "UniaxialStrainRotate_PrescribedDeformation.inp";
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
    xmlNewChild(mpm, nullptr, BAD_CAST "prescribed_deformation_file", BAD_CAST prescribed.c_str());
    xmlNewChild(mpm, nullptr, BAD_CAST "minimum_subcycles_for_F", BAD_CAST "-2");
    auto ero = xmlNewChild(mpm, nullptr, BAD_CAST "erosion", BAD_CAST "");
    xmlNewProp(ero, BAD_CAST "algorithm", BAD_CAST "none");

    // Physical constants
    auto pc = xmlNewChild(rootNode, nullptr, BAD_CAST "PhysicalConstants", BAD_CAST "");
    xmlNewChild(pc, nullptr, BAD_CAST "gravity", BAD_CAST "[0,0,0]");

    // Material properties
    auto matProp = 
      xmlNewChild(rootNode, nullptr, BAD_CAST "MaterialProperties", BAD_CAST "");
    mpm = xmlNewChild(matProp, nullptr, BAD_CAST "MPM", BAD_CAST "");
    auto mat = xmlNewChild(mpm, nullptr, BAD_CAST "material", BAD_CAST "");
    xmlNewProp(mat, BAD_CAST "name", BAD_CAST "TabularPlasticCap");

    // General properties
    xmlNewChild(mat, nullptr, BAD_CAST "density", BAD_CAST "1050");
    xmlNewChild(mat, nullptr, BAD_CAST "melt_temp", BAD_CAST "3695.0");
    xmlNewChild(mat, nullptr, BAD_CAST "room_temp", BAD_CAST "294.0");
    xmlNewChild(mat, nullptr, BAD_CAST "thermal_conductivity", BAD_CAST "174.0e-7");
    xmlNewChild(mat, nullptr, BAD_CAST "specific_heat", BAD_CAST "134.0e-8");
    auto cm = xmlNewChild(mat, nullptr, BAD_CAST "constitutive_model", BAD_CAST "");
    xmlNewProp(cm, BAD_CAST "type", BAD_CAST "tabular_plasticity_cap");

    // Elastic properties
    std::string table_elastic = 
      std::string(currPath) + "/" + "tabular_linear_elastic.json";
    auto elastic = xmlNewChild(cm, nullptr, BAD_CAST "elastic_moduli_model", 
                               BAD_CAST "");
    xmlNewProp(elastic, BAD_CAST "type", BAD_CAST "tabular");
    xmlNewChild(elastic, nullptr, BAD_CAST "filename", 
                BAD_CAST table_elastic.c_str());
    xmlNewChild(elastic, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "PlasticStrainVol, TotalStrainVol");
    xmlNewChild(elastic, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "Pressure");
    auto interp_elastic = xmlNewChild(elastic, nullptr, BAD_CAST "interpolation",
                                      BAD_CAST "");
    xmlNewProp(interp_elastic, BAD_CAST "type", BAD_CAST "linear");
    xmlNewChild(elastic, nullptr, BAD_CAST "G0", 
                BAD_CAST "1.0e4");
    xmlNewChild(elastic, nullptr, BAD_CAST "nu", 
              BAD_CAST "0.2");

    // Yield criterion
    std::string table_yield = 
      std::string(currPath) + "/" + "tabular_drucker_prager.json";
    auto yield = xmlNewChild(cm, nullptr, BAD_CAST "plastic_yield_condition", 
                             BAD_CAST "");
    xmlNewProp(yield, BAD_CAST "type", BAD_CAST "tabular_cap");
    xmlNewChild(yield, nullptr, BAD_CAST "filename", 
                BAD_CAST table_yield.c_str());
    xmlNewChild(yield, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "Pressure");
    xmlNewChild(yield, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "SqrtJ2");
    auto yield_interp = xmlNewChild(yield, nullptr, BAD_CAST "interpolation",
                              BAD_CAST "");
    xmlNewProp(yield_interp, BAD_CAST "type", BAD_CAST "linear");
    xmlNewChild(yield, nullptr, BAD_CAST "cap_ellipticity_ratio", 
                BAD_CAST "0.7");

    // Hydrostat
    std::string table_hydrostat = 
      std::string(currPath) + "/" + "DrySand_HydrostatData.json";
    xmlNewChild(cm, nullptr, BAD_CAST "filename", 
                BAD_CAST table_hydrostat.c_str());
    xmlNewChild(cm, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "TotalStrainVol");
    xmlNewChild(cm, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "Pressure");
    auto cm_interp = xmlNewChild(cm, nullptr, BAD_CAST "interpolation",
                              BAD_CAST "");
    xmlNewProp(cm_interp, BAD_CAST "type", BAD_CAST "linear");
    
    // Cap evolution
    std::string table_cap = 
      std::string(currPath) + "/" + "tabular_cap.json";
    auto cap = xmlNewChild(cm, nullptr, BAD_CAST "internal_variable_model", 
                           BAD_CAST "");
    xmlNewProp(cap, BAD_CAST "type", BAD_CAST "tabular_cap");
    xmlNewChild(cap, nullptr, BAD_CAST "filename", 
                BAD_CAST table_cap.c_str());
    xmlNewChild(cap, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "PlasticStrainVol");
    xmlNewChild(cap, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "Pressure");
    auto cap_interp = xmlNewChild(cap, nullptr, BAD_CAST "interpolation",
                                  BAD_CAST "");
    xmlNewProp(cap_interp, BAD_CAST "type", BAD_CAST "linear");

    // Geometry
    auto geom = xmlNewChild(mat, nullptr, BAD_CAST "geom_object", BAD_CAST "");
    auto box = xmlNewChild(geom, nullptr, BAD_CAST "box", BAD_CAST "");
    xmlNewProp(box, BAD_CAST "label", BAD_CAST "Plate1");
    xmlNewChild(box, nullptr, BAD_CAST "min", BAD_CAST "[0.0,0.0,0.0]");
    xmlNewChild(box, nullptr, BAD_CAST "max", BAD_CAST "[1.0,1.0,1.0]");
    xmlNewChild(geom, nullptr, BAD_CAST "res", BAD_CAST "[1,1,1]");
    xmlNewChild(geom, nullptr, BAD_CAST "velocity", BAD_CAST "[0,0,0]");
    xmlNewChild(geom, nullptr, BAD_CAST "temperature", BAD_CAST "294");
    xmlNewChild(geom, nullptr, BAD_CAST "color", BAD_CAST "0");

    // Contact
    auto contact = xmlNewChild(mpm, nullptr, BAD_CAST "contact", BAD_CAST "");
    xmlNewChild(contact, nullptr, BAD_CAST "type", BAD_CAST "null");
    xmlNewChild(contact, nullptr, BAD_CAST "materials", BAD_CAST "[0]");
    xmlNewChild(contact, nullptr, BAD_CAST "mu", BAD_CAST "0.1");

    // Grid
    auto grid = xmlNewChild(rootNode, nullptr, BAD_CAST "Grid", BAD_CAST "");
    xmlNewChild(grid, nullptr, BAD_CAST "BoundaryConditions", BAD_CAST "");
    auto level = xmlNewChild(grid, nullptr, BAD_CAST "Level", BAD_CAST "");
    box = xmlNewChild(level, nullptr, BAD_CAST "Box", BAD_CAST "");
    xmlNewProp(box, BAD_CAST "label", BAD_CAST "1");
    xmlNewChild(box, nullptr, BAD_CAST "lower", BAD_CAST "[-2,-2,-2]");
    xmlNewChild(box, nullptr, BAD_CAST "upper", BAD_CAST "[3,3,3]");
    xmlNewChild(box, nullptr, BAD_CAST "resolution", BAD_CAST "[5,5,5]");
    xmlNewChild(box, nullptr, BAD_CAST "extraCells", BAD_CAST "[0,0,0]");
    xmlNewChild(box, nullptr, BAD_CAST "patches", BAD_CAST "[1,1,1]");

    // Print the document to stdout
    //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);
    std::string ups_file = std::string(currPath) + "/" + 
                           "UniaxialStrainRotateTabularPlasticityCap.ups";
    xmlSaveFormatFileEnc(ups_file.c_str(), doc, "ISO-8859-1", 1);

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

int main(int argc, char** argv, char* env[]) {

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new VaangoEnv(argc, argv, env));
  return RUN_ALL_TESTS();
}

TEST(TabularPlasticityCapTest, singleParticleTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string ups_file = std::string(currPath) + "/" + 
                         "UniaxialStrainRotateTabularPlasticityCap.ups";
  //std::cout << "Created Filename = " << ups_file << std::endl;

  // Remove existing uda
  std::string uda_file = std::string(currPath) + "/" + 
                         "UniaxialStrainRotateTabularPlasticityCap.uda.000";
  Dir::removeDir(uda_file.c_str());

  char * start_addr = (char*)sbrk(0);
  bool thrownException = false;

  try {
    ProblemSpecP ups = VaangoEnv::createInput();
    ups->getNode()->_private = (void *) ups_file.c_str();
    //std::cout << "Filename = " << static_cast<char*>(ups->getNode()->_private) << std::endl;

    const ProcessorGroup* world = Uintah::Parallel::getRootProcessorGroup();
    SimulationController* ctl = scinew AMRSimulationController(world, false, ups);
    
    RegridderCommon* reg = 0;
    SolverInterface* solve = SolverFactory::create(ups, world, "");

    UintahParallelComponent* comp = ComponentFactory::create(ups, world, false, "");
    SimulationInterface* sim = dynamic_cast<SimulationInterface*>(comp);
    ctl->attachPort("sim", sim);
    comp->attachPort("solver", solve);
    comp->attachPort("regridder", reg);

    LoadBalancerCommon* lbc = LoadBalancerFactory::create(ups, world);
    lbc->attachPort("sim", sim);

    DataArchiver* dataarchiver = scinew DataArchiver(world, -1);
    Output* output = dataarchiver;
    ctl->attachPort("output", dataarchiver);
    dataarchiver->attachPort("load balancer", lbc);
    comp->attachPort("output", dataarchiver);
    dataarchiver->attachPort("sim", sim);

    SchedulerCommon* sched = SchedulerFactory::create(ups, world, output);
    sched->attachPort("load balancer", lbc);
    ctl->attachPort("scheduler", sched);
    lbc->attachPort("scheduler", sched);
    comp->attachPort("scheduler", sched);
    sched->setStartAddr( start_addr );
    sched->addReference();

    ctl->run();
    delete ctl;

    sched->removeReference();
    delete sched;
    delete lbc;
    delete sim;
    delete solve;
    delete output; 

  } catch (ProblemSetupException& e) {
    std::cout << e.message() << std::endl;
    thrownException = true;
  } catch (Exception& e) {
    std::cout << e.message() << std::endl;
    thrownException = true;
  } catch (...) {
    std::cout << "**ERROR** Unknown exception" << std::endl;
    thrownException = true;
  }

}
