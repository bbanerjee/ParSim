#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>

#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>
#include <CCA/Components/LoadBalancers/LoadBalancerFactory.h>
#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/Regridder/RegridderCommon.h>
#include <CCA/Components/Schedulers/SchedulerCommon.h>
#include <CCA/Components/Schedulers/SchedulerFactory.h>
#include <CCA/Components/SimulationController/AMRSimulationController.h>
#include <CCA/Components/Solvers/SolverFactory.h>

#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/SolverInterface.h>

#include <Core/Exceptions/Exception.h>
#include <Core/Malloc/Allocator.h>
#include <Core/OS/Dir.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/Environment.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Uintah;

class DEMEnv : public ::testing::Environment
{
public:
  int d_argc;
  char** d_argv;
  char** d_env;

  explicit DEMEnv(int argc, char** argv, char* env[])
  {
    d_argc = argc;
    d_argv = argv;
    d_env  = env;
  }

  virtual ~DEMEnv() {}

  virtual void SetUp()
  {
    Uintah::Parallel::initializeManager(d_argc, d_argv);
    Uintah::create_sci_environment(d_env, 0, true);
  }

  virtual void TearDown()
  {
    Uintah::Parallel::finalizeManager();
  }

  static ProblemSpecP createInput()
  {
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
    xmlDocSetRootElement(doc, rootNode);

    auto meta = xmlNewChild(rootNode, nullptr, BAD_CAST "Meta", BAD_CAST "");
    xmlNewChild(meta, nullptr, BAD_CAST "title", BAD_CAST "Unit test Hybrid DEM-MPM");

    auto simComp = xmlNewChild(rootNode, nullptr, BAD_CAST "SimulationComponent", BAD_CAST "");
    xmlNewProp(simComp, BAD_CAST "type", BAD_CAST "mpm");

    auto time = xmlNewChild(rootNode, nullptr, BAD_CAST "Time", BAD_CAST "");
    xmlNewChild(time, nullptr, BAD_CAST "maxTime", BAD_CAST "0.1");
    xmlNewChild(time, nullptr, BAD_CAST "initTime", BAD_CAST "0.0");
    xmlNewChild(time, nullptr, BAD_CAST "delt_min", BAD_CAST "1.0e-6");
    xmlNewChild(time, nullptr, BAD_CAST "delt_max", BAD_CAST "0.01");
    xmlNewChild(time, nullptr, BAD_CAST "timestep_multiplier", BAD_CAST "0.1");
    xmlNewChild(time, nullptr, BAD_CAST "max_timesteps", BAD_CAST "2");

    auto da = xmlNewChild(rootNode, nullptr, BAD_CAST "DataArchiver", BAD_CAST "");
    xmlNewChild(da, nullptr, BAD_CAST "filebase", BAD_CAST "DEMTest.uda");
    xmlNewChild(da, nullptr, BAD_CAST "outputTimestepInterval", BAD_CAST "1");
    auto save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST "");
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.x");
    save = xmlNewChild(da, nullptr, BAD_CAST "save", BAD_CAST "");
    xmlNewProp(save, BAD_CAST "label", BAD_CAST "p.externalforce");

    auto mpm = xmlNewChild(rootNode, nullptr, BAD_CAST "MPM", BAD_CAST "");
    xmlNewChild(mpm, nullptr, BAD_CAST "time_integrator", BAD_CAST "explicit");
    xmlNewChild(mpm, nullptr, BAD_CAST "enable_dem", BAD_CAST "true");

    auto pc = xmlNewChild(rootNode, nullptr, BAD_CAST "PhysicalConstants", BAD_CAST "");
    xmlNewChild(pc, nullptr, BAD_CAST "gravity", BAD_CAST "[0,0,0]");

    auto matProp = xmlNewChild(rootNode, nullptr, BAD_CAST "MaterialProperties", BAD_CAST "");
    auto mpm_mat = xmlNewChild(matProp, nullptr, BAD_CAST "MPM", BAD_CAST "");
    auto mat = xmlNewChild(mpm_mat, nullptr, BAD_CAST "material", BAD_CAST "");
    xmlNewChild(mat, nullptr, BAD_CAST "density", BAD_CAST "1000");
    xmlNewChild(mat, nullptr, BAD_CAST "thermal_conductivity", BAD_CAST "1.0");
    xmlNewChild(mat, nullptr, BAD_CAST "specific_heat", BAD_CAST "1.0");
    xmlNewChild(mat, nullptr, BAD_CAST "is_dem_material", BAD_CAST "true");
    xmlNewChild(mat, nullptr, BAD_CAST "radius", BAD_CAST "0.5");
    xmlNewChild(mat, nullptr, BAD_CAST "kn", BAD_CAST "1.0e6");

    auto cm = xmlNewChild(mat, nullptr, BAD_CAST "constitutive_model", BAD_CAST "");
    xmlNewProp(cm, BAD_CAST "type", BAD_CAST "comp_neo_hook_ideal_gas");
    xmlNewChild(cm, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e6");
    xmlNewChild(cm, nullptr, BAD_CAST "shear_modulus", BAD_CAST "1.0e5");

    // Two spheres close to each other
    auto geom = xmlNewChild(mat, nullptr, BAD_CAST "geom_object", BAD_CAST "");
    auto sphere1 = xmlNewChild(geom, nullptr, BAD_CAST "sphere", BAD_CAST "");
    xmlNewChild(sphere1, nullptr, BAD_CAST "origin", BAD_CAST "[0,0,0]");
    xmlNewChild(sphere1, nullptr, BAD_CAST "radius", BAD_CAST "0.5");
    xmlNewChild(geom, nullptr, BAD_CAST "res", BAD_CAST "[1,1,1]");
    xmlNewChild(geom, nullptr, BAD_CAST "velocity", BAD_CAST "[0,0,0]");
    xmlNewChild(geom, nullptr, BAD_CAST "temperature", BAD_CAST "300");

    geom = xmlNewChild(mat, nullptr, BAD_CAST "geom_object", BAD_CAST "");
    auto sphere2 = xmlNewChild(geom, nullptr, BAD_CAST "sphere", BAD_CAST "");
    xmlNewChild(sphere2, nullptr, BAD_CAST "origin", BAD_CAST "[0.8,0,0]"); // Overlap by 0.2
    xmlNewChild(sphere2, nullptr, BAD_CAST "radius", BAD_CAST "0.5");
    xmlNewChild(geom, nullptr, BAD_CAST "res", BAD_CAST "[1,1,1]");
    xmlNewChild(geom, nullptr, BAD_CAST "velocity", BAD_CAST "[0,0,0]");
    xmlNewChild(geom, nullptr, BAD_CAST "temperature", BAD_CAST "300");

    auto grid = xmlNewChild(rootNode, nullptr, BAD_CAST "Grid", BAD_CAST "");
    xmlNewChild(grid, nullptr, BAD_CAST "BoundaryConditions", BAD_CAST "");
    auto level = xmlNewChild(grid, nullptr, BAD_CAST "Level", BAD_CAST "");
    auto box = xmlNewChild(level, nullptr, BAD_CAST "Box", BAD_CAST "");
    xmlNewChild(box, nullptr, BAD_CAST "lower", BAD_CAST "[-2,-2,-2]");
    xmlNewChild(box, nullptr, BAD_CAST "upper", BAD_CAST "[2,2,2]");
    xmlNewChild(box, nullptr, BAD_CAST "resolution", BAD_CAST "[10,10,10]");
    xmlNewChild(box, nullptr, BAD_CAST "patches", BAD_CAST "[1,1,1]");

    std::string ups_file = "DEMTest.ups";
    xmlSaveFormatFileEnc(ups_file.c_str(), doc, "ISO-8859-1", 1);

    ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    return ps;
  }
};

TEST(DEMTest, CollisionTest)
{
  char* start_addr = (char*)sbrk(0);
  try {
    ProblemSpecP ups = DEMEnv::createInput();
    std::string ups_file = "DEMTest.ups";
    ups->getNode()->_private = (void*)ups_file.c_str();

    const ProcessorGroup* world = Uintah::Parallel::getRootProcessorGroup();
    std::unique_ptr<SimulationController> ctl = std::make_unique<AMRSimulationController>(world, ups);
    std::unique_ptr<UintahParallelComponent> simComp = ComponentFactory::create(ups, world, nullptr, "");
    SimulationInterface* simulator = dynamic_cast<SimulationInterface*>(simComp.get());
    simulator->problemSetup(ups);
    ctl->attachPort("simulator", simulator);

    std::shared_ptr<SolverInterface> solver = SolverFactory::create(ups, world, "");
    UintahParallelComponent* solverComp = dynamic_cast<UintahParallelComponent*>(solver.get());
    simComp->attachPort("solver", solver.get());
    solverComp->attachPort("simulator", simulator);

    std::unique_ptr<LoadBalancerCommon> lbc = LoadBalancerFactory::create(ups, world);
    lbc->attachPort("simulator", simulator);
    ctl->attachPort("load balancer", lbc.get());
    simComp->attachPort("load balancer", lbc.get());

    SchedulerCommon* sched = SchedulerFactory::create(ups, world);
    sched->attachPort("load balancer", lbc.get());
    sched->attachPort("simulator", simulator);
    ctl->attachPort("scheduler", sched);
    simComp->attachPort("scheduler", sched);
    lbc->attachPort("scheduler", sched);

    sched->setStartAddr(start_addr);
    sched->addReference();

    std::unique_ptr<DataArchiver> dataarchiver = std::make_unique<DataArchiver>(world, -1);
    dataarchiver->attachPort("simulator", simulator);
    dataarchiver->attachPort("load balancer", lbc.get());
    ctl->attachPort("output", dataarchiver.get());
    simComp->attachPort("output", dataarchiver.get());
    sched->attachPort("output", dataarchiver.get());

    sched->getComponents();
    lbc->getComponents();
    solverComp->getComponents();
    dataarchiver->getComponents();
    simComp->getComponents();
    ctl->getComponents();

    ctl->run();

    dataarchiver->releaseComponents();
    sched->releaseComponents();
    lbc->releaseComponents();
    solverComp->releaseComponents();
    simComp->releaseComponents();
    ctl->releaseComponents();
    sched->removeReference();
    delete sched;

  } catch (Exception& e) {
    std::cout << e.message() << std::endl;
    throw;
  } catch (...) {
    std::cout << "**ERROR** Unknown exception" << std::endl;
    throw;
  }
}

int main(int argc, char** argv, char* env[])
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new DEMEnv(argc, argv, env));
  return RUN_ALL_TESTS();
}
