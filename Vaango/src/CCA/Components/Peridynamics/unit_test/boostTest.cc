/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

//#include <CCA/Components/Peridynamics/unit_test/utBlank.h>
#include <CCA/Components/Peridynamics/unit_test/utParticleTest.h>

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/Grid/MaterialManagerP.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Components/MPM/Contact/Contact.h>

#include <Core/Parallel/UintahParallelComponent.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InvalidGrid.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Tracker/TrackerClient.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/SimulationController/AMRSimulationController.h>
#include <CCA/Components/Models/ModelFactory.h>
#include <CCA/Components/Solvers/CGSolver.h>
#include <CCA/Components/Solvers/DirectSolve.h>

#include <CCA/Components/ReduceUda/UdaReducer.h>
#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <CCA/Components/Solvers/SolverFactory.h>
#include <CCA/Components/Regridder/RegridderFactory.h>
#include <CCA/Components/LoadBalancers/LoadBalancerFactory.h>
#include <CCA/Components/Schedulers/SchedulerFactory.h>
#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/Exception.h>
#include <Core/Exceptions/InternalError.h>

#include <Core/Thread/Time.h>

#include <Core/Util/DebugStream.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>

#include <sci_defs/hypre_defs.h>
#include <sci_defs/malloc_defs.h>
#include <sci_defs/mpi_defs.h>
#include <sci_defs/uintah_defs.h>
#include <sci_defs/cuda_defs.h>

#include <Core/Malloc/Allocator.h>

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>
#include <sys/stat.h>

#include <time.h>
#include <iomanip>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

using namespace boost::unit_test;
using namespace Vaango;

extern Uintah::Mutex cerrLock;
static Uintah::DebugStream stackDebug("ExceptionStack", true);

void abortCleanupFunc() {
  Uintah::Parallel::finalizeManager( Uintah::Parallel::Abort );
}

void runTest(const string filename) {
//    string filename = "inputs/utParticleTest.ups";
    char * start_addr = (char*)sbrk(0);
    bool do_AMR=false;
    int  udaSuffix = -1;
    string udaDir; 
    bool validateUps = true;
    
    // Read input file
    Uintah::ProblemSpecP ups = Uintah::ProblemSpecReader().readInputFile( filename, validateUps );
    const Uintah::ProcessorGroup* world = Uintah::Parallel::getRootProcessorGroup();
    Uintah::SimulationController* ctl = scinew Uintah::AMRSimulationController(world, do_AMR, ups);
    Uintah::RegridderCommon* reg = 0;
    Uintah::SolverInterface* solve = 0;

    // Create the components
    Uintah::UintahParallelComponent* comp = scinew utParticleTest(world);
    Uintah::SimulationInterface* sim = dynamic_cast<Uintah::SimulationInterface*>(comp);

    ctl->attachPort("sim", sim);
    comp->attachPort("solver", solve);
    comp->attachPort("regridder", reg);

    // Load balancer
    Uintah::LoadBalancerCommon* lbc = Uintah::LoadBalancerFactory::create(ups, world);
    lbc->attachPort("sim", sim);
    
    // Output
    Uintah::DataArchiver* dataarchiver = scinew Uintah::DataArchiver(world, udaSuffix);
    Uintah::Output* output = dataarchiver;
    ctl->attachPort("output", dataarchiver);
    dataarchiver->attachPort("load balancer", lbc);
    comp->attachPort("output", dataarchiver);
    dataarchiver->attachPort("sim", sim);
    
    // Scheduler
    Uintah::SchedulerCommon* sched = Uintah::SchedulerFactory::create(ups, world, output);
    sched->attachPort("load balancer", lbc);
    ctl->attachPort("scheduler", sched);
    lbc->attachPort("scheduler", sched);
    comp->attachPort("scheduler", sched);

    sched->setStartAddr( start_addr );
    sched->addReference();
    
    ups = 0;

    // Run
    ctl->run();
    delete ctl;

    sched->removeReference();
    delete sched;
    delete lbc;
    delete sim;
    delete output;
}

bool init_unit_test() {
    bool success = true;
    string params[] = {"inputs/utParticleTest.ups"};
    test_suite* proc = BOOST_TEST_SUITE("boostTest0");
    proc->add(BOOST_PARAM_TEST_CASE(&runTest, params, params+1));    
    framework::master_test_suite().add( proc );
    return success;
} 

int main(int argc, char *argv[], char *env[]) {
    int returnCode = -1;
    Uintah::Thread::setDefaultAbortMode("exit");
    Uintah::Thread::self()->setCleanupFunction( &abortCleanupFunc );
    //Uintah::Parallel::determineIfRunningUnderMPI( argc, argv );
    Uintah::create_sci_environment( env, 0, true );
    bool thrownException = false;

    try {
        Uintah::Parallel::initializeManager( argc, argv );
        returnCode = unit_test_main(init_unit_test, argc, argv);
    } 
    catch (Uintah::Exception& e) {
        cerrLock.lock();
        std::cout << "\n\n" << Uintah::Parallel::getMPIRank() << " Caught exception: " << e.message() << "\n\n";
        if(e.stackTrace())
            stackDebug << "Stack trace: " << e.stackTrace() << '\n';
        cerrLock.unlock();
        thrownException = true;
    }
  
    Uintah::TypeDescription::deleteAll();
    Uintah::Parallel::finalizeManager(thrownException ? Uintah::Parallel::Abort : Uintah::Parallel::NormalShutdown);

    if (thrownException) {
        if( Uintah::Parallel::getMPIRank() == 0 ) {
            std::cout << "\n\nAN EXCEPTION WAS THROWN... Goodbye.\n\n";
        }
      Uintah::Thread::exitAll(1);
    }
  
    if( Uintah::Parallel::getMPIRank() == 0 ) {
        std::cout << "Sus: going down successfully\n";
    }

    Uintah::Thread::exitAll(0);
    return(returnCode);
} // end main()
