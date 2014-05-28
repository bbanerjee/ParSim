/*! \file utTest0.cpp */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <vector>

using namespace boost::unit_test;

void Test0() {
    BOOST_CHECK_CLOSE(1.0, 1.0, 1.0e-12);
    BOOST_CHECK_CLOSE(0.0, 1.0, 1.0e-12);
}


bool init_unit_test_suite() {
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utTest0");
  proc->add(BOOST_TEST_CASE(&Test0));
  framework::master_test_suite().add(proc);

  return success;
}

bool init_unit_test() {
  return init_unit_test_suite();
}

int main (int argc, char* argv[]) {
    int numProcs = 1;
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

    int returnCode = -1;
    returnCode = unit_test_main(init_unit_test, argc, argv);
  
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  
    return returnCode;
}
