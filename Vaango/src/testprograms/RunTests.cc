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

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

/*********************************************
 * RunTests.cc
 * Wayne Witzel
 * SCI Institute,
 * University of Utah
 *
 * Program for running tests and reporting them
 * in a tree of test suites.
 *
 * Add methods to populate this suite tree with more
 * suites below where it says:
 * ADD MORE POPULATING METHODS ABOVE FOR OTHER TEST SUITES
 */

#include <testprograms/TestSuite/SuiteTree.h>
#include <testprograms/TestMatrix3/testmatrix3.h>
#include <testprograms/TestConsecutiveRangeSet/TestConsecutiveRangeSet.h>
#include <testprograms/TestRangeTree/TestRangeTree.h>
#include <testprograms/TestBoxGrouper/TestBoxGrouper.h>

#include <cstdlib>
#include <iostream>

#ifndef _WIN32
#  include <unistd.h>
#else
#  include <process.h>
#endif


using namespace Uintah;

void
usage(char* prog_name)
{
  std::cerr <<  "usage: " << prog_name << " [-e|-a|-h]\n";
  std::cerr <<  "\t-e:  expands test suite tree even where all tests have passed\n";
  std::cerr <<  "\t-a:  reports all suites (not just failed ones)\n";
  std::cerr <<  "\t-v:  verbose mode\n";
  std::cerr <<  "\t-h:  lists this help information\n";
}

int
main(int argc, char* argv[])
{
  bool expandAll = false;
  bool reportAll = false;
  bool verbose = false;
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      usage(argv[0]);
      return 1;
    }
    else {
      int j = 1;
      while (argv[i][j] != '\0') {
	switch (argv[i][j]) {
	case 'a':
	  reportAll = true;
	  break;
	case 'e':
	  expandAll = true;
	  break;
	case 'h':
	  usage(argv[0]);
	  return 0;
	case 'v':
	  verbose = true;
	  break;
	default:
	  std::cerr <<  "unkown option: " << argv[i][j] << std::endl;
	  usage(argv[0]);
	  return 1;
	}
	j++;
      }
    }
  }
  
  srand(getpid());
  SuiteTreeNode* suites = new SuiteTreeNode("All Tests");

  // populate the suites tree
  suites->addSubTree(matrix3TestTree());
  suites->addSubTree(ConsecutiveRangeSetTestTree());
  suites->addSubTree(RangeTreeTestTree(verbose, 20000));
  suites->addSubTree(BoxGrouperTestTree(verbose));

  /* ADD MORE POPULATING METHODS ABOVE FOR OTHER TEST SUITES */

  // show suite tree summary
  bool allPassed;
  std::cout << suites->composeSubSummary("", expandAll, allPassed) << std::endl;
  
  // report failed suites to itemized the tests the failed
  if (reportAll)
    suites->reportAllSuites();
  else if (!allPassed) {
     std::list<Suite*> failedSuites;
    suites->appendFailedSuites(failedSuites);
    
    std::cout << "Failed Suite Reports:" << std::endl << std::endl;
    for (list<Suite*>::iterator it = failedSuites.begin();
	 it != failedSuites.end(); it++) {
      (*it)->report();
      std::cout << std::endl;
    }

    if (failedSuites.size() == 1)
      std::cout << "\n1 suite failed.\n";
    else
      std::cout << std::endl << failedSuites.size() << " suites failed.\n";
  }

  if (allPassed) {
    std::cout << "All tests passed!\n";
    return 1;
  }

  return 0;
}

