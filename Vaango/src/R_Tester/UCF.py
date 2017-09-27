#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2017 Parresia Research Limited, New Zealand
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

#!/usr/bin/env python

from os import chdir,getcwd,mkdir,system,environ
from sys import argv,exit,platform
from helpers.runVaangoTests import runVaangoTests, inputs_root
from helpers.modUPS import modUPS

#______________________________________________________________________
#Performance and other UCF related tests

#  Test syntax: ( "folder name", "input file", # processors, "OS", ["flags1","flag2"])
#  flags: 
#       gpu:                    - run test if machine is gpu enabled
#       no_uda_comparison:      - skip the uda comparisons
#       no_memoryTest:          - skip all memory checks
#       no_restart:             - skip the restart tests
#       no_dbg:                 - skip all debug compilation tests
#       no_opt:                 - skip all optimized compilation tests
#       do_performance_test:    - Run the performance test, log and plot simulation runtime.
#                                 (You cannot perform uda comparsions with this flag set)
#       doesTestRun:            - Checks if a test successfully runs
#       abs_tolerance=[double]  - absolute tolerance used in comparisons
#       rel_tolerance=[double]  - relative tolerance used in comparisons
#       exactComparison         - set absolute/relative tolerance = 0  for uda comparisons
#       startFromCheckpoint     - start test from checkpoint. (/home/csafe-tester/CheckPoints/..../testname.uda.000)
#       vaango_options="string"    - Additional command line options for vaango command
#
#  Notes: 
#  1) The "folder name" must be the same as input file without the extension.
#  2) If the processors is > 1.0 then an mpirun command will be used
#  3) Performance_tests are not run on a debug build.
#______________________________________________________________________

UNUSED = [ ("dlb_mpm",        "mpmDLBTest.ups",     8, "NONE"), \
           ("dlb_mpmice",     "mpmiceDLBTest.ups",  8, "NONE") 
]

NIGHTLYTESTS = [ ("ice_perf_test",          "icePerformanceTest.ups",             1, "Linux", ["do_performance_test"]),  \
                 ("mpmice_perf_test",       "mpmicePerformanceTest.ups",          1, "Linux", ["do_performance_test"]), \
                 ("LBwoRegrid",             "LBwoRegrid.ups",                     2, "Linux", ["exactComparison"]), \
                 ("switchExample_impm_mpm", "Switcher/switchExample_impm_mpm.ups",1, "Linux", ["no_memoryTest"]), \
                 ("switchExample3",         "Switcher/switchExample3.ups",        1, "Linux", ["no_restart","no_memoryTest"]), \
               ]

LOCALTESTS = [ #("switchExample_impm_mpm", "Switcher/switchExample_impm_mpm.ups",1, "Linux", ["no_memoryTest"]), \
               #("switchExample3",         "Switcher/switchExample3.ups",        1, "Linux", ["no_restart","no_memoryTest"]), \
               #("ice_perf_test",          "icePerformanceTest.ups",             1, "Linux", ["do_performance_test"]),  \
               #("mpmice_perf_test",       "mpmicePerformanceTest.ups",          1, "Linux", ["do_performance_test"]), \
               #("LBwoRegrid",             "LBwoRegrid.ups",                     2, "Linux", ["exactComparison"])
             ]

#__________________________________

def getNightlyTests() :
  return NIGHTLYTESTS

def getLocalTests() :
  return LOCALTESTS

#__________________________________

if __name__ == "__main__":

  print "starting UCF tests"
  if environ['WHICH_TESTS'] == "local":
    TESTS = LOCALTESTS
  else:
    TESTS = NIGHTLYTESTS

  result = runVaangoTests(argv, TESTS, "UCF")
  exit( result )
  
