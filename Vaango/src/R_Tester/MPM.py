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

from sys import argv,exit
from os import environ
from helpers.runVaangoTests import runVaangoTests

#______________________________________________________________________
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

NIGHTLYTESTS = [  
                  #("Centrifuge_2D_AGR_SimPBC_dense_layer_lores_loforce", "SoilPlasticityDamage/RotatingCoords/Centrifuge_2D_AGR_SimPBC_dense_layer_lores_loforce.ups",                       4,  "Linux", ["exactComparison"] ), \
                  #("Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained_delay_offset_nobucket_initstress", "SoilPlasticityDamage/Centrifuge/Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained_delay_offset_nobucket_initstress.ups",                       4,  "Linux", ["exactComparison"] ), \
#                  ("DropBunny2Frags", "SoilPlasticityDamage/BunnyHummer/DropBunny2Frags.ups",                       1,  "Linux", ["exactComparison"] ), \
  ("const_test_viscoelastic_fortran.ups", 
   "ViscoElastic/const_test_viscoelastic_fortran.ups",                       
    1,  "Linux", ["exactComparison"] ), \
                  ("const_test_brittle_damage", "SoilPlasticityDamage/const_test_brittle_damage.ups",                       1,  "Linux", ["exactComparison"] ), \
                  ("AreniscaTest_01a_UniaxialStrainCompressExtend", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_01a_UniaxialStrainCompressExtend.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_01_UniaxialStrainRotate", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_01_UniaxialStrainRotate.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_02_VertexTreatment", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_02_VertexTreatment.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_03_UniaxialStrainNoHardening", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_03_UniaxialStrainNoHardening.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_04_CurvedYieldSurface", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_04_CurvedYieldSurface.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_05_HydrostaticCompressionFixedCap", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_05_HydrostaticCompressionFixedCap.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_06_UniaxialStrainCapEvolution", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_06_UniaxialStrainCapEvolution.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_07_HydrostaticCompressionCapEvolution", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_07_HydrostaticCompressionCapEvolution.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_08_LoadingUnloading", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_08_LoadingUnloading.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_09_FluidFilledPoreSpace", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_09_FluidFilledPoreSpace.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_10_PureIsochoricStrainRates", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_10_PureIsochoricStrainRates.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_10_TransientStressEigenvaluesConstVectors", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_10_TransientStressEigenvaluesConstVectors.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_11_UniaxialStrainJ2plasticity", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_11_UniaxialStrainJ2plasticity.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_12_NonlinearElasticity", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_12_NonlinearElasticity.ups", 1, "Linux", ["exactComparison"] ), \
                  #("AreniscaTest_13_UniaxialStrainRateDependence", "SoilPlasticityDamage/Arenisca/Arenisca3/AreniscaTest_13_UniaxialStrainRateDependence.ups", 1, "Linux", ["exactComparison"] ), \
                  ("disks_complex",                       "disks_complex.ups",                       4,  "Linux", ["exactComparison"] ), \
                  ("heatcond2mat",                        "heatcond2mat.ups",                        1,  "Linux", ["exactComparison"] ),  \
                  ("inclined_plane_sphere",               "inclined_plane_sphere.ups",               1,  "Linux", ["exactComparison"] ),  \
                  #("foam_crush",                          "foam_crush.ups",                          4,  "Linux", ["exactComparison"] ),  \
                  #("periodic_disks",                      "periodic_disks.ups",                      1,  "Linux", ["exactComparison"] ),  \
                  #("periodic_spheres3D",                  "periodic_spheres3D.ups",                  8,  "Linux", ["exactComparison"] ),  \
                  ("const_test_hypo",                     "const_test_hypo.ups",                     1,  "Linux", ["exactComparison"] ),  \
                  ("const_test_cmr",                      "const_test_cmr.ups",                      1,  "Linux", ["exactComparison"] ),  \
                  ("const_test_nhp",                      "const_test_nhp.ups",                      1,  "Linux", ["exactComparison"] ),  \
                  ("const_test_vs",                       "const_test_vs.ups",                       1,  "Linux", ["exactComparison"] ),  \
                  ("adiCuJC4000s696K",                    "adiCuJC4000s696K.ups",                    1,  "Linux", ["exactComparison"] ),  \
                  #("adiCuMTS4000s696K",                   "adiCuMTS4000s696K.ups",                   1,  "Linux", ["exactComparison"] ),  \
                  #("adiCuPTW4000s696K",                   "adiCuPTW4000s696K.ups",                   1,  "Linux", ["exactComparison"] ),  \
                  #("adiCuSCG4000s696K",                   "adiCuSCG4000s696K.ups",                   1,  "Linux", ["exactComparison"] ),  \
                  #("adiCuZA4000s696K",                    "adiCuZA4000s696K.ups",                    1,  "Linux", ["exactComparison"] ),  \
                  #("test_corrug_plate",                   "test_corrug_plate.ups",                   1,  "Linux", ["exactComparison"] ),  \
                  #("test_cyl_pene_no_ero",                "test_cyl_pene_no_ero.ups",                1,  "Linux", ["exactComparison"] ),  \
                  #("test_gurson_beckerdrucker_mts",       "test_gurson_beckerdrucker_mts.ups",       1,  "Linux", ["exactComparison"] ),  \
                  #("test_hypoviscoelastic_radial_return", "test_hypoviscoelastic_radial_return.ups", 1,  "Linux", ["exactComparison"] ),  \
                  #("Charpy",                              "Charpy.ups",                              8,  "Linux", ["exactComparison"] ),  \
                  #("disks_complex",                       "disks_complex.ups",                       4,  "Darwin", ["doesTestRun"]    ),     \
                  #("heatcond2mat",                        "heatcond2mat.ups",                        1,  "Darwin", ["doesTestRun"]    ),     \
                  #("inclined_plane_sphere",               "inclined_plane_sphere.ups",               1,  "Darwin", ["doesTestRun"]    ),     \
                  #("const_test_cmr",                      "const_test_cmr.ups",                      1,  "Darwin", ["doesTestRun"]    ),     \
                  #("const_test_nhp",                      "const_test_nhp.ups",                      1,  "Darwin", ["doesTestRun"]    ),     \
                  #("adiCuJC4000s696K",                    "adiCuJC4000s696K.ups",                    1,  "Darwin", ["doesTestRun"]    ),     \
                  #("adiCuMTS4000s696K",                   "adiCuMTS4000s696K.ups",                   1,  "Darwin", ["doesTestRun"]    ),     \
                  #("adiCuPTW4000s696K",                   "adiCuPTW4000s696K.ups",                   1,  "Darwin", ["doesTestRun"]    ),     \
                  #("adiCuSCG4000s696K",                   "adiCuSCG4000s696K.ups",                   1,  "Darwin", ["doesTestRun"]    ),     \
                  #("adiCuZA4000s696K",                    "adiCuZA4000s696K.ups",                    1,  "Darwin", ["doesTestRun"]    ),     \
                  #("test_corrug_plate",                   "test_corrug_plate.ups",                   1,  "Darwin", ["doesTestRun"]    ),     \
                  #("test_cyl_pene_no_ero",                "test_cyl_pene_no_ero.ups",                1,  "Darwin", ["doesTestRun"]    ),     \
                  #("test_gurson_beckerdrucker_mts",       "test_gurson_beckerdrucker_mts.ups",       1,  "Darwin", ["doesTestRun"]    ) 
            ]
              
# Tests that are run during local regression testing              
LOCALTESTS = NIGHTLYTESTS

#__________________________________

def getNightlyTests() :
  return NIGHTLYTESTS

def getLocalTests() :
  return LOCALTESTS

#__________________________________

if __name__ == "__main__":

  if environ['WHICH_TESTS'] == "local":
    TESTS = LOCALTESTS
  else:
    TESTS = NIGHTLYTESTS

  result = runVaangoTests(argv, TESTS, "MPM")
  exit( result )

