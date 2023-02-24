#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import multiprocessing
import subprocess as sub_proc

from CamClayTestSuite_PostProc import *

# Make sure python 3 is being used
if sys.version_info < (2, 8, 0):
  print("The test suite can only be run using Python 3.6 and newer.")

# Automated script for running and post processing CamClay verification tetsts.
# Requires that ../opt/StandAlone/ and ../opt/StandAlone/tools/extractors/ be on working PATH
# If you would like to run tests with more than one processor than mpirun must also be on PATH
# defaults to max number of processors or number of patches whichever is lower

# List for tests requiring restart with damping off and new end time
# ups end_time
RESTART_LIST = [[], []]

# Post processing list, ie tests that have a method to do their post processing
POST_PROCESS_LIST = [
    'CamClay_01_UniaxialStrainRotate.ups',
    'CamClay_02_VertexTreatment.ups',
    'CamClay_03_UniaxialStrainNoHardening.ups',
    'CamClay_04_CurvedYieldSurface.ups',
    'CamClay_05_HydrostaticCompressionFixedCap.ups',
    'CamClay_06_UniaxialStrainCapEvolution.ups',
    'CamClay_07_HydrostaticCompressionCapEvolution.ups',
    'CamClay_08_LoadingUnloading.ups',
    'CamClay_09_MultiaxialStrainLoadUnload.ups',
]

# get the current directory (will be working directory)
cur_dir = os.getcwd()

# get vaango/src path as environmental variable
vaango_src_path = os.path.abspath(os.environ['VAANGO_SRC'])

# The vaango executable, e.g., /home/banerjee/ParSim/Vaango/runs/vaango
# Typically a link to the executable created in the build directory is placed
# in the runs directory using, e.g.,
# ln -s /home/banerjee/ParSim/Vaango/dbg/StandAlone/vaango \
#       /home/banerjee/ParSim/Vaango/runs/vaango
# A link to partextract is also kept at the same place
vaango_exe = os.path.abspath(os.environ['VAANGO_EXE'])
partextract_exe = os.path.abspath(os.environ['PARTEXTRACT_EXE'])
print("vaango_exe", vaango_exe)

# construct default paths based on location of vaango_src_path
# print("cur_dir = ", cur_dir)
# print("vaango_src_path = ", vaango_src_path)
default_inputs_path = os.path.join(
    vaango_src_path, 'StandAlone/inputs/MPM/SoilPlasticityDamage/CamClay/')
default_working_dir = os.path.join(cur_dir, 'CamClay_test_runs')

# If the working directory does not exist then make it.
if not os.path.exists(default_working_dir):
  os.makedirs(default_working_dir)

# Make plots directory
default_plot_dir = os.path.join(default_working_dir, 'Plots')
if not os.path.exists(default_plot_dir):
  os.makedirs(default_plot_dir)

# Build list of tets in the default inputs directory
# TEST_LIST = []
# tmp_test_list = os.listdir(default_inputs_path)
# for item in tmp_test_list:
#  if os.path.isfile(item):
#    if '.ups' in item:
#      TEST_LIST.append(os.path.abspath(item))
# TEST_LIST.sort()

# for test in TEST_LIST:
#  print test

TEST_LIST = []
for test in POST_PROCESS_LIST:
  TEST_LIST.append(os.path.join(default_inputs_path, test))

### COMMENT ME OUT!!!!!!! ###
TEST_LIST = [
      #TEST_LIST[0], #Test 01
      #TEST_LIST[1], #Test 02
      #TEST_LIST[2], #Test 03
      #TEST_LIST[3], #Test 04
      #TEST_LIST[4], #Test 05
      #TEST_LIST[5], #Test 06
      #TEST_LIST[6], #Test 07
      #TEST_LIST[7], #Test 08
      TEST_LIST[8], #Test 08
]
### --------------------- ###

for test in TEST_LIST:
  print(test)


def copy_test_to(from_file, to_file):
  '''Reads in test ups file at from_file and writes to to_file while replacing
  absolute references to prescribed deformation files with relative ones. Also
  copys deformation file to same root folder.'''
  # Read in the from file and close
  F_from_file = open(from_file, "r")
  from_file_lines = F_from_file.read().split('\n')
  F_from_file.close()

  # Open/create the to file
  F_to_file = open(to_file, "w")
  to_file_root = os.path.split(to_file)[0]

  # Copy the ups but change the Prescribed def filebase also copy this file
  start_tag = '<prescribed_deformation_file>'
  end_tag = '</prescribed_deformation_file>'
  for line in from_file_lines:

    if start_tag in line and end_tag in line:
      def_file = line.split(start_tag)[1].split(end_tag)[0].strip()
      print("def_file_path = " + def_file)
      print("def_file_name = " + os.path.basename(def_file))
      def_file = os.path.basename(def_file)
      line = line.split(start_tag)[0] + start_tag + def_file + end_tag
      shutil.copyfile(os.path.join(default_inputs_path, def_file),
                      os.path.join(to_file_root, def_file))

    F_to_file.write(line + '\n')

  F_to_file.close()
  return os.path.abspath(to_file)


def setup_restart(uda_path, new_end_time):
  # Fix input file
  input_file = os.path.join(uda_path, 'input.xml')
  F = open(input_file, "r+")
  all_lines = F.read().split('\n')
  F.seek(0)
  for line in all_lines:
    if '<maxTime>' in line:
      line = '    <maxTime>' + format(new_end_time, '1.6e') + '</maxTime>'
    F.write(line + '\n')
  F.close()

  # Find checkpoint dirs and change damping
  for item in os.listdir(os.path.join(uda_path, "checkpoints")):
    tmp_file = os.path.join(uda_path, "checkpoints", item, "timestep.xml")
    if os.path.isfile(tmp_file):
      F = open(tmp_file, "r+")
      all_lines = F.read().split('\n')
      F.seek(0)
      for line in all_lines:
        if '<artificial_damping_coeff>' in line:
          line = '    <artificial_damping_coeff>0.0</artificial_damping_coeff>'
        F.write(line + '\n')
      F.close()


def run_test(ups_path,
             WITH_MPI=False,
             NUM_PROCS=1,
             RESTART=False,
             DAMPING_OFF_NEW_END_TIME=False,
             POST_PROC_ONLY=False):
  ''' '''
  print('\nRunning test:\t', os.path.split(ups_path)[1])

  # Determine root path
  root_path = os.path.split(os.path.abspath(ups_path))[0]

  # Determine uda path
  print("Root path = ", root_path)
  print("Ups path = ", ups_path)
  F_ups = open(ups_path, "r")
  ups_lines = F_ups.read()
  uda_path = os.path.join(
      root_path,
      ups_lines.split('<filebase>')[1].split('</filebase>')[0].strip())
  F_ups.close()
  print("UDA path = ", uda_path)

  # Change current working directory to root path
  os.chdir(root_path)
  # Open runlog
  F_log = open(
      os.path.join(root_path, 'TEST_RUNLOG_'+
                   os.path.split(ups_path)[1]), "w")

  # Construct the argument list for subprocess to use.
  if not (WITH_MPI) or int(NUM_PROCS) <= 1:
    print("vaango_exe", vaango_exe)
    args = [vaango_exe, os.path.split(ups_path)[1]]
  else:
    args = [
        'mpirun', '-np',
        str(int(NUM_PROCS)), vaango_exe, 
        os.path.split(ups_path)[1]
    ]

  if POST_PROC_ONLY:
    uda_path = uda_path + '.000'
  else:
    # Run the test and wait for it to complete
    print('Running: ', args)
    tmp = sub_proc.Popen(args, stdout=F_log, stderr=sub_proc.PIPE)
    dummy = tmp.wait()
    print('Run complete: ')

    F_log.close()
    # If test calls for retstart
    if RESTART:
      # If turn damping off and run to new end time
      if DAMPING_OFF_NEW_END_TIME:
        # Setup the restart by setting damping to zero and modifying end time
        print(
            'Setting <artificial_damping_coeff> to zero and restarting with new end time of ',
            format(NEW_END_TIME, '1.4e'))
        setup_restart(uda_path, DAMPING_OFF_NEW_END_TIME)
        print('Done.\nRestarting...')
        # Open new runlog
        F_log = open(
            os.path.join(root_path, 'TEST_RUNLOG_RESTART_',
                         os.path.split(ups_path)[1]), "w")
        # Construct the argument list
        if not (WITH_MPI) or NUM_PROCS <= 1:
          args = [vaango_exe, '-restart', '-move', uda_path + '.000']
        else:
          args = [
              'mpirun', '-np',
              str(int(NUM_PROCS)), vaango_exe, '-restart', '-move',
              uda_path + '.000'
          ]
        # Run the test and wait for it to complete
        tmp = sub_proc.Popen(args, stdout=F_log, stderr=sub_proc.PIPE)
        dummy = tmp.wait()
        F_log.close()
        uda_path = uda_path + '.001'
    else:
      uda_path = uda_path + '.000'

  print('Test done.')
  print(os.path.dirname(os.path.dirname(uda_path)))
  os.chdir(os.path.dirname(os.path.dirname(uda_path)))
  return uda_path


def clear_uda(uda_path):
  print('Deleting uda...')
  tmp = sub_proc.Popen(['rm', '-rf', uda_path],
                       stdout=sub_proc.PIPE,
                       stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  tmp = sub_proc.Popen(['rm', uda_path.split('.uda')[0] + '.uda'],
                       stdout=sub_proc.PIPE,
                       stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done')


def run_all_tests(TEST_METHODS=False, CLEAR_UDA=False, POST_PROC_ONLY=False):
  global default_working_dir, TEST_LIST, RESTART_LIST, MPI_FLAG, NUM_CPUS, POST_PROCESS_LIST
  print('#-- Running All Tests --#')
  for test in TEST_LIST:
    # Copy to working directory
    ups_path = copy_test_to(
        test, os.path.join(default_working_dir,
                           os.path.split(test)[1]))
    # Run
    if test not in RESTART_LIST[0]:
      uda_path = run_test(ups_path,
                          WITH_MPI=MPI_FLAG,
                          NUM_PROCS=NUM_CPUS,
                          POST_PROC_ONLY=POST_PROC_ONLY)
    else:
      new_end_time = RESTART_LIST[1][RESTART_LIST[0].index(test)]
      uda_path = run_test(ups_path,
                          WITH_MPI=MPI_FLAG,
                          NUM_PROCS=NUM_CPUS,
                          DAMPING_OFF_NEW_END_TIME=new_end_time,
                          POST_PROC_ONLY=POST_PROC_ONLY)
    # Post process if called for
    if TEST_METHODS:
      test_yield_surface(uda_path)
    else:
      post_proc(test, uda_path, default_plot_dir)
    # Clean up the uda
    if CLEAR_UDA:
      clear_uda(uda_path)


def post_proc(test, uda_path, save_path):
  global POST_PROCESS_LIST
  test_name = os.path.split(test)[1]
  if test_name in POST_PROCESS_LIST:
    if test_name == 'CamClay_01_UniaxialStrainRotate.ups':
      test01_postProc(uda_path, save_path)
    if test_name == 'CamClay_02_VertexTreatment.ups':
      test02_postProc(uda_path, save_path)
    if test_name == 'CamClay_03_UniaxialStrainNoHardening.ups':
      test03_postProc(uda_path, save_path)
    if test_name == 'CamClay_04_CurvedYieldSurface.ups':
      test04_postProc(uda_path, save_path)
    if test_name == 'CamClay_05_HydrostaticCompressionFixedCap.ups':
      test05_postProc(uda_path, save_path)
    if test_name == 'CamClay_06_UniaxialStrainCapEvolution.ups':
      test06_postProc(uda_path, save_path)
    if test_name == 'CamClay_07_HydrostaticCompressionCapEvolution.ups':
      test07_postProc(uda_path, save_path)
    if test_name == 'CamClay_08_LoadingUnloading.ups':
      test08_postProc(uda_path, save_path)
    if test_name == 'CamClay_09_MultiaxialStrainLoadUnload.ups':
      test09_postProc(uda_path, save_path)
  else:
    print('\nERROR: test: ', test, '\n\tNot on post processing list.\n')


if __name__ == "__main__":
  # Determine number of CPU cores
  NUM_CPUS = multiprocessing.cpu_count()
  # Set MPI Flag intially to false
  MPI_FLAG = False

  # Not setup yet. Need to scan ups to determine # patches
  if False:
    ABORT = False
    if len(sys.argv) == 3:
      if sys.argcv[1] == '-mpirun':
        MPI_FLAG = True
        try:
          NUM_CPUS = int(sys.argv[2])
        except:
          NUM_CPUS = 1
          print(
              '\nError: invalid number of processors entered with -mpirun flag')
      else:
        print(
            '\nInvalid Arguments entered!\n\tuse: CamClaySuite.py -mpirun <# processor cores>'
        )
        ABORT = True
    else:
      not_done = True
      while not_done:
        mpi_check = raw_input(
            "Would you like to run using mpirun? (Y/N/(A)bort)\n").lower()
        if mpi_check == 'y':
          not_done_2 = True
          MPI_FLAG = True
          while not_done_2:
            try:
              num_cores = int(
                  raw_input("Enter the number of cores to run on:\n"))
              NUM_CPUS = num_cores
              not_done_2 = False
              not_done = False
            except:
              print('Invalid entry. Please try again.')
        elif mpi_check == 'n':
          not_done = False
        elif mpi_check == 'a':
          not_done = False
          ABORT = True
      if not (ABORT):
        if multiprocessing.cpu_count() < NUM_CPUS:
          NUM_CPUS = multiprocessing.cpu_count()
          print(
              '\nWarning: number of cores requested more than are available locally.\n\t# cores set to: ',
              NUM_CPUS)
        print(' ')
        run_all_tests()

  TEST_METHODS = False
  POST_PROC_ONLY = True
  #POST_PROC_ONLY = False
  #CLEAR_UDA = True
  CLEAR_UDA = False
  run_all_tests(TEST_METHODS, CLEAR_UDA, POST_PROC_ONLY)
