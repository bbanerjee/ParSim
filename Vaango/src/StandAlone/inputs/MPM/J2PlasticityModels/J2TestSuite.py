#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import multiprocessing
import subprocess as sub_proc

from J2TestSuite_PostProc import *

#Automated script for running and post processing Arenisca verification tetsts.
#Post processing list, ie tests that have a method to do their post processing
POST_PROCESS_LIST = [
    'J2Test_01_UniaxialStrainRotateJ2Lin.ups',
    'J2Test_01a_UniaxialStressJ2Lin.ups',
    'J2Test_01b_UniaxialStrainLoadUnloadJ2Lin.ups',
    'J2Test_02_UniaxialStrainLoadUnloadJ2NonLin.ups',
    'J2Test_03_UniaxialStrainRotateDPLin.ups',
    'J2Test_03a_UniaxialStrainLoadUnloadDPLin.ups',
    'J2Test_04_UniaxialStrainLoadUnloadDPNonLin.ups',
    'J2Test_05_UniaxialStrainLoadUnloadNonLinDPNonLin.ups',
    'J2Test_06_HydrostaticCompression.ups',
    'J2Test_07_HydrostaticLoadUnload.ups',
    'J2Test_08_UniaxialStrainCompresson.ups',
    'J2Test_09_UniaxialStrainTension.ups',
    'J2Test_10_UniaxialStrainRotate.ups',
    'J2Test_11_TriaxialStrainTension.ups',
    'J2Test_12_UniaxialStrainLoadUnload.ups',
    'J2Test_13_MultiaxialStrainLoadUnload.ups',
]

#get uintah/src path as enviornmental variable
#uintah_src_path = os.path.abspath(os.environ['UINTAH_SRC'])
uintah_src_path = os.path.abspath(".")

# The uintah executable, e.g., /home/banerjee/ParSim/Vaango/runs/vaango
# Typically a link to the executable created in the build directory is placed
# in the runs directory using, e.g.,
# ln -s ~/ParSim/Vaango/dbg/StandAlone/vaango \
#       ~/ParSim/Vaango/runs/vaango
# A link to partextract is also kept at the same place
#vaango_exe = os.path.abspath(os.environ['UINTAH_EXE'])
#partextract_exe = os.path.abspath(os.environ['PARTEXTRACT_EXE'])
vaango_exe = os.path.abspath("../../vaango_opt")
partextract_exe = os.path.abspath("../../partextract")

#construct default paths based on location of uintah_src_path
default_inputs_path = uintah_src_path
default_working_dir = uintah_src_path + '/test_run'
#If the working directory does not exist then make it.
if not os.path.exists(default_working_dir):
  os.makedirs(default_working_dir)

#Make plots directory
default_plot_dir = default_working_dir + '/Plots'
if not os.path.exists(default_plot_dir):
  os.makedirs(default_plot_dir)

TEST_LIST = []
for test in POST_PROCESS_LIST:
  TEST_LIST.append(default_inputs_path + '/' + test)

#
# Uncomment tests yo wish to run
#
TEST_LIST = [
    #  TEST_LIST[0], #Test 01
    #  TEST_LIST[1], #Test 01a
    #  TEST_LIST[2], #Test 01b
    #  TEST_LIST[3], #Test 02
    #  TEST_LIST[4], #Test 03
    #  TEST_LIST[5], #Test 03a
    #  TEST_LIST[6], #Test 04
    #  TEST_LIST[7], #Test 05
    #  TEST_LIST[8], #Test 06
    #  TEST_LIST[9], #Test 07
    #  TEST_LIST[10], #Test 08
    #  TEST_LIST[11], #Test 09
    #  TEST_LIST[12], #Test 10
    #  TEST_LIST[13], #Test 11
    #  TEST_LIST[14], #Test 12
    #  TEST_LIST[15], #Test 13
]
### --------------------- ###

for test in TEST_LIST:
  print(test)


def copy_test_to(from_file, to_file):
  '''Reads in test ups file at from_file and writes to to_file while replacing
absolute references to prescribed deformation files with relative ones. Also
copys deformation file to same root folder.'''
  #Read in the from file and close
  F_from_file = open(from_file, "r")
  from_file_lines = F_from_file.read().split('\n')
  F_from_file.close()

  #Open/create the to file
  F_to_file = open(to_file, "w")
  to_file_root = os.path.split(to_file)[0]

  #Copy the ups but change the Prescribed def filebase also copy this file
  start_tag = '<prescribed_deformation_file>'
  end_tag = '</prescribed_deformation_file>'
  for line in from_file_lines:

    if start_tag in line and end_tag in line:
      def_file = line.split(start_tag)[1].split(end_tag)[0].strip()
      def_file = os.path.basename(def_file)
      line = line.split(start_tag)[0] + start_tag + def_file + end_tag
      shutil.copyfile(os.path.join(default_inputs_path, def_file),
                      os.path.join(to_file_root, def_file))

    F_to_file.write(line + '\n')

  F_to_file.close()
  return os.path.abspath(to_file)


def setup_restart(uda_path, new_end_time):
  #Fix input file
  input_file = uda_path + '/input.xml'
  F = open(input_file, "r+")
  all_lines = F.read().split('\n')
  F.seek(0)
  for line in all_lines:
    if '<maxTime>' in line:
      line = '    <maxTime>' + format(new_end_time, '1.6e') + '</maxTime>'
    F.write(line + '\n')
  F.close()

  #Find checkpoint dirs and change damping
  for item in os.listdir(uda_path + "/checkpoints/"):
    tmp_file = uda_path + "/checkpoints/" + item + "/timestep.xml"
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

  #Determine root path
  root_path = os.path.split(os.path.abspath(ups_path))[0]

  #Determine uda path
  print("Root path = ", root_path)
  print("Ups path = ", ups_path)
  F_ups = open(ups_path, "r")
  ups_lines = F_ups.read()
  uda_path = os.path.join(
      root_path,
      ups_lines.split('<filebase>')[1].split('</filebase>')[0].strip())
  F_ups.close()
  print("UDA path = ", uda_path)

  if POST_PROC_ONLY:
    uda_path = uda_path
  else:
    #Change current working directory to root path
    os.chdir(root_path)

    #Open runlog
    F_log = open(
        os.path.join(root_path, 'TEST_RUNLOG_' + os.path.split(ups_path)[1]),
        "w")

    #Construct the argument list for subprocess to use.
    if not (WITH_MPI) or int(NUM_PROCS) <= 1:
      args = [vaango_exe, os.path.split(ups_path)[1]]
    else:
      args = [
          'mpirun', '-np',
          str(int(NUM_PROCS)), vaango_exe,
          os.path.split(ups_path)[1]
      ]

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
  global default_working_dir, TEST_LIST, MPI_FLAG, NUM_CPUS, POST_PROCESS_LIST
  print('#-- Running All Tests --#')
  for test in TEST_LIST:

    # Copy to working directory
    ups_path = copy_test_to(
        test, os.path.join(default_working_dir,
                           os.path.split(test)[1]))

    #Run
    uda_path = run_test(ups_path,
                        WITH_MPI=MPI_FLAG,
                        NUM_PROCS=NUM_CPUS,
                        POST_PROC_ONLY=POST_PROC_ONLY)

    #Post process if called for
    if TEST_METHODS:
      test_yield_surface(uda_path)
    else:
      post_proc(test, uda_path, default_plot_dir, POST_PROCESS_LIST)

    #Clean up the uda
    if CLEAR_UDA:
      clear_uda(uda_path)


if __name__ == "__main__":

  MPI_FLAG = False
  NUM_CPUS = 1
  TEST_METHODS = False
  POST_PROC_ONLY = False
  if len(sys.argv) == 2:
    if sys.argv[1] == '-post':
      POST_PROC_ONLY = True

  #CLEAR_UDA = True
  CLEAR_UDA = False
  run_all_tests(TEST_METHODS, CLEAR_UDA, POST_PROC_ONLY)
