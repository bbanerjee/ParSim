#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import multiprocessing
import subprocess as sub_proc  
import tempfile
import argparse
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib import ticker

#--------------------------------------------------------------------------
# Create a Box class and associated methods
# Alternatively: Use a named tuple : from collections import namedtuple
#--------------------------------------------------------------------------
class Box():

  def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax):
    self.xmin = xmin
    self.ymin = ymin
    self.zmin = zmin
    self.xmax = xmax
    self.ymax = ymax
    self.zmax = zmax

  def __repr__(self):
    return "Box()"

  def __str__(self):
    return ("["+str(self.xmin)+","+str(self.ymin)+","+str(self.zmin)+","
      +str(self.xmax)+","+str(self.ymax)+","+str(self.zmax)+"]")

  def invalid(self):
    if (self.xmax <= self.xmin or self.ymax <= self.ymin or self.zmax <= self.zmin):
      return True
    return False

  def update(self, xmin, ymin, zmin, xmax, ymax, zmax):
    self.xmin = xmin
    self.ymin = ymin
    self.zmin = zmin
    self.xmax = xmax
    self.ymax = ymax
    self.zmax = zmax

#--------------------------------------------------------------------------
# Main function (parse arguments)
#--------------------------------------------------------------------------
def main(argv):

  uda_name = ''
  selection_box = Box(0,0,0,0,0,0) 
  material_id = 0 
  timestep_lo = 0
  timestep_hi = 0

  # Parse comman line arguments
  parser = argparse.ArgumentParser(
    usage="extract_centrifuge_pos_data --box xmin ymin zmin xmax ymax zmax --mat <mat id> --timestep <int> --uda <uda file_name_prefix>")
  parser.add_argument('--box', nargs=6, help = 'xmin ymin zmin xmax ymax zmax', type=float, required=True)
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--timestep', default=0, type=int, required=False)
  parser.add_argument('--uda', required=True)
  args = parser.parse_args()

  # Parse box
  vals = args.box
  xmin = float(vals[0])
  ymin = float(vals[1])
  zmin = float(vals[2])
  xmax = float(vals[3])
  ymax = float(vals[4])
  zmax = float(vals[5])
  selection_box.update(xmin, ymin, zmin, xmax, ymax, zmax)
  if selection_box.invalid():
    print "**ERROR** Invalid selction box."
    sys.exit()

  # Parse mat
  material_id = int(args.mat)

  # Parse time step lo/hi
  timestep = int(args.timestep)

  # Parse uda
  uda_name = args.uda
  
  return uda_name, selection_box, material_id, timestep

#--------------------------------------------------------------------------
# Define the executables to be used and the working directories 
# Two executables are needed - selectpart and extractPosVelMasVol
# selectpart:  selects a list of particleIDs inside a selection box
# extractPosVelMasVol: extracts particle positions and velocities for each particleID 
#              in the selection box
# In bash:
#   export SELECTPART_EXE=/home/..../selectpart
#   export EXTRACTPOSVEL_EXE=/home/..../extractPosVelMasVol
#--------------------------------------------------------------------------
def set_execs_and_dirs():

  # Define the executables
  if (os.getenv('SELECTPART_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export SELECTPART_EXE=/path/to/selectpart"
    sys.exit()
  select_exe = os.path.abspath(os.environ['SELECTPART_EXE'])

  if (os.getenv('EXTRACTPOSVEL_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export EXTRACTPOSVEL_EXE=/path/to/extractPosVelMasVol"
    sys.exit()
  extract_exe = os.path.abspath(os.environ['EXTRACTPOSVEL_EXE'])

  # Define the working and output directories
  # The uda files will be read from the working directory
  if (os.getenv('WORKING_DIR') == None):
    print "**ERROR** Missing environment var. In bash: export WORKING_DIR=/path/to/uda/directory"
    sys.exit()
  work_dir = os.path.abspath(os.environ['WORKING_DIR'])

  # Make a directory for saving the plots and extracted data
  plot_dir = os.path.join(work_dir, "DataAndPlots")
  if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

  return select_exe, extract_exe, work_dir, plot_dir

#--------------------------------------------------------------------------
# Define the output file names
#--------------------------------------------------------------------------
def set_output_files(plot_dir, uda_name):

  # Define the selected particle list file + error file
  part_list_file = os.path.join(plot_dir, uda_name + ".partlist")
  error_file = os.path.join(plot_dir, uda_name + ".error")

  return part_list_file, error_file

#--------------------------------------------------------------------------
# Run selectpart
#--------------------------------------------------------------------------
def selectpart(select_exe, work_dir, uda_name, 
               selection_box, material_id, part_list_file, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Set the time steps to the first one
  timestep_low = 0
  timestep_high = 0

  # Create the command + arguments
  box = selection_box
  args = [select_exe, "-box", 
          str(box.xmin), str(box.ymin), str(box.zmin), 
          str(box.xmax), str(box.ymin), str(box.zmin), 
          str(box.xmin), str(box.ymax), str(box.zmin), 
          str(box.xmin), str(box.ymin), str(box.zmax), 
          "-mat", str(material_id), 
          "-timesteplow", str(timestep_low), "-timestephigh", str(timestep_high),
          "-uda", uda_dir]
  print args

  # Run the command and save output in part_list_file
  part_list = open(part_list_file, "w+")
  err_file = open(error_file+".select", "w+")
  tmp = sub_proc.Popen(args, stdout=part_list, stderr=err_file)
  dummy = tmp.wait()

  # Read the part list file to get a list of particle ids
  part_list.seek(0)
  particle_ids = []
  x = []
  y = []
  z = []
  for line in part_list:
    line = line.strip().split()
    particle_ids.append(line[3])
    x.append(line[4])
    y.append(line[5])
    z.append(line[6])

  # close and return
  part_list.close()
  err_file.close()
  print "Num. particle selected = ", len(particle_ids)
  print "X-range = ", min(x), ",", max(x)
  print "Y-range = ", min(y), ",", max(y)
  print "Z-range = ", min(z), ",", max(z)
  return zip(particle_ids, x, y, z)

#--------------------------------------------------------------------------
# Run extractPosVelMasVol for each particle id in the selected list
#--------------------------------------------------------------------------
def extractPosVelMasVol(extract_exe, work_dir, uda_name, plot_dir, material_id, time_id,
               part_list_file, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create the command + arguments
  output_file = os.path.join(plot_dir, uda_name + ".pos_" + str(material_id))
  args_pos = [extract_exe, "-m", str(material_id),
    "-p", part_list_file, "-timestep", str(time_id), "-uda", uda_dir, "-o", output_file]
  print args_pos

  # Run the command and save output in part_data_files
  pos_data = tempfile.TemporaryFile()
  #stress_data = open(stress_file, "w+")
  err_file = open(error_file+".posextract", "w+")
  #tmp = sub_proc.Popen(args_pos, stdout=pos_data, stderr=sub_proc.PIPE)
  tmp = sub_proc.Popen(args_pos, stdout=pos_data, stderr=err_file)
  dummy = tmp.wait()

  return 

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  uda_name, selection_box, material_id, timestep = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Selection box : ', selection_box
  print '  Material id : ', material_id
  print '  Timesteps  : ', timestep

  select_exe, extract_exe, work_dir, plot_dir = set_execs_and_dirs()
  print select_exe, extract_exe, work_dir, plot_dir

  part_list_file, error_file = set_output_files(plot_dir, uda_name)
  print part_list_file, error_file

  part_list_file = part_list_file + ".mat" + str(material_id)
  particle_ids = selectpart(select_exe, work_dir, uda_name, 
    selection_box, material_id, part_list_file, error_file)
  #print particle_ids

  extractPosVelMasVol(extract_exe, work_dir, uda_name, plot_dir, material_id, timestep, part_list_file, error_file)

