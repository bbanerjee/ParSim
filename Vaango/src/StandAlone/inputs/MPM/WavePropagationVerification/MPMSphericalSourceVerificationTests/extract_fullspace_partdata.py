#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import multiprocessing
import subprocess as sub_proc  
import tempfile
import argparse
import math
import numpy as np

#--------------------------------------------------------------------------
# Create a Box class and associated methods
# Alternatively: Use a named tuple : from collections import namedtuple
#--------------------------------------------------------------------------
class Box():

  def __init__(self, x, y, z, dx, dy, dz):
    self.xmin = x - dx
    self.ymin = y - dy
    self.zmin = z - dz
    self.xmax = x + dx
    self.ymax = y + dy
    self.zmax = z + dz

  def __repr__(self):
    return "Box()"

  def __str__(self):
    return ("["+str(self.xmin)+","+str(self.ymin)+","+str(self.zmin)+","
      +str(self.xmax)+","+str(self.ymax)+","+str(self.zmax)+"]")

  def invalid(self):
    if (self.xmax <= self.xmin or self.ymax <= self.ymin or self.zmax <= self.zmin):
      return True
    return False

  def update(self, x, y, z, dx, dy, dz):
    self.xmin = x - dx
    self.ymin = y - dy
    self.zmin = z - dz
    self.xmax = x + dx
    self.ymax = y + dy
    self.zmax = z + dz

#--------------------------------------------------------------------------
# Main function (parse arguments)
#--------------------------------------------------------------------------
def main(argv):

  uda_name = ''
  material_id = 0 

  # Parse comman line arguments
  parser = argparse.ArgumentParser(
    usage="extract_fullspace_partdata --point <x> <y> <z> --theta <value> --phi <value> --dx <value> --mat <mat id> --uda <uda file_name_prefix>")
  parser.add_argument('--point', nargs=3, help = 'x y z', type=float, required=True)
  parser.add_argument('--theta', default = 180, type=float, required=True, choices=xrange(0,360))
  parser.add_argument('--phi', default = 0, type=float, required=True, choices=xrange(-90,90))
  parser.add_argument('--dx', default=0.01, type=float, required=True)
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--uda', required=True)
  args = parser.parse_args()

  vals = args.point
  if (vals):
    x = float(vals[0])
    y = float(vals[1])
    z = float(vals[2])

  theta = math.radians(float(args.theta))
  phi = math.radians(float(args.phi))

  dx = float(args.dx)

  # Parse mat
  material_id = int(args.mat)

  # Parse uda
  uda_name = args.uda
  
  return uda_name, x, y, z, theta, phi, dx, material_id

#--------------------------------------------------------------------------
# Define the executables to be used and the working directories 
# Two executables are needed - selectpart and partextract
# selectpart:  selects a list of particleIDs inside a selection box
# partextract: extracts particle variables corresponding to each particleID 
#              in the selection box
# In bash:
#   export SELECTPART_EXE=/home/..../selectpart
#   export PARTEXTRACT_EXE=/home/..../partextract
#--------------------------------------------------------------------------
def set_execs_and_dirs():

  # Define the executables
  if (os.getenv('SELECTPART_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export SELECTPART_EXE=/path/to/selectpart"
    sys.exit()
  select_exe = os.path.abspath(os.environ['SELECTPART_EXE'])

  if (os.getenv('PARTEXTRACT_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export PARTEXTRACT_EXE=/path/to/partextract"
    sys.exit()
  extract_exe = os.path.abspath(os.environ['PARTEXTRACT_EXE'])

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
               selection_box, material_id, 
               timestep, part_list_file, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file doe snot exist", uda_dir
    sys.exit()

  # Create the command + arguments
  box = selection_box
  args = [select_exe, "-box", str(box.xmin), str(box.ymin), str(box.zmin), 
    str(box.xmax), str(box.ymax), str(box.zmax), "-mat", str(material_id), 
    "-timesteplow", str(timestep), "-timestephigh", str(timestep),
    uda_dir]
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
  return zip(particle_ids, x, y, z)

#--------------------------------------------------------------------------
# Run partextract for each particle id in the selected list
#--------------------------------------------------------------------------
def partextract(extract_exe, work_dir, uda_name, particle_id, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create the command + arguments (for the position)
  args_pos = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.x", uda_dir]
  print args_pos

  # Run the command and save output in part_data_files
  pos_data = tempfile.TemporaryFile()
  err_file = open(error_file+".pos_extract", "w+")
  tmp = sub_proc.Popen(args_pos, stdout=pos_data, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of positions as a function of time
  pos_data.seek(0)
  times = []
  x = []
  y = []
  z = []
  for line in pos_data:
    line = line.strip().split()
    times.append(float(line[0]))
    x.append(np.float64(line[4]))
    y.append(np.float64(line[5]))
    z.append(np.float64(line[6]))

  # close and return
  pos_data.close()
  err_file.close()

  # Create the command + arguments (for the velocity)
  args_vel = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.velocity", uda_dir]
  print args_vel

  # Run the command and save output in part_data_files
  vel_data = tempfile.TemporaryFile()
  err_file = open(error_file+".vel_extract", "w+")
  tmp = sub_proc.Popen(args_vel, stdout=vel_data, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of velocities as a function of time
  vel_data.seek(0)
  times = []
  vx = []
  vy = []
  vz = []
  for line in vel_data:
    line = line.strip().split()
    times.append(float(line[0]))
    vx.append(np.float64(line[4]))
    vy.append(np.float64(line[5]))
    vz.append(np.float64(line[6]))

  # close and return
  vel_data.close()

  # Create the command + arguments (for the stress)
  args_stress = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.stress", uda_dir]
  print args_stress

  # Run the command and save output in part_data_files
  stress_data = tempfile.TemporaryFile()
  err_file = open(error_file+".stress_extract", "w+")
  tmp = sub_proc.Popen(args_stress, stdout=stress_data, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of stresse as a function of time
  stress_data.seek(0)
  times = []
  sigx = []
  sigy = []
  sigz = []
  for line in stress_data:
    line = line.strip().split()
    times.append(float(line[0]))
    sigx.append(np.float64(line[4]))
    sigy.append(np.float64(line[8]))
    sigz.append(np.float64(line[12]))

  # close and return
  stress_data.close()

  # Join all the data into one data frame with columns
  particle_data = zip(times, x, y, z, vx, vy, vz, sigx, sigy, sigz)

  return particle_data

#-----------------------------------------------------------------------
# Detector positions
#  -- These will have to be changed as needed
#-----------------------------------------------------------------------
def getPositions(xcen, ycen, zcen, radii, theta, phi):

  # Generate detector positions
  x = []
  y = []
  z = []
  for radius in radii:
    x.append(xcen + radius*math.cos(theta)*math.cos(phi))
    y.append(ycen + radius*math.sin(theta)*math.cos(phi))
    z.append(zcen + radius*math.sin(phi))
  print x
  print y
  print z
  return zip(x, y, z)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  # These are the radial distances from the input point
  # at which "detector" particles will be selected
  radii = [0.1, 0.152, 0.203, 0.254, 0.305, 0.356]

  # The spherical angles of the detector
  # theta = azimuthal angle
  # phi = polar angle
  theta = math.radians(180)
  phi = math.radians(0)

  uda_name, x, y, z, theta, phi, dx, material_id = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Point : ', x, y, z
  print '  Detector orientation : ', theta, phi
  print '  dx : ', dx
  print '  Material id : ', material_id

  select_exe, extract_exe, work_dir, plot_dir = set_execs_and_dirs()
  print select_exe, extract_exe, work_dir, plot_dir

  part_list_file, error_file = set_output_files(plot_dir, uda_name)
  print part_list_file, error_file

  # Get the locations of the detector (change if needed)
  positions = getPositions(x, y, z, radii, theta, phi)

  # Select the particles from the first timestep
  particle_ids = []
  for idx, pos in enumerate(positions):
    # Create selection box
    xpos = pos[0]
    ypos = pos[1]
    zpos = pos[2]
    selection_box = Box(float(xpos), float(ypos), float(zpos), 
                        float(dx), float(dx), float(dx))

    # Get the particle ids inside the box
    timestep = 0
    part_select = selectpart(select_exe, work_dir, uda_name, selection_box, material_id, 
                             timestep, part_list_file, error_file)
    #print particle_ids

    # part_select may contain id, x, y, z for more than one particle
    # save the first id
    particle_ids.append(part_select[0])
    
  # Loop thru particles extract data
  # Write all data to output file
  data_file = os.path.join(plot_dir, uda_name+".partdata")
  fid = open(data_file, 'a+')
  count = 0
  for particle_id in particle_ids:
    count = count + 1
    data = partextract(extract_exe, work_dir, uda_name, particle_id[0], error_file)
    for time_data in data:
      fid.write('%s %s %s %s %s %s %s %s %s %s\n' % (time_data[0], time_data[1], time_data[2], 
        time_data[3], time_data[4], time_data[5], time_data[6],
        time_data[7], time_data[8], time_data[9]))
    print particle_id, ":", count, " of ", len(particle_ids)

  fid.close()
