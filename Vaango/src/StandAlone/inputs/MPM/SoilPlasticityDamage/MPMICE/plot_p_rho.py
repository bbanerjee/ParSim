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

  def __init__(self, xmin, ymin, zmin, dx, dy, dz):
    self.xmin = xmin
    self.ymin = ymin
    self.zmin = zmin
    self.xmax = xmin + dx
    self.ymax = ymin + dy
    self.zmax = zmin + dz

  def __repr__(self):
    return "Box()"

  def __str__(self):
    return ("["+str(self.xmin)+","+str(self.ymin)+","+str(self.zmin)+","
      +str(self.xmax)+","+str(self.ymax)+","+str(self.zmax)+"]")

  def invalid(self):
    if (self.xmax <= self.xmin or self.ymax <= self.ymin or self.zmax <= self.zmin):
      return True
    return False

  def update(self, xmin, ymin, zmin, dx, dy, dz):
    self.xmin = xmin
    self.ymin = ymin
    self.zmin = zmin
    self.xmax = xmin + dx
    self.ymax = ymin + dy
    self.zmax = zmin + dz

#--------------------------------------------------------------------------
# Main function (parse arguments)
#--------------------------------------------------------------------------
def main(argv):

  uda_name = ''
  selection_box = Box(0,0,0,0,0,0) 
  material_id = 0 

  # Parse comman line arguments
  parser = argparse.ArgumentParser(
    usage="plot_p_rho --box xmin ymin zmin dx dy dz --mat <mat id> --uda <uda file_name_prefix>")
  parser.add_argument('--box', nargs=6, help = 'xmin ymin zmin dx dy dz', type=float, required=True)
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--uda', required=True)
  args = parser.parse_args()

  # Parse box
  vals = args.box
  xmin = float(vals[0])
  ymin = float(vals[1])
  zmin = float(vals[2])
  dx = float(vals[3])
  dy = float(vals[4])
  dz = float(vals[5])
  selection_box.update(xmin, ymin, zmin, dx, dy, dz)
  if selection_box.invalid():
    print "**ERROR** Invalid selction box."
    sys.exit()

  # Parse mat
  material_id = int(args.mat)

  # Parse uda
  uda_name = args.uda
  
  return uda_name, selection_box, material_id  

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
               selection_box, material_id, part_list_file, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file doe snot exist", uda_dir
    sys.exit()

  # Set the time steps to the first one
  timestep_low = 0
  timestep_high = 0

  # Create the command + arguments
  box = selection_box
  args = [select_exe, "-box", str(box.xmin), str(box.ymin), str(box.zmin), 
    str(box.xmax), str(box.ymax), str(box.zmax), "-mat", str(material_id), 
    "-timesteplow", str(timestep_low), "-timestephigh", str(timestep_high),
    uda_dir]
  print args

  # Run the command and save output in part_list_file
  part_list = open(part_list_file, "w+")
  err_file = open(error_file, "w+")
  tmp = sub_proc.Popen(args, stdout=part_list, stderr=err_file)
  dummy = tmp.wait()

  # Read the part list file to get a list of particle ids
  part_list.seek(0)
  particle_ids = []
  for line in part_list:
    line = line.strip().split()
    particle_ids.append(line[3])

  # close and return
  part_list.close()
  err_file.close()
  return particle_ids

#--------------------------------------------------------------------------
# Run partextract for each particle id in the selected list
#--------------------------------------------------------------------------
def partextract(extract_exe, work_dir, uda_name, plot_dir, particle_id):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create stress, velocity, mass, and volume data files
  #stress_file = os.path.join(plot_dir, uda_name+".stress."+str(particle_id))
  #velocity_file = os.path.join(plot_dir, uda_name+".velocity."+str(particle_id))
  #mass_file = os.path.join(plot_dir, uda_name+".mass."+str(particle_id))
  #volume_file = os.path.join(plot_dir, uda_name+".volume."+str(particle_id))

  # Create the command + arguments
  args_stress = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.stress", uda_dir]
  args_velocity = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.velocity", uda_dir]
  args_mass = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.mass", uda_dir]
  args_volume = [extract_exe, "-partid", str(particle_id),
    "-partvar", "p.volume", uda_dir]
  print args_velocity

  # Run the command and save output in part_data_files
  stress_data = tempfile.TemporaryFile()
  #stress_data = open(stress_file, "w+")
  tmp = sub_proc.Popen(args_stress, stdout=stress_data, stderr=sub_proc.PIPE)
  dummy = tmp.wait()

  # Read the part data file to get a list of stresses as a function of time
  stress_data.seek(0)
  times = []
  pressures = []
  for line in stress_data:
    line = line.strip().split()
    times.append(float(line[0]))
    S11 = np.float64(line[4])
    S22 = np.float64(line[8])
    S33 = np.float64(line[12])
    p = (S11+S22+S33)/3.0
    pressures.append(p)

  # close and return
  stress_data.close()

  # Run the command and save output in part_data_files
  velocity_data = tempfile.TemporaryFile()
  #velocity_data = open(velocity_file, "w+")
  tmp = sub_proc.Popen(args_velocity, stdout=velocity_data, stderr=sub_proc.PIPE)
  dummy = tmp.wait()

  # Read the part data file to get a list of velocities as a function of time
  velocity_data.seek(0)
  velocities = []
  for line in velocity_data:
    line = line.strip().split()
    v1 = np.float64(line[4])
    v2 = np.float64(line[5])
    v3 = np.float64(line[6])
    #velocities.append(np.array([v1, v2, v3]))
    velocities.append([v1, v2, v3])

  # close and return
  velocity_data.close()

  # Run the command and save output in part_data_files
  mass_data = tempfile.TemporaryFile()
  #mass_data = open(mass_file, "w+")
  tmp = sub_proc.Popen(args_mass, stdout=mass_data, stderr=sub_proc.PIPE)
  dummy = tmp.wait()

  # Read the part data file to get a list of masses as a function of time
  mass_data.seek(0)
  masses = []
  for line in mass_data:
    line = line.strip().split()
    m = np.float64(line[4])
    masses.append(m)

  # close and return
  mass_data.close()

  # Run the command and save output in part_data_files
  volume_data = tempfile.TemporaryFile()
  #volume_data = open(volume_file, "w+")
  tmp = sub_proc.Popen(args_volume, stdout=volume_data, stderr=sub_proc.PIPE)
  dummy = tmp.wait()

  # Read the part data file to get a list of volumes as a function of time
  volume_data.seek(0)
  volumes = []
  for line in volume_data:
    line = line.strip().split()
    v = np.float64(line[4])
    volumes.append(v)

  # close and return
  volume_data.close()

  # Compute density
  density = np.array(masses, dtype=np.float)/np.array(volumes, dtype=np.float)

  # Join all the data into one data frame with columns
  # time, pressure,, vx, vy, vz, mass, volume, density
  particle_data = zip(times, pressures, velocities, masses, volumes, density)

  # Write all data to output file
  data_file = os.path.join(plot_dir, uda_name+".data."+str(particle_id))
  fid = open(data_file, 'w')
  for data in particle_data:
    #data_list = [data[0], data[1], data[2][0], data[2][1], data[2][2], data[3], data[4], data[5]]
    #print data_list
    fid.write('%s %s %s %s %s %s %s %s\n' % (data[0], data[1], 
      data[2][0], data[2][1], data[2][2], data[3], data[4], data[5]))

  fid.close()
    
  return particle_data

#--------------------------------------------------------------------------
# Set up plot parameters
#--------------------------------------------------------------------------
def set_plot_parameters():

  #Set matplotlib defaults to desired values
  #S#Set the legend to best fit
  fontSize = 16
  markers = None
  plt.rcParams['legend.loc']='best'
  #Set font size
  plt.rcParams['mathtext.it'] = 'serif:bold'
  plt.rcParams['mathtext.rm'] = 'serif:bold'
  plt.rcParams['mathtext.sf'] = 'serif:bold'
  plt.rcParams['font.size']=fontSize
  plt.rcParams['font.weight']='bold'
  plt.rcParams['axes.labelsize']='medium'
  #plt.rcParams['axes.labelweight']='bold'
  plt.rcParams['legend.fontsize']='medium'
  #Set linewidth
  lineWidth = 2
  plt.rcParams['lines.linewidth']=lineWidth
  #Set markersize
  plt.rcParams['lines.markersize'] = 8
  #Set padding for tick labels and size
  plt.rcParams['xtick.major.pad']  = 12
  plt.rcParams['ytick.major.pad']  = 8
  plt.rcParams['xtick.major.size'] = 6
  plt.rcParams['xtick.minor.size'] = 3
  plt.rcParams['ytick.major.size'] = 6
  plt.rcParams['ytick.minor.size'] = 3
  #resolution
  plt.rcParams['figure.dpi']=120
  
  font = {'family' : 'serif',
          'weight' : 'bold',
          'size'   : fontSize}
  rc('font', **font)
  rc('text', usetex=True)

#--------------------------------------------------------------------------
# Save as pdf
#--------------------------------------------------------------------------
def savePDF(fig, name, size='1920x1080'):
  res = float(plt.rcParams['figure.dpi'])
  #Add Check for file already existing as name.png
  if size == '640x480':
    size = [640/res,480/res]
  if size == '1080x768':
    size = [1080/res,768/res]
  if size == '1152x768':
    size = [1152/res,768/res]
  if size == '1280x854':
    size = [1280/res,854/res]
  if size == '1280x960':
    size = [1280/res,960/res]
  if size == '1920x1080':
    size = [1920/res,1080/res]
  #set the figure size for saving
  fig.set_size_inches(size[0],size[1])
  #save at speciified resolution
  fig.savefig(name+'.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

#--------------------------------------------------------------------------
# For bold math fonts
#--------------------------------------------------------------------------
def str_to_mathbf(string):
  #Only works with single spaces no leading space
  string = string.split()
  return_string = ''
  for elem in string:
    elem = r'$\mathbf{'+elem+'}$'
    return_string+=elem+'  '
  return return_string[0:-1]

#--------------------------------------------------------------------------
# For exponential format
#--------------------------------------------------------------------------
def exp_fmt(x,loc):
  tmp = format(x,'1.2e').split('e')
  lead = tmp[0]
  exp = str(int(tmp[1]))
  if exp=='0' and lead=='0.00':
    return r'$\mathbf{0.00}$'
  else:
    if int(exp)<10 and int(exp)>0:
      exp = '+0'+exp
    elif int(exp)>-10 and int(exp)<0:
      exp = '-0'+exp.split('-')[1]
    elif int(exp)>10:
      exp = '+'+exp
    return r'$\mathbf{'+lead+r'\cdot{}10^{'+exp+'}}$'

#--------------------------------------------------------------------------
# Plot the pressure as a function of time
#--------------------------------------------------------------------------
def plot_pressure(particle_data, particle_id, plt_color):

  # Get the dat ainto plottable form
  times = []
  pressures = []
  for data in particle_data:
    times.append(data[0])
    pressures.append(data[1])

  # Plot p-t
  plt.figure("pressure")
  ax = plt.subplot(111)
  plt.plot(np.array(times), np.array(pressures), 'b-', 
    linewidth=2, color=plt_color, label=str(particle_id))

  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('Pressure (Pa)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax.xaxis.set_major_formatter(formatter_exp)
  #ax.yaxis.set_major_formatter(formatter_exp)    
  ax.xaxis.set_major_formatter(formatter_int)
  ax.yaxis.set_major_formatter(formatter_int)

  #ax.set_xlim(Xlims[0],Xlims[1])
  #ax.set_ylim(Ylims[0],Ylims[1])
  
  return ax

#--------------------------------------------------------------------------
# Plot the density as a function of time
#--------------------------------------------------------------------------
def plot_density(particle_data, particle_id, plt_color):

  # Get the dat ainto plottable form
  times = []
  densities = []
  for data in particle_data:
    times.append(data[0])
    densities.append(data[5])

  # Plot p-t
  plt.figure("density")
  ax = plt.subplot(111)
  plt.plot(np.array(times), np.array(densities), 'b-', 
    linewidth=2, color=plt_color, label=str(particle_id))

  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('Density (kg/m^3)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax.xaxis.set_major_formatter(formatter_exp)
  #ax.yaxis.set_major_formatter(formatter_exp)    
  ax.xaxis.set_major_formatter(formatter_int)
  ax.yaxis.set_major_formatter(formatter_int)

  #ax.set_xlim(Xlims[0],Xlims[1])
  #ax.set_ylim(Ylims[0],Ylims[1])
  
  return ax




#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  uda_name, selection_box, material_id = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Selection box : ', selection_box
  print '  Material id : ', material_id

  select_exe, extract_exe, work_dir, plot_dir = set_execs_and_dirs()
  print select_exe, extract_exe, work_dir, plot_dir

  part_list_file, error_file = set_output_files(plot_dir, uda_name)
  print part_list_file, error_file

  particle_ids = selectpart(select_exe, work_dir, uda_name, 
    selection_box, material_id, part_list_file, error_file)
  print particle_ids

  # Set up figures
  set_plot_parameters()
  fig1 = plt.figure("pressure")
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70, param_text,ha='left',va='top',size='x-small')

  fig2 = plt.figure("density")
  plt.clf()

  # Loop thru particles
  count = 0
  for particle_id in particle_ids:
    count = count + 1
    particle_data = partextract(extract_exe, 
      work_dir, uda_name, plot_dir, particle_id)
   
    plt_color = cm.BrBG(float(count)/float(len(particle_ids)))
    plot_pressure(particle_data, particle_id, plt_color)
    plot_density(particle_data, particle_id, plt_color)

  plt.figure("pressure")
  plt.grid(True)
  plt.legend(loc='best', prop={'size':10})
  savePDF(fig1, 'pressure_vs_time', size='1280x960')

  plt.figure("density")
  plt.grid(True)
  #plt.legend(loc='best', prop={'size':10})
  savePDF(fig2, 'density_vs_time', size='1280x960')
  plt.show()

