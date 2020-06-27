#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import shlex
import multiprocessing
import subprocess as sub_proc  
import tempfile
import argparse
import numpy as np
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib import ticker

#--------------------------------------------------------------------------
# Main function (parse arguments)
#--------------------------------------------------------------------------
def main(argv):

  uda_name = ''
  material_id = 0 

  # Parse comman line arguments
  parser = argparse.ArgumentParser(
    usage="extract_centrifuge_mom --box <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> --mat <mat id> --tlo <int timestep low> --thi <int timestep hi> --tstep <int stepsize> --uda <uda file_name_prefix>")
  parser.add_argument('--box', nargs=6, help = 'xmin ymin zmin xmax ymax zmax', type=float, required=True)
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--uda', required=True)
  parser.add_argument('--tlo', default=0, type=int, required=False)
  parser.add_argument('--thi', default=0, type=int, required=True)
  parser.add_argument('--tstep', default=1, type=int, required=False)
  args = parser.parse_args()

  vals = args.box
  if (vals):
    xmin = float(vals[0])
    ymin = float(vals[1])
    zmin = float(vals[2])
    xmax = float(vals[3])
    ymax = float(vals[4])
    zmax = float(vals[5])

  # Parse mat
  material_id = int(args.mat)

  # Parse time step lo/hi
  timestep_lo = int(args.tlo)
  timestep_hi = int(args.thi)
  timestep_inc = int(args.tstep)

  # Parse uda
  uda_name = args.uda
  
  ptmin = [xmin, ymin, zmin]
  ptmax = [xmax, ymax, zmax]

  return uda_name, ptmin, ptmax, material_id, timestep_lo, timestep_hi, timestep_inc

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

  if (os.getenv('EXTRACTPVEC_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export EXTRACTPVEC_EXE=/path/to/extractPVec"
    sys.exit()
  extract_vec_exe = os.path.abspath(os.environ['EXTRACTPVEC_EXE'])

  if (os.getenv('EXTRACTPSCA_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export EXTRACTPSCA_EXE=/path/to/extractPscalar"
    sys.exit()
  extract_sca_exe = os.path.abspath(os.environ['EXTRACTPSCA_EXE'])

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

  return select_exe, extract_vec_exe, extract_sca_exe, work_dir, plot_dir

#--------------------------------------------------------------------------
# Run selectpart
#--------------------------------------------------------------------------
def selectpart(select_exe, work_dir, uda_name, plot_dir,
               ptmin, ptmax, material_id, tlo, thi, tinc, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create the command + arguments
  pt0 = ptmin;
  pt1 = [ptmax[0], ptmin[1], ptmin[2]]
  pt2 = [ptmin[0], ptmax[1], ptmin[2]]
  pt3 = [ptmin[0], ptmin[1], ptmax[2]]
  command = (select_exe + " -box " +  
             str(pt0[0]) + " " +  str(pt0[1]) + " " +  str(pt0[2]) + " " +  
             str(pt1[0]) + " " +  str(pt1[1]) + " " +  str(pt1[2]) + " " +  
             str(pt2[0]) + " " +  str(pt2[1]) + " " +  str(pt2[2]) + " " +
             str(pt3[0]) + " " +  str(pt3[1]) + " " +  str(pt3[2]) + " " +
             " -mat " +  str(material_id) +  " -timesteplow 0 " + " -timestephigh 0 " +
             " -uda " +  uda_dir)
  print(command)
  args = shlex.split(command)
  print(args)

  # Run the command and save output in partlist_file
  partlist_file_name = uda_name + ".mom.partlist"
  partlist_file_dir = os.path.join(plot_dir, partlist_file_name)
  partlist_file = open(partlist_file_dir, "w+")
  err_file = open(error_file+".select", "w+")
  tmp = sub_proc.Popen(args, stdout=partlist_file, stderr=err_file)
  dummy = tmp.wait()

  # Read the part list file to get a list of particle ids
  partlist_file.seek(0)
  particle_ids = []
  x = []
  y = []
  z = []
  for line in partlist_file:
    line = line.strip().split()

    particle_ids.append(line[3])
    x.append(line[4])
    y.append(line[5])
    z.append(line[6])

  # close and return
  partlist_file.close()
  err_file.close()
  return partlist_file_dir, zip(particle_ids, x, y, z)


#--------------------------------------------------------------------------
# Run partextract for each particle id in the selected list
#--------------------------------------------------------------------------
def partextract(extract_vec_exe, extract_sca_exe, work_dir, uda_name, plot_dir, partlist,
                partlist_file, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create velocity file
  velocity_file_name = uda_name + ".velocity"
  velocity_file_dir = os.path.join(plot_dir, velocity_file_name)

  # Create the command + arguments (for the velocity)
  command = (extract_vec_exe + 
             " -partvar " + " p.velocity " +
             " -m " + str(material_id) +
             " -p " + partlist_file + 
             " -uda " + uda_dir +
             " -o " + velocity_file_dir + 
             " -timefiles")
  print(command)
  args_velocity = shlex.split(command)
  print(args_velocity)

  # Run the command 
  velocity_stdout = tempfile.TemporaryFile()
  err_file = open(error_file+".velocity_extract", "w+")
  tmp = sub_proc.Popen(args_velocity, stdout=velocity_stdout, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of velocitys as a function of time
  velocity_data = open(velocity_file_dir, "r+")
  velocity_data.seek(0)
  particles = defaultdict(list)
  times = defaultdict(list)
  velocities = defaultdict(list)
  for line in velocity_data:
    line = line.strip().split()

    time = float(line[0])
    particle = line[3]
    times[particle].append(time)
    particles[time].append(particle)
    v1 = np.float64(line[4])
    v2 = np.float64(line[5])
    v3 = np.float64(line[6])
    vel = np.array([v1, v2, v3])
    velocities[time].append(vel)

  # close and return
  velocity_data.close()
  err_file.close()

  # Create mass file
  mass_file_name = uda_name + ".mass"
  mass_file_dir = os.path.join(plot_dir, mass_file_name)

  # Create the command + arguments (for the mass)
  command = (extract_sca_exe + 
             " -partvar " + " p.mass " +
             " -m " + str(material_id) +
             " -p " + partlist_file + 
             " -uda " + uda_dir +
             " -o " + mass_file_dir + 
             " -timefiles")
  print(command)
  args_mass = shlex.split(command)
  print(args_mass)

  # Run the command 
  mass_stdout = tempfile.TemporaryFile()
  err_file = open(error_file+".mass_extract", "w+")
  tmp = sub_proc.Popen(args_mass, stdout=mass_stdout, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of masses as a function of time
  mass_data = open(mass_file_dir, "r+")
  mass_data.seek(0)
  masses = defaultdict(list)
  for line in mass_data:
    line = line.strip().split()

    time = float(line[0])
    particle = line[3]
    mm = np.float64(line[4])
    masses[time].append(mm)

  # close and return
  mass_data.close()
  err_file.close()

  
  # Compute momentum
  def calc_mom(vel, mass):
    mom = np.array([vel[0]*mass, vel[1]*mass, vel[2]*mass])
    return(mom)

  part0 = partlist[0][0]
  times_list = times[part0]
  total_momentum_list = []
  for time in times_list:
    momenta = map(lambda vel, mass: calc_mom(vel, mass),
                  velocities[time], masses[time])
    total_momentum = reduce(lambda mom1, mom2: np.array([mom1[0]+mom2[0],
                                                         mom1[1]+mom2[1],
                                                         mom1[2]+mom2[2]]),
                            momenta)
    total_momentum_mag = np.sqrt(total_momentum[0]**2 +
                                 total_momentum[1]**2 +
                                 total_momentum[2]**2)
    total_momentum_list.append(total_momentum_mag)


  # Join all the data into one data frame with columns
  # time, sum of normal momentum component
  particle_data = zip(times_list, total_momentum_list)
  #print particle_data

  return particle_data

#--------------------------------------------------------------------------
# Set up plot parameters
#--------------------------------------------------------------------------
def set_plot_parameters():

  #Set matplotlib defaults to desired values
  #S#Set the legend to best fit
  fontSize = 16
  markers = None
  plt.rcParams['legend.loc']='right'
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

def savePNG(fig, name, size='1920x1080'):
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
  fig.savefig(name+'.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

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
# Plot the total momentum as a function of time
#--------------------------------------------------------------------------
def plot_momentum(fig, mom_data, plt_color, line_style):

  # Function that returns t
  # Get the data into plottable form
  times = []
  mom = []
  for data in mom_data:
    times.append(data[0]*1.0e3)
    mom.append(data[1])

  # Plot x-t
  plt.figure(fig)
  ax = plt.subplot(111)
  plt.plot(np.array(times), np.array(mom), linestyle=line_style,
    linewidth=2, color=plt_color)

  plt.xlabel(str_to_mathbf('Time (ms)'))
  plt.ylabel(str_to_mathbf('Total momentum (kg-m/s)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax.xaxis.set_major_formatter(formatter_exp)
  #ax.yaxis.set_major_formatter(formatter_exp)    
  ax.xaxis.set_major_formatter(formatter_int)
  ax.yaxis.set_major_formatter(formatter_int)

  #ax.set_xlim(0, 10)
  #ax.set_ylim(-10, 50)
  #ax.set_ylim(Ylims[0],Ylims[1])
  ax.grid(True)
  
  return ax

#--------------------------------------------------------------------------
# Save the total momentum as a function of time
#--------------------------------------------------------------------------
def save_momentum(mom_data, plt_file):

  save_file = plt_file + ".data"
  f = open(save_file, 'w')
  f.write("Time (s)" + " " + "Total momentum (kg-m/s)" + "\n")

  # Get the data into writeable form
  length = len(list(mom_data))
  index = 0
  for data in mom_data:
    time = data[0]
    mom = data[1]
    if (index == length - 1):
      f.write(str(time) + " " + str(mom))
    else:
      f.write(str(time) + " " + str(mom) + "\n")
    index = index + 1

  f.close()

  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  # Read the command line
  uda_name, ptmin, ptmax, material_id, tlo, thi, tinc = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Box : '
  print '    Pmin =  [', ptmin[0], ',', ptmin[1], ',', ptmin[2], ']; '
  print '    Pmax =  [', ptmax[0], ',', ptmax[1], ',', ptmax[2], ']; '
  print '  Material id : ', material_id
  print '  t_low = ', tlo, ' t_high = ', thi, ' t_inc = ', tinc

  # Set up executables and working directories
  select_exe, extract_vec_exe, extract_sca_exe, work_dir, plot_dir = set_execs_and_dirs()
  print select_exe, extract_vec_exe, extract_sca_exe, work_dir, plot_dir

  # Create the particle list
  error_file = os.path.join(plot_dir, uda_name + ".error")
  partlist_file, part_select = selectpart(select_exe, work_dir, uda_name, plot_dir,
                                          ptmin, ptmax, material_id, 
                                          tlo, thi, tinc, error_file)
  #print part_select

  # Extract the particle tractions assuming normals remain fixed
  data = partextract(extract_vec_exe, extract_sca_exe, work_dir, uda_name, plot_dir, 
                     part_select, partlist_file, 
                     error_file)

  # Set up figures
  set_plot_parameters()
  fig1 = plt.figure("momentum")
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70, param_text,ha='left',va='top',size='x-small')

  #colors = ['Red', 'Green', 'Blue', 'Magenta', 'Violet']
  #plt_color = colors[idx]
  plt_color = cm.PiYG(float(2)/float(len(data[0])))
  line_style = '-'
  plot_momentum('momentum', data, plt_color, line_style)

  plt.figure('momentum')
  #plt.legend(bbox_to_anchor=(1.05,1), loc=2, prop={'size':10})
  plt_file = os.path.join(plot_dir, uda_name+".mom")
  savePNG(fig1, plt_file, size='1280x960')
  plt.show()

  # Save the momentum data
  save_momentum(data, plt_file)

