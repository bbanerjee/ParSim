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
    usage="extract_centrifuge_fext --plane <x0> <y0> <z0> <x1> <y1> <z1> <x2> <y2> <z2> --surf_area <area> --mat <mat id> --tlo <int timestep low> --thi <int timestep hi> --tstep <int stepSize> --uda <uda file_name_prefix>")
  parser.add_argument('--plane', nargs=9, help = 'x0 y0 z0 x1 y1 z1 x2 y2 z2', type=float, required=True)
  parser.add_argument('--surf_area', default=0, type=float, required=True)
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--uda', required=True)
  parser.add_argument('--tlo', default=0, type=int, required=False)
  parser.add_argument('--thi', default=0, type=int, required=True)
  parser.add_argument('--tstep', default=1, type=int, required=False)
  args = parser.parse_args()

  vals = args.plane
  if (vals):
    x0 = float(vals[0])
    y0 = float(vals[1])
    z0 = float(vals[2])
    x1 = float(vals[3])
    y1 = float(vals[4])
    z1 = float(vals[5])
    x2 = float(vals[6])
    y2 = float(vals[7])
    z2 = float(vals[8])

  surf_area = float(args.surf_area)

  # Parse mat
  material_id = int(args.mat)

  # Parse time step lo/hi
  timestep_lo = int(args.tlo)
  timestep_hi = int(args.thi)
  timestep_inc = int(args.tstep)

  # Parse uda
  uda_name = args.uda
  
  pt0 = [x0, y0, z0]
  pt1 = [x1, y1, z1]
  pt2 = [x2, y2, z2]

  return uda_name, pt0, pt1, pt2, surf_area, material_id, timestep_lo, timestep_hi, timestep_inc

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
  extract_exe = os.path.abspath(os.environ['EXTRACTPMAT_EXE'])

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
# Run selectpart
#--------------------------------------------------------------------------
def selectpart(select_exe, work_dir, uda_name, plot_dir,
               pt0, pt1, pt2, material_id, tlo, thi, tinc, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create the command + arguments
  command = (select_exe + " -plane " +  
             str(pt0[0]) + " " +  str(pt0[1]) + " " +  str(pt0[2]) + " " +  
             str(pt1[0]) + " " +  str(pt1[1]) + " " +  str(pt1[2]) + " " +  
             str(pt2[0]) + " " +  str(pt2[1]) + " " +  str(pt2[2]) + 
             " -mat " +  str(material_id) +  " -timesteplow 0 " + " -timestephigh 0 " +
             " -uda " +  uda_dir)
  print(command)
  args = shlex.split(command)
  print(args)

  # Run the command and save output in partlist_file
  partlist_file_name = uda_name + ".stress.partlist"
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

    # Special for abaqus geometry pieces where particles are
    # not distributed accrording to the grid
    if (abs(np.float64(line[6]) - pt0[2]) < 1.0e-6):
      particle_ids.append(line[3])
      x.append(line[4])
      y.append(line[5])
      z.append(line[6])
      print(line[4], line[5], line[6])

  # close and return
  partlist_file.close()
  err_file.close()
  return partlist_file_dir, zip(particle_ids, x, y, z)


#--------------------------------------------------------------------------
# Run partextract for each particle id in the selected list
#--------------------------------------------------------------------------
def partextract(extract_exe, work_dir, uda_name, plot_dir, partlist,
                partlist_file, surf_area, surf_normal, pt0, error_file):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create stress file
  stress_file_name = uda_name + ".stress"
  stress_file_dir = os.path.join(plot_dir, stress_file_name)

  # Create the command + arguments (for the stress)
  command = (extract_exe + 
             " -partvar " + " p.stress " +
             " -m " + str(material_id) +
             " -p " + partlist_file + 
             " -uda " + uda_dir +
             " -o " + stress_file_dir + 
             " -timefiles")
  print(command)
  args_stress = shlex.split(command)
  print(args_stress)

  # Run the command 
  stress_stdout = tempfile.TemporaryFile()
  err_file = open(error_file+".stress_extract", "w+")
  tmp = sub_proc.Popen(args_stress, stdout=stress_stdout, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of stresss as a function of time
  stress_data = open(stress_file_dir, "r+")
  stress_data.seek(0)
  particles = defaultdict(list)
  times = defaultdict(list)
  traction_n = defaultdict(list)
  for line in stress_data:
    line = line.strip().split()

    # Special for abaqus geometry pieces where particles are
    # not distributed accrording to the grid
    if (abs(np.float64(line[15]) - pt0[2]) < 1.0e-6):
      time = float(line[0])
      particle = line[3]
      times[particle].append(time)
      particles[time].append(particle)
      S11 = np.float64(line[4])
      S12 = np.float64(line[5])
      S13 = np.float64(line[6])
      S22 = np.float64(line[8])
      S23 = np.float64(line[9])
      S33 = np.float64(line[12])
      sigma = np.array([[S11, S12, S13],[S12, S22, S23],[S13, S23, S33]])
      traction = sigma.dot(surf_normal)
      traction_n[time].append(traction.dot(surf_normal))

  # close and return
  stress_data.close()
  err_file.close()

  # Filter the particle data
  part0 = partlist[0][0]
  times_list = times[part0]
  #print particles[times_list[0]]
  #print len(particles[times_list[0]])
  numPart = float(len(particles[times_list[0]]))
  traction_list = []
  for time in times_list:
    print sum(traction_n[time])
    traction_list.append(sum(traction_n[time])*surf_area/numPart)

  # Join all the data into one data frame with columns
  # time, sum of normal traction component
  particle_data = zip(times_list, traction_list)
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
# Plot the pressure as a function of time
#--------------------------------------------------------------------------
def plot_traction(fig, press_data, plt_color, line_style):

  # Function that returns t
  # Get the dat ainto plottable form
  times = []
  press = []
  for data in press_data:
    times.append(data[0]*1.0e3)
    press.append(-data[1]*1.0e-3)
    #press.append(-data[1]*0.000145037738)

  # Plot x-t
  plt.figure(fig)
  ax = plt.subplot(111)
  plt.plot(np.array(times), np.array(press), linestyle=line_style,
    linewidth=2, color=plt_color)

  plt.xlabel(str_to_mathbf('Time (ms)'))
  plt.ylabel(str_to_mathbf('Impact force (kN)'))

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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  # Read the command line
  uda_name, pt0, pt1, pt2, surf_area, material_id, tlo, thi, tinc = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Plane : '
  print '    P0 =  [', pt0[0], ',', pt0[1], ',', pt0[2], ']; '
  print '    P1 =  [', pt1[0], ',', pt1[1], ',', pt1[2], ']; '
  print '    P2 =  [', pt2[0], ',', pt2[1], ',', pt2[2], ']; '
  print '  Material id : ', material_id
  print '  t_low = ', tlo, ' t_high = ', thi, ' t_inc = ', tinc

  # Set up executables and working directories
  select_exe, extract_exe, work_dir, plot_dir = set_execs_and_dirs()
  print select_exe, extract_exe, work_dir, plot_dir

  # Create the particle list
  error_file = os.path.join(plot_dir, uda_name + ".error")
  partlist_file, part_select = selectpart(select_exe, work_dir, uda_name, plot_dir,
                                           pt0, pt1, pt2, material_id, 
                                           tlo, thi, tinc, error_file)
  #print part_select

  # Compute normal to the plane
  p0 = np.array(pt0)
  p1 = np.array(pt1)
  p2 = np.array(pt2)
  v1 = p1 - p0
  v2 = p2 - p0
  surf_normal = np.cross(v1, v2)
  surf_normal = surf_normal/np.linalg.norm(surf_normal)
  print "Normal = ", surf_normal

  # Extract the particle tractions assuming normals remain fixed
  data = partextract(extract_exe, work_dir, uda_name, plot_dir, 
                     part_select, partlist_file, 
                     surf_area, surf_normal, pt0, error_file)

  # Set up figures
  set_plot_parameters()
  fig1 = plt.figure("traction")
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70, param_text,ha='left',va='top',size='x-small')

  #colors = ['Red', 'Green', 'Blue', 'Magenta', 'Violet']
  #plt_color = colors[idx]
  plt_color = cm.PiYG(float(2)/float(len(data[0])))
  line_style = '-'
  plot_traction('traction', data, plt_color, line_style)

  plt.figure('traction')
  #plt.legend(bbox_to_anchor=(1.05,1), loc=2, prop={'size':10})
  plt_file = os.path.join(plot_dir, uda_name+".fext")
  savePNG(fig1, plt_file, size='1280x960')
  plt.show()


